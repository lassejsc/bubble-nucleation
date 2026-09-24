module percolation_tools
    USE,intrinsic :: ISO_FORTRAN_ENV, only : qp => real64
    USE,intrinsic :: ieee_arithmetic
  implicit none


    !***********************************************************************!
    !
    ! Some derived types used
    ! 
    !_______________________________________________________________________!
    type bubble
        real(qp) :: coordinates(3)
        real(qp) :: radius
        integer :: cell_coordinates(3)
        logical :: in_stack=.false.
        logical :: visited_global=.false.
    end type

    type bstack
        type(bubble), pointer :: bubble
    end type bstack

    type cells
        integer :: count=0
        type(bub_array),allocatable:: bubblearray(:)
    end type cells

    type bub_array
        type(bubble),pointer :: bubble
    end type bub_array
contains

    subroutine add_new_bubble(bubble_array_deb,bubble,cell,N_bub,r_max,N_f,N_cells)
      !   Subroutine for adding bubbles to bubble_array_deb
        !   Inputs:
        !       bubble_array_deb, type(bstack), array for the bubbles
        !       bubble, real(8) 4-Array, the bubble to be added
        !       cell, type(cells) 3-Dimensional array, the cells array to be updated
        !       N_bub, integer, number of bubbles to be updated
        !       r_max, float, maximal size of a bubble
        !       N_f, integer, final mean number of bubbles
        !       N_cells, integer, number of cells in one direction.
        !_______________________________________________________________________!
        real(qp) :: bubble(4)
        type(cells) :: cell(:,:,:)
        integer :: cell_coords(3)
        real(qp) :: r_max
        type(bstack),allocatable :: bubble_array_deb(:)
        integer :: N_bub,N_f,N_cells
        N_bub=N_bub+1

        allocate(bubble_array_deb(N_bub)%bubble)
       
        bubble_array_deb(N_bub)%bubble%coordinates = bubble(:3)
        bubble_array_deb(N_bub)%bubble%radius=bubble(4)        
        bubble_array_deb(N_bub)%bubble%in_stack=.false.
        bubble_array_deb(N_bub)%bubble%visited_global=.false.
        cell_coords=floor(bubble(:3)/(2*r_max))+1

        bubble_array_deb(N_bub)%bubble%cell_coordinates=cell_coords
        !Allocating the cell array if not allocated
        if (.not. allocated(cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray)) then
            allocate(cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray(ceiling(3.0_qp*N_f/N_cells)))  !Consumes ton of time
        endif
        cell(cell_coords(1),cell_coords(2),cell_coords(3))%count=cell(cell_coords(1),cell_coords(2),cell_coords(3))%count+1
 
        !Update the pointer in the cell%bubblearray which is array of pointers :)
        cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray(cell(cell_coords(1),cell_coords(2),cell_coords(3))%count)%bubble => bubble_array_deb(N_bub)%bubble
        
    end subroutine
 subroutine find_clusters(bubble_array_deb,boundary,r_max,cell,N_f,N_cells,successes)
        !***********************************************************************!
        ! Cluster finding algorithm.
        ! Inputs are the array of bubbles as type bstack, etc,
        ! 
        ! Updates the integer successes if system percolates.
        !_______________________________________________________________________!

        real(qp) volume
        integer:: successes
        real(qp) pi,r
        integer boundary,i,x,y,z,n_cells,cell_coords(3),j,k,l,N_f,t,cluster_count
        real(qp) :: random_loc(4),r_max
        type(cells), allocatable :: cell(:,:,:)
        type(cells),allocatable :: clusters(:)
        !The resulting cluster NOTE! the cells type is used here too, but misleading naming but it is just type(bstack) and integer count for number of bubbles
        type(cells), allocatable :: visited_local  
        integer :: limits(3,3)
        type(bstack),allocatable:: bubble_array_deb(:)
        type(bstack) ,allocatable :: stack(:)
        type(bubble), pointer :: current_bubble,neighbor
        integer :: stack_counter


        volume = boundary**3

        cluster_count=0
        if (.not. allocated(stack)) then
            allocate(stack(3*N_f),clusters(3*N_f))
            allocate(visited_local)
            allocate(visited_local%bubblearray(3*N_f))
        endif


       
        do i=1,N_f
            if (bubble_array_deb(i)%bubble%visited_global .eqv. .false.) then
                stack_counter=1
                stack(stack_counter)%bubble=>bubble_array_deb(i)%bubble !Add first bubble to stack
                visited_local%count=0

                do while(stack_counter /= 0)
                    
                    !Remove from the top of the stack and set as current bubble 
                    current_bubble => stack(stack_counter)%bubble 
                    stack_counter=stack_counter-1

                    current_bubble%in_stack = .false.
                    current_bubble%visited_global=.true.

                    visited_local%count=visited_local%count+1
                    !Adds current bubble to the visited bubbles which, in the end, is the cluster
                    visited_local%bubblearray(visited_local%count)%bubble => current_bubble  


                    cell_coords=current_bubble%cell_coordinates
                    limits = search_range(cell_coords,n_cells,.false.) !Limit the search to be within the volume
                    do x=limits(1,1),limits(1,3)
                        do y=limits(2,1),limits(2,3)
                            do z=limits(3,1),limits(3,3)
                                j=cell_coords(1)+x
                                k=cell_coords(2)+y
                                l=cell_coords(3)+z

                                do t=1,cell(j,k,l)%count
                                    !We check each bubble in each surrounding cell

                                    neighbor=> cell(j,k,l)%bubblearray(t)%bubble    
                                    if ((.not. neighbor%visited_global) .and. (.not. neighbor%in_stack)) then
                                        if (norm2(current_bubble%coordinates-neighbor%coordinates) <= current_bubble%radius+neighbor%radius) then
                                            !If bubbles overlap this executes
            
                                            neighbor%in_stack=.true.
                                            neighbor%visited_global=.true.

                                            stack_counter=stack_counter+1
                                            stack(stack_counter)%bubble => neighbor

                                        end if
                                    end if
                                    
                                enddo

                            enddo
                        enddo
                    enddo
                enddo 
    
                if (percolation_checker(visited_local,boundary)) then
                    successes=successes+1
                    !!If one wants to get the percolating cluster, the following loop can be uncommented.

                    !do t=1,visited_local%count
                    !    print*, visited_local%bubblearray(t)%bubble%coordinates,visited_local%bubblearray(t)%bubble%radius
                    !end do
                    
                    exit
                endif

            end if
        enddo
   

    end subroutine

    function search_range(cell_coords,n_cells,wrap) result(limits)
        !***********************************************************************!
        ! function used for limiting the range of cell searched to only the 
        ! nearby ones
        !
        !   Inputs:
        !       cell_coords, an integer 3-array of the cell coordinates
        !       n_cells, integer, number of cells in one direction.
        !_______________________________________________________________________!        
        integer :: cell_coords(3)
        integer :: n_cells
        logical :: wrap
        integer :: limits(3,3)
        integer i
        limits(:,1)=-1
        limits(:,2)=0
        limits(:,3)=1
        do i=1,3
            if (cell_coords(i)==1) then
              if (.not. wrap) then
                limits(i,1)=0
              else 
                limits(i,1)=n_cells-cell_coords(i) 
              endif
            endif
            if (cell_coords(i)==n_cells) then
              if (.not. wrap) then
                limits(i,3)=0
              else 
                limits(i,3)=1-cell_coords(i) 
              endif
            endif
        enddo 
    end function

    real(qp) function get_shortest_distance(current_bubble,neighbor,boundary) result(dist)
        !***********************************************************************!
        !   Function for getting the shortest distance in periodic boundary conditions between two bubbles
        !   Inputs:
        !       current_bubble, real(8) 4-Dim array, bubble1
        !       neighbor, real(8) 4-Dim array, bubble2
        !       boundary, integer, size of the boundary
        !       
        !   Output:
        !       dist, real(8), distance between current_bubble and neighbor
        !_______________________________________________________________________!
        implicit none
        real(qp) current_bubble(:),neighbor(:),distance(3)
        integer boundary,i

        distance=[( abs(current_bubble(i)-neighbor(i)), i=1,3 )]

        if (distance(1) > boundary/2.0_qp) then 
            distance(1) = boundary - distance(1)
        end if
        if (distance(2) > boundary/2.0_qp) then
            distance(2) = boundary - distance(2)
        end if    
        if (distance(3) > boundary/2.0_qp) then
            distance(3) = boundary - distance(3)
        end if

        dist = sqrt(distance(1)**2.0_qp+distance(2)**2.0_qp+distance(3)**2.0_qp)

    end function

    function percolation_checker(cluster,boundary) result(bool)
        !***********************************************************************!
        ! function for checking for percolationg when given a cluster of bubbles
        ! 
        !   Inputs:
        !       cluster, type(cells) (1-Dimensional !!), array for the cluster
        !        boundary, integer, size of the boundary
        !   Output:
        !        bool, boolean value for whether percolation occured
        !_______________________________________________________________________!        
        type(cells),allocatable :: cluster
        type(bubble), pointer :: current, opposite
        integer i,j,L,boundary
        logical bool,xminperc,xmaxperc,yminperc,ymaxperc,zminperc,zmaxperc
        bool=.false.
        L=boundary
        xminperc=.false.
        xmaxperc=.false.
        yminperc=.false.
        ymaxperc=.false.
        zminperc=.false.
        zmaxperc=.false. 
        


        do i=1,cluster%count
            current =>cluster%bubblearray(i)%bubble

            !If overlap happens at cooridnate zero, make the boundary at L smaller by the amount of overlap
            if (current%coordinates(1)-current%radius <= 0 ) then !TODO: add not here, and test whether it makes things faster (no other effects)
                xminperc=.true.
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(1)+opposite%radius >= L-(current%radius-current%coordinates(1)) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates, L)<current%radius+opposite%radius) then

                        xmaxperc=.true.
                        exit
                    end if
                enddo
            end if

            !Same expect opposite, if we intersect boundary at L, move the boundary at zero right by the amount of overlap
            if (current%coordinates(1)+current%radius >= L) then
                xmaxperc=.true.
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(1)-opposite%radius <= (current%radius+current%coordinates(1)-L) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates,L)<current%radius+opposite%radius) then

                        xminperc=.true. 
                        exit
                    end if
                enddo
            end if

            !It is mandatory due to periodic boundaries to do the check both ways in this manner, as the check is not symmetric , the first bubble can overlap either boundary
            ! while no bubble without periodic boundaries overlap the boundary on the opposite side

            !Same repeated for all axes. 

            if (current%coordinates(2)-current%radius <= 0 ) then
                yminperc=.true.
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(2)+opposite%radius >= L-(current%radius-current%coordinates(2)) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates,L)<current%radius+opposite%radius) then
                        ymaxperc=.true.
                        exit
                    end if
                enddo
            end if

            if (current%coordinates(2)+current%radius >= L ) then
                ymaxperc=.true.
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(2)-opposite%radius <= (current%radius+current%coordinates(2)-L) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates,L)<current%radius+opposite%radius) then

                        yminperc=.true. 
                        exit
                    end if
                enddo
            end if



            if (current%coordinates(3)-current%radius <= 0 ) then
                zminperc=.true.
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(3)+opposite%radius >= L-(current%radius-current%coordinates(3)) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates,L)<current%radius+opposite%radius) then
                        zmaxperc=.true.
                        exit
                    end if
                enddo
            end if

            if (current%coordinates(3)+current%radius >= L ) then
                zmaxperc=.true.
                
                do j=1,cluster%count
                    opposite => cluster%bubblearray(j)%bubble
                    if (opposite%coordinates(3)-opposite%radius <= (current%radius+current%coordinates(3)-L) &
                    .and. get_shortest_distance(opposite%coordinates,current%coordinates,L)<current%radius+opposite%radius) then
                        
                        zminperc=.true. 
                        exit
                    end if
                enddo
            end if


            !If cluster wraps around in at least one coordinate direction.
            if ((xmaxperc .and. xminperc) .or. (ymaxperc .and. yminperc) .or. (zmaxperc .and. zminperc) ) then
                 bool=.true.
                 exit
            endif

        enddo

    end function
end module percolation_tools



