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
        integer :: image(3)=0  
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

    subroutine add_new_bubble(bubble_array_deb,bubble,cell,largest_bubbles,N_bub,r_max,N_f,N_cells,r_final)
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
        type(cells) :: cell(:,:,:),largest_bubbles
        integer :: cell_coords(3)
        real(qp) :: r_max,r_refine,r_final
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

        if (r_final>=r_max) then 
          if (.not. allocated(largest_bubbles%bubblearray)) then
            allocate(largest_bubbles%bubblearray(N_f))
          endif

          largest_bubbles%count=largest_bubbles%count+1 
          largest_bubbles%bubblearray(largest_bubbles%count)%bubble=>bubble_array_deb(N_bub)%bubble
          return
        endif
        
        !Allocating the cell array if not allocated
        if (.not. allocated(cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray)) then
            allocate(cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray(ceiling(3.0_qp*N_f/N_cells)))  !Consumes ton of time
        endif
        cell(cell_coords(1),cell_coords(2),cell_coords(3))%count=cell(cell_coords(1),cell_coords(2),cell_coords(3))%count+1
 
        !Update the pointer in the cell%bubblearray which is array of pointers :)
        cell(cell_coords(1),cell_coords(2),cell_coords(3))%bubblearray(cell(cell_coords(1),cell_coords(2),cell_coords(3))%count)%bubble => bubble_array_deb(N_bub)%bubble
        
    end subroutine

 subroutine find_clusters(bubble_array_deb,largest_bubbles,boundary,r_max,cell,N_f,N_cells,successes)
        !***********************************************************************!
        ! Cluster finding algorithm.
        ! Inputs are the array of bubbles as type bstack, etc,
        ! 
        ! Updates the integer successes if system percolates.
        !_______________________________________________________________________!

        real(qp) volume
        integer:: successes
        integer :: shift(3)
        real(qp) pi,r
        integer boundary,i,x,y,z,n_cells,cell_coords(3),j,k,l,N_f,t,cluster_count
        real(qp) :: random_loc(4),r_max
        type(cells), allocatable :: cell(:,:,:),largest_bubbles
        type(cells),allocatable :: clusters(:)
        !The resulting cluster NOTE! the cells type is used here too, but misleading naming but it is just type(bstack) and integer count for number of bubbles
        type(cells), allocatable :: visited_local  
        integer :: limits(3,3),cell_bounds(3,2)
        type(bstack),allocatable:: bubble_array_deb(:)
        type(bstack) ,allocatable :: stack(:)
        type(bubble), pointer :: current_bubble,neighbor
        integer :: stack_counter,root_indx

        cell_bounds(:,1)=1
        cell_bounds(:,2)=N_cells
!        print *,n_cells,boundary
        volume = boundary**3

        cluster_count=0
        if (.not. allocated(stack)) then
            allocate(stack(3*N_f),clusters(3*N_f))
            allocate(visited_local)
            allocate(visited_local%bubblearray(3*N_f))
        endif


       
        do i=1,N_f
                if (bubble_array_deb(i)%bubble%visited_global) cycle
                stack_counter=1
                stack(stack_counter)%bubble=>bubble_array_deb(i)%bubble !Add first bubble to stack
                visited_local%count=0
                bubble_array_deb(i)%bubble%image=0
                do while(stack_counter /= 0)
                    
                    !Remove from the top of the stack and set as current bubble 
                    current_bubble => stack(stack_counter)%bubble 
                    stack_counter=stack_counter-1

                    ! current_bubble%in_stack = .false.
                    current_bubble%visited_global=.true.

                    !visited_local%count=visited_local%count+1
                    !Adds current bubble to the visited bubbles which, in the end, is the cluster
                    !visited_local%bubblearray(visited_local%count)%bubble => current_bubble  


                    cell_coords=current_bubble%cell_coordinates
                    !limits = search_range(cell_coords,n_cells,.true.) !Limit the search to be within the volume
                    !
                    do t=1,largest_bubbles%count
                    
                    neighbor=> largest_bubbles%bubblearray(t)%bubble    
                    shift=nint((current_bubble%coordinates-neighbor%coordinates)/real(boundary,qp)) 
                    if (associated(neighbor,current_bubble)) cycle
                    if (get_shortest_distance(current_bubble%coordinates,neighbor%coordinates,boundary) <= current_bubble%radius+neighbor%radius) then
                      if ((.not. neighbor%visited_global)) then !.and. (.not. neighbor%in_stack)) then
                            neighbor%image=current_bubble%image+shift  
                            neighbor%visited_global=.true.
                            stack_counter=stack_counter+1
                            stack(stack_counter)%bubble => neighbor
                      else if (any(current_bubble%image+shift/= neighbor%image)) then
                            !$omp atomic
                            successes=successes+1
                            return
                        endif
                    end if
                    
                    enddo


                    do x=-1,1
                      do y=-1,1
                        do z=-1,1
                                j=modulo(cell_coords(1)-1+x,N_cells)+1
                                k=modulo(cell_coords(2)-1+y,N_cells)+1
                                l=modulo(cell_coords(3)-1+z,N_cells)+1
                                 
                                do t=1,cell(j,k,l)%count
                                    
                                    !We check each bubble in each surrounding cell

                                    neighbor=> cell(j,k,l)%bubblearray(t)%bubble    
                                    shift=nint((current_bubble%coordinates-neighbor%coordinates)/real(boundary,qp)) 
                                    ! shift(2)=floor((cell_coords(2)+y)/(1.0_qp*N_cells))
                                    ! shift(3)=floor((cell_coords(3)+z)/(1.0_qp*N_cells))
                                    if (associated(neighbor,current_bubble)) cycle
                                    
                                    if (get_shortest_distance(current_bubble%coordinates,neighbor%coordinates,boundary) <= current_bubble%radius+neighbor%radius) then
                                 
                                    !if (norm2(neighbor%coordinates-current_bubble%coordinates-shift*real(boundary,qp)) <= current_bubble%radius+neighbor%radius) then

                                      if ((.not. neighbor%visited_global)) then !.and. (.not. neighbor%in_stack)) then
                                        
                                            !If bubbles overlap this executes
                                            neighbor%image=current_bubble%image+shift  
                                            

                                            ! neighbor%in_stack=.true.
                                            neighbor%visited_global=.true.
                                            ! neighbor%visited_local=.true.

                                            stack_counter=stack_counter+1
                                            stack(stack_counter)%bubble => neighbor

                                      else if (any(current_bubble%image+shift/= neighbor%image)) then
                                            !$omp atomic
                                            successes=successes+1
                                            return
                                        endif
                                    end if
                                    
                                enddo

                            enddo
                        enddo
                    enddo
                enddo 
    
                if (percolation_checker(visited_local,boundary,cell_bounds)) then
                    print*,"Why are we here?"
                    CALL ABORT
                    !successes=successes+1
                    !!If one wants to get the percolating cluster, the following loop can be uncommented.

                    !do t=1,visited_local%count
                    !    print*, visited_local%bubblearray(t)%bubble%coordinates,visited_local%bubblearray(t)%bubble%radius
                    !end do
                    
                    exit
                endif

        enddo
   

    end subroutine

    ! function search_range(cell_coords,n_cells,wrap) result(limits)
    !     !***********************************************************************!
    !     ! function used for limiting the range of cell searched to only the 
    !     ! nearby ones
    !     !
    !     !   Inputs:
    !     !       cell_coords, an integer 3-array of the cell coordinates
    !     !       n_cells, integer, number of cells in one direction.
    !     !_______________________________________________________________________!        
    !     integer :: cell_coords(3)
    !     integer :: n_cells
    !     logical :: wrap
    !     integer :: limits(3,3)
    !     integer i
    !     limits(:,1)=-1
    !     limits(:,2)=0
    !     limits(:,3)=1
    !     do i=1,3
    !         if (cell_coords(i)==1) then
    !           if (.not. wrap) then
    !             limits(i,1)=0
    !           else 
    !             limits(i,1)=n_cells-cell_coords(i) 
    !           endif
    !         endif
    !         if (cell_coords(i)==n_cells) then
    !           if (.not. wrap) then
    !             limits(i,3)=0
    !           else 
    !             limits(i,3)=1-cell_coords(i) 
    !           endif
    !         endif
    !     enddo 
    ! end function

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
      function percolation_checker(cluster,boundary,cell_bounds) result(bool)!,cell_min,cell_max) result(bool)
        type(cells),allocatable :: cluster
        integer :: cell_min,cell_max !Min and max of possible cell coordinates
        type(bubble), pointer :: current, opposite
        integer i,j,k,L,boundary,indx
        logical bool!,perc(3)       
        integer :: cell_bounds(3,2)
        bool=.false.
        ! perc(:)=.false.
        do i=1,cluster%count
          current =>cluster%bubblearray(i)%bubble
          do k=1,3 
            indx=0
            if (current%cell_coordinates(k)==cell_bounds(k,1)) then !lower bound
              indx=2
            endif
            if (current%cell_coordinates(k)==cell_bounds(k,2)) then !upper bound
              indx=1
            endif

            if (indx==0) then
                
              cycle
            endif

            do j=1,cluster%count
              if (i .ne. j) then
                opposite=>cluster%bubblearray(j)%bubble
                  
            !    print *,current%cell_coordinates(k),opposite%cell_coordinates(k)
                  if (opposite%cell_coordinates(k)==cell_bounds(k,indx)&
                    .and. get_shortest_distance(current%coordinates,opposite%coordinates,boundary)<=current%radius+opposite%radius &
                    ) then
                      ! perc(k)=.true.
                      bool=.true.
                      exit
                  endif 
                
              end if
            enddo

          enddo 
        enddo

    end function

  end module percolation_tools



