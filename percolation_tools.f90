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
            if (current%coordinates(1)-current%radius <= 0 ) then
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



