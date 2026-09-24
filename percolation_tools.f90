module percolation_tools
    USE,intrinsic :: ISO_FORTRAN_ENV, only : qp => real64
    USE,intrinsic :: ieee_arithmetic
  implicit none



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

end module percolation_tools



