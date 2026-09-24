program simu
    USE omp_lib
    USE percolation_tools
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
    
    !Stuff pulled from the mtfort90 code which was originally made by Hiroshi Takano and edited by Richard Woloshyn.
    !All credit of the original Mersenne-Twister goes to them and Makoto Matsumoto and Takuji Nishimura.

    !I removed the save functionality of the MT code as it was not central to functionality, where the property was used has been commented
    !In case something was broken due to removing it.

    ! Default seed
    integer, parameter :: defaultsd = 4357
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1
    integer:: test
! the array for the state vector
    integer, dimension(0:N-1) :: mt !Had save
    integer                :: mti = N1 !Had save

    !***********************************************************************!
    !
    ! Type declarations
    ! 
    !_______________________________________________________________________!
    type(bstack),allocatable:: bubble_array_deb(:)

    integer :: istat,seed,un
    integer :: nargs,i
    character(len=10), allocatable :: args(:)
    real(qp) :: pi = 3.141592653589793_qp
    integer boundary,j
    real(qp) :: tau_f,p_f
    
    type(cells),allocatable :: cell(:,:,:)

    real(qp) :: p_0,beta,volume,N_f

    real(qp) :: r_max,T1,T2
    integer :: N_bub,N_cells,successes,repeats
    logical :: inside
    integer threads
    character(len=20) arg
    integer :: testing=10
    integer :: test_var
    !***********************************************************************!
    ! Get cmd arguments
    ! tau_f          the time to simulate to
    ! p_f            nucleation rate
    ! boundary       volume boundary
    ! repeats        number of repeats
    ! threads        number of threads to run the simulation on
    !_______________________________________________________________________!

    nargs = command_argument_count()
    if (nargs /= 5) then
        call get_command_argument(0,arg)
        print*,"Usage:",trim(arg)," tau_f p_f boundary repeats threads"
        stop
    end if
    allocate(args(nargs))

    do i=1,nargs
        call get_command_argument(i,args(i))
    enddo
    read(args(1),*) tau_f
    read(args(2),*) p_f
    read(args(3),*) boundary
    read(args(4),*) repeats
    read(args(5),*) threads
    !***********************************************************************!
    ! Generate the bubbles, check for percolation
    ! Parallelized such that each thread runs repeats number of simulations
    ! so final probability is successes/(repeats*threads)
    !_______________________________________________________________________!


    !Parameters for simulation, based on fixing nucleation rate p_f and boundary size
    beta=(8.0_qp*pi*p_f)**(1.0_qp/4.0_qp) 
    volume=boundary**3 
    N_f=p_f*volume/beta !Final mean number of bubbles


    successes=0

    allocate(bubble_array_deb(3*int(N_f)))    
    call omp_set_num_threads(threads)

    !$OMP PARALLEL PRIVATE(bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,seed,mt,mti,test,testing) SHARED(boundary,tau_f,p_f,repeats,successes)

  !Move to first inside loop

    !For whatever reason seed = 0 breaks everything. 
    do i=1,repeats
        call sed(seed,mt,mti,test_var)
        !Generates the bubbles into the volume
        call bubbles(tau_f,p_f,boundary,bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,mt,mti)
        !Finds clusters and checks for percolation
        call find_clusters(bubble_array_deb,boundary,r_max,cell,N_bub,N_cells,successes)

        !Deallocate for next repeat
        do j=1,N_bub
            deallocate(bubble_array_deb(j)%bubble)
        enddo
        cell(:,:,:)%count=0
    enddo

    !$OMP END PARALLEL

    print *,tau_f,successes!/(1.0_qp*repeats*omp_get_max_threads())!,omp_get_thread_num(),omp_get_max_threads()

    contains
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

    subroutine bubbles(tau_f,p_0,boundary,bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,mt,mti)
        !***********************************************************************!
        ! subroutine for generating the bubbles in a fast transition.
        !
        !   Inputs(not in order):
        !       tau_f, float, final value of tau to simulate to
        !       p_0, float, should be set to 1, the nucleation rate at t_f
        !       boundary,integer, size of the volume boundary
        !       N_f, integer, final mean number of bubbles.
        !
        !   The following should only be supplied to the subroutine but not set beforehand.
        !       bubble_array_deb, type(bstack), array for the bubbles to be placed in, details near beginning of file
        !       r_max, integer, used for setting the size of the cells
        !       cell, type(cells) 3-Dimensional array for representing the cells
        !       N_bub, integer, used for number of bubbles
        !       N_cells, integer, used for number of cells
        !       mt, see Mersenne-Twister 
        !       mti, see Mersenne-Twister
        !       
        !       
        !_______________________________________________________________________!
        type(cells),allocatable :: cell(:,:,:)
        type(bstack),allocatable:: bubble_array_deb(:)
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        real(qp) :: tau_f,p_0,beta,volume,N_f
        real(qp) :: tau_initial,tau1,tau2
        type(bubble), pointer ::neighbor
        real(qp) :: r_max
        real(qp) :: location(4),distance
        integer boundary
        integer :: i
        integer :: N_bub,N_cells
        logical :: inside
        integer :: N_cells_prev
        integer cell_coord(3),x,y,z,j,k,l,h,limits(3,3),t
        integer:: seed
        beta=(8.0_qp*pi*p_0)**(1.0_qp/4.0_qp) 

        N_cells=0
        volume=boundary**3 
        N_bub=0

        N_f=p_0*volume/beta 
    
        tau_initial=ieee_value(1.0_qp,ieee_negative_inf)

        !if (.not. allocated(bubble_array_deb)) allocate(bubble_array_deb(2*int(N_f)))

        tau1 = next_nucleation_time(tau_initial,N_f,mt,mti) !First nucleation time  

        r_max=1/beta*(tau_f-tau1)

        !Get the boundary and r_max to be divisible with each other so we get even number of cells        
        !Not so trivial, we try to minimize the N_cells, this is done by increasing r_max until boundary is divisible by 2*r_max
        !As r_max can be a float, we look at it in order of 10^-3, increasing r_max by one until boundary/(2*r_max) is even.

        r_max=ceiling((r_max-floor(r_max))*100)+100*floor(r_max)


        do while(mod(boundary*100,int(2.0_qp*r_max))/=0) 
            r_max = r_max + 1
            if (r_max > boundary*100) then
                r_max=boundary*100
                exit
            end if 
        enddo

        r_max=r_max/100.0_qp

        N_cells_prev=N_cells
        N_cells=int(boundary/(2.0_qp*r_max))
        if (N_cells==0) N_cells=1 !In the case things fail, for example if the boundary is too small

        !The following is for allocating more cells in the case we have more cells in subsequent repetitions
        if (.not. allocated(cell) .or. N_cells_prev < N_cells) then
            if (allocated(cell)) then
                deallocate(cell)
            endif
            allocate(cell(N_cells,N_cells,N_cells))
        endif

        !Add the first bubble
        location = get_random_location(boundary,mt,mti) 
        call add_new_bubble(bubble_array_deb,location,cell,N_bub,r_max,int(N_f),N_cells)


        do while(tau1 < tau_f)
            inside = .false.
            tau2=next_nucleation_time(tau1,N_f,mt,mti) 
            location=get_random_location(boundary,mt,mti) !New nucleation location

            !Update the size of all bubbles
            do i=1,N_bub 
                bubble_array_deb(i)%bubble%radius=bubble_array_deb(i)%bubble%radius + 1.0_qp/beta*(tau2-tau1)
            enddo 
            cell_coord = floor(location(:3)/(2*r_max))+1

            !Check if new location is inside an existing bubble or not
            limits = search_range(cell_coord,n_cells,.true.) !Limit the search to be within the volume
            do x=1,3
                do y=1,3
                    do z=1,3
                        j=cell_coord(1)+limits(1,x)
                        k=cell_coord(2)+limits(2,y)
                        l=cell_coord(3)+limits(3,z)
                            do t=1,cell(j,k,l)%count
                                neighbor=> cell(j,k,l)%bubblearray(t)%bubble    
                                distance = get_shortest_distance(neighbor%coordinates,location,boundary)
                                if (distance <= neighbor%radius) then
                                    inside = .true.
                                    exit
                                end if
                            enddo

                    enddo
                enddo
            enddo

            
            if (.not. inside) then
                call add_new_bubble(bubble_array_deb,location,cell,N_bub,r_max,int(N_f),int(N_cells))

                !!Following can be used to get number of bubbles as function of tau, useful for comparing against theory.
                !print*,tau1,N_bub               

            end if
            tau1=tau2
        end do
        !!If one wants to get all the locations and radii of the bubbles at time tau, the following can be outputted

        !do i=1,N_bub
        !    print*,bubble_array_deb(i)%bubble%coordinates, bubble_array_deb(i)%bubble%radius
        !enddo

    end subroutine

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
    subroutine sed(seed,mt,mti,test_var)
        integer seed,un,istat,test_var
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        !***********************************************************************!
        !
        ! Seeding from /dev/urandom
        !   seed, integer, initial seed
        !   mt,mti, see Mersenne Twister
        !_______________________________________________________________________!

        open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            print*,"Could not access /dev/urandom/"
        endif 
        seed=test_var
        seed=1
        call sgrnd(seed,mt,mti)

    end subroutine

    subroutine add_new_bubble(bubble_array_deb,bubble,cell,N_bub,r_max,N_f,N_cells)
        !***********************************************************************!
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


    function next_nucleation_time(tau1,N_f,mt,mti) result(tau2)
        !***********************************************************************!
        !   Function for getting the next nucleation time tau2, based on poisson statistics
        !   Inputs:
        !       tau1, real(8), the initial time tau1 to start from 
        !       N_f, integer, final mean number of bubbles
        !       mt,mti, see Mersenne-Twister
        !   Output:
        !       tau2,real(8), the time when next bubble could be nucleated.
        !_______________________________________________________________________!
        real(qp) :: tau1,tau2,N_f,r
        integer :: seed
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  

        r=grnd(mt,mti)
        tau2=log(exp(tau1)-1.0_qp/N_f*log(r))  

    end function

    function get_random_location(boundary,mt,mti) result(random_location)
        !***********************************************************************!
        !   Function for getting the next nucleation time tau2, based on poisson statistics
        !   Inputs:
        !       boundary, integer, size of the boundary
        !       mt,mti, see Mersenne-Twister
        !   Output:
        !       random_location, real(8) 4-Dimensional array, a random location
        !_______________________________________________________________________!
        implicit none
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        real(qp) random_location(4),r
        integer boundary
        integer i

        do i=1,3
            r=grnd(mt,mti)
            random_location(i)=min(boundary*r,real(boundary-1.0e-10_qp,kind=qp))
        end do
        random_location(4)=0.0_qp

    end function
  
    !Rest is the MT code crammed into this sloppily to get it to work parallel, see the original MT code, I have no idea about this myself
    subroutine sgrnd(seed,mt,mti)
        implicit none
    !   
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  !Had save
    !      setting initial seeds to mt[N] using
    !      the generator Line 25 of Table 1 in
    !      [KNUTH 1981, The Art of Computer Programming
    !         Vol. 2 (2nd Ed.), pp102]
    !
        integer, intent(in) :: seed
        mt(0) = iand(seed,-1)
        do mti=1,N-1
          mt(mti) = iand(69069 * mt(mti-1),-1)
        enddo
    !
        return
      end subroutine sgrnd
    real(8) function grnd(mt,mti)
    implicit integer(a-z)
    integer, dimension(0:N-1) :: mt !Had save
    integer                :: mti  
! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                        if sgrnd() has not been called,


        print*,"Error or something I don't know, I just crammed entierity of the mersenne twister into here because"
        print*, "I couldn't get the module variables to be private for each thread, so something might've gone wrong here if this prints out :D"

!                              a default initial seed is used
      endif

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    endif

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      grnd=dble(y)/(2.0d0**32-1.0d0)
    endif

    return
  end function grnd

end program simu
