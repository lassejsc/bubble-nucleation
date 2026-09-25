
program simu
    USE omp_lib
    USE percolation_tools
    USE,intrinsic :: ISO_FORTRAN_ENV, only : qp => real64
    USE,intrinsic :: ieee_arithmetic
    implicit none
       
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
    integer boundary
    real(qp) :: tau_f,p_f
    
    type(cells),allocatable :: cell(:,:,:)

    real(qp) :: p_0,beta,volume,N_f

    real(qp) :: r_max,T1,T2
    integer :: N_bub,N_cells,successes,repeats,j
    logical :: inside
    real(qp) :: phi_c
    integer threads
    character(len=20) arg


    !***********************************************************************!
    ! Get cmd arguments
    ! phi            target fractional volume
    ! r              sphere radius
    ! boundary       volume boundary
    ! repeats        number of repeats
    ! threads        number of threads to run the simulation on
    !_______________________________________________________________________!

    nargs = command_argument_count()
    if (nargs /= 5) then
        call get_command_argument(0,arg)
        print*,"Usage:",trim(arg)," phi r boundary repeats threads"
        stop
    end if
    allocate(args(nargs))

    do i=1,nargs
        call get_command_argument(i,args(i))
    enddo
    read(args(1),*) phi_c
    read(args(2),*) r_max
    read(args(3),*) boundary
    read(args(4),*) repeats
    read(args(5),*) threads
    !***********************************************************************!
    ! Generate the bubbles, check for percolation
    ! Parallelized so each thread runs repeats number of simulations
    ! so final prob is successes/repeats*threads
    !_______________________________________________________________________!


!    repeats=50
    beta=(8.0_qp*pi*p_f)**(1.0_qp/4.0_qp) !ok
    volume=boundary**3 !ok
    successes=0


    call omp_set_num_threads(threads)
    N_f=ceiling(-log(1.0_qp-phi_c)/(4.0_qp/3.0_qp*pi*r_max**(3.0_qp)*boundary**(-3.0_qp)))
    !print*,"init",N_f
 !   allocate(bubble_array_deb(2*int(N_f)))    


    !$OMP PARALLEL PRIVATE(bubble_array_deb,cell,N_bub,N_cells,boundary,seed,mt,mti) SHARED(repeats,successes,r_max,N_f)
    allocate(bubble_array_deb(int(N_f)))    
    read(args(3),*) boundary

    call sed(seed,mt,mti)
    do i=1,repeats
        call bubbles_eq(boundary,bubble_array_deb,cell,N_f,N_bub,N_cells,r_max,mt,mti)
        call find_clusters(bubble_array_deb,boundary,r_max,cell,N_bub,N_cells,successes)
        do j=1,N_bub
            deallocate(bubble_array_deb(j)%bubble)    
        enddo
        cell(:,:,:)%count=0
    enddo
    !$OMP END PARALLEL
    print *,phi_c,successes!,omp_get_thread_num(),omp_get_max_threads()
    contains

    
    subroutine bubbles_eq(boundary,bubble_array_deb,cell,N_f,N_bub,N_cells,r_max,mt,mti)
        type(cells),allocatable :: cell(:,:,:)
        type(bstack),allocatable:: bubble_array_deb(:)
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        real(qp) :: tau_f,p_0,beta,volume,N_f
        real(qp) :: tau_initial,tau1,tau2
        real(qp) :: r_max,rcell
        real(qp) :: location(4),distance
        integer boundary
        integer :: i
        integer :: N_bub,N_cells
        logical :: inside
        integer :: N_cells_prev


        N_cells=0
        volume=boundary**3 !ok
        N_bub=0

        

        rcell=ceiling((r_max-floor(r_max))*100)+100*floor(r_max)
        
        do while(mod(boundary*100,int(2.0_qp*rcell))/=0) 
            rcell = rcell + 1
        
            if (rcell > boundary*100) then
                rcell=boundary*100
                exit
            end if 
        enddo

        rcell=rcell/100.0_qp
        N_cells_prev=N_cells
        N_cells=boundary/(2*int(rcell))

        if (N_cells==0) N_cells=1
        if (.not. allocated(cell) .or. N_cells_prev < N_cells) then
            if (allocated(cell)) then
                deallocate(cell)
            endif
            allocate(cell(N_cells,N_cells,N_cells))
        endif

        do i=1,int(N_f)
            location = get_random_location(boundary,mt,mti) 
            location(4)=r_max
            call add_new_bubble(bubble_array_deb,location,cell,N_bub,rcell,int(N_f),N_cells)
        enddo

        !do i=1,N_bub
        !    print*,bubble_array_deb(i)%bubble%coordinates, bubble_array_deb(i)%bubble%radius
        !enddo
    end subroutine

    subroutine sed(seed,mt,mti)
        integer seed,un,istat
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        !***********************************************************************!
        !
        ! Seeding from /dev/urandom
        ! 
        !_______________________________________________________________________!

        open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            print*,"Could not access /dev/urandom/"
        endif 
        call sgrnd(seed,mt,mti)
    end subroutine


    function get_random_location(boundary,mt,mti) result(random_location)
        implicit none
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        real(qp) random_location(4),r
        integer boundary
        integer i
        integer seed

        !print*,seed,omp_get_thread_num(),grnd()
        
        do i=1,3
            r=grnd(mt,mti)
            random_location(i)=min(boundary*r,real(boundary-1.0e-10_qp,kind=qp))
        end do
        random_location(4)=0.0_qp

    end function


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
