program simu
    USE percolation_tools
    USE omp_lib
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
    integer boundary,j
    real(qp) :: tau_f,p_f
    
    type(cells),allocatable :: cell(:,:,:)

    real(qp) :: p_0,beta,volume,N_f

    real(qp) :: r_max,T1,T2
    integer :: N_bub,N_cells,successes,repeats,N_ff
    logical :: inside
    integer threads
    character(len=20) arg
    integer :: testing=10
    real(qp) :: hubble
    !***********************************************************************!
    ! Get cmd arguments
    ! tau_f          the time to simulate to
    ! p_f            nucleation rate
    ! boundary       volume boundary
    ! repeats        number of repeats
    ! threads        number of threads to run the simulation on
    !_______________________________________________________________________!

    nargs = command_argument_count()
    if (nargs /= 7) then
        call get_command_argument(0,arg)
        print*,"Usage:",trim(arg)," tau_f p_f boundary repeats threads hubble N_final"
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
    read(args(6),*) hubble
    read(args(7),*) N_f
    !***********************************************************************!
    ! Generate the bubbles, check for percolation
    ! Parallelized so each thread runs repeats number of simulations
    ! so final prob is successes/repeats*threads
    !_______________________________________________________________________!


!    repeats=50
    !beta=(8.0_qp*pi*p_f)**(1.0_qp/4.0_qp) !ok
    volume=boundary**3 !ok
 !   N_f=20000
    N_ff=N_f
!    N_ff=20000

    successes=0

    allocate(bubble_array_deb(3*int(N_f)))    
    call omp_set_num_threads(threads)

    !$OMP PARALLEL PRIVATE(bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,seed,mt,mti,test,testing) SHARED(N_ff,boundary,tau_f,p_f,repeats,successes)

    !seed=0
    !hubble=0.1
    
    call sed(seed,mt,mti)
    do i=1,repeats
        !print*,seed,omp_get_thread_num()
        call bubbles(tau_f,p_f,boundary,bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,mt,mti)
        call find_clusters(bubble_array_deb,boundary,r_max,cell,N_bub,N_cells,successes)

        do j=1,N_bub
            deallocate(bubble_array_deb(j)%bubble)
        enddo
        cell(:,:,:)%count=0
    enddo

    !$OMP END PARALLEL

    print *,tau_f,successes!/(1.0_qp*repeats*omp_get_max_threads())!,omp_get_thread_num(),omp_get_max_threads()

    contains


    subroutine bubbles(tau_f,p_0,boundary,bubble_array_deb,r_max,cell,N_f,N_bub,N_cells,mt,mti)
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
        integer:: seed,io
        real(qp) :: V_overlap
        !beta=(8.0_qp*pi*p_0)**(1.0_qp/4.0_qp) !ok

        N_cells=0
        volume=boundary**3 !ok
        N_bub=0

        !N_f=p_0*volume/beta 
        N_f=3*hubble**4/(p_0*volume)
      
        !tau_initial=ieee_value(1.0_qp,ieee_negative_inf)
        tau_initial= -1/hubble
        !if (.not. allocated(bubble_array_deb)) allocate(bubble_array_deb(2*int(N_f)))

        tau1 = next_nucleation_time(tau_initial,N_f,mt,mti) !First nucleation time  
        !print *,tau1
        !r_max=ceiling(2.0_qp/beta*(tau_f-tau1))/2.0_qp !
        r_max=-tau1
        !print*,r_max
        r_max=ceiling((r_max-floor(r_max))*100)+100*floor(r_max)
        !print*,r_max,"after ceilin"
        !Get the boundary and r_max to be divisible with each other so we get even number of cells

        !print*,mod(boundary*100,int(2.0_qp*r_max)),"mod_preloo"
        do while(mod(boundary*100,int(2.0_qp*r_max))/=0) 
            r_max = r_max + 1

            !r_max=ceiling((r_max-floor(r_max))*100)+floor(r_max)*100
           ! print*,r_max, int(2.0_qp*r_max),"r_max,int"
           ! print*,mod(boundary*100,int(2.0_qp*r_max)),r_max,"mod"
            if (r_max > boundary*100) then
                r_max=boundary*100
                exit
            end if 
        enddo
        !print*,"end",r_max,int(2.0_qp*r_max),mod(boundary*100,int(2.0_qp*r_max))
        r_max=r_max/100.0_qp
        N_cells_prev=N_cells
        N_cells=int(boundary/(2.0_qp*r_max))
        if (N_cells==0) N_cells=1
        if (.not. allocated(cell) .or. N_cells_prev < N_cells) then
            if (allocated(cell)) then
                deallocate(cell)
            endif
            allocate(cell(N_cells,N_cells,N_cells))
        endif

!        print*,r_max,N_cells,1/beta*(tau_f-tau1),2*r_max,boundary/(2*r_max)
        !print*,boundary,int(2.0_qp*r_max),r_max,1/beta*(tau_f-tau1),(ceiling(2/beta*(tau_f-tau1))),ceiling(2/beta*(tau_f-tau1))/2.0_qp,N_cells,omp_get_thread_num()
        !First bubble
        location = get_random_location(boundary,mt,mti) 
        call add_new_bubble(bubble_array_deb,location,cell,N_bub,r_max,int(N_ff),N_cells)

        do while(tau1 < tau_f)
            inside = .false.
            tau2=next_nucleation_time(tau1,N_f,mt,mti) !ok
            !print*,tau2,N_bub


            location=get_random_location(boundary,mt,mti)
            V_overlap=0
            do i=1,N_bub 

                bubble_array_deb(i)%bubble%radius=bubble_array_deb(i)%bubble%radius + (tau2-tau1)
                V_overlap=V_overlap+4*pi/3*(bubble_array_deb(i)%bubble%radius)**3
            enddo 
            !print *,tau1,exp(-V_overlap/volume),N_bub!,V_overlap
            if (exp(-V_overlap/volume)<=1.0e-4_qp) then
                exit
            end if
            !print*,bubble_array_deb(N_bub)%bubble%radius

    !            print*,N_bub
            cell_coord = floor(location(:3)/(2*r_max))+1
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
                call add_new_bubble(bubble_array_deb,location,cell,N_bub,r_max,int(N_ff),int(N_cells))
                !print*,tau1,N_bub               

            end if
            tau1=tau2
        end do
        !open(newunit=io,file="coords_expand.dat",status="new",action="write")
        !do i=1,N_bub
            !print*,bubble_array_deb(i)%bubble%coordinates, bubble_array_deb(i)%bubble%radius
        !    write(io,*) bubble_array_deb(i)%bubble%coordinates, bubble_array_deb(i)%bubble%radius
        !enddo
        !close(io)
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
        seed=1
        call sgrnd(seed,mt,mti)
       ! print*,grnd(),seed,omp_get_thread_num()
    end subroutine

    

    function next_nucleation_time(tau1,N_f,mt,mti) result(tau2)
        real(qp) :: tau1,tau2,N_f,r
        integer :: seed
        integer, dimension(0:N-1) :: mt !Had save
        integer                :: mti  
        !print*,seed,omp_get_thread_num(),grnd(),"nucleation"
        
        r=grnd(mt,mti)
        !print*,tau1,N_f,log(r)
        !print*,N_f*log(r)
        !print *,tau1**(-3.0_qp)+N_f*log(r)
        tau2=-(abs(tau1**(-3.0_qp)+N_f*log(r)))**(-1/3.0_qp) !Wrong???
        !print *,"TAU2",tau2
        if (tau2>0.0_qp) then
            print*,"Something wrong tau2>0"
            stop
        end if
    end function

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
