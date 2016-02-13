module random_functions
    !Module for random number generation and such
    !possible improvements: https://www.ucl.ac.uk/~ucakarc/work/randgen.html
    implicit none

contains

    real(8) function mean_array(A,N) result(M)
        implicit none
        integer, intent(in) ::  N
        real(8), dimension(N), intent(in)   ::  A
        M=sum(A)/N
    end function mean_array

    real(8) function variance_array(A,N) result(V)
        implicit none
        integer, intent(in) ::  N
        real(8), dimension(N), intent(in)   ::  A
        integer :: i
        real(8) :: T=0, M
        M=mean_array(A,N)

        do i=1,N
        T=T+(A(i)-M)**2
        enddo
        V=sqrt(T/(N-1))
    end function variance_array

    real(8) function correlation_length(g,N) result(Xi)
        implicit none
        integer, intent(in) :: N
        real(8), dimension(N), intent(in)   :: g
        !"We define the correlation or the connectivity length by $\xi$ as some average distance of two sites belonging to the same
        !cluster: $\xi^2=\dfrac{\sum _r r^2 g(r)}{\sum _r g(r)}$ "
    end function        

    real(8) function drand() result(r)
        !   Here the function drand() generates a random double precision floating point number between 0 and 1.
        implicit none
        !the precision problem seems solved.
        call random_number(r)
    end function drand

    subroutine init_random_seed()
        !copied from https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
        use iso_fortran_env, only: int64
        !implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t

        call random_seed(size = n)
        allocate(seed(n))
        open(newunit=un,file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
        if(istat == 0) then
            !First try if the OS provides a random number generator
            read(un) seed
            close(un)
        else
            ! Fallback to XOR:ing the current time and pid. The PID is useful in case one launches multiple instances of the same program in parallel.
            call system_clock(t)
            if(t == 0) then
                call date_and_time(values=dt)    
                t=((dt(1)-1970)* 365_int64*24*60*60*1000+dt(2)* 31_int64*24*60*60*1000 +dt(3)* 24_int64*60*60*1000 &
                    +dt(5)* 60*60*1000+dt(6)* 60*1000 +dt(7)* 1000 +dt(8))
            endif
            pid = getpid()
            t=ieor(t,int(pid,kind(t)))
            do i=1,n
            seed(i)=lcg(t)
            enddo
        endif
        call random_seed(put=seed)

    contains
        !This simple PRNG might not be good enough for real work, but is sufficient for seeding a better PRNG.
        function lcg(s)
            integer :: lcg
            integer(int64) :: s
            if (s == 0) then
                s = 104729
            else
                s = mod(s, 4294967296_int64)
            end if
            s = mod(s*279470273_int64,4294967291_int64)
            lcg=int(mod(s,int(huge(0),int64)),kind(0))
        end function lcg
    end subroutine init_random_seed

end module random_functions

!program test
!    use random_functions
!end program test
