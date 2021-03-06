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

    real(8) function gammln(x) result(res)
        implicit none
        real(8), intent(in) :: x
        res=log(gamma(x))
        end function

    FUNCTION fact(n)
        integer n,i
        real(8) fact
        fact=1
        do i=1,n
        fact=fact*i
        enddo
        end function

    real(8) function factln(n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(8)::temp
        !temp = PRODUCT((/(i, i=1,n)/))
        res=log(fact(n))
    end function factln
    
    real(8) function factln2(n) result(res) !bug ici
        !the functions factln and binomial are inspired by numerical recipies hosted at
        !http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c6-1.pdf
        !Returns    $\ln(n!)$.
        implicit none
        integer, intent(in) :: n
        real(8), dimension(101) :: a=0        !A static array is automatically initialized to zero.
        if (n < 0) then 
            write(*,*) "Negative factorial in routine factln"
        else if (n==0 .or. n==1) then 
            res=0.0
        !else if (n <= 100) then return a[n] ? a[n] : (a[n]=gammln(n+1.0));    !In range of table.
        else if (n <= 100) then 
            if(a(n)==0) then
                res=a(n)
            else
                res=gammln(n+1.0_8)
            endif
        else
                res=gammln(n+1.0_8);  !Out of range of table.
        endif
    end function

    ! bindist() is the function that returns the values of non-normalized binomial distribution probabilities
    ! the method used to compute them is the one suggested by Newman and Ziff in their paper
    ! which uses a recursive relation
    !
    ! inputs: n,k,p
    ! output: the (non-normalized) k-th coefficient B(n,k,p) of binomial distribution of parameters n and p
    !
    ! the way it works:
    !  - compute k0 = the value of k for which B(n,k,p) is maximal ( k0 is the integer between p*(n+1)-1 and p*(n+1) )
    !  - set the value of B(n,k0,p) to 1	(this is why we need later to properly normalize the coefficients)
    !  - check if k is greater or smaller than k0
    !  - if k < k0 THEN
    !	- compute B(n,k,p) from B(n,k+1,p)
    !    ELSE
    !	- compute B(n,k,p) from B(n,k-1,p)
    real(8) recursive function bindist(n,k,p)  result(b)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: k
        real(8), intent(in) :: p
        integer	            :: k0

        k0=floor(p*(n+1))

        if(k.eq.k0) then
            b=1.0
        else
            if(k.lt.k0) then
                b=bindist(n,k+1,p)*(k+1)/(n-k)*(1-p)/p
            else
                b=bindist(n,k-1,p)*(n-k+1)/k*p/(1-p)
            end if
        end if
        return
    end function bindist


    function bindist_norecu(n,k,p) result(b)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: k
        real(8), intent(in) :: p
        real(8)		    :: b
        integer	            :: k0
        integer		    :: j

        k0=floor(p*(n+1))
        b=1.0

        j=k
        if(k.lt.k0) then
            do while(j.lt.k0)
            b=b*(j+1)/(n-j)*(1-p)/p
            j=j+1
            end do
        else
            do while(j.gt.k0)
            b=b*(n-j+1)/j*p/(1-p)
            j=j-1
            end do
        end if

        return
    end function bindist_norecu

    integer(8) function binomial(n,k) result(res)
        implicit none
        integer, intent(in) :: n,k
        !Returns the binomial coefficient $\binom{n,k}$  as a floating-point number.
        res=floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
        !The floor function cleans up roundoff error for smaller values of n and k
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
!    write(*,*) gammln(1.0_8), 0 !ln(0!))=0
!    write(*,*) gammln(2.0_8), 0 !ln(1!)=0
!    write(*,*) gammln(3._8), 1.4 !ln(2)=1.4
!    write(*,*) gammln(4._8), 0 !ln(6)=1.4
!    write(*,*) factln(10)  !ln(6)=1.4
!    write(*,*) binomial(5,3)
!end program test
