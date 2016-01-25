! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. It
! calculates, with algorithmic complexity O(n) wether or not there is
! percolation.

!The version of the program we have given in Ap-
!pendix A measures the largest cluster size as a function
!of number of occupied sites.

!Below is a simple translation of the functions in fortran
!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:

module constants_mcp
    integer, parameter :: L=128   !Linear dimension
    integer, parameter :: N=L*L
    integer :: EMPTY=(-N-1) !?
    integer, dimension(N) :: ptr, order   !Array of pointers, Nearest neighbors
    integer, dimension(N,4) :: nn      !Occupation order
end module

subroutine  boundaries()
    !Next we set up the array nn()() which contains a list of the nearest neighbors of each site. Only this array need be changed in order for the program to work with a lattice of different topology.
    use constants_mcp
    integer :: i

    do i=0, N
    nn(i,1) = mod((i+1),N)
    nn(i,2) = mod((i+N-1),N)
    nn(i,3) = mod((i+L),N)
    nn(i,4) = mod((i+N-L),N)
    if (mod(i,L)==0) then
        nn(i,2) = i+L-1
    endif
    if (mod((i+1),L)==0) then
        nn(i,1) = i-L+1 
    endif
    enddo
end subroutine

function drand() result(r)
    !Here the function drand() generates a random double precision floating point number between 0 and 1. Used in permutation
    implicit none
    real(4) :: r !bug real4/real8 to be investigated
    call random_number(r)
end function drand

subroutine permutation
    !Now we generate the random order in which the sites will be occupied, by randomly permuting the integers from 0 to N − 1:
    use constants_mcp
    integer :: i, j, temp
    do i=0, N
    order(i) = i
    enddo
    do i=0, N
    j = i + (N-i)*drand()
    temp = order(i)
    order(i) = order(j)
    order(j) = temp
    enddo
end subroutine

recursive function findroot(i) result(res)
!We also define a function which performs the “find" operation, returning the label of the root site of a cluster, as well as path compression.
    use constants_mcp
    integer, intent(in) ::  i
    integer             ::  res
    if(ptr(i)<0) then
        res=i
    else
        ptr(i) = findroot(ptr(i))
    endif
end function

!subroutine percolate
!    integer :: i,j,s1,s2,r1,r2,big=0
!    do i=0, N
!        ptr(i) = EMPTY
!    enddo
!    do i=0, N
!        r1 = s1 = order(i)
!        ptr(s1) = -1
!        do j=0, 4
!            s2 = nn(s1)(j)
!            if (ptr(s2)<>EMPTY) then
!                r2 = findroot(s2)
!                if (r2<>r1) then
!                    if (ptr(r1)>ptr(r2)) then
!                        ptr(r2) += ptr(r1)
!                        ptr(r1) = r2
!                        r1 = r2
!                    else
!                        ptr(r1) += ptr(r2)
!                        ptr(r2) = r1
!                    endif
!                    if (-ptr(r1)>big) then
!                        big = -ptr(r1)
!                    endif
!                endif
!            endif
!        write(*,*), i, i, i+1, big
!        enddo
!    enddo
!end subroutine
!
program main
    !allocate somewhere
    !Consiter using function subroutine init_random_seed() defined here:
    !https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
    !call random_seed(size = n) to init the seed for random number generation

    call boundaries
    call permutation
    !    call percolate
end program
