! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. 
! algorithmic complexity O(n)

!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:
!Below is a simple translation of the functions in fortran

module constants_mcp
    implicit none
    integer, parameter :: L=128   !Linear dimension
    integer, parameter :: N=L*L
    integer, parameter :: EMPTY=(-N-1) !?
    integer, dimension(N) :: ptr, order   !Array of pointers, Nearest neighbors
    integer, dimension(N,4) :: nn      !Occupation order
contains

    subroutine boundaries()
        !Next we set up the array nn()() which contains a list of the nearest neighbors of each site. Only this array need be changed in order for the program to work with a lattice of different topology.
        implicit none
        integer :: i

        do i=1, N
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

    subroutine permutation
        !Now we generate the random order in which the sites will be occupied, by randomly permuting the integers from 0 to N − 1:
        use random_functions
        implicit none
        integer :: i, j, temp
        do i=1, N
        order(i) = i    !Create a list of all the bonds in any convenient order. Positions in this list are numbered from 1 to M .
        enddo
        do i=1, N
        j = i + (N-i)*drand()   !Choose a number j uniformly at random in the range i ≤ j ≤ M .
        temp = order(i)         !Exchange the bonds in positions i and j. (If i = j then nothing happens.)
        order(i) = order(j)
        order(j) = temp
        enddo
        !The first step in our algorithm is to decide an order in which the bonds will be occupied. We wish to choose this order uniformly at random from all possible such orders, i.e., we wish to choose a random permutation of the bonds.
    end subroutine

    recursive integer(8) function findroot(i) result(res)
        !We also define a function which performs the “find" operation, returning the label of the root site of a cluster, as well as path compression.
        implicit none
        integer, intent(in) ::  i
        !integer             ::  res
        if(ptr(i)<0) then
            res=i   !If ptr[i] < 0, “i” is the root of the cluster, and |ptr[i]| gives the number of sites belonging to the cluster.
        else
            ptr(i) = findroot(ptr(i))
        endif
    end function

    subroutine percolate
        implicit none
        !when a bond is added to the lattice we must find which clusters the sites at either end belong to. If the two sites belong to different clusters, those clusters must be amalgamated into a single cluster; otherwise, if the two belong to the same cluster, we need do nothing. Algorithms which achieve these steps are known as “union/find” algorithms.
        integer :: i,j,s1,s2,r1,r2,big=0
        do i=1, N
        ptr(i) = EMPTY
        enddo
        do i=1, N
        r1 = order(i)
        s1 = r1
        ptr(s1) = -1
        do j=1, 4
        s2 = nn(s1,j)
        if (ptr(s2) /= EMPTY) then
            r2 = findroot(s2)
            if (r2 /= r1) then
                if (ptr(r1)>ptr(r2)) then
                    ptr(r2) = ptr(r2)+ ptr(r1)
                    ptr(r1) = r2
                    r1 = r2
                else
                    ptr(r1) = ptr(r1)+ ptr(r2)
                    ptr(r2) = r1
                endif
                if (-ptr(r1)>big) then
                    big = -ptr(r1)
                endif
            endif
        endif
        write(*,*), i, i, i+1, big
        enddo
        enddo
    end subroutine

end module constants_mcp

program main
    use constants_mcp
    use random_functions
    !allocate somewhere for dynamical arrays
    call random_seed() !to init the seed for random number generation

    call boundaries
    call permutation
    call percolate
end program
