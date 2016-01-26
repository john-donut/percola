! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. 
! algorithmic complexity O(n)

!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:
!Below is a simple translation of the functions in fortran

!Site percolation can be handled by only a very slight modification of the algorithm. In this case, sites are occupied in random
!order on the lattice, and each site occupied is declared to be a new cluster of size 1. Then one adds, one by one, all bonds between this site and the occupied sites, if any, adjacent to it, using the algorithm described above.-> how?

module constants_mcp
    implicit none
    integer, parameter :: L=10   !Linear dimension
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
        nn(i,1) = mod(i,N)
        nn(i,2) = mod((i+N-2),N)
        nn(i,3) = mod((i+L-1),N)
        nn(i,4) = mod((i+N-L-1),N)
        if (mod(i,L)==0) then
            nn(i,2) = i+L-1
        endif
        if (mod((i+1),L)==0) then   !conditions aux limites périodiques
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
        do i=1, N   !<- bug might come from here
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
        integer, intent(in) ::  i   !This function takes an integer argument, which is the label of a site, and returns the label of the root site of the cluster to which that site belongs.
        
        if(ptr(i)<0) then   !is negative for all root nodes and contains the label of the site pointed to otherwise.
            res=i   !If ptr[i] < 0, “i” is the root of the cluster, and |ptr[i]| gives the number of sites belonging to the cluster.
        else
            !ptr(i) = findroot(ptr(i))
            res = findroot(ptr(i))
        endif
        !path compression: In the find part of the algorithm, trees are traversed to find their root sites. If two initial sites lead to the same root, then they belong to the same cluster. In addition, after the traversal is completed, all pointers along the path traversed are changed to point directly the root of their tree.
    end function

    subroutine percolate
        implicit none
        !when a bond is added to the lattice we must find which clusters the sites at either end belong to. If the two sites belong to different clusters, those clusters must be amalgamated into a single cluster; otherwise, if the two belong to the same cluster, we need do nothing. Algorithms which achieve these steps are known as “union/find” algorithms.
        !1. Initially all sites are clusters in their own right. Each is its own root site, and contains a record of its own size, which is 1.
        !2. Bonds are occupied in random order on the lattice.
        !3. Each bond added joins together two sites. We follow pointers from each of these sites separately until we reach the root sites of the clusters to which they belong. Then we go back along the paths we followed through each tree and adjust all pointers along those paths to point directly to the corresponding root sites.
        !4. If the two root sites are the same site, we need do nothing further.
        !5. If the two root nodes are different, we examine the cluster sizes stored in them, and add a pointer from the root of the smaller cluster to the root of the larger, thereby making the smaller tree a subtree of the larger one. If the two are the same size, we may choose whichever tree we like to be the subtree of the other. We also update the size of the larger cluster by adding the size of the smaller one to it.
        integer :: i,j,s1,s2,r1,r2,big=0
        do i=1, N
        ptr(i) = EMPTY
        enddo
        do i=1, N   !Sites are occupied in the order specified by the array order[]
        r1 = order(i)
        s1 = r1
        ptr(s1) = -1    !1. Initially all sites are clusters in their own right. Each is its own root site, and contains a record of its own size, which is 1.
        do j=1, 4
        s2 = nn(s1,j)   !for each occupied neighbour
        if (ptr(s2) /= EMPTY) then
            r2 = findroot(s2)   !The function findroot() is called to find the roots of each of the adjacent sites.
            if (r2 /= r1) then      
                if (ptr(r1)>ptr(r2)) then
                !5. If the two roots nodes are different, we examine the cluster sizes stored in them, and add a pointer from the root of the smaller cluster to the root of the larger, thereby making the smaller tree a subtree of the larger one. If the two are the same size, we may choose whichever tree we like to be the subtree of the other. We also update the size of the larger cluster by adding the size of the smaller one to it.
                    ptr(r2) = ptr(r2)+ ptr(r1)
                    ptr(r1) = r2
                    r1 = r2
                else        !4. If the two root sites are the same site, we need do nothing further.
                    ptr(r1) = ptr(r1)+ ptr(r2)
                    ptr(r2) = r1
                endif
                if (-ptr(r1)>big) then  !vérifie le plus grand cluster
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
