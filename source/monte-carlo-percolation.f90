! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. 
! algorithmic complexity O(n)

!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:
!Below is a simple translation of the functions in fortran

!Site percolation can be handled by only a very slight modification of the algorithm. In this case, sites are occupied in random
!order on the lattice, and each site occupied is declared to be a new cluster of size 1. Then one adds, one by one, all bonds between this site and the occupied sites, if any, adjacent to it, using the algorithm described above.-> how?

module constants_mcp
    implicit none
    integer, parameter :: L=1500   !Linear dimension
    integer, parameter :: nrepet=100 !nombre de répétition
    integer, parameter :: N=L*L
    real, dimension(L**2) :: perc_prob_n=0 ! probability that there exists a percolating cluster as function of n=number of occupied sites
    integer, parameter :: EMPTY=(-N-1) !?
    integer, dimension(N) :: ptr, order   !Array of pointers, Nearest neighbors
    integer, dimension(N,4) :: nn      !Occupation order
    integer, dimension(L**2) :: pp 
contains

    integer(8) function dwhere(x,y) result(i)
        integer, intent(in) :: x,y
        i=x+(y-1)*L
    end function

    subroutine boundaries()
        !Next we set up the array nn()() which contains a list of the nearest neighbors of each site. Only this array need be changed in order for the program to work with a lattice of different topology.
        implicit none
        integer :: i,j, x, y

        do i=1, N  
        x=mod(i-1,L)+1
        y=(i-x)/L+1
        nn(i,1)=dwhere((mod(i-2+L,L)+1),y) !gauche
        nn(i,2)=dwhere((mod(i,L)+1),y) !droite
        nn(i,3)=dwhere(x,mod(y,L)+1)
        nn(i,4)=dwhere(x,mod(y-2+L,L)+1)
        !write(*,*), "matrice",i, " " , (nn(i,j),j=1,4)
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
        integer :: i,j,k,s1,s2,r1,r2,nb_fusion,nb_cluster=N,big=0
        integer	:: crx,cry
        integer, dimension(L**2,4) :: touch_border ! array used to determine if a cluster is crossing the system along one of the two directions
        integer, dimension(N) :: ns
        do i=1,L
        touch_border((i-1)*L+1,1)=1 ! sites on right border
        touch_border(i*L,2)=1   ! sites on left border
        touch_border(i,3)=1     ! sites on top border
        touch_border(L*(L-1)+i,4)=1 ! sites on bottom border
        end do
        crx=0
        cry=0
        do i=1, N
        ptr(i) = EMPTY
        enddo
        do i=1, N   !Sites are occupied in the order specified by the array order[]
        nb_fusion = 0
        ns(1) = ns(1)+1
        r1 = order(i)
        s1 = r1
        ptr(s1) = -1    !1. Initially all sites are clusters in their own right. Each is its own root site, and contains a record of its own size, which is 1.
        pp(i)=0
        do j=1, 4
        s2 = nn(s1,j)   !for each occupied neighbour
        if (ptr(s2) /= EMPTY) then

            r2 = findroot(s2)   !The function findroot() is called to find the roots of each of the adjacent sites.
            if (r2 /= r1) then
                nb_fusion = nb_fusion+1      
                if (ptr(r1)>ptr(r2)) then
                    !5. If the two roots nodes are different, we examine the cluster sizes stored in them, and add a pointer from the root of the smaller cluster to the root of the larger, thereby making the smaller tree a subtree of the larger one. If the two are the same size, we may choose whichever tree we like to be the subtree of the other. We also update the size of the larger cluster by adding the size of the smaller one to it.
                    do k=1, 4
                    touch_border(r2,k)=ior(touch_border(r1,k),touch_border(r2,k))
                    end do
                    ns(abs(ptr(r1))) = ns(abs(ptr(r1)))-1
                    ns(abs(ptr(r2))) = ns(abs(ptr(r2)))-1
                    ns(abs(ptr(r1)+ptr(r2))) = ns(abs(ptr(r1)+ptr(r2)))+1
                    ptr(r2) = ptr(r2)+ ptr(r1)
                    ptr(r1) = r2
                    r1 = r2
                else        !4. If the two root sites are the same site, we need do nothing further.
                    ns(abs(ptr(r1))) = ns(abs(ptr(r1)))-1
                    ns(abs(ptr(r2))) = ns(abs(ptr(r2)))-1
                    ns(abs(ptr(r1)+ptr(r2))) = ns(abs(ptr(r1)+ptr(r2)))+1
                    ptr(r1) = ptr(r1)+ ptr(r2)
                    do k=1,4
                    touch_border(r1,k)=ior(touch_border(r1,k),touch_border(r2,k))
                    end do
                    ptr(r2) = r1
                endif
                if (-ptr(r1)>big) then  !vérifie le plus grand cluster
                    big = -ptr(r1)
                endif
            endif
        endif
        !write(10,*), i, i, i+1, big
        enddo
        if(touch_border(r1,1)+touch_border(r1,2).gt.1) then
            crx=1   ! the cluster is crossing the system along direction x
        end if
        if(touch_border(r1,3)+touch_border(r1,4).gt.1) then
            cry=1   ! the cluster is crossing the system along direction y
        end if
        if(crx+cry.ge.1) then
            pp(i)=1
        end if
        pp(i+1)=pp(i)
        !call susceptibilite(i, nb_fusion, nb_cluster)
        enddo
        pp(L**2)=1
        perc_prob_n=perc_prob_n+pp
    end subroutine

    subroutine susceptibilite(i, nb_fusion, nb_cluster)

        integer :: i, nb_fusion, nb_cluster
        real(8) :: moyenne_microcanonique

        nb_cluster = nb_cluster-nb_fusion
        moyenne_microcanonique = 1.*N/nb_cluster
        write(11,*) i, moyenne_microcanonique

    end subroutine susceptibilite

end module constants_mcp

program main
    use constants_mcp
    use random_functions
    integer :: m
    character(len=20) :: filename
    !allocate somewhere for dynamical arrays
    call random_seed() !to init the seed for random number generation

    write (filename, "('clusters',I4.4,'.dat')") L
    open (unit=10,file=filename)
    open(unit = 11, file = 'susceptibilite.res')
    open(unit = 14, file = 'percolation.res')
    moyenne_observable:do m=1,nrepet
    call boundaries
    call permutation
    call percolate
    enddo moyenne_observable
    do m=1, N
    write(14,*) m,perc_prob_n(m)/nrepet
    enddo
    close(14)
    close(10)
    close(11)
end program main
