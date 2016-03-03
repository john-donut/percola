! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. 
! algorithmic complexity O(n)

!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:
!Below is a simple translation of the functions in fortran

!Site percolation can be handled by only a very slight modification of the algorithm. In this case, sites are occupied in random
!order on the lattice, and each site occupied is declared to be a new cluster of size 1. Then one adds, one by one, all bonds between this site and the occupied sites, if any, adjacent to it, using the algorithm described above.-> how?

! This program contains some additional lines that let you compute the size of the largest non-percolating cluster and the susceptibility

module constants_mcp
    implicit none
    integer, parameter :: L=50	  !Linear dimension
    integer, parameter :: nrepet=100 !nombre de répétition
    integer, parameter :: N=L*L
    real, dimension(L**2) :: perc_prob_n=0.0 ! probability that there exists a percolating cluster as function of n=number of occupied sites
    real, dimension(L**2) :: p_infinite_n=0.0	! probability that a site belongs to the percolating cluster (order parameter) as function of n=number of occupied sites
    real(8), dimension(L**2)	:: lcn=0 	    ! size of the largest non-percolating cluster as function of n
    real(8), dimension(L**2)	:: susceptibility=0 ! susceptibility (average size of non-percolating clusters) 
    integer, parameter :: EMPTY=(-N-1)
    integer, dimension(N) :: ptr, order   
    integer, dimension(N,4) :: nn      		!Array of pointers, Nearest neighbors
    integer, dimension(N,2) :: disp		! displacements
!	disp(i,1) is the displacement, along direction x, of site i from the root site of the cluster he belongs to
!	disp(i,2) is the displacement, along direction y, of site i from the root site of the cluster he belongs to
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
        nn(i,1)=dwhere((mod(i-2+L,L)+1),y) 	!gauche
        nn(i,2)=dwhere((mod(i,L)+1),y) 		!droite
        nn(i,3)=dwhere(x,mod(y,L)+1)		!bottom
        nn(i,4)=dwhere(x,mod(y-2+L,L)+1)	!top
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
        j = i + (N-i)*rand()   !Choose a number j uniformly at random in the range i ≤ j ≤ M .
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


!!! below there is the routine to compute the displacement of a site s from its root, given pointers and displacements of all intermediate sites
recursive function compute_disp(s) result(d)
  implicit none
	integer, intent(in)			:: s
	integer, dimension(2)			:: d
	integer, dimension(2)			:: dparent	
	integer :: parent

	d=0

	if(ptr(s).gt.0) then
		parent=ptr(s)			! 'parent' site of s
		dparent(1)=disp(s,1)	! this quantity is the x-displacement of s from its parent
		dparent(2)=disp(s,2)	! this quantity is the y-displacement of s from its parent
		d=compute_disp(parent)+dparent	! recursive relation
	end if

	return
end function compute_disp

    subroutine PBCpercolate
        implicit none
        !when a bond is added to the lattice we must find which clusters the sites at either end belong to. If the two sites belong to different clusters, those clusters must be amalgamated into a single cluster; otherwise, if the two belong to the same cluster, we need do nothing. Algorithms which achieve these steps are known as “union/find” algorithms.
        !1. Initially all sites are clusters in their own right. Each is its own root site, and contains a record of its own size, which is 1.
        !2. Bonds are occupied in random order on the lattice.
        !3. Each bond added joins together two sites. We follow pointers from each of these sites separately until we reach the root sites of the clusters to which they belong. Then we go back along the paths we followed through each tree and adjust all pointers along those paths to point directly to the corresponding root sites.
        !4. If the two root sites are the same site, we need do nothing further.
        !5. If the two root nodes are different, we examine the cluster sizes stored in them, and add a pointer from the root of the smaller cluster to the root of the larger, thereby making the smaller tree a subtree of the larger one. If the two are the same size, we may choose whichever tree we like to be the subtree of the other. We also update the size of ethe larger cluster by adding the size of the smaller one to it.
        integer :: i,j,k,s1,s2,r1,r2,big=0
!	integer :: wrapping
	integer :: rmin,r_perc,rlc		! rlc is the root of the largest non-percolating cluster, r_perc the root of the percolating cluster
	integer :: dx,dy,s,r
!	integer, dimension(4)    :: b,r
	integer, dimension(L**2) :: wrapping
   	integer, dimension(L**2) :: pp
        integer, dimension(L**2) :: psites		! psites(n) is the size of the percolating cluster for the configuration with n occupied sites
	integer(kind=8), dimension(L**2) :: lc		! largest non-percolating cluster size
	integer(kind=8), dimension(L**2) :: susc	! susc(n) is the sum of the squares of sizes of non-percolating cluster, in the config. with n occupied sites

        do i=1, N
          ptr(i)=empty
	  disp(i,:)=0
	  wrapping(i)=0	! wrapping(s)=1 if s is a root site of a wrapping cluster, wrapping(s)=0 in all other cases
	  pp(i)=0
          psites(i)=0
	  lc(i)=0
	  susc(i)=0
        enddo

	r_perc=empty		! r_perc will contain the position of root site of the percolating cluster (that we choose to be the largest wrapping cluster)
	rlc=empty		! rlc will contain the position of the root of the largest non-percolating cluster
    
        do i=1, N-1  !Sites are occupied in the order specified by the array order[]
        r1 = order(i)
        s1 = r1
        ptr(s1) = -1    !1. Initially all sites are clusters in their own right. Each is its own root site, and contains a record of its own size, which is 1.
        do j=1,4
        s2 = nn(s1,j)   !for each occupied neighbour
        if (ptr(s2) /= EMPTY) then 
            r2 = findroot(s2)   !The function findroot() is called to find the roots of each of the adjacent sites.
	    disp(s2,:)=compute_disp(s2)
	    if(s2.ne.r2) ptr(s2)=r2	! assign r2 (the root) to ptr(s2), so that now s2 points directly towards its root (path compression)

            if (r2 /= r1) then
		if(r2.ne.r_perc) then
			susc(i)=susc(i)-ptr(r2)**2	! <--- here we remove the contribution given by the cluster that s2 belongs to from the summation susc, beacuse this cluster will merge with the cluster that s1 belongs to (this operation must be done only if the cluster that s2 belongs to is not the percolating cluster!)
		end if
                if (ptr(r1)>ptr(r2)) then
		    !5. If the two roots nodes are different, we examine the cluster sizes stored in them, and add a pointer from the root of the smaller cluster to the root of the larger, thereby making the smaller tree a subtree of the larger one. If the two are the same size, we may choose whichever tree we like to be the subtree of the other. We also update the size of the larger cluster by adding the size of the smaller one to it.
                    ptr(r2) = ptr(r2)+ ptr(r1)
                    ptr(r1) = r2
!		    we need also to update the displacement of the root of the smaller cluster
!		    this is easily done if we know the displacement of s1 from its root and the displacement of s2 from its root
!		    in general it is given by  displacement(s2) -displacement(s1) -e
!		    with e the unit vector connecting s1 to s2
		    disp(r1,1)=-1*disp(s1,1)+disp(s2,1)
		    disp(r1,2)=-1*disp(s1,2)+disp(s2,2)
		    select case (j)
			case(1)
			     disp(r1,1)=disp(r1,1)+1 ! s2 is on the left of s1, so we must add a spacing (+1,0) to disp(r1)
			case(2)
			     disp(r1,1)=disp(r1,1)-1 ! s2 is on the right of s1, so we must add a spacing (-1,0) to disp(r1)
			case(3)
			     disp(r1,2)=disp(r1,2)-1 ! s2 is below s1, so we must add a spacing (0,-1) to disp(r1)
			case(4)
			     disp(r1,2)=disp(r1,2)+1 ! s2 is above s1, so we must add a spacing (0,+1) to disp(r1)
		    end select
		    rmin=r1		! we also retain memory of which was the smaller cluster
                    r1 = r2
                else
!		    do the same things as above, but now the roles of r1 and r2 are reversed
                    ptr(r1) = ptr(r1)+ ptr(r2)
                    ptr(r2) = r1
		    disp(r2,1)=-1*disp(s2,1)+disp(s1,1)
		    disp(r2,2)=-1*disp(s2,2)+disp(s1,2)
		    select case (j)
			case(1)
			     disp(r2,1)=disp(r2,1)-1 ! s2 is on the left of s1, so we must add a spacing (-1,0) to disp(r2)
			case(2)
			     disp(r2,1)=disp(r2,1)+1 ! s2 is on the right of s1, so we must add a spacing (+1,0) to disp(r2)
			case(3)
			     disp(r2,2)=disp(r2,2)+1 ! s2 is below s1, so we must add a spacing (0,+1) to disp(r2)
			case(4)
			     disp(r2,2)=disp(r2,2)-1 ! s2 is above s1, so we must add a spacing (0,-1) to disp(r2)
		    end select
		    rmin=r2	
                endif
!		this is where we update the displacements and pointers of sites belonging to the smaller cluster (so we need to know rmin)
		do s=1,L**2
			if(ptr(s).eq.rmin) then
				disp(s,:)=compute_disp(s)
				ptr(s)=r1
			end if
		end do

	    else	! <--- enter this ELSE when the r1=r2, i.e. when the roots of s1 and s2 coincide (which means we are adding a bond between sites belonging to the same cluster)
	   ! we need to check for wrapping condition
		dx=disp(s2,1)-disp(s1,1)
		dy=disp(s2,2)-disp(s1,2)
		select case (j)
			case(1)
			     if(dx.ne.-1) wrapping(r1)=1	
			case(2)
			     if(dx.ne.1) wrapping(r1)=1
			case(3)
			     if(dy.ne.1) wrapping(r1)=1
			case(4)
			     if(dy.ne.-1) wrapping(r1)=1
		end select
            endif
        endif

        enddo
	
	if(wrapping(r1).eq.1) then
		pp(i)=1		! <--- the first time that a wrapping cluster is encountered pp is set to 1 and all values pp(n) with larger n will be equal to 1
		if(-ptr(r1).gt.psites(i)) then	! <--- check if the wrapping cluster is the largest one
!			before we update r_perc, by giving to it the value r1
!			we need to re-introduce the contribution of the cluster r_perc to the summation susc.
!			i.e. we encountered a wrapping cluster that is bigger than the one that we had at the previous step of the cycle
!			so we take the new wrapping cluster as the 'percolating' cluster
!			but also we need to re-include the contribution of the 2nd largest wrapping cluster into summation susc
			if((r_perc.ne.EMPTY)) then
				if((r1.ne.r_perc).and.(ptr(r_perc).ne.r1)) then ! <--- this is done only if the two cluster are in fact distinct (and if there was already a perc. cluster)
				susc(i)=susc(i)+ptr(r_perc)**2
				endif
			end if
			r_perc=r1
		end if
	end if

	if(r_perc.ne.empty) psites(i)=-ptr(r_perc)	! <--- update psites if there is a percolating cluster

!	now we update the size of the largest non-percolating cluster

	if(ptr(s1).ne.r_perc) then	! check if the new site s1 that was added has been attached to the percolating cluster, if yes do nothing
		if(-ptr(r1).gt.lc(i)) then
			rlc=r1
			lc(i)=-ptr(r1)
		end if
		if(s1.ne.r_perc) susc(i)=susc(i)+ptr(r1)**2
	end if

!	but it may happen that the new site s1 that was added has connected the percolating cluster with the largest non-percolating cluster, causing them to merge
!	if this happened then we need to recompute lc from scratch !
!	by exploring all the non-percolating clusters and choosing the largest
	if((ptr(rlc).eq.r_perc).or.(rlc.eq.r_perc)) then
		lc(i)=0
		do s=1,L**2
			if((ptr(s).ne.EMPTY).and.(s.ne.r_perc).and.(-ptr(s).gt.lc(i))) then
				rlc=s
				lc(i)=-ptr(s)
			end if
		end do
	end if



	lc(i+1)=lc(i)
        pp(i+1)=pp(i)
        psites(i+1)=psites(i)
	susc(i+1)=susc(i)

        enddo

!	trivial values
        pp(L**2)=1
        psites(L**2)=L**2
	susc(L**2)=0
	lc(L**2)=0	

        perc_prob_n=perc_prob_n+1.0*pp
        p_infinite_n=p_infinite_n+1.0*psites/L**2
	susceptibility=susceptibility+1.0*susc
	lcn=lcn+1.0*lc

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

    call boundaries
    open(unit = 14, file = 'susc.res')
    moyenne_observable:do m=1,nrepet
!    call random_seed() !to init the seed for random number generation
    call permutation
    call PBCpercolate
!	print*,"sample n. = ",m
    enddo moyenne_observable
    p_infinite_n=p_infinite_n/nrepet
    perc_prob_n=perc_prob_n/nrepet
    lcn=lcn/nrepet
    susceptibility=susceptibility/nrepet
    do m=1, N
    write(14,*) m,perc_prob_n(m),p_infinite_n(m),lcn(m),susceptibility(m)
    enddo
    close(14)
    print*, 'done!'
end program main
