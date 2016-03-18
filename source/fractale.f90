! Algorithm inspired from M. E. J. Newman & R. M. Ziff 2001 paper. 
! algorithmic complexity O(n)

!In this appendix we give a complete program in C for our algorithm for site percolation on a square lattice of N = L × L sites with periodic boundary conditions. This program prints out the size of the largest cluster on the lattice as a function of number of occupied sites n for values of n from 1 to N . The entire program consists of 73 lines of code. First we set up some constants and global variables:
!Below is a simple translation of the functions in fortran

!Site percolation can be handled by only a very slight modification of the algorithm. In this case, sites are occupied in random
!order on the lattice, and each site occupied is declared to be a new cluster of size 1. Then one adds, one by one, all bonds between this site and the occupied sites, if any, adjacent to it, using the algorithm described above.-> how?

! This program contains some additional lines that let you compute the size of the largest non-percolating cluster and the susceptibility

module constants_mcp
    implicit none
    integer ::leon,eric
    integer, parameter 		   :: dp=100		! number of different values of p in interval (0,1) that are considered
	!	other variables
    real(8)		   		   :: delta_p		! distance between consecutive values of p, we choose M equidistant points in the interval (0,1)
  !	n-dependent averages
    real(8), dimension(dp)		   :: perc_prob_p	! probability that there exists a percolating cluster as function of p
    real(8), dimension(dp)		   :: p_infinite_p	! probability that a site belongs to the percolating cluster as function of p
    real(8), dimension(dp)		   :: lc_p		! size of the largest non-percolating cluster as function of p
    real(8), dimension(dp)		   :: susc_p		! susceptibility as function of p
    real(8), dimension(10)		   :: Masse
    integer, parameter :: L=50	  !Linear dimension
    integer, parameter :: nrepet=100 !nombre de répétition
    integer, parameter :: N=L*L
    real, dimension(L**2) :: perc_prob_n=0.0 ! probability that there exists a percolating cluster as function of n=number of occupied sites
    real, dimension(L**2) :: p_infinite_n=0.0	! probability that a site belongs to the percolating cluster (order parameter) as function of n=number of occupied sites
    real(8), dimension(L**2)	:: lcn=0 	    ! size of the largest non-percolating cluster as function of n
    real(8), dimension(L**2)	:: susceptibility=0 ! susceptibility (average size of non-percolating clusters) 
    real(8), dimension(L,L)	:: mass
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

subroutine convolution
implicit none
	real(8)		          :: p
	real(8)			  :: norm
!	integer		   	  :: realization,i,n
	integer		   	  :: i,n
!	integer, dimension(L**2)  :: pp,perc_sites,lc
	real(8), dimension(:), allocatable  :: coeff
	real(8), dimension(:), allocatable :: perc_prob_n1	! probability that there exists a percolating cluster as function of n=number of occupied sites
	real(8), dimension(:), allocatable :: p_infinite_n1	! probability that a site belongs to the percolating cluster (order parameter) as function of n
	real(8), dimension(:), allocatable :: lc_n1		! size of the largest non-percolating cluster as function of n
	real(8), dimension(:), allocatable :: susc_n1		! susceptibility as function of n
!	p-dependent averages

	delta_p=1.0/(dp+1)	! with this choise, the values of p are p_i = i * delta_p , i=1,2,3,..., M
				! you can however decide to choose the values of p differently, not forced to choose them uniformly in (0,1)

	perc_prob_p=0.0
	p_infinite_p=0.0
	lc_p=0.0
	susc_p=0.0

	allocate(perc_prob_n1(L**2))	
	allocate(p_infinite_n1(L**2))
	allocate(lc_n1(L**2))
	allocate(susc_n1(L**2))
	allocate(coeff(L**2))

	open(111,file='susc.res',status='old',action='read')


!		Read values from input file
		do i=1,L**2
			read(111,*) n,perc_prob_n1(i),p_infinite_n1(i),lc_n1(i),susc_n1(i)
		end do

		close(111)
!		Now compute the p-dependent averages Q(p) by convolving the sequence Q(n) with the binomial distribution
		do i=1,dp
			p=i*delta_p
!			compute the coefficients of the summation
!			here we use the method suggested by Newman and Ziff in their paper to compute the binomial distribution in a faster way
!			first we compute the non-normalized binomial coefficients, then we compute the normalization constant
			do n=1,L**2
!				coeff(n)=bindist(L**2,n,p) ! the function bindist returns the values of the (non-normalized) binomial distribution probabilities
				coeff(n)=bindist_norecu(L**2,n,p) ! the function bindist returns the values of the (non-normalized) binomial distribution probabilities
			end do

			do n=1,L**2
				perc_prob_p(i)=perc_prob_p(i)+coeff(n)*perc_prob_n1(n)
				p_infinite_p(i)=p_infinite_p(i)+coeff(n)*p_infinite_n1(n)
				lc_p(i)=lc_p(i)+coeff(n)*lc_n1(n)
				susc_p(i)=susc_p(i)+coeff(n)*susc_n1(n)
			end do

!			compute normalization factor for the probabilities
			norm=0.0
			do n=1,L**2
				norm=norm+coeff(n)
			end do

!			divide the averages by the correct normalization factor
			
			perc_prob_p(i)=perc_prob_p(i)/norm
			p_infinite_p(i)=p_infinite_p(i)/norm
			lc_p(i)=lc_p(i)/norm
			susc_p(i)=susc_p(i)/norm
			susc_p(i)=susc_p(i)/L**2/(p-p_infinite_p(i))
		end do


!		write results into file
		open(222,file='out.dat',status='replace',action='write')

		write(222,*) "# Results for site percolation on square lattice"
		write(222,*) "# L= ", L
		write(222,*) "# number of p points = ", dp
		write(222,*)
		write(222,'(x,f9.6,xx,f9.6,xx,f9.6,xx,f14.3,xx,f14.3)') 0.0,0.0,0.0,0.0,0.0
		do i=1,dp
			write(222,'(x,f9.6,xx,f9.6,xx,f9.6,xx,f14.3,xx,f14.3)') i*delta_p,perc_prob_p(i),p_infinite_p(i),lc_p(i),susc_p(i)
		end do
		write(222,'(x,f9.6,xx,f9.6,xx,f9.6,xx,f14.3,xx,f14.3)') 1.0,1.0,1.0,0.0,0.0
		close(222)

	print*, 'done!'
end subroutine

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

    subroutine permutation
        !Now we generate the random order in which the sites will be occupied, by randomly permuting the integers from 0 to N − 1:
        use random_functions
        implicit none
        integer :: i, j, temp
        do i=1, N
        order(i) = i    !Create a list of all the bonds in any convenient order. Positions in this list are numbered from 1 to M .
        enddo
        do i=1, N   !<- bug might come from here
        j = i + floor((N-i)*rand())   !Choose a number j uniformly at random in the range i ≤ j ≤ M .
	if(j.ne.N) j=j+1	! this operation must be done, because it may happen that the random number generated is 0, then j would be equal to i
!		exchange the positions of i-th and j-th elements of the array order
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
        integer :: i,j,k,a,z,s,s1,s2,r1,r2,big=0
!	integer :: wrapping
	integer :: rmin,r_perc,rlc,perco		! rlc is the root of the largest non-percolating cluster, r_perc the root of the percolating cluster
	integer :: dx,dy,r,c,kk
!	integer, dimension(4)    :: b,r
	integer, dimension(L**2) :: wrapping
   	integer, dimension(L**2) :: pp
        integer, dimension(L**2) :: psites		! psites(n) is the size of the percolating cluster for the configuration with n occupied sites
	integer(kind=8), dimension(L**2) :: lc		! largest non-percolating cluster size
	integer(kind=8), dimension(L**2) :: susc	! susc(n) is the sum of the squares of sizes of non-percolating cluster, in the config. with n occupied sites
	kk=0
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
    if(eric==1 .and. leon==1) then
    eric=eric+1
    mass=0.0
    do a=1,L
		do z=1,L
			s=(a-1)*L+z
			if(ptr(s).ne.empty) then
				r=findroot(s)
				if(wrapping(r)==1) mass(a,z)=1
			endif
		end do
	end do
    call fractal
    endif
  
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

	if(leon==1) then
!	here we print the system configuration
	open(15,file='table.dat')
	write(15,*) "# n= ", i	! <--- we print the number of occupied sites at which percolation occured
	do a=1,L
		do z=1,L
			s=(a-1)*L+z
			if(ptr(s).eq.empty) then
				c=0	! <--- empty sites are represented by zeros
			else
			 	c=1	! <--- occupied sites are represente by 1s
				r=findroot(s)
				if(wrapping(r)==1) c=2 ! <--- if an occupied site belongs to a crossing cluster, we represent it by a 2
			endif
			write(15,advance='no',fmt='(x,i1)') c	
			if(z.eq.L) write(15,*)			! <--- go to next line when we printed the entire row
		end do
	end do
! 	at the end we printed a matrix of 0s, 1s and 2s
    open(16, file='anim.plot')
    write(16,*) 'set autoscale'
    write(16,*) 'set size ratio 1'
    write(16,*) 'set xtics 0.5,1'
    write(16,*) 'set ytics 0.5,1'
    write(16,*) 'set grid front linetype -1'
    write(16,*) 'set format x ""'
    write(16,*) 'set format y ""'
    write(16,*) 'set xrange[-0.5:49.5]'
    write(16,*) 'set yrange[-0.5:49.5]'
    write(16,*) 'unset key'
    write(16,*) 'unset label'
    write(16,*) 'set cbrange [0.0:2.0]'
    write(16,*) 'set palette defined ( 0 ''white'', 1 ''red'', 2 ''blue'')'
    write(16,*) 'unset colorbox'
    write(16,*) 'set term png'
    write(16,'(a,i5.5,a)') "set output 'a",kk,".png'"
    write(16,*) 'plot ''table.dat'' matrix with image'
    close(16)
    call system('gnuplot anim.plot')
	close(15)
	endif

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

    subroutine fractal
        real(8) ::d
        integer :: rayon,i,j

    Masse = 0
    do i=15,35
		do j=15,35
            		d=(i-25)**2+(j-25)**2
			rayon=ceiling (sqrt(d))
            if(mass(i,j)==1 .and. rayon<=10 .and. rayon>0) then
            Masse(rayon)=Masse(rayon)+1
            endif
		end do
	end do
    do i=1,9
    Masse(i+1)=Masse(i)+Masse(i+1)
    enddo
    end subroutine fractal

end module constants_mcp

program main
    use constants_mcp
    use random_functions
    integer :: m
    leon=1
   
    open(unit = 14, file = 'susc.res')
    open(unit = 56, file = 'masse.res')
    moyenne_observable:do m=1,nrepet
    eric=1
    call random_seed() !to init the seed for random number generation
    open(unit = 15, file = 'perc-config.dat',status='replace',action='write') ! <--- file that will contain the percolating configuration
    call boundaries
    call permutation
    call PBCpercolate
!	print*,"sample n. = ",m
    close(15)
    leon=leon+1
    enddo moyenne_observable
    p_infinite_n=p_infinite_n/nrepet
    perc_prob_n=perc_prob_n/nrepet
    lcn=lcn/nrepet
    susceptibility=susceptibility/nrepet
    Masse=Masse/nrepet
    do m=1, N
    write(14,*) m,perc_prob_n(m),p_infinite_n(m),lcn(m),susceptibility(m)
    enddo
    do m=1,10
    write(56,*) m, Masse(m)
    enddo
    close(14)
    close(56)
    print*, 'done!'
    call convolution
    !call system('avconv -r 25    -i a%05d.png -vcodec mjpeg -qscale 1 -y test.avi')
    call system('"rm" *.png')
end program main
