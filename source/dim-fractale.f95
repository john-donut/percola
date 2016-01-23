subroutine dim_fractale(nom,n,result_dimension)
!Calcule la dimension box-counting (ou Minkowski–Bouligand) des diagrammes de phases obtenus précédemment.

    implicit none
    real, dimension(:), allocatable     :: a,b,c,E1,E2,E12
    character (len=40), intent(in)      :: nom 
    integer                             :: i, j, k, P 
    integer, intent(in)                 :: n
    real(8)                             :: max_x, max_y, min_x, min_y,dx, dy, eps
    real(8), intent(inout)              :: result_dimension
    integer, dimension(:,:), allocatable:: quadrillage

allocate(a(n))!contiennent les donées, temps
allocate(b(n))!theta
allocate(c(n))!d/dt theta
allocate(E1(n))!Energies, inutiles mais par compatilibité avec les fichiers
allocate(E2(n))
allocate(E12(n))

open(1,file=nom)
read(1,*) a,b,c,E1,E2,E12
close(1)
deallocate(E1,E2,E12)!inutiles ici

max_x = maxval(b) !theta
max_y = maxval(c) !dtheta/dt
min_x = minval(b)
min_y = minval(c)

dx = ((max_x-min_x)/n)
dy = ((max_y-min_y)/n)

allocate(quadrillage(n,n))
quadrillage=0
P=0

do i=1, n
    j = floor((b(i)-min_x)/dx)+1 !selon x, la case du tableau
    k = floor((c(i)-min_y)/dx)+1 !selon y, la case du tableau
    if (j < n) then 
        if (k < n) then
            if (quadrillage(j, k)==0) then
                quadrillage(j,k)=1;   P=P+1
            endif
        endif
    endif
enddo

eps = dx*dy !surface élémentaire rectangulaire
result_dimension = log(1.*P)/log(1/eps)
write(*,*) "Dimension fractale estimee=", result_dimension

deallocate(a,b,c,quadrillage)
end subroutine dim_fractale

program test
    real(8) :: dimension_test
    dimension_test = 0
    call dim_fractale("attracteurs/attracteur-300K.dat", 10000, dimension_test)
end program test

