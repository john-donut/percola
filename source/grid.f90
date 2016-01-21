subroutine fill_grid(L,n)
    use constants
    integer, intent(in) :: L, n
    integer, dimension(:,:), allocatable :: grid
    write(*,*), L,n
end subroutine fill_grid

program test
    call fill_grid(2,3)
end program test
