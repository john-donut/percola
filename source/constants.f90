module constants
    integer:: dim=2     !pour commencer
    real(8):: pi=3.14

end module

program test
    use constants
    write(*,*), pi, dim
end program test
