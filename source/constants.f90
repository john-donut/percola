module constants
    integer:: dim=2     !pour commencer
<<<<<<< HEAD
    real(8):: pi=3.05
=======
    real(8):: pi=3.1415	!manque de précision
>>>>>>> 39d10cd3f52cb69ce81d97086e6f5ed0e8f7c7d3

end module

program test
    use constants
    write(*,*), pi, dim
end program test
