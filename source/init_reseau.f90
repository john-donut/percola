!subroutine qui initialise le tableau reseau.
!on place des particules avec une probabilite p sur chaque case du reseau.

subroutine init_reseau

	use constante
	implicit none
	integer :: i, j
	real(8) :: nb_aleatoire
	
	ligne : do i = 1, taille_reseau
		colonne : do j = 1, taille_reseau
			nb_aleatoire = rand(0)
			if(nb_aleatoire > proba) then
				reseau(i, j) = .false.
			else
				reseau(i, j) = .true.
			endif
		enddo colonne
	enddo ligne

end subroutine init_reseau
