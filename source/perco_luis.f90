!subroutine verifiant s'il existe un chemin percolant
!utilise l'algorithme de luis

subroutine perco_luis(percolation_ouinon)

	use constante
	implicit none
	integer, dimension(taille_reseau, taille_reseau) :: label
	integer :: i, j, i_verification, max_verification
	logical, intent(inout) :: percolation_ouinon
	
	percolation_ouinon = .false.
	label = 0
	!s'il y a une particule sur la case (i,j) alors on l'etiquette 1 sur la case sinon on met 0
	do i = 1, taille_reseau
		do j = 1, taille_reseau
			if(reseau(i, j) == .true.) then
				label(i, j) = 1
			else
				label(i, j) = 0
			endif
		enddo
	enddo
	
	!on remplace tous les 1 de la premiere colonne par des 2
	colonne_1 : do i = 1, taille_reseau
		if(label(i, 1) == 1) label(i, 1) = 2
	enddo colonne_1
	
	!on parcourt le reseau plusieurs fois pour construire le chemin percolant
	colonne : do i = 2, taille_reseau
		colonne_verification : do i_verification = 2, max_verification !cette boucle en plus permet de ne pas oublier certain chemin complique (retour en arriere, etc)
			ligne : do j = 1, taille_reseau
				if(label(i, j) == 2) then !si la case(i, j) admet une "bonne" particule
					
					if(i == taille_reseau) then !si on est sur le bord droit du reseau
						if(label(i, j-1) == 1) r(i, j-1) = 2
						if(label(i, j+1) == 1) r(i, j+1) = 2
						if(label(i-1, j) == 1) r(i-1, j) = 2
					
					elseif(j == 1) then !si on est sur le bord superieur du reseau
						if(label(i, j+1) == 1) r(i, j+1) = 2
						if(label(i-1, j) == 1) r(i-1, j) = 2
						if(label(i+1, j) == 1) r(i+1, j) = 2
					
					elseif(j == taille_reseau) then !si on est sur le bord inferieur du reseau
						if(label(i, j-1) == 1) r(i, j-1) = 2
						if(label(i-1, j) == 1) r(i-1, j) = 2
						if(label(i+1, j) == 1) r(i+1, j) = 2
					
					elseif(i == taille_reseau && j = 1) then !si on est sur le coin superieur droit
						if(label(i-1, j) == 1) r(i-1, j) = 2
						if(label(i, j+1) == 1) r(i, j+1) = 2
					
					elseif(i == taille_reseau && j = taille_reseau) then !si on est sur le coin inferieur droit
						if(label(i, j-1) == 1) r(i, j-1) = 2
						if(label(i-1, j) == 1) r(i-1, j) = 2
					
					else !si on est pas sur les bords
						if(label(i, j-1) == 1) r(i, j-1) = 2 !et si la case du haut admet une particule alors c'est une "bonne" particule
						if(label(i, j+1) == 1) r(i, j+1) = 2 !et si la case du bad admet une particule alors c'est une "bonne" particule
						if(label(i-1, j) == 1) r(i-1, j) = 2 !et si la case de gauche admet une particule alors c'est une "bonne" particule
						if(label(i+1, j) == 1) r(i+1, j) = 2 !et si la case de droite admet une particule alors c'est une "bonne" particule
					endif
					if(label(i, taille_reseau) == 2) percolation_ouinon = .true. !il y a percolation quand on a une "bonne" particule sur la derniere colonne		
				endif
			enddo ligne
		enddo colonne_verification
		max_verification = max_verification+1
	enddo colonne
	
end subroutine perco_luis
