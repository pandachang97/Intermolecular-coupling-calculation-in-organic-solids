!**********************************************************************!
!   Calculate the charge transfer integrals from the output of a       !
!   Gaussian calculation.                                              !
!**********************************************************************!
program g09_tq
    use g09_commonvar
    implicit none

    character(200) :: fch, logf, fout, buff
    integer, parameter :: fno = 103
    real(8), allocatable :: tq(:), tq_sc(:)
    integer p, b, c, i, j, estate, num_of_homo, num_of_lumo
    integer, allocatable :: basis_to_atom_map(:)
    logical warn, dimer, intra_tq
    real(8) tq_mu(3), tq_sum, Cou_coupling
    real(8), allocatable :: Gross_orb_pop(:), Gross_orb_pop_sc(:)
    integer  ::  g1_zeta,  g2_zeta, N_mole, pos


    ! initialize the program, read options and input file
    dimer = .false.
    N_mole = 1


    call tq_init(fch, logf, fout, estate, dimer, intra_tq, g1_zeta, g2_zeta )
    
    if ( dimer ) N_mole = 2
    ! read in info from gaussian files
    do j = 1 , N_mole
    !! working on name of second molecule
        if ( j == 2 ) then
        !!add _m2 in fch file
        buff = fch
        buff = trim( adjustl( buff ) )
        pos = scan( buff , '.' )   
        buff = buff( 1: pos - 1 )
        buff = trim( adjustl( buff ) ) //'_m2.fch'
        fch = buff
        !! add _m2 in log file
        buff = logf
        buff = trim( adjustl( buff ) )
        pos = scan( buff , '.' )
        buff = buff( 1: pos - 1 )
        buff = trim( adjustl( buff ) ) //'_m2.log'
        logf = buff
        !! add _m2 in out file
        buff = fout
        buff = trim( adjustl( buff ) )
        pos = scan( buff , '.' )
        buff = buff( 1: pos - 1 )
        buff = trim( adjustl( buff ) ) //'_m2.out'
        fout = buff
        endif
    print*, '>> Getting info about the molecule'
    call g09_calcInfo( fch )
!   call g09_mocoeff( fch )
    call g09_overlap( logf )
!   call g09_cicoeff( logf, estate )
!   call g09_mulliken_pop( logf )       
    call g09_tdm( logf )       
 
!   write(*,*) 'The overlap matrix'
!   do i = 1, g09_task_numBasisFunctions
!     write(*,'(5f10.5)')  overlap( i, : ) 
!   end do

!write(*,'(a,(f10.6) )') 'The cic(20,24) is', cic(20,24) 
!write(*,'(a,(f10.6) )') 'The cic(25,19) is', cic(25,19) 

    if ( g09_task_numBasisFunctions .ne. g09_task_numBasisFunctionsUsed ) then
        print*, "*****************************************************************"
        print*, "WARNING: NumBasisFunctions != NumBasisFunctions Used."
        print*, "This is because you are using a pure Cartesian d/f shells basis set."
        print*, "*****************************************************************"
    end if

    ! map the ao basis functions to the atoms
        if ( allocated( basis_to_atom_map ) ) deallocate( basis_to_atom_map )
    allocate( basis_to_atom_map(g09_task_numBasisFunctions) )
    c = 0
    do p = 1, g09_task_numAtoms
        do b = 1, g09_atom(p)%basisNum
        c = c + 1
        basis_to_atom_map(c) = p
        end do
    end do


    !! Second, sum over the matrix elements in Transition Mulliken Population
    !! matrix and store them in the Gross orbital populations vector
    !! (Gross_orb_pop)

    
    if ( allocated ( Gross_orb_pop ) ) deallocate( Gross_orb_pop )
    allocate(Gross_orb_pop( g09_task_numBasisFunctions ) )

    !!Initialize the vector
        Gross_orb_pop = 0.d0

    ! loop over all basis functions
!   warn = .true.
    do b = 1, g09_task_numBasisFunctions
        p = basis_to_atom_map(b)
        do c = 1, g09_task_numBasisFunctions

        ! truncate if the overlap is zero
        if ( overlap(b,c) == 0.d0 .or. tdm(b,c) == 0.d0 ) cycle
        if ( b /= c ) then
        Gross_orb_pop(b) = Gross_orb_pop(b) + ( tdm(b,c) + tdm(c,b) ) * overlap(b,c) * 5.d-1
        else
        Gross_orb_pop(b) = Gross_orb_pop(b) + tdm(b,c) 
        endif
        end do
    end do

        Gross_orb_pop = Gross_orb_pop * (-1.d0)
    ! calculate the transition charges with transition Mulliken population 

    if ( allocated ( tq ) ) deallocate( tq )
    allocate(tq(g09_task_numAtoms))
    tq = 0.d0

    do i = 1 , g09_task_numBasisFunctions
        p = basis_to_atom_map(i) 

        tq(p) = tq(p) + Gross_orb_pop(i)

    enddo

    ! normalize the transition charges to account for the molecular
    ! orbitals being doubly occupied, the other factor of root 2 comes
    ! when the CI coefficients are normalized to 1
!   tq = tq * dsqrt(2.d0)

    ! calcualte the tansition dipole moment and sum of transition charges
    tq_mu = 0.d0
    tq_sum = 0.d0
    do p = 1, g09_task_numAtoms
        tq_sum = tq_sum + tq(p) 
        tq_mu(1) = tq_mu(1) + g09_atom(p)%x*tq(p)
        tq_mu(2) = tq_mu(2) + g09_atom(p)%y*tq(p)
        tq_mu(3) = tq_mu(3) + g09_atom(p)%z*tq(p)
    end do

    open( unit = fno, file = trim(fout), status = 'new', action = 'write')
    write(fno,*) "The Mulliken transition charges from transition density &
matrix in Gaussian 16 "
    write( fno, '(a10, i4)' ) 'NumAtoms ', g09_task_numAtoms
    write( fno, * ) 'Atomic Number, X (bohr), Y(bohr), Z(bohr), TQ (au)'
    do p = 1, g09_task_numAtoms 
        write( fno, '(i4,",",4(f14.7,","))' ) g09_atom(p)%atomicNum, &
            g09_atom(p)%x, g09_atom(p)%y, g09_atom(p)%z, tq(p)
    end do
    write( fno, * ) 'tq sum,', tq_sum
    write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (au)', tq_mu(1), &
                                         tq_mu(2), tq_mu(3)

    tq_mu = tq_mu * au_to_debye
    write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (debye)', tq_mu(1), &
                                         tq_mu(2), tq_mu(3)
    write ( fno, '(a, (f14.7) )' ) 'The strength of tq dipole (debye)' , &
          sqrt ( sum ( tq_mu ** 2 ) )

    print*, '>> Wrote output to file ', trim(adjustl(fout))
    print*, ' '

!!!Now, have a try to calculate the transition density matrix by CI coefficients
!!! and the molecular orbital coefficients.
!   if ( allocated ( tdm_sc ) ) deallocate( tdm_sc )
!   allocate(tdm_sc( g09_task_numBasisFunctions, g09_task_numBasisFunctions ))
!!!Find the occupied orbitals and unoccupied orbitals
!   num_of_homo = g09_task_numElectrons / 2 + mod( g09_task_numElectrons, 2 )
!   num_of_lumo = num_of_homo + 1 
!!! calculate the transition density matrix by cic and moc
!!! In the TDDFT case, the formula of TDM is described in this article below
!!! doi: 10.20944/preprints202106.0498.v1 
!!! Molecules 2021, 26, 4245
!   do b = 1 ,  g09_task_numBasisFunctions
!       do c = 1 , g09_task_numBasisFunctions
!           do i = 1 , num_of_homo
!               do j = num_of_lumo , g09_task_numBasisFunctionsUsed
!               tdm_sc(b,c) = cic(i,j) * moc(b,i) * moc(c,j)
!               enddo
!           enddo
!       enddo
!   enddo
!!!debug of tdm_sc
!!do i = 1 , g09_task_numBasisFunctions
!!write (*, '(5f10.5)')  tdm_sc(i,:)
!!enddo
!!!Now, we have a new transition density matrix

!   if ( allocated ( Gross_orb_pop_sc ) ) deallocate( Gross_orb_pop_sc )
!   allocate(Gross_orb_pop_sc( g09_task_numBasisFunctions ) )

!   Gross_orb_pop_sc = 0.d0

!   do b = 1, g09_task_numBasisFunctions
!       p = basis_to_atom_map(b)
!       do c = 1, g09_task_numBasisFunctions

!       ! truncate if the overlap is zero
!       if ( overlap(b,c) == 0.d0 .or. tdm_sc(b,c) == 0.d0 ) cycle
!       if ( b /= c ) then
!       Gross_orb_pop_sc(b) = Gross_orb_pop_sc(b) + tdm_sc(b,c) * overlap(b,c) * 5.d-1
!       else
!       Gross_orb_pop_sc(b) = Gross_orb_pop_sc(b) + tdm_sc(b,c) * overlap(b,c) 
!       endif
!       end do
!   end do

!       Gross_orb_pop_sc = Gross_orb_pop_sc * (-1.d0)

!   ! calculate the transition charges with transition Mulliken population 

!   if ( allocated ( tq_sc ) ) deallocate( tq_sc )
!   allocate(tq_sc( g09_task_numAtoms ))
!   tq_sc = 0.d0

!   do i = 1 , g09_task_numBasisFunctions
!       p = basis_to_atom_map(i) 

!       tq_sc(p) = tq_sc(p) + Gross_orb_pop_sc(i)

!   enddo


!   ! normalize the transition charges to account for the molecular
!   ! orbitals being doubly occupied, the other factor of root 2 comes
!   ! when the CI coefficients are normalized to 1
!   tq = tq * dsqrt(2.d0)

!   ! calcualte the tansition dipole moment and sum of transition charges
!   tq_mu = 0.d0
!   tq_sum = 0.d0
!   do p = 1, g09_task_numAtoms
!       tq_sum = tq_sum + tq(p) 
!       tq_mu(1) = tq_mu(1) + g09_atom(p)%x*tq_sc(p)
!       tq_mu(2) = tq_mu(2) + g09_atom(p)%y*tq_sc(p)
!       tq_mu(3) = tq_mu(3) + g09_atom(p)%z*tq_sc(p)
!   end do

!   write(fno,*) "The Mulliken transition charges from transition density &
! matrix by CI coefficient and molecular orbitals "
!   write( fno, '(a10, i4)' ) 'NumAtoms ', g09_task_numAtoms
!   write( fno, * ) 'Atomic Number, X (bohr), Y(bohr), Z(bohr), TQ (au)'
!   do p = 1, g09_task_numAtoms 
!       write( fno, '(i4,",",4(f14.7,","))' ) g09_atom(p)%atomicNum, &
!           g09_atom(p)%x, g09_atom(p)%y, g09_atom(p)%z, tq_sc(p)
!   end do
!   write( fno, * ) 'tq sum,', tq_sum
!   write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (au)', tq_mu(1), &
!                                        tq_mu(2), tq_mu(3)
!   tq_mu = tq_mu * au_to_debye
!   write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (debye)', tq_mu(1), &
!                                        tq_mu(2), tq_mu(3)
!   write ( fno, '(a, (f14.7) )' ) 'The strength of tq dipole (debye)' , &
!         sqrt ( sum ( tq_mu ** 2 ) )

!   print*, '>> Wrote output to file ', trim(adjustl(fout))
!   print*, ' '

    !!Now, calculating the coulomb coupling between group1 (g1_alpha to g1_zeta)
    !! and group2 (g2_alpha to g2_zeta)
!   Cou_coupling = 0.d0
!   do i = g1_alpha , g1_zeta
!       do j = g2_alpha, g2_zeta
!        Cou_coupling = Cou_coupling + tq(i)*tq(j)/ &
!             dsqrt( ( g09_atom(i)%x - g09_atom(j)%x )**2 + &
!                    ( g09_atom(i)%y - g09_atom(j)%y )**2 + &
!                    ( g09_atom(i)%z - g09_atom(j)%z )**2 )
!       enddo
!   enddo
!   !!Now, convert it from hatree to eV and output
!   Cou_coupling = Cou_coupling * hartree_to_ev
!   write(fno,'(a,i4,a,i4,a)') 'The Coulomb Coupling between group 1( atom ', g1_alpha ,' to ', g1_zeta, ' ) '
!   write(fno,'(a,i4,a,i4,a)') 'and the group 2( atom ', g2_alpha ,' to ', g2_zeta, ' )  is shown as below'
!   write(fno, '(f12.8,3x, a)')  Cou_coupling, 'eV'
    close ( fno )
    enddo

 end program


!**********************************************************************!
!   Calculate the coupling from two tq files                           ! 
!**********************************************************************!
subroutine calc_cpl_from_tq( tqf1, tqf2 , intra_tq, g1_zeta,  g2_zeta )
    use g09_commonvar
    implicit none
    
    character(100), intent(in) :: tqf1, tqf2
    character(200)   :: buff
    integer :: fno = 101
    type( g09_atom_type ), allocatable :: g09_atom1(:), g09_atom2(:)
    integer p1, p2, numAtoms1, numAtoms2, pos
    real(8), allocatable :: tq1(:), tq2(:)
    real(8) cpl, mu1(3), mu2(3), com1(3), com2(3), r(3), rmag   
    logical, intent(in)    :: intra_tq
    integer, intent(in)    :: g1_zeta, g2_zeta
    real*8      :: transe_cpl, pointd_cpl
    character(100) :: output_dat

    ! read the first tq file
    open( unit = fno, file = trim(tqf1), status = 'old', action = 'read' )
    !!ignore the first line
    read( fno, * ) buff
    read( fno, '(10X, i4)' ) numAtoms1
    if ( intra_tq ) then
        numAtoms1 = g1_zeta
    endif
    if ( allocated( g09_atom1 )) deallocate( g09_atom1 )
    if ( allocated( tq1 )) deallocate( tq1 )
    allocate( g09_atom1( numAtoms1 ), tq1( numAtoms1 ) )
    read( fno, '(X)' )
    do p1 = 1, numAtoms1
        read( fno, '(5X, 4(f14.7, X))' ) g09_atom1(p1)%x, &
            g09_atom1(p1)%y, g09_atom1(p1)%z, tq1(p1)
    end do
    close( fno )

    ! read the second tq file
    open( unit = fno, file = trim(tqf2), status = 'old', action = 'read' )
    !!ignore the first line
    read( fno, * ) buff
    read( fno, '(10X, i4)' ) numAtoms2
    if ( intra_tq ) then
        numAtoms2 = g2_zeta
    endif
    if ( allocated( g09_atom2 )) deallocate( g09_atom2 )
    if ( allocated( tq2 )) deallocate( tq2 )
    allocate( g09_atom2( numAtoms2 ), tq2( numAtoms2 ) )
    read( fno, '(X)' )
    do p2 = 1, numAtoms2
        read( fno, '(5X, 4(f14.7, X))' ) g09_atom2(p2)%x, &
            g09_atom2(p2)%y, g09_atom2(p2)%z, tq2(p2)
    end do
    close( fno )

    ! calculate the coupling (in Hartree)
    cpl = 0.d0
    do p1 = 1, numAtoms1
    do p2 = 1, numAtoms2
        cpl = cpl + tq1(p1)*tq2(p2)/ &
              dsqrt( ( g09_atom1(p1)%x - g09_atom2(p2)%x )**2 + &
                     ( g09_atom1(p1)%y - g09_atom2(p2)%y )**2 + &
                     ( g09_atom1(p1)%z - g09_atom2(p2)%z )**2 )
    end do
    end do
   
    ! convert to wavenumber and print out
    cpl = cpl * hartree_to_ev
    print'(a,f14.4)', ' Transition charge coupling (eV): ', cpl
    !! write the transition charge coupling into a dat file
!   buff = trim( adjustl( tqf1 ) )
!   pos =  scan( buff , '.')
!   buff = buff( 1 : pos - 1 )
!   buff = trim( adjustl( buff ) ) // '.dat'
!   
!   open( unit = fno, file = trim(adjustl( buff ) ), action = 'write' ) 
!   write( fno, '( a, f14.4)') 'Coulomb coupling(eV): ' , cpl
!   close(fno)
        transe_cpl = cpl
 
    ! also calculate according to the point dipole approximation
    ! for comparison

    ! first calculate the transition dipole moments and center of mass
    ! really this isnt the center of mass since the coordinates arent
    ! weighted by mass, but...
    mu1 = 0.d0
    com1 = 0.d0
    do p1 = 1, numAtoms1
        mu1(1) = mu1(1) + g09_atom1(p1)%x*tq1(p1)
        mu1(2) = mu1(2) + g09_atom1(p1)%y*tq1(p1)
        mu1(3) = mu1(3) + g09_atom1(p1)%z*tq1(p1)
        com1(1) = com1(1) + g09_atom1(p1)%x
        com1(2) = com1(2) + g09_atom1(p1)%y
        com1(3) = com1(3) + g09_atom1(p1)%z
    end do
    com1 = com1 / (1.d0 * numAtoms1 )

    mu2 = 0.d0
    com2 = 0.d0
    do p2 = 1, numAtoms2
        mu2(1) = mu2(1) + g09_atom2(p2)%x*tq2(p2)
        mu2(2) = mu2(2) + g09_atom2(p2)%y*tq2(p2)
        mu2(3) = mu2(3) + g09_atom2(p2)%z*tq2(p2)
        com2(1) = com2(1) + g09_atom2(p2)%x
        com2(2) = com2(2) + g09_atom2(p2)%y
        com2(3) = com2(3) + g09_atom2(p2)%z
    end do
    com2 = com2 / (1.d0 * numAtoms2 )
     
    ! the displacement unit vector and the magnitude
    r = com2 - com1
    rmag = dsqrt(sum(r**2))
    r = r/rmag
    
    ! the transition dipole coupling
    cpl = ( mu1(1)*mu2(1) + mu1(2)*mu2(2) + mu1(3)*mu2(3) - &
          3.d0 * ( mu1(1) * r(1) + mu1(2) * r(2) + mu1(3) * r(3) ) * &
                 ( mu2(1) * r(1) + mu2(2) * r(2) + mu2(3) * r(3) ) ) / &
                 rmag**3
    cpl = cpl * hartree_to_ev
    print'(a,f14.4)', ' Transition dipole coupling (eV): ', cpl
        pointd_cpl = cpl
!! working on the name of output file
   output_dat = trim( adjustl( 'cpl_result.dat' ) )

 !! working on printing out result
  open ( unit = fno, file = trim( output_dat ), status = 'new', action = 'write' )
        write ( fno, '( a, 2f10.5 )') "The Coulomb coupling by Trans_charg and Dipole appro. is: " , &
  transe_cpl, pointd_cpl
  
        close( fno )

end subroutine


!**********************************************************************!
!   Initialize the program. Read command line options and input file   ! 
!**********************************************************************!
subroutine tq_init(fch, logf, fout, estate, dimer, intra_tq, g1_zeta,  g2_zeta )
    implicit none

    character(200), intent(out) :: fch, logf, fout
    integer, intent(out) :: estate, g1_zeta, g2_zeta
    integer nargs, narg, ios, line, pos, N_mole , i
    integer, parameter :: fno = 67, fno2 = 68
    character(32) arg, fin, label, task
    character(100) buff, emethod, tqf1, tqf2, fxyz, fxyz2
    logical exists, makeinput, dimer, exists2, intra_tq, caluclating_Coulomb
    character(1) g09_estate 

    makeinput = .false.
    dimer = .false.
    exists2 = .false.
    intra_tq = .false.
    caluclating_Coulomb = .false.
    fin = ''
    fout = ''
    fxyz = ''
    fxyz2 = ''
    N_mole = 1
    buff = ''

    ! check if any command line arguments are found
    nargs = command_argument_count()

    if ( nargs > 0 ) then
        narg = 1
        do while ( narg <= nargs )
            call get_command_argument( narg, arg )
            arg = adjustl(arg)
            select case ( arg )
                case('--help', '-h')
                    call print_help()
                    stop
                case('-i')
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No input file given with option -i'
                        call print_help()
                        stop
                    end if
                    ! get the name of the input file
                    call get_command_argument(narg, fin)
                    fin = adjustl(fin)
                    ! check if it exists
                    inquire( file = trim( fin ), exist = exists )
                    if ( .not. exists ) then
                        print'(3a)', ' Input file ', trim(adjustl(fin)), &
                                     ' does not exist.'
                        call print_help()
                        stop
                    end if
                case('-o')
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No output file given with option -o'
                        call print_help()
                        stop
                    end if
                    ! get the name of the input file
                    call get_command_argument(narg, fout)
                    fout = adjustl(fout)
                case('-m')
                    makeinput = .true. 
                case('-c')
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No tq file given with option -c'
                        call print_help()
                        stop
                    end if
                    ! get the name of the tq file
                    call get_command_argument(narg, tqf1)
                    tqf1 = adjustl(tqf1)
                    ! get the name of the second tq file
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No second tq file given with option -c'
                        call print_help()
                        stop
                    end if
                    ! get the name of the tq file
                    call get_command_argument(narg, tqf2)
                    tqf2 = adjustl(tqf2)
                    caluclating_Coulomb = .true.
                case default
                    print'(3a)', ' Option ', adjustl(arg), ' unknown'
                    call print_help()
                    stop
            end select
            narg = narg + 1
        end do
    else
        call print_help()
        stop
    end if

    if ( fin == '' ) then
        print*, ' No input file was given. Cannot proceed'
        call print_help()
        stop
    else
    ! read the input file
    print*,  '  '
        print*, 'Reading input file: ', trim(adjustl(fin))
        open( unit = fno, file = fin, action = 'read' )
        line = 0
        ios = 0
        do while ( ios == 0 )
            ! read the file line by line
            read( fno, '(a)', iostat = ios ) buff
            buff = trim( adjustl( buff ) )
            ! check in out status is ok
            if ( ios == 0 ) then
                ! increase line number
                line = line + 1
                ! find the position of the separator
                pos = scan( buff, ' ' )
                ! store the label
                label = buff( 1:pos )
                ! store the parameter
                buff = buff( pos + 1 : )
                ! ignore comment lines
                if ( label(1:1) ==  '#' ) cycle

                ! set the parameters
                select case (label)
                    case('emethod')
                        read( buff, *, iostat = ios ) emethod
                        print*, 'emethod : ', trim(adjustl(emethod))
                    case('task')
                        read( buff, *, iostat = ios ) task
                        print*, 'task : ', trim(adjustl(task))
                    case('fch')
                        read( buff, *, iostat = ios ) fch
                        print*, 'fch: ', trim(adjustl(fch))
                    case('log')
                        read( buff, *, iostat = ios ) logf
                        print*, 'log: ', trim(adjustl(logf))
                    case('xyz_file_m1')
                        read( buff, *, iostat = ios ) fxyz
                        print*, 'fxyz: ', trim(adjustl(fxyz))
                    case('xyz_file_m2')
                        read( buff, *, iostat = ios ) fxyz2
                        print*, 'fxyz2: ', trim(adjustl(fxyz2))
                    case('dimer_pack')
                        read( buff, *, iostat = ios ) dimer
                        print*, 'dimer: ', dimer
                    case('estate')
                        read( buff, *, iostat = ios ) estate
                        print*, 'estate: ', estate
                        write(g09_estate,'(I1)') estate 
!                       write(*,*) "The estate is ", g09_estate
                    case('intra_Cou')
                        read( buff, *, iostat = ios ) intra_tq
                        print*, 'intra Coulomb coupling: ', intra_tq
                    case('group1_last') 
                        read( buff, *, iostat = ios ) g1_zeta
                        print*, 'group1_last: ', g1_zeta
                    case('group2_last') 
                        read( buff, *, iostat = ios ) g2_zeta
                        print*, 'group2_last: ', g2_zeta
                    case default
                        print*, 'Label ', trim(adjustl(label)), &
                                ' unknown'
                end select
            end if
        end do
        close( fno )
        print*, ' '
    end if

    if ( makeinput ) then
        ! set default method and basis if they are not set in the input file
        if ( emethod == '' ) emethod = 'CIS(Nstates=1,Root=1)/cc-pVTZ'
        if ( dimer ) N_mole = 2 
!       print*, 'I am working on two molecules', N_mole
       do i = 1 , N_mole 
            if ( i == 2 ) then
            !! change xyz file name for second molecule
            buff = trim( adjustl( fxyz ) )
            pos = scan( buff , '1' )
            if( pos == 0 ) stop 'wrong xyz name used, please use ***1.xyz instead'
            buff = buff( 1 : pos - 1 )
            buff =trim( adjustl( buff ) )//'2.xyz'
            fxyz = buff 
!           print*, fxyz , N_mole
            !! change com file name for second molecule
            buff = trim( adjustl( task ) )
            buff = trim( adjustl( buff ) )//'_m2'
            task = buff 
!           print*, task, N_mole
            endif
        ! make sure the xyz files exists
        inquire( file = trim(fxyz), exist = exists )
        if ( .not. exists ) then
            print*, 'The file ', trim(adjustl(fxyz)), &
                    ' does not exist. Aborting...'
            stop
        end if
 
        ! setup the gaussian calculation
        open( unit = fno, file = trim(adjustl(task))//'.com', &
              action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))
        write( fno, * ) '%NProc=24'
        write( fno, * ) '%mem=32GB'
        write( fno, * ) '#P '//trim(emethod)//' SP IOP(9/40=5)  IOP(3/33=1)' &
                        //' Pop=Full NoSymm IOP(6/22=-4) IOP(6/29='//trim(adjustl(g09_estate))//') '
        write( fno, * )
        write( fno, * ) trim(adjustl(task))
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyz), action = 'read' )
!       read( fno2, * )
!       read( fno2, * )
        do
            read( fno2, '(a)', end=105) buff
            write( fno, * ) buff
        end do
105     continue
        close( fno2 )
        !!Here, if we want to calculate the dimer


        ! write a blank line so g09 doesnt crash
        write( fno, * )

        close( fno )
    enddo
    
        print*, 'Successfuly made the G16 input files'
        stop
    end if
 !! Calculating the Coulomb coupling with transtion charges in the output file
    if ( caluclating_Coulomb ) then
    call calc_cpl_from_tq( tqf1, tqf2, intra_tq, g1_zeta, g2_zeta )
    stop
    endif

end subroutine


!**********************************************************************!
!   Read in the Mulliken population from a Gaussian output file        !
!**********************************************************************!
subroutine g09_mulliken_pop( logf )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: logf
    character(100)   :: format_string, buff_basis , buff_atom, buff_atom_number
    integer, parameter :: fno = 509 
    logical :: exists, special
    character(100) read_line
    integer  :: col(5), row, i, j, m, n, pos, nimei
    

    ! check that the log file exists
    inquire( file = logf, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Log file file ', trim(adjustl(logf)), ' not found.'
        stop
    end if

    ! open the log file and read the Full Mulliken population analysis
    open( unit = fno, file = logf, status = 'old', &
          action = 'read' )


    ! allocate space for the Full Mulliken population analysis
    if ( allocated ( mulliken_pop_matrix ) ) deallocate( mulliken_pop_matrix )
    allocate( mulliken_pop_matrix ( g09_task_numBasisFunctions,          &
                       g09_task_numBasisFunctions ) )

    ! read the Full Mulliken population analysis
    do
        read( fno, '(a)', end = 1001 ) read_line
        if ( index( read_line , 'Full Mulliken population analysis' ) .ne. 0 ) then
            do 
                ! first read the column numbers. The formatting in the
                ! log file requires a little extra work to read these
                ! values
                col = 0
                read( fno, *, iostat=nimei, end=102 ) col(:)
                if (nimei<0)   stop   
!               write(*,*) col(:)
                102 continue
                ! read the info in each row
                do
                    read ( fno, '(i4, 17x, *(f10.5))', end=105, iostat=nimei) row,  &
                     ( mulliken_pop_matrix( row, col(i) ), i = 1, &
                             maxval( (col(:) ) - col(1) ) + 1 )
                    105 continue
        
                    if ( row == g09_task_numBasisFunctions ) exit
                end do
                if ( maxval( col ) == g09_task_numBasisFunctions ) exit
            end do
            exit
        end if
    end do


    ! close the file
    close( fno )

    ! for convience fill in the upper diagonal too
    do i = 1, g09_task_numBasisFunctions
    do j = i+1, g09_task_numBasisFunctions
        mulliken_pop_matrix( i, j ) = mulliken_pop_matrix( j, i )
    end do
    end do

    
!   do i = 1, g09_task_numBasisFunctions
!     write(*,'(5f10.5)')  mulliken_pop_matrix( i, : ) 
!   end do
    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', ' Full Mulliken population analysis not found in file ', trim(adjustl(logf))
    print'(2a)', ' Try running Gausion 09 with IOP(6/22=-4) and IOP(6/29=1) '
    stop
  
end subroutine

!**********************************************************************!
!   Read the Transition Density Matrix(TDM) from a Gaussian output file!
!**********************************************************************!
subroutine g09_tdm( logf )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: logf
    character(100)   :: format_string
    integer, parameter :: fno = 519 
    logical :: exists, special
    character(100) read_line
    integer  :: col(5), row, i, j, m, n, pos, nimei
    

    ! check that the log file exists
    inquire( file = logf, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Log file file ', trim(adjustl(logf)), ' not found.'
        stop
    end if

    ! open the log file and read the Full Mulliken population analysis
    open( unit = fno, file = logf, status = 'old', &
          action = 'read' )


    ! allocate space for the Full Mulliken population analysis
    if ( allocated ( tdm ) ) deallocate( tdm )
    allocate( tdm( g09_task_numBasisFunctions,          &
                       g09_task_numBasisFunctions ) )

    ! read the Full Mulliken population analysis
    do
        read( fno, '(a)', end = 1001 ) read_line
        if ( index( read_line , 'Density Matrix' ) .ne. 0 ) then
            do 
                ! first read the column numbers. The formatting in the
                ! log file requires a little extra work to read these
                ! values
                col = 0
                read( fno, '(18x,(5i10) )', iostat=nimei, end=112 ) col(:)
                if (nimei /= 0)   stop 'Error reading column of TDM. Aborting'
!               write(*,*) col(:)
                112 continue
                ! read the info in each row
                do
                    read ( fno, '(i4, 17x, 5(f10.5))', end=115, iostat=nimei) row,  &
                     ( tdm( row, col(i) ), i = 1, &
                             maxval( (col(:) ) - col(1) ) + 1 )
                if (nimei /= 0)   stop 'Error reading elements of TDM. Aborting'
                    
                    115 continue
        
                    if ( row == g09_task_numBasisFunctions ) exit
                end do
                if ( maxval( col ) == g09_task_numBasisFunctions ) exit
            end do
            exit
        end if
    end do


    ! close the file
    close( fno )

    ! for convience fill in the upper diagonal too
    do i = 1, g09_task_numBasisFunctions
    do j = i+1, g09_task_numBasisFunctions
        tdm( i, j ) = tdm( j, i )
    end do
    end do

!   write (*,*) 'The transition density matrix is'
!   do i = 1, g09_task_numBasisFunctions
!     write(*,'(5f10.5)')  tdm( i, : ) 
!   end do
    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', ' Transition Density Matrix not found in file ', trim(adjustl(logf))
    print'(2a)', ' Try running Gausion 09 with IOP(6/22=-4) and IOP(6/29=1) '
    stop
  
end subroutine


!**********************************************************************!
!   Print help info for user                                           ! 
!**********************************************************************!
subroutine print_help()
    implicit none

    print*, ' '
    print*, ' g09_tq version 1.0'
    print*, ' '
    print*, ' This  program  calculates  the transition  charges  for  a '
    print*, ' specified excited state  of a  molecule.  The  calculation '
    print*, ' uses output from a Gaussian 09 calculation.'
    print*,  '  '
    print*, ' Usage:'
    print*,  '  '
    print*, ' -- help, -h: print this message'
    print*, ' -i file    : specifies the input file containing the names '
    print*, '              of the formatted checkpoint and log files. If '
    print*, '              the -m option is also given, the  input  file '
    print*, '              is read for information about the Gaussian 09 '
    print*, '              calculation to add to the input files.'
    print*, ' -o file    : specifies the output file'
    print*, ' -m         : makes one Gaussian 09 input files from  a xyz '
    print*, '              files of atomic coordinates.'
    print*, ' '
    print*, ' To  use  this  program to calculate the transition charges '
    print*, ' you  must  first  run  an  excited  state  calculation  in '
    print*, ' Gaussian 09 of the molecule of interest. The -m  option of '
    print*, ' this program can  be  used  to  generate  the  appropriate '
    print*, ' Gaussian 09 input files if xyz  coordinate files  for  the '
    print*, ' molecule of interest is provided.'
    print*, ' '
    print*, ' The  resulting  checkpoint  file   must  be  converted  to '
    print*, ' formatted checkpoint files using the Gaussian 09 formcheck '
    print*, ' utility before  it  can be read by this program.'
    print*,  '  '

end subroutine
