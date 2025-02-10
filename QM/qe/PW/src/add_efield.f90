!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by J. Tobik
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
!          25/06/2009 (Riccardo Sabatini)
!               reformulation using a unique saw(x) function (included in 
!               cell_base) in all e-field related routines and inclusion of 
!               a macroscopic electronic dipole contribution in the mixing 
!               scheme. 
!

!
!--------------------------------------------------------------------------
SUBROUTINE add_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds an electric field to the local potential. The
  !   field is made artificially periodic by introducing a saw-tooth
  !   potential. The field is parallel to a reciprocal lattice vector bg, 
  !   according to the index edir.
  !
  !   if dipfield is false the electric field correction is added to the
  !   potential given as input (the bare local potential) only
  !   at the first call to this routine. In the following calls
  !   the routine exit.
  !
  !   if dipfield is true the dipole moment per unit surface is calculated
  !   and used to cancel the electric field due to periodic boundary
  !   conditions. This potential is added to the Hartree and xc potential
  !   in v_of_rho. NB: in this case the electric field contribution to the 
  !   band energy is subtracted by deband.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, saw, &
                            eopreg, forcefield, el_dipole, ion_dipole, tot_dipole
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity
! DFT-CES start
  USE scf,           ONLY : ext
  USE input_parameters, ONLY : dft_ces
  USE scatter_mod, ONLY : gather_grid
! DFT-CES end
  
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr)! ef is added to this potential
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  REAL(DP),INTENT(IN)    :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  LOGICAL,INTENT(IN)     :: iflag ! set to true to force recalculation of field
  !
  ! local variables
  !
  INTEGER :: idx,  i, j, k, j0, k0
  INTEGER :: ir, na, ipol
  REAL(DP) :: length, vamp, value, sawarg, bmod
! DFT-CES start
  REAL(DP) :: ces_dipole
! DFT-CES end

  LOGICAL :: first=.TRUE.
  SAVE first
  !
  REAL(DP), ALLOCATABLE :: v_ext_r(:,:), ext_write_r(:), v_xyz(:,:,:)
  INTEGER :: i1, i2, i3, ounit
  INTEGER, EXTERNAL :: find_free_unit
  !
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT.tefield) RETURN
  ! efield only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.
  IF ((.NOT.dipfield).AND.(.NOT.first) .AND..NOT. iflag) RETURN
  first=.FALSE.

  IF ((edir.lt.1).or.(edir.gt.3)) THEN
     CALL errore('add_efield',' wrong edir',1)
  ENDIF

  !---------------------
  !  Variable initialization
  !---------------------

  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)

  tot_dipole = 0._dp
  el_dipole  = 0._dp
  ion_dipole = 0._dp
! DFT-CES start
  ces_dipole = 0._dp
! DFT-CES end
  
  !---------------------
  !  Calculate dipole
  !---------------------
  
  if (dipfield) then
  !
  ! dipole correction is active 
  !
     CALL compute_el_dip(emaxpos, eopreg, edir, rho, el_dipole)
     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)
! DFT-CES start
     if (dft_ces) then
       CALL compute_el_dip(emaxpos, eopreg, edir, ext%of_r, ces_dipole)
     endif
     tot_dipole  = -el_dipole + ion_dipole -ces_dipole
! DFT-CES end
     CALL mp_bcast(tot_dipole, 0, intra_image_comm)
  !  
  !  E_{TOT} = -e^{2} \left( eamp - dip \right) dip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*(eamp-tot_dipole/2.d0)*tot_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2} \left( eamp - dip \right) z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *(eamp - tot_dipole) &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  else
  !
  ! dipole correction is not active
  !

     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)

  !  
  !  E_{TOT} = -e^{2} eamp * iondip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*eamp*ion_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2}  eamp z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *eamp &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  end if

  !
  !  Calculate potential and print values 
  !   
  
  length=(1._dp-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
  
  vamp=e2*(eamp-tot_dipole)*length

  IF (ionode) THEN
       !
       ! Output data
       !
       WRITE( stdout,*)
       WRITE( stdout,'(5x,"Adding external electric field":)')

       IF (dipfield) then
          WRITE( stdout,'(/5x,"Computed dipole along edir(",i1,") : ")' ) edir
! DFT-CES start
          if (dft_ces) then
            WRITE( stdout, '(5X,"System dipole  ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
     (-el_dipole+ion_dipole)*(omega/fpi), (-el_dipole+ion_dipole)*(omega/fpi)*au_debye
            WRITE( stdout, '(8X,"CES dipole  ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                           -ces_dipole*(omega/fpi), -(ces_dipole*(omega/fpi)*au_debye)
          endif
! DFT-CES end

          !
          !  If verbose prints also the different components
          !
          IF ( iverbosity > 0 ) THEN
              WRITE( stdout, '(8X,"Elec. dipole ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            el_dipole, (el_dipole*au_debye)
              WRITE( stdout, '(8X,"Ion. dipole  ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                          ion_dipole, (ion_dipole*au_debye)
          ENDIF

          WRITE( stdout, '(8X,"Dipole       ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            (tot_dipole* (omega/fpi)),   &
                                            ((tot_dipole* (omega/fpi))*au_debye)  

          WRITE( stdout, '(8x,"Dipole field ", 1F15.4," Ry au, ")') &
                                             tot_dipole
          WRITE( stdout,*)

       ENDIF

       IF (abs(eamp)>0._dp) WRITE( stdout, &
          '(8x,"E field amplitude [Ha a.u.]: ", es11.4)') eamp 
        
       WRITE( stdout,'(8x,"Potential amp.   ", f11.4," Ry")') vamp 
       WRITE( stdout,'(8x,"Total length     ", f11.4," bohr")') length
       WRITE( stdout,*)     
  ENDIF


  !
  !------------------------------
  !  Add potential
  !  
  !  V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v} 
  !          Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod} 
  !          
  !---------------------

  !
  ! Loop in the charge array
  !
  ALLOCATE( v_ext_r(dfftp%nnr,nspin)) ! ### added for external potential write
  j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
  DO ir = 1, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
     !
     ! ... three dimensional indexes
     !
     idx = ir -1
     k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
     idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
     k   = k + k0
     j   = idx / dfftp%nr1x
     idx = idx - dfftp%nr1x * j
     j   = j + j0
     i   = idx

     ! ... do not include points outside the physical range

     IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
 
     if (edir.eq.1) sawarg = DBLE(i)/DBLE(dfftp%nr1)
     if (edir.eq.2) sawarg = DBLE(j)/DBLE(dfftp%nr2)
     if (edir.eq.3) sawarg = DBLE(k)/DBLE(dfftp%nr3)
     
     value = e2*(eamp - tot_dipole)*saw(emaxpos,eopreg,sawarg) * (alat/bmod)

     vpoten(ir) = vpoten(ir) + value
     v_ext_r(ir,1) = value ! ### added for external potential write

  END DO
  
  ALLOCATE( ext_write_r(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  ALLOCATE( v_xyz(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x) )
#if defined(__MPI)
  CALL gather_grid(dfftp, v_ext_r(:,1), ext_write_r)
#else
     ext_write_r = v_ext_r(:,1)
#endif
  !
  if (ionode) then
     ounit = find_free_unit()
     OPEN (unit=ounit, file='v_saw.cube',form='formatted',status='unknown')
     !
     i=0
     do i3=1,dfftp%nr3x
        do i2=1,dfftp%nr2x
           do i1=1,dfftp%nr1x
              i=i+1
              v_xyz(i1,i2,i3) = ext_write_r(i)
           end do
        end do
     end do
     CALL write_cubefile_diag (v_xyz, ounit)
     close(ounit)
     write( stdout,*) "                                  "
     WRITE( stdout,*)'  ## Saw potential for dipc (v_saw.cube) has been written.## '
  end if
  !
  DEALLOCATE(v_ext_r)
  DEALLOCATE(ext_write_r)
  DEALLOCATE(v_xyz)
  !
  
  RETURN

END SUBROUTINE add_efield
! 
