! DFT-CES core subroutines. Written by H.-K. Lim & H. Kim
! Copyright (C) 2016 M-design group @ KAIST
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_ces
  !----------------------------------------------------------------------------
  !
  ! ... This routine reads and computes the MD rho and MD pot
  !
  USE scf,                  ONLY : ext, ext_rep
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dfftp
  USE input_parameters,     ONLY : rho_ces, pauli_rep_ces
  USE wavefunctions_module, ONLY : psic
  USE lsda_mod,             ONLY : nspin
  USE fft_interfaces,       ONLY : fwfft
  USE scatter_mod,          ONLY : scatter_grid, gather_grid
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : prefix, tmp_dir, seqopn 
  USE ions_base,            ONLY : nat, tau
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: rho_ext_read(:,:), rho_xyz(:,:,:)
  REAL(DP), ALLOCATABLE :: ext_write_r(:), v_xyz(:,:,:)
  LOGICAL :: file_exst
  REAL(DP) :: real_temp(3), glen(3), ehart_ext, charge_ext
  !
  INTEGER :: iunextpot, nat_ext, nr1_ext, nr2_ext, nr3_ext, i, j, i1, i2, i3, ounit, forc_unit
  INTEGER :: multiplier
  LOGICAL :: exst
  CHARACTER(len=80) :: readfile
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  allocate(rho_ext_read(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3))
  !
  if (ionode) then
    do j=1, 2
    !
    multiplier = -1
    if(j==1) readfile = rho_ces
    if(j==2) readfile = pauli_rep_ces
    if(j==2) multiplier = 1
    !
    allocate(rho_xyz(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x))
    !
    inquire(file= "./" // trim(readfile), exist=exst)
    if(.not.exst) call errore ('init_ces', 'rho_ces does not exist', 1)
    !
    iunextpot = find_free_unit()
    open(unit=iunextpot, file="./"//trim(readfile), status='old')
    !
    read(iunextpot,*)
    read(iunextpot,*)
    ! Skipping header in rho_ces
    read(iunextpot,'(I5,3F12.6)') nat_ext, real_temp(3)
    read(iunextpot,'(I5,3F12.6)') nr1_ext, real_temp(3)
    if(nr1_ext.ne.dfftp%nr1x) call errore('init_ces','Mismatched nr1',1)
    read(iunextpot,'(I5,3F12.6)') nr2_ext, real_temp(3)
    if(nr2_ext.ne.dfftp%nr2x) call errore('init_ces','Mismatched nr2',1)
    read(iunextpot,'(I5,3F12.6)') nr3_ext, real_temp(3)
    if(nr3_ext.ne.dfftp%nr3x) call errore('init_ces','Mismatched nr3',1)
    ! Checking number of grid points in each direction
    do i=1,nat_ext
      read(iunextpot,*)
    enddo
    ! Passing atoms section
    read(iunextpot,*)(((rho_xyz(i1,i2,i3),i3=1,nr3_ext),i2=1,nr2_ext),i1=1,nr1_ext)
    ! Parsing MD rho
    i=0
    do i3=1,nr3_ext
      do i2=1,nr2_ext
        do i1=1,nr1_ext
          i=i+1
          rho_ext_read(i,j) = multiplier*rho_xyz(i1,i2,i3) 
          ! due to MDrho.cube has opposite sign of charge
        enddo
      enddo
    enddo
    ! transforming row-major to column-major order
    close(iunextpot)
    deallocate(rho_xyz)
    !
    enddo
  endif
  ! ionode end
#if defined (__MPI)
  call scatter_grid(dfftp, rho_ext_read(:,1), ext%of_r(:,1))   ! MD charge density
  call scatter_grid(dfftp, rho_ext_read(:,2), ext_rep%of_r(:,1)) ! rep potential
#else
  ext%of_r(:,1) = rho_ext_read(:,1)
  ext_rep%of_r(:,1) = rho_ext_read(:,2)
#endif
  !
  psic(:)=ext%of_r(:,1)
  !
  call fwfft('Rho', psic, dfftp)
  !
  ext%of_g(:,1)=psic(dfftp%nl(:))
  !
  if(nspin==2) then
    ext%of_g(:,1)=ext%of_g(:,1)/2.d0
    ext%of_g(:,2)=ext%of_g(:,1)
  endif
  !
  ! MD rho r-space to g-space
  !
  ext%kin_r(:,:)=0.d0
  call v_h(ext%of_g, ehart_ext, charge_ext, ext%kin_r)
  !
  ALLOCATE( ext_write_r(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  ALLOCATE( v_xyz(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x) )
  !
  do j=1,2
  !
  if(j==1) readfile = 'v_md.cube'
  !
#ifdef __MPI
  if(j==1)  CALL gather_grid(dfftp, ext%kin_r(:,1), ext_write_r)
#else
  if(j==1) ext_write_r = ext%kin_r(:,1)
#endif
  !
  if (ionode) then
     ounit = find_free_unit()
     OPEN (unit = ounit, file = readfile, form = 'formatted', status = 'unknown')
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
!     CALL write_cubefile_diag (v_xyz, ounit)
     close(ounit)
!     write( stdout,*) "                                  "
!     WRITE( stdout, '(a,1x,a)') '  ## MD potential file has been written:', readfile
  end if
  !
  enddo
  !
  DEALLOCATE(ext_write_r)
  DEALLOCATE(v_xyz)
  !
  write(stdout,*) "                                  "
  write(stdout,*) "    DFT-CES: Initiation completed."
  !
  deallocate(rho_ext_read)
  !
  RETURN
  !
END SUBROUTINE init_ces

!----------------------------------------------------------------------------
SUBROUTINE v_ces( v )
  !----------------------------------------------------------------------------
  !
  ! ... This routine adds MD potential to Hatree potential
  !
  USE scf,                  ONLY : rho, ext, ext_rep, scf_type, vltot
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm
  USE ions_base,            ONLY : nat, tau, extfor, ityp, extfor, atm
  USE cell_base,            ONLY : alat, at, bg
  USE input_parameters,     ONLY : dft_ces
  USE lsda_mod,             ONLY : nspin
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : prefix, tmp_dir
  USE cell_base,            ONLY : at
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE ener,                 ONLY : E_ces
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !
  LOGICAL :: exst
  REAL(DP) :: glen(3), ehart_ext, charge_ext, e_ext_ion, e_ext_disp, ne_ext_disp 
  REAL(DP) :: e_ext_el, e_ext_rep_el
  REAL(DP) :: tpos(3), inpos(3), saveinpos(3,nat)
  INTEGER :: is, i, nt, at_num, forc_unit
  INTEGER, EXTERNAL :: atomic_number
  !
  do is=1, nspin
    v(:,is) = v(:,is) + ext%kin_r(:,is)+ext_rep%of_r(:,1)
  enddo
  ! Adding MD pot to Hartree pot
!  DO i = 1,nat
!    write(stdout,'("Hi check ",5F12.6," ")') tau(:,i)
!    WRITE (unit = forc_unit, FMT = *) tau(1,i)
!  ENDDO

  !
  E_ces=0.d0
  glen(1)=alat*at(1,1)/dble(dfftp%nr1x)
  glen(2)=alat*at(2,2)/dble(dfftp%nr2x)
  glen(3)=alat*at(3,3)/dble(dfftp%nr3x)
  !
  e_ext_ion=0.d0
  e_ext_ion=sum(ext%of_r(:,1)*vltot)*glen(1)*glen(2)*glen(3)
  e_ext_disp=0.d0
  ne_ext_disp=0.d0
#ifdef __MPI
  call mp_sum(e_ext_ion, intra_bgrp_comm)
#endif
  ! MD rho * ion pot energy calculation
  !
  e_ext_el=0.d0
  e_ext_rep_el=0.d0
  do is=1, nspin
    e_ext_el = e_ext_el + sum(rho%of_r(:,is)*ext%kin_r(:,is)) &
               *glen(1)*glen(2)*glen(3)
    e_ext_rep_el = e_ext_rep_el + sum(rho%of_r(:,is)*ext_rep%of_r(:,is), MASK=rho%of_r(:,is) .GT. 0.0000001) &
               *glen(1)*glen(2)*glen(3)

  enddo
#ifdef __MPI
  call mp_sum(e_ext_el, intra_bgrp_comm)
  call mp_sum(e_ext_rep_el, intra_bgrp_comm)
#endif
  ! QE rho * MD pot energy calculation
  !
  if (ionode) then
    !      
    inquire(file = 'meanforce', exist=exst)
    if(.not.exst) THEN
      !      
      OPEN (unit = forc_unit, file = 'meanforce', status = 'unknown')
      !OPEN (unit = forc_unit, file = trim(prefix)//trim( prefix )// '.init', status = 'unknown')
      !      
      DO i=1,nat
        nt = ityp(i)
        at_num = atomic_number(trim(atm(nt)))
        tpos = matmul( transpose(bg), tau(:,i) )
        tpos = tpos - nint(tpos - 0.5d0)
        inpos = alat * matmul( at, tpos )
        saveinpos(1,i)=inpos(1)
        saveinpos(2,i)=inpos(2)
        saveinpos(3,i)=inpos(3)
        WRITE (unit = forc_unit, FMT = *) saveinpos(:,i)
!     WRITE(ounit,'(I5,5F12.6)') at_num, inpos
!        write(stdout,'(" DFT-CES: CHECK ",I5," ",5F12.6," ",F15.6," ",F15.6," ",F15.6,"")') at_num, inpos, extfor(1,i), &
!        extfor(2,i), extfor(3,i)
      ENDDO
      !
    else
      OPEN (unit = forc_unit, file = 'meanforce', status = 'OLD',action='READ')
      !      
      do i=1,nat
        !
        READ ( unit = forc_unit, FMT = *) saveinpos(:,i)
!        write(stdout,'(" DFT-CES: CHECK ",I5," ",F15.6," ",F15.6," ",F15.6,"")') at_num, saveinpos(1,i), saveinpos(2,i), &
!        saveinpos(3,i)
        !
      enddo
      !
    endif
    CLOSE (unit = forc_unit)
    !
  endif

!  write(stdout, '("before the show")')


  DO i=1,nat
     nt = ityp(i)
     at_num = atomic_number(trim(atm(nt)))
     tpos = matmul( transpose(bg), tau(:,i) )
     tpos = tpos - nint(tpos - 0.5d0)
     inpos = alat * matmul( at, tpos ) 
     ne_ext_disp = ne_ext_disp + (inpos(1) - saveinpos(1,i)) * extfor(1,i) + (inpos(2) - saveinpos(2,i)) * extfor(2,i) + &
                  (inpos(3) - saveinpos(3,i)) * extfor(3,i)  
!...extfor in Ry/bohr 
!...inpos in bohr
!...summation of projection onto the vector of forces of displacement vector.
!     write(stdout,'(" DFT-CES: CHECK_ini ",I5," ",5F12.6," ",F15.6," ",F15.6," ",F15.6,"")') at_num, inpos, extfor(1,i), &
!     extfor(2,i), extfor(3,i)
!     write(stdout,'(" DFT-CES: check aft ",F15.6," ")') saveinpos(1,i)
!     write(stdout,'(" DFT-CES: e_ext_disp ",F15.6," ")') e_ext_disp 
  ENDDO
  e_ext_disp = - ne_ext_disp
!  write(stdout,'(" DFT-CES: final e_ext_disp ",F15.6," ")') e_ext_disp
  !
  E_ces = e_ext_ion + e_ext_el + e_ext_rep_el + e_ext_disp
  ! Total external interaction energy
  write(stdout,*) "        "
  write(stdout,'("  DFT-CES: QM_ion*rho_ces = ",F17.8," Ry")') e_ext_ion
  write(stdout,'("  DFT-CES: QM_rho*rho_ces = ",F17.8," Ry")') e_ext_el
  write(stdout,'("  DFT-CES: QM_rho*pauli_rep_ces = ",F17.8," Ry")') e_ext_rep_el
  write(stdout,'("  DFT-CES: Displace*meanforce = ",F17.8," Ry")') e_ext_disp
  write(stdout,'("  DFT-CES: MD pot added. E_ces= ",F17.8," Ry")') E_ces
  !
  RETURN
  !
END SUBROUTINE v_ces

!----------------------------------------------------------------------------
SUBROUTINE ces_force( force_ces )
  !----------------------------------------------------------------------------
  !
  ! ... This routine calculate external forces on ions
  !
  USE scf,                  ONLY : ext
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dfftp
  USE scatter_mod,          ONLY : gather_grid
  USE io_global,            ONLY : ionode, stdout
  USE cell_base,            ONLY : at, bg, alat
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, zv, tau
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: force_ces(3,nat)
  REAL(DP), ALLOCATABLE :: v_ext_r(:), v_grad_r(:,:), v_ext_tot(:)
  REAL(DP) :: tpos(3), inpos(3), gridValue, triInter
  INTEGER :: ipol, na
  !
  ALLOCATE( v_ext_r(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  ALLOCATE( v_ext_tot(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  ALLOCATE( v_grad_r(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3) )
  !
#ifdef __MPI
  call gather_grid(dfftp, ext%kin_r(:,1), v_ext_r(:))
#else
  v_ext_r(:) = ext%kin_r(:,1)
#endif
  ! 
  v_ext_tot = v_ext_r 
  !
  if (ionode) then
    call potgrad(v_ext_tot, v_grad_r)
    !
    do na=1, nat
      tpos = matmul( transpose(bg), tau(:,na) )
      tpos = tpos - nint(tpos - 0.5d0)
      inpos = alat * matmul( at, tpos )
      do ipol=1, 3
        force_ces(ipol, na) = -1 * zv(ityp(na)) &
                  * triInter(v_grad_r(:,ipol),inpos(1),inpos(2),inpos(3))
      enddo
      !
      write(stdout,'(''### Nuclear force on ion ''I4'' by external field''3F12.6)') na, force_ces(:,na)
      write(stdout,'(''### ZV on ion ''I4'' by external field''3F12.6)') na, zv(ityp(na))
      !
    enddo
    !
  endif
  !
  DEALLOCATE( v_ext_r )
  DEALLOCATE( v_ext_tot )
  DEALLOCATE( v_grad_r )
  !
  RETURN
  !
END SUBROUTINE ces_force

!----------------------------------------------------------------------------
SUBROUTINE potgrad( vin, vout)
  !----------------------------------------------------------------------------
  !
  ! ... This routine calculate gradients of MD potential
  !
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE fft_base,      ONLY : dfftp
  USE cell_base,     ONLY : celldm, at, alat
  !
  IMPLICIT NONE
  !
  REAL(DP)       :: vin(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
  REAL(DP)       :: vout(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3)
  INTEGER :: i, j, k, im, jm, km
  REAL(DP) :: glen(3)
  !
  glen(1) = alat*at(1,1)/dble(dfftp%nr1x)
  glen(2) = alat*at(2,2)/dble(dfftp%nr2x)
  glen(3) = alat*at(3,3)/dble(dfftp%nr3x)
  !
  do k=0, dfftp%nr3x - 1
    do j=0, dfftp%nr2x - 1
      do i=0, dfftp%nr1x - 1
        !
        if (i==0) then
          im = dfftp%nr1x - 1
        else
          im = i - 1
        end if
        if (j==0) then
          jm = dfftp%nr2x - 1
        else
          jm = j - 1
        end if
        if (k==0) then
          km = dfftp%nr3x - 1
        else
          km = k - 1
        end if
        ! x-dir potential gradient
        vout((i+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1,1) =                   &
        -1*(vin((MOD((i+1),dfftp%nr1x)+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1)&
        - vin((im+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1))/(2*glen(1))
        ! y-dir potential gradient
        vout((i+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1,2) =                   &
        -1*(vin((i+MOD((j+1),dfftp%nr2x)*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1)&
        - vin((i+jm*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1))/(2*glen(2))
        ! z-dir potential gradient
        vout((i+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1,3) =                   &
        -1*(vin((i+j*dfftp%nr1x+MOD((k+1),dfftp%nr3x)*dfftp%nr1x*dfftp%nr2x)+1)&
        - vin((i+j*dfftp%nr1x+km*dfftp%nr1x*dfftp%nr2x)+1))/(2*glen(3))
        !
        end do
    end do
  end do
  !
  return
  !
END SUBROUTINE potgrad

SUBROUTINE write_cubefile_diag(grid_data, ounit)

  USE cell_base,     ONLY : alat, at, bg
  USE ions_base,     ONLY : nat, ityp, tau, atm
  USE fft_base,      ONLY : dfftp
  USE kinds,         ONLY : DP

  IMPLICIT NONE
  INTEGER, INTENT(IN):: ounit
  REAL(DP), INTENT(IN) :: grid_data(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  INTEGER          :: i, nt, i1, i2, i3, at_num
  INTEGER, EXTERNAL:: atomic_number
  real(DP)    :: at_chrg, tpos(3), inpos(3)

  WRITE(ounit,*) 'Cubfile created from PWScf calculation'
  WRITE(ounit,*) 'Total SCF Density'
  WRITE(ounit,'(I5,3F12.6)') nat, 0.0d0, 0.0d0, 0.0d0
  WRITE(ounit,'(I5,3F12.6)') dfftp%nr1, (alat*at(i,1)/dble(dfftp%nr1),i=1,3)
  WRITE(ounit,'(I5,3F12.6)') dfftp%nr2, (alat*at(i,2)/dble(dfftp%nr2),i=1,3)
  WRITE(ounit,'(I5,3F12.6)') dfftp%nr3, (alat*at(i,3)/dble(dfftp%nr3),i=1,3)

  DO i=1,nat
     nt = ityp(i)
     at_num = atomic_number(trim(atm(nt)))
     at_chrg= dble(at_num)
     tpos = matmul( transpose(bg), tau(:,i) )
     tpos = tpos - nint(tpos - 0.5d0)
     inpos = alat * matmul( at, tpos )
     WRITE(ounit,'(I5,5F12.6)') at_num, at_chrg, inpos
  ENDDO

  DO i1=1,dfftp%nr1
     DO i2=1,dfftp%nr2
        WRITE(ounit,'(6E13.5)') (grid_data(i1,i2,i3),i3=1,dfftp%nr3)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE write_cubefile_diag

!----------------------------------------------------------------------------
FUNCTION triInter( vin, x, y, z)
  !----------------------------------------------------------------------------
  !
  ! ... This function returns interpolated grid value for arbitrary position
  !
  USE kinds,         ONLY : DP
  USE fft_base,      ONLY : dfftp
  USE cell_base,     ONLY : celldm, at, alat
  USE io_global,     ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP) :: vin(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x), triInter, gridValue
  REAL(DP) :: x, y, z, px, py, pz, xd, yd, zd, c00, c10, c01, c11, c0, c1
  INTEGER  :: i, j, k
  REAL(DP) :: glen(3)
  !
  glen(1) = alat*at(1,1)/dble(dfftp%nr1x)
  glen(2) = alat*at(2,2)/dble(dfftp%nr2x)
  glen(3) = alat*at(3,3)/dble(dfftp%nr3x)
  !
  px = DMOD(x,alat*at(1,1))/glen(1)
  py = DMOD(y,alat*at(2,2))/glen(2)
  pz = DMOD(z,alat*at(3,3))/glen(3)
  !
  xd = px - INT(px)
  yd = py - INT(py)
  zd = pz - INT(pz)
  !
  c00 = gridValue(vin,INT(px),INT(py),INT(pz))*(1-xd)     &
       +gridValue(vin,1+INT(px),INT(py),INT(pz))*xd
  c10 = gridValue(vin,INT(px),1+INT(py),INT(pz))*(1-xd)   &
       +gridValue(vin,1+INT(px),1+INT(py),INT(pz))*xd
  c01 = gridValue(vin,INT(px),INT(py),1+INT(pz))*(1-xd)   &
       +gridValue(vin,1+INT(px),INT(py),1+INT(pz))*xd
  c11 = gridValue(vin,INT(px),1+INT(py),1+INT(pz))*(1-xd) &
       +gridValue(vin,1+INT(px),1+INT(py),1+INT(pz))*xd
  !
  c0 = c00*(1-yd) + c10*yd
  c1 = c01*(1-yd) + c11*yd
  !
  triInter = c0*(1-zd) + c1*zd
  !
END FUNCTION

FUNCTION gridValue (vin, i, j, k)
  !
  USE kinds,         ONLY : DP
  USE fft_base,      ONLY : dfftp
  ! 
  IMPLICIT NONE
  !
  REAL(DP) :: vin(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x), gridValue
  INTEGER  :: i, j, k
  !
  gridValue = vin((i+j*dfftp%nr1x+k*dfftp%nr1x*dfftp%nr2x)+1)
  !
END FUNCTION
