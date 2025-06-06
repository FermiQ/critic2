! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Calculations using molecular wavefunctions.
submodule (molcalc) proc
  implicit none

  !xx! private procedures
  ! subroutine molcalc_nelec()
  ! subroutine molcalc_peach()
  ! subroutine molcalc_hfenergy()
  ! subroutine molcalc_expression(expr)

contains

  !> Driver for molecular calculations
  module subroutine molcalc_driver(line)
    use tools_io, only: ferror, faterr, uout, isexpression_or_word, lower, equal,&
       lgetword, getword
    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word, expr, savevar
    integer :: lp, lpo, imode
    integer, parameter :: imode_none = 0
    integer, parameter :: imode_nelec = 1
    integer, parameter :: imode_peach = 2
    integer, parameter :: imode_hf = 3
    integer, parameter :: imode_expr = 4

    write (uout,'("* MOLCALC: calculations using meshes ")')

    imode = imode_none
    savevar = ""
    expr = ""

    lp = 1
    do while (.true.)
       lpo = lp
       word = lgetword(line,lp)
       if (equal(word,'assign')) then
          savevar = getword(line,lp)
          if (len_trim(savevar) == 0) then
             call ferror('molcalc_driver','Zero-length variable name in MOLCALC/ASSIGN',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'peach')) then
          imode = imode_peach
       elseif (equal(word,'hf')) then
          imode = imode_hf
       else if (len_trim(word) == 0 .and. imode == imode_none) then
          imode = imode_nelec
          exit
       else if (len_trim(word) > 0) then
          if (isexpression_or_word(expr,line,lpo)) then
             lp = lpo
             imode = imode_expr
          else
             call ferror('molcalc_driver','Wrong syntax in MOLCALC',faterr,syntax=.true.)
             return
          end if
       else
          exit
       end if
    end do

    if (imode == imode_nelec) then
       call molcalc_nelec(savevar)
    elseif (imode == imode_peach) then
       call molcalc_peach(savevar)
    elseif (imode == imode_hf) then
       call molcalc_hfenergy(savevar)
    elseif (imode == imode_expr) then
       call molcalc_expression(expr,savevar)
    end if

  end subroutine molcalc_driver

  !xx! private procedures

  subroutine molcalc_nelec(savevar)
    use arithmetic, only: setvariable
    use systemmod, only: sy
    use meshmod, only: mesh
    use global, only: mesh_type, mesh_level
    use tools_io, only: string, uout
    use param, only: im_rho
    character*(*), intent(in) :: savevar

    type(mesh) :: m
    integer :: prop(1)
    real*8 :: nelec

    call m%gen(sy%c,mesh_type,mesh_level)

    write (uout,'("+ Mesh integral of the reference field")')
    call m%report()

    prop(1) = im_rho
    call m%fill(sy%f(sy%iref),prop,.not.sy%c%ismolecule)

    nelec = sum(m%f(:,1) * m%w)
    write (uout,'("+ Volume (bohr^3) = ",A)') string(sum(m%w),'f',14,8)
    write (uout,'("+ Field integral = ",A)') string(nelec,'f',14,8)
    write (uout,*)
    if (len_trim(savevar) > 0) call setvariable(trim(savevar),nelec)

  end subroutine molcalc_nelec

  subroutine molcalc_peach(savevar)
    use arithmetic, only: setvariable
    use systemmod, only: sy
    use meshmod, only: mesh
    use fieldmod, only: type_wfn
    use global, only: mesh_type, mesh_level
    use tools_io, only: ferror, faterr, getline, uin, ucopy, string, isinteger, isreal,&
       lgetword, equal, uout
    use types, only: realloc
    character*(*), intent(in) :: savevar

    type(mesh) :: m
    integer :: i, n, lp
    logical :: ok
    real*8 :: lam, dden, oia
    character(len=:), allocatable :: line, word
    integer, allocatable :: imo1(:), imo2(:), prop(:)
    real*8, allocatable :: kk(:)

    if (.not.sy%c%ismolecule) then
       call ferror("molcalc_driver","MOLCALC can not be used with crystals",faterr,syntax=.true.)
       return
    end if
    if (sy%f(sy%iref)%type /= type_wfn) then
       call ferror("molcalc_driver","PEACH can be used with molecular wavefunctions only",faterr,syntax=.true.)
       return
    end if

    write (uout,'("+ Measure of overlap between orbitals in an excitation (PEACH). ")')
    write (uout,'("  Please cite: Peach et al., J. Chem. Phys. 128 (2008) 044118. (10.1063/1.2831900)")')
    call m%report()
    allocate(imo1(10),imo2(10),kk(10))
    n = 0
    do while (getline(uin,line,ucopy=ucopy))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'endmolcalc').or.equal(word,'end')) then
          exit
       end if
       n = n + 1
       if (n > ubound(imo1,1)) then
          call realloc(imo1,2*n)
          call realloc(imo2,2*n)
          call realloc(kk,2*n)
       end if

       ! first orbital
       ok = isinteger(imo1(n),word)
       if (.not.ok) goto 999

       ! second orbital
       word = lgetword(line,lp)
       if (equal(word,"->")) &
          word = lgetword(line,lp)
       ok = isinteger(imo2(n),word)
       if (.not.ok) goto 999

       ! coefficient
       word = lgetword(line,lp)
       ok = isreal(kk(n),word)
       if (.not.ok) goto 999
    enddo
    if (n == 0) then
       call ferror("molcalc_driver","No MOs in PEACH",faterr,syntax=.true.)
       return
    end if

    ! generate and fill the mesh
    call m%gen(sy%c,mesh_type,mesh_level)
    allocate(prop(2*n))
    do i = 1, n
       prop(i) = 100 + imo1(i)
       prop(n+i) = 100 + imo2(i)
    end do
    call m%fill(sy%f(sy%iref),prop,.false.)
    deallocate(prop)

    lam = 0d0
    dden = 0d0
    do i = 1, n
       oia = sum(abs(m%f(:,i)) * abs(m%f(:,n+i)) * m%w)
       lam = lam + kk(i) * kk(i) * oia
       dden = dden + kk(i) * kk(i)
    end do
    lam = lam / dden

    write (uout,'("+ PEACH = ",A)') string(lam,'f',8,3)
    write (uout,*)
    if (len_trim(savevar) > 0) call setvariable(trim(savevar),lam)

    deallocate(imo1,imo2,kk)
    return
999 continue
    call ferror("molcalc_peach","error reading line: " // string(line),faterr)

  end subroutine molcalc_peach

  !> Compute an expression in the molecular mesh. Save the result in variable
  !> savevar.
  subroutine molcalc_expression(expr,savevar)
    use systemmod, only: sy
    use meshmod, only: mesh
    use global, only: mesh_type, mesh_level
    use arithmetic, only: setvariable
    use tools_io, only: string, uout
    use types, only: scalar_value
    character*(*), intent(in) :: expr, savevar

    type(mesh) :: m
    real*8, allocatable :: ff(:)
    real*8 :: fval, fsum
    integer :: i
    logical :: ok
    character(len=:), allocatable :: errmsg

    call m%gen(sy%c,mesh_type,mesh_level)

    write (uout,'("+ Molecular mesh integral calculation")')
    write (uout,'("  Expression: ",A)') trim(expr)
    call m%report()

    allocate(ff(m%n))
    !$omp parallel do private(fval,ok,errmsg)
    do i = 1, m%n
       fval = sy%eval(expr,errmsg,m%x(:,i))
       !$omp critical (save)
       ff(i) = fval
       !$omp end critical (save)
    end do
    !$omp end parallel do
    fsum = sum(ff * m%w)
    write (uout,'("+ Integral(",A,") = ",A/)') string(expr), string(fsum,'f',14,8)
    deallocate(ff)
    if (len_trim(savevar) > 0) call setvariable(trim(savevar),fsum)

  end subroutine molcalc_expression

  subroutine molcalc_hfenergy(savevar)
    use tools_io, only: ferror, faterr
#ifdef HAVE_CINT
    use tools_io, only: uout, string
    use arithmetic, only: setvariable
    use systemmod, only: sy
    use fieldmod, only: type_wfn
    use wfn_private, only: wfn_rhf
#endif
    character*(*), intent(in) :: savevar

#ifdef HAVE_CINT
    integer :: ioff, joff, koff, loff
    integer :: nbas, nbast
    integer :: i, j, k, l, di, dj, dk, dl, is0, shls(4), i1, j1
    real*8, allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:)
    real*8, external :: CINTgto_norm
    integer, external :: CINTcgto_cart, CINT1e_kin_cart, CINT1e_ovlp_cart, CINT1e_nuc_cart
    integer, external :: CINTcgto_spheric, CINT1e_kin_sph, CINT1e_ovlp_sph, CINT1e_nuc_sph
    integer, external :: CINT2e_cart, CINT2e_sph
    real*8, allocatable :: hmn(:,:), pmn(:,:), smn(:,:), jmn(:,:), kmn(:,:)
    real*8, allocatable :: vmn(:,:)
    real*8 :: ee, enuc, etot, dij

    if (.not.sy%goodfield(id=sy%iref,type=type_wfn)) then
       call ferror("molcalc_hfenergy","HF requires a molecular wavefunction",faterr)
    end if
    if (.not.allocated(sy%f(sy%iref)%wfn%cint)) then
       call ferror("molcalc_hfenergy","No basis set information present in this field",faterr)
    end if
    if (sy%f(sy%iref)%wfn%wfntyp /= wfn_rhf) then
       call ferror("molcalc_hfenergy","HF only implemented for restricted wavefunctions",faterr)
    end if

    associate(cint => sy%f(sy%iref)%wfn%cint)

      nbas = cint%nbas
      nbast = cint%nbast
      allocate(hmn(nbast,nbast),smn(nbast,nbast),pmn(nbast,nbast))
      allocate(jmn(nbast,nbast),kmn(nbast,nbast),vmn(nbast,nbast))
      ! allocate(eri(nbast,nbast,nbast,nbast)) ! in-core - cannot be done for most systems

      ioff = 0
      hmn = 0d0
      smn = 0d0
      do i = 1, nbas
         di = CINTcgto(i-1)
         joff = 0
         do j = 1, nbas
            dj = CINTcgto(j-1)
            if (j >= i) then
               allocate(buf1e(di,dj,1))
               shls(1) = i-1
               shls(2) = j-1

               ! kinetic energy
               if (cint%lsph) then
                  is0 = CINT1e_kin_sph(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               else
                  is0 = CINT1e_kin_cart(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               end if
               hmn(ioff+1:ioff+di,joff+1:joff+dj) = hmn(ioff+1:ioff+di,joff+1:joff+dj) + buf1e(:,:,1)

               ! nuclear attraction
               if (cint%lsph) then
                  is0 = CINT1e_nuc_sph(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               else
                  is0 = CINT1e_nuc_cart(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               end if
               hmn(ioff+1:ioff+di,joff+1:joff+dj) = hmn(ioff+1:ioff+di,joff+1:joff+dj) + buf1e(:,:,1)

               ! overlap
               if (cint%lsph) then
                  is0 = CINT1e_ovlp_sph(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               else
                  is0 = CINT1e_ovlp_cart(buf1e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env)
               end if
               smn(ioff+1:ioff+di,joff+1:joff+dj) = buf1e(:,:,1)

               ! propagate to the upper half
               hmn(joff+1:joff+dj,ioff+1:ioff+di) = transpose(hmn(ioff+1:ioff+di,joff+1:joff+dj))
               smn(joff+1:joff+dj,ioff+1:ioff+di) = transpose(smn(ioff+1:ioff+di,joff+1:joff+dj))

               deallocate(buf1e)
            end if
            joff = joff + dj
         end do
         ioff = ioff + di
      end do

      ! make the 1-dm
      pmn = matmul(transpose(cint%moc),cint%moc) * 2d0

      ! fixme: use the ERI symmetries
      jmn = 0d0
      kmn = 0d0
      ioff = 0
      do i = 1, nbas
         di = CINTcgto(i-1)
         joff = 0
         do j = 1, nbas
            dj = CINTcgto(j-1)
            koff = 0
            do k = 1, nbas
               dk = CINTcgto(k-1)
               loff = 0
               do l = 1, nbas
                  dl = CINTcgto(l-1)

                  allocate(buf2e(di,dj,dk,dl,1))
                  shls = (/i-1,j-1,k-1,l-1/)
                  if (cint%lsph) then
                     is0 = CINT2e_sph(buf2e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env,0_8)
                  else
                     is0 = CINT2e_cart(buf2e,shls,cint%atm,cint%natm,cint%bas,cint%nbas,cint%env,0_8)
                  end if

                  ! eri(ioff+1:ioff+di,joff+1:joff+dj,koff+1:koff+dk,loff+1:loff+dl) = buf2e(:,:,:,:,1)
                  do j1 = loff+1, loff+dl
                     do i1 = joff+1, joff+dj
                        kmn(i1,j1) = kmn(i1,j1) + sum(pmn(ioff+1:ioff+di,koff+1:koff+dk) * buf2e(:,i1-joff,:,j1-loff,1))
                     end do
                     do i1 = koff+1, koff+dk
                        jmn(i1,j1) = jmn(i1,j1) + sum(pmn(ioff+1:ioff+di,joff+1:joff+dj) * buf2e(:,:,i1-koff,j1-loff,1))
                     end do
                  end do

                  deallocate(buf2e)

                  loff = loff + dl
               end do
               koff = koff + dk
            end do
            joff = joff + dj
         end do
         ioff = ioff + di
      end do

      ! calculate V
      vmn = jmn - 0.5d0 * kmn

      ! calculate energies
      ee = sum(pmn * (hmn + 0.5d0 * vmn))
      enuc = 0d0
      do i = 1, sy%c%ncel
         do j = i+1, sy%c%ncel
            dij = norm2(sy%c%at(i)%r - sy%c%at(j)%r)
            enuc = enuc + sy%c%spc(sy%c%at(i)%is)%z * sy%c%spc(sy%c%at(j)%is)%z / dij
         end do
      end do
      etot = enuc + ee

      ! energies
      write (uout,'("+ Total energy = ",A," Hartree")') string(etot,'f',decimal=10)
      write (uout,'("  Number of electrons = ",A)') string(sum(pmn * smn),'f',decimal=10)
      write (uout,*)
      if (len_trim(savevar) > 0) call setvariable(trim(savevar),etot)

    end associate

  contains
    integer function CINTcgto(i)
      integer :: i
      if (sy%f(sy%iref)%wfn%cint%lsph) then
         CINTcgto = CINTcgto_spheric(i, sy%f(sy%iref)%wfn%cint%bas)
      else
         CINTcgto = CINTcgto_cart(i, sy%f(sy%iref)%wfn%cint%bas)
      end if
    end function CINTcgto
#else
    call ferror("molcalc_hfenergy","HF requires the CINT library",faterr)
#endif
  end subroutine molcalc_hfenergy

end submodule proc
