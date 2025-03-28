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

!> Integration by recursive division of the IWS.
submodule (qtree) proc
  implicit none

  !xx! private procedures
  ! subroutine qtree_setsph1(lvl,verbose)
  ! subroutine qtree_setsph2(verbose)

  integer, parameter :: INT_spherequad_type = INT_lebedev !< type of angular quadrature
  integer, parameter :: INT_spherequad_ntheta = 30 !< number of nodes in theta (polar)
  integer, parameter :: INT_spherequad_nphi = 30 !< number of nodes in phi (azimuthal)
  integer, parameter :: INT_spherequad_nleb = 170 !< number of lebedev nodes

contains

  !> Main driver for the QTREE integration.
  module subroutine qtree_driver(line)
    use systemmod, only: sy
    use qtree_basic, only: minl, sphfactor, sphintfactor, gradient_mode,&
       qtree_ode_mode, stepsize, docontacts, killext,&
       ode_abserr, integ_mode, integ_scheme, keastnum, plot_mode, prop_mode, mpstep,&
       qtreefac, cub_abs, cub_rel, cub_mpts, ws_origin, ws_scale, autosph, checkbeta,&
       plotsticks, color_allocate, setsph_lvl, vcutoff, r_betaint, r_betagp
    use global, only: eval_next, dunit0, iunit
    use tools_io, only: isinteger, lgetword, ferror, faterr, equal
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word
    integer :: lp, level, plevel, idum, idum2
    real*8 :: rdum
    logical :: ok

    ! default values
    level = 6
    plevel = 0
    minl = 4
    if (allocated(sphfactor)) deallocate(sphfactor)
    allocate(sphfactor(sy%c%nneq))
    sphfactor = 0d0
    if (allocated(sphintfactor)) deallocate(sphintfactor)
    allocate(sphintfactor(sy%c%nneq))
    sphintfactor = 1d0
    gradient_mode = 1
    qtree_ode_mode = 8
    stepsize = 0.1d0
    ode_abserr = -1d0
    integ_mode = 0
    integ_scheme = 1
    keastnum = 5
    plot_mode = 0
    prop_mode = 1
    mpstep = 0
    qtreefac = 2d0
    cub_abs = 1d0
    cub_rel = 1d-6
    cub_mpts = 1000
    docontacts = .false.
    ws_origin = (/0d0, 0d0, 0d0/)
    ws_scale = -1d0
    killext = .true.
    autosph = 2
    checkbeta = .false.
    plotsticks = .true.
    color_allocate = -1
    setsph_lvl = 6
    vcutoff = 0d0
    if (allocated(r_betaint)) deallocate(r_betaint)
    allocate(r_betaint(sy%c%nneq))
    r_betaint = 0d0
    if (allocated(r_betagp)) deallocate(r_betagp)
    allocate(r_betagp(sy%c%nneq))
    r_betagp = 0d0

    ! read input
    lp = 1
    ok = isinteger(level,line,lp)
    ok = ok .and. isinteger(plevel,line,lp)
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"minl")) then
          ok = eval_next(minl,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in MINL keyword',faterr,line)
       elseif (equal(word,"gradient_mode")) then
          ok = eval_next(gradient_mode,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in GRADIENT_MODE keyword',faterr,line)
       elseif (equal(word,"qtree_ode_mode")) then
          ok = eval_next(qtree_ode_mode,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in QTREE_ODE_MODE keyword',faterr,line)
       elseif (equal(word,"stepsize")) then
          ok = eval_next(stepsize,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in STEPSIZE keyword',faterr,line)
          stepsize = stepsize / dunit0(iunit)
       elseif (equal(word,"ode_abserr")) then
          ok = eval_next(ode_abserr,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in ODE_ABSERR keyword',faterr,line)
          stepsize = stepsize / dunit0(iunit)
       elseif (equal(word,"sphfactor")) then
          do while (isinteger(idum,line,lp))
             ok = eval_next(rdum,line,lp)
             if (.not.ok) &
                call ferror('qtree_driver','error in SPHFACTOR keyword',faterr,line)
             if (idum < 0 .or. idum > sy%c%nneq) &
                call ferror('qtree_driver','error in SPHFACTOR keyword: atom id out of range',faterr,line)
             sphfactor(idum) = rdum
          end do
       elseif (equal(word,"sphintfactor")) then
          do while (isinteger(idum,line,lp))
             ok = eval_next(rdum,line,lp)
             if (.not.ok) &
                call ferror('qtree_driver','error in SPHINTFACTOR keyword',faterr,line)
             if (idum < 0 .or. idum > sy%c%nneq) &
                call ferror('qtree_driver','error in SPHINTFACTOR keyword: atom id out of range',faterr,line)
             sphintfactor(idum) = rdum
          end do
       elseif (equal(word,"integ_mode")) then
          ok = eval_next(idum,line,lp)
          ok = ok .and. eval_next(idum2,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in INTEG_MODE keyword',faterr,line)
          if (idum < 0 .or. idum > 20) &
             call ferror('qtree_driver','error in INTEG_MODE level: out of range',faterr,line)
          integ_mode(idum) = idum2
       elseif (equal(word,"integ_scheme")) then
          ok = eval_next(integ_scheme,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in INTEG_SCHEME keyword',faterr,line)
       elseif (equal(word,"keastnum")) then
          ok = eval_next(keastnum,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in KEASTNUM keyword',faterr,line)
       elseif (equal(word,"plot_mode")) then
          ok = eval_next(plot_mode,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in PLOT_MODE keyword',faterr,line)
       elseif (equal(word,"prop_mode")) then
          ok = eval_next(prop_mode,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in PROP_MODE keyword',faterr,line)
       elseif (equal(word,"mpstep")) then
          ok = eval_next(mpstep,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in MPSTEP keyword',faterr,line)
       elseif (equal(word,"qtreefac")) then
          ok = eval_next(qtreefac,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in MPSTEP keyword',faterr,line)
       elseif (equal(word,"cub_abs")) then
          ok = eval_next(cub_abs,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in CUB_ABS keyword',faterr,line)
       elseif (equal(word,"cub_rel")) then
          ok = eval_next(cub_rel,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in CUB_REL keyword',faterr,line)
       elseif (equal(word,"cub_mpts")) then
          ok = eval_next(cub_mpts,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in CUB_MPTS keyword',faterr,line)
       elseif (equal(word,"docontacts")) then
          docontacts = .true.
       elseif (equal(word,"ws_origin")) then
          ok = eval_next(ws_origin(1),line,lp)
          ok = ok .and. eval_next(ws_origin(2),line,lp)
          ok = ok .and. eval_next(ws_origin(3),line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in WS_ORIGIN keyword',faterr,line)
       elseif (equal(word,"ws_scale")) then
          ok = eval_next(ws_scale,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in WS_SCALE keyword',faterr,line)
       elseif (equal(word,"nokillext")) then
          killext = .false.
       elseif (equal(word,"autosph")) then
          ok = eval_next(autosph,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in AUTOSPH keyword',faterr,line)
       elseif (equal(word,"checkbeta")) then
          checkbeta = .true.
       elseif (equal(word,'noplotsticks')) then
          plotsticks = .false.
       elseif (equal(word,"color_allocate")) then
          ok = eval_next(color_allocate,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in COLOR_ALLOCATE keyword',faterr,line)
       elseif (equal(word,'setsph_lvl')) then
          ok = eval_next(setsph_lvl,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in SETSPH_LVL keyword',faterr,line)
       elseif (equal(word,'vcutoff')) then
          ok = eval_next(vcutoff,line,lp)
          if (.not.ok) &
             call ferror('qtree_driver','error in VCUTOFF keyword',faterr,line)
       elseif (len_trim(word) == 0) then
          exit
       else
          call ferror('qtree_driver','unknown keyword in QTREE',faterr,line,syntax=.true.)
       end if
    end do

    ! run qtree
    call qtree_integration(level,plevel)

  end subroutine qtree_driver

  !> Main routine for the QTREE integration.
  module subroutine qtree_integration(lvl, plvl)
    use systemmod, only: sy
    use integration, only: int_output_header, int_output_fields
    use qtree_tetrawork, only: paint_inside_spheres, tetrah_subdivide, integ_corner_sum
    use qtree_utils, only: open_difftess, close_difftess, small_writetess, getkeast
    use qtree_basic, only: qtreeidx, qtreei, qtreer, nterm, ngrd_term, ngrd_int, nlocate,&
       nlocate_sloppy, korder, maxl, nnuc, intcorner_deferred, kxyz, kw, nt_orig, minlen,&
       maxlen, periodic, leqv, tcontact, tcontact_void, r_betaint, lustick, atprop, &
       tvol, savefgr, savelapgr, ndiff, ngrd1, ngrd2, bvec, perm3, cmat, torig, lrotm,&
       eps_tetrah_contact, cindex, leqvf, qtree_initialize, qtree_checksymmetry,&
       qtree_cleanup, minl, sphfactor, gradient_mode, qtree_ode_mode, stepsize, ode_abserr,&
       integ_mode, integ_scheme, plot_mode, prop_mode, qtreefac, docontacts, ws_scale,&
       autosph, plotsticks, color_allocate, setsph_lvl, vcutoff
    use CUI, only: cubpack_info
    use keast, only: keast_rule, keast_order_num
    use global, only: quiet, &
       int_radquad_errprop, int_gauleg, int_qags, int_radquad_type,&
       fileroot
    use tools_io, only: uout, faterr, ferror, warning, string, fopen_write, tictac, fclose
    use bisect, only: sphereintegrals_lebedev, sphereintegrals_gauleg
    use types, only: basindat, int_result, out_field

    integer, intent(in) :: lvl
    integer, intent(in) :: plvl

    integer :: i, j, k, tt, tto, ff, ffo, pp
    integer(qtreeidx) :: idx, idxo
    integer :: iv0(3,4), idum, al1, al2, al3, al4, opid, cid
    logical :: invop
    real*8 :: xp(3)
    integer :: l2
    logical :: ok
    real*8 :: xx(3)
    character*50 :: roottess
    real*8, allocatable :: sphereprop(:,:)
    real*8 :: abserr
    integer :: neval, meaneval, vin(3), vino(3)
    character(10) :: pname
    real*8 :: xface_orig(3,3), xface_end(3,3), vtotal, memsum
    real*8 :: rface_orig(3,3), rface_end(3,3), face_diff(3)
    integer :: nalloc, ralloc
    integer(qtreeidx) :: siz
    ! data arrays
    integer(qtreei), allocatable :: trm(:,:)
    real(qtreer), allocatable :: fgr(:,:), lapgr(:,:), vgr(:)
    real*8, allocatable :: acum_atprop(:,:)
    ! for the output
    type(basindat) :: bas
    type(int_result), allocatable :: res(:)

    real*8, parameter :: eps = 1d-6

    if (.not.quiet) call tictac("Start QTREE")

    ! header
    write (uout,'("* QTREE integration ")')
    write (uout,'("  Please cite: ")')
    write (uout,'("    A. Otero-de-la-Roza et al., Comput. Phys. Commun. 180 (2009) 157. (10.1002/jcc.21620)"/)')
    if (sy%c%ismolecule) then
       call ferror("qtree","QTREE only available for crystals",faterr)
    end if

    ! consistency warnings
    if (plot_mode > 0 .and. color_allocate == 0) then
       call ferror('qtree','Non zero PLOT_MODE is not compatible with a single color array',warning)
       write (uout,'("* PLOT mode is set to 0.")')
       write (uout,*)
       plot_mode = 0
    end if
    !$ if (docontacts) then
    !$    call ferror('qtree','DOCONTACTS is not compatible with multi-threaded calculations.',warning)
    !$    write (uout,'("* The contacts are NOT going to be used.")')
    !$    write (uout,*)
    !$    docontacts = .false.
    !$ end if
    if (docontacts .and. color_allocate == 0) then
       call ferror('qtree','DOCONTACTS is not compatible with a single color array',warning)
       write (uout,'("* The contacts are NOT going to be used.")')
       write (uout,*)
       docontacts = .false.
    end if

    ! allocate
    allocate(sphereprop(sy%npropi,sy%c%nneq))

    ! Determine sphere radii
    write (uout,'("+ Pre-calculating the beta-sphere radii at lvl : ",A)') string(setsph_lvl)

    if (autosph == 1) then
       call qtree_setsph1(min(setsph_lvl,lvl),.false.)
    else
       call qtree_setsph2(.false.)
    end if

    write (uout,'("+ Initializing QTREE")')
    call qtree_initialize(lvl,plvl,acum_atprop,trm,fgr,lapgr,vgr,.false.)
    l2 = 2**maxl
    nterm = 0
    ngrd_term = 0
    ngrd_int = 0
    nlocate = 0
    nlocate_sloppy = 0
    korder = 0

    ! mask for the property output
    allocate(res(sy%npropi))
    do i = 1, sy%npropi
       res(i)%reason = ""
       res(i)%done = .false.
       allocate(res(i)%psum(nnuc))
       res(i)%psum = 0d0
       res(i)%outmode = out_field
    end do
    if (prop_mode == 0) then
       res(1)%done = .true.
       do i = 2, sy%npropi
          res(i)%reason = "because prop_mode = 0"
       end do
    else if (prop_mode == 1) then
       res(1)%done = .true.
       res(2)%done = .true.
       do i = 3, sy%npropi
          res(i)%reason = "because prop_mode = 1"
       end do
    else if (prop_mode == 2) then
       res(1)%done = .true.
       res(2)%done = .true.
       res(3)%done = .true.
       do i = 4, sy%npropi
          res(i)%reason = "because prop_mode = 2"
       end do
    else
       res(:)%done = .true.
    end if

    ! header
    write (uout,'("  Maximum subdivision level: ",A)') string(maxl)
    write (uout,'("  Tetrahedra pre-split level: ",A)') string(plvl)
    write (uout,'("  Number of attractors: ",A)') string(nnuc)
    write (uout,'("  GRADIENT mode: ",A)') string(gradient_mode)
    write (uout,'("  ODE mode: ",A)') string(qtree_ode_mode)
    write (uout,'("  STEP size in ODE integration: ",A)') string(stepsize,'e',decimal=6)
    write (uout,'("  Abs. error in ODE integration: ",A)') string(ode_abserr,'e',decimal=6)
    write (uout,'("  PROPERTY calculation level: ",A)') string(prop_mode)
    write (uout,'("  INTEGRATION scheme: ",A)') string(integ_scheme)
    write (uout,'("  INTEGRATION modes, per level: ")')
    do i = minl+1, maxl
       write (uout,'("  Level ",A,": ",A)') string(i), string(integ_mode(i))
    end do
    write (uout,'("  CORNER INTEGRATION deferred flag: ",L1)') intcorner_deferred
    write (uout,'("  PLOT mode : ",A)') string(plot_mode)

    if (any(integ_mode(minl+1:maxl) >= 1 .and. integ_mode(minl+1:maxl) <= 10)) then
       write (uout,'("+ Using the KEAST library ")')
       write (uout,'("  P. Keast, Comput. Meth. Appl. Mech. 55 (1986) 339-348.")')
       call getkeast()
       do i = 1, 10
          call keast_order_num(i,korder(i))
          call keast_rule(i,korder(i),kxyz(:,:,i),kw(:,i))
       end do
    else if (any(integ_mode(minl+1:maxl) == 12)) then
       write (uout,'("+ Using the CUBPACK library ")')
       write (uout,'("  R. Cools and A. Haegemans, ACM Trans. Math. Softw. 29 (2003) 287-296.")')
       call cubpack_info(uout)
    end if
    write (uout,*)

    ! Integration spehres
    if (INT_radquad_errprop > 0 .and. INT_radquad_errprop <= sy%npropi) then
       pname = sy%propi(INT_radquad_errprop)%prop_name
    else
       pname = "max       "
    end if
    write (uout,'("+ BETA sphere integration details")')
    if (INT_spherequad_type == INT_gauleg) then
       write (uout,'("  Method : Gauss-Legendre, non-adaptive quadrature ")')
       write (uout,'("  Polar angle (theta) num. of nodes: ",A)') string(INT_spherequad_ntheta)
       write (uout,'("  Azimuthal angle (phi) num. of nodes: ",A)') string(INT_spherequad_nphi)
    else if (INT_spherequad_type == INT_lebedev) then
       write (uout,'("  Method : Lebedev, non-adaptive quadrature ")')
       write (uout,'("  Number of nodes: ",A)') string(INT_spherequad_nleb)
    end if
    if (any(sphfactor(1:nnuc) < -1d-12)) then
       write (uout,'("  Using Rodriguez et al. strategy for beta-spehre radius (JCC 30(2009)1082).")')
    end if
    if (INT_radquad_type == INT_qags .or. &
        INT_radquad_type == INT_qags .or. &
        INT_radquad_type == INT_qags) then
       write (uout,'("  Using the QUADPACK library ")')
       write (uout,'("  R. Piessens, E. deDoncker-Kapenga, C. Uberhuber and D. Kahaner,")')
       write (uout,'("  Quadpack: a subroutine package for automatic integration, Springer-Verlag 1983.")')
    end if

    write (uout,'("+ Initial number of tetrahedra in IWS: ",A)') string(nt_orig)

    if (ws_scale > 0d0) then
       vtotal = sy%c%omega / ws_scale**3
    else
       vtotal = sy%c%omega
    end if
    write (uout,'("  Volume of the primitive unit cell: ",A)') string(sy%c%omega,'f',decimal=6)
    write (uout,'("  Volume of the integration region: ",A)') string(vtotal,'f',decimal=6)
    write (uout,'("  Shortest tetrahedron side: ",A)') string(minlen,'f',decimal=6)
    write (uout,'("  Longest tetrahedron side (estimated): ",A)') string(maxlen,'f',decimal=6)
    write (uout,'("  QTREE proj. distance to grid points: ",A)') &
       string(min(minlen / 2**maxl / qtreefac,0.1d0),'e',decimal=7)

    ! Check consistency of the symmetry
    call qtree_checksymmetry()

    write (uout,'("  Calculating the tetrahedra contacts?: ",L)') docontacts
    write (uout,'("  Assuming periodic cell?: ",L)') periodic

    ! contact info
    if (docontacts) then
       ! write (uout,'("* CONTACTS of the tetrahedra faces")')
       ! write (uout,'("* Permutations: ")')
       ! write (uout,'("  1(123), 2(231), 3(312), 4(132), 5(321), 6(213)")')
       do i = 1, nt_orig
          ! write (uout,'("* Tetrahedron: ", A)') string(i)
          do j = 1, 4
             al1 = 24 * leqv
             al2 = 6 * leqv
             al3 = leqv
             al4 = leqv
             idum = tcontact(i,j)
             if (idum == tcontact_void) cycle
             invop = (idum < 0)
             idum = abs(idum) - 1
             tt = idum / al1 + 1
             idum = idum - (tt-1)*al1
             ff = idum / al2 + 1
             idum = idum - (ff-1)*al2
             pp = idum / al3 + 1
             idum = idum - (pp-1)*al3
             cid = idum / al4 + 1
             idum = idum - (cid-1)*al4
             opid = idum + 1
             ! write (uout,'("*  Face ",A," contact: t=",A," f=",A," p=",A," cv=",A," op=",A," inv? ",L1)') &
             !    string(j), string(tt), string(ff), string(pp), string(cid), string(opid), invop
          end do
       end do
       ! write (uout,*)
    end if

    ! max level output
    write (uout,'("+ Integrating to max. level: ",A)') string(maxl)
    write (uout,'("  Tetrahedra pre-split level: ",A)') string(plvl)
    ! allocate output
    nalloc = size(trm,2)
    siz = size(trm,1)
    ! write (uout,'("* The COLOR_ALLOCATE flag is: ",A)') string(color_allocate)
    ! write (uout,'("* Color vectors (trm):")')
    ! write (uout,'("  + Number allocated: ",A)') string(nalloc)
    ! write (uout,'("  + Elements per vector: ",A)') string(siz)
    ! write (uout,'("  + Integer kind of elements: ",A)') string(qtreei)
    ! write (uout,*)
    memsum = nalloc * siz * qtreei

    ralloc = 0
    if (allocated(fgr)) ralloc = ralloc + size(fgr,2)
    if (allocated(lapgr)) ralloc = ralloc + size(lapgr,2)
    if (allocated(vgr)) ralloc = ralloc + 1
    siz = 0
    if (allocated(fgr)) then
       siz = size(fgr,1)
    else if (allocated(vgr)) then
       siz = size(vgr)
    end if
    ! write (uout,'("* F vector in memory?: ",L)') savefgr
    ! write (uout,'("* DEL2F vector in memory?: ",L)') savelapgr
    ! write (uout,'("* CORNER VOLUME vector in memory?: ",L)') allocated(vgr)
    ! write (uout,'("* Real vectors (fgr, lapgr, vgr): ")')
    ! write (uout,'("  + Number allocated: ",A)') string(ralloc)
    ! write (uout,'("  + Elements per vector: ",A)') string(siz)
    ! write (uout,'("  + Real kind of elements: ",A)') string(qtreer)
    ! write (uout,*)
    memsum = memsum + ralloc * siz * qtreer

    ! write (uout,'("* Estimated memory requirement : ",A," MB")') string(memsum / 1024**2,'f',decimal=2)
    ! write (uout,*)

    !!!!! THE TETRAHEDRA INTEGRATION STARTS HERE !!!!!

    ! Integrate the sphere properties
    ! write (uout,'("* BETA sphere integration")')
    sphereprop = 0d0
    do i = 1, nnuc
       xx = sy%c%x2c(sy%f(sy%iref)%cp(i)%x)

       sphereprop(:,i) = 0d0
       if (all(integ_mode(minl+1:maxl) /= 0) .and. r_betaint(i) > eps) then
          if (INT_spherequad_type == INT_gauleg) then
             call sphereintegrals_gauleg(xx,r_betaint(i), &
                INT_spherequad_ntheta,INT_spherequad_nphi,sphereprop(:,i),&
                abserr,neval,meaneval)
          else
             call sphereintegrals_lebedev(xx,r_betaint(i),INT_spherequad_nleb,&
                sphereprop(:,i),abserr,neval,meaneval)
          end if
       end if
       write (uout,'("+ Integrating the beta-sphere for ncp: ",A)') string(i)
       ! write (uout,'(" Non-equivalent ncp: ",A)') string(i)
       ! write (uout,'(" NCP at: ",3(A,"  "))') (string(cp(i)%x(j),'f',decimal=6),j=1,3)
       ! write (uout,'(" Sphere factor : ",A)') string(sphfactor(i),'f',decimal=6)
       ! write (uout,'(" Beta-sphere radius for gradpaths : ",A)') string(r_betagp(i),'f',decimal=6)
       ! write (uout,'(" Beta-sphere radius for integrations : ",A)') string(r_betaint(i),'f',decimal=6)
       ! write (uout,'(" Checking beta-sphere in basin?: ",L)') checkbeta
       ! if (all(integ_mode(minl+1:maxl) /= 0) .and. qtree_active(i) .and. r_betaint(i) > eps) then
       !    write (uout,'(" Integrated radial error (",A,"): ",A)') string(pname), string(abserr,'e',decimal=6)
       !    if (abserr > 1d-1 .and. INT_radquad_errprop == 2) then
       !       write (uout,'(/"! The adaptive radial integration error is too large.")')
       !       write (uout,'("! This is usually caused by: i) discontinuities in the rmt, and ")')
       !       write (uout,'("! ii) beta-spheres that are too large and include part of the interstitial  ")')
       !       write (uout,'("! To avoid this problem, decrease beta-sphere size with sphfactor."/)')
       !       call ferror('qtree_integration','Radial integration error is too large',warning)
       !    end if
       !    write (uout,'(" Number of evaluations: ",A)') string(neval)
       !    write (uout,'(" Avg. evaluations per ray: ",A)') string(meaneval)
       !    write (uout,'("id    property    Integral (sph.)")')
       !    write (uout,'(35("-"))')
       !    do j = 1, sy%npropi
       !       write (uout,'(99(A," "))') string(j,length=3,justify=ioj_left),&
       !          string(sy%propi(j)%prop_name,length=10,justify=ioj_center),&
       !          string(sphereprop(j,i),'e',decimal=10,length=18,justify=6)
       !    end do
       ! end if
       ! write (uout,*)
    end do

    ! mark the grid point inside the beta-spheres, if
    ! all the color arrays have been allocated
    if (color_allocate == 1) then
       write (uout,'("+ Marking inside-spheres grid points")')
       write (uout,*)
       do tt = 1, nt_orig
          call paint_inside_spheres(tt,tt,trm)
       end do
    end if

    ! open the differences tess file
    if (gradient_mode < 0) then
       write (roottess,'(A,A,I2.2)') trim(fileroot),"_diffterm",maxl
       call open_difftess(roottess)
       write (uout,'(A,A)') "* DIFF. TESS file opened : ", trim(roottess) // ".tess"
       write (uout,*)
    end if

    ! open the sticks files
    if (plot_mode > 0 .and. plotsticks) then
       do i = 0, maxl
          write (roottess,'(A,A,I2.2,A,I2.2)') trim(fileroot),"_level",maxl,".",i
          lustick(i) = fopen_write(trim(roottess) // ".stick")
          write (uout,'(A,A)') "* STICK file opened : ", trim(roottess) // ".stick"
       end do
       write (uout,*)
    end if

    ! zero the atomic properties array
    atprop = 0d0

    ! header
    write (uout,'("th | Task")')
    write (uout,'("-------------------------------------------------------")')

    ! integrate recursively each IWST
    ! allocatable arrays in private -> openmp 3.0 specs
    !$omp parallel do private(acum_atprop,iv0) firstprivate(trm,fgr,lapgr,vgr)
    do tt = 1, nt_orig

       !$omp critical (IO)
       write (uout,'(I3,"|",A)')  tt, " Starting QTREE "
       !$omp end critical (IO)
       if (tvol(tt) < vcutoff) then
          !$omp critical (IO)
          write (uout,'(I3,"|",A)')  tt, " The volume of the tetrahedron is lower than the cutoff, skipped. "
          !$omp end critical (IO)
          cycle
       end if

       ! nullify the accumulator
       acum_atprop = 0d0

       ! mark the grid point inside the beta-spheres, if
       ! only one trm has been allocated
       if (color_allocate == 0) then
          trm(:,1) = 0
          if (savefgr) fgr(:,1) = -1d0
          if (savelapgr) lapgr(:,1) = 0d0
          !$omp critical (IO)
          write (uout,'(I3,"|",A)')  tt, " Marking inside-sphere grid points. "
          !$omp end critical (IO)
          call paint_inside_spheres(tt,1,trm)
       end if
       ! nullify the corner volume array (vgr)
       if (allocated(vgr)) vgr(:) = 0d0

       !$omp critical (IO)
       write (uout,'(I3,"|",A)')  tt, " Finding grid vertex colors."
       !$omp end critical (IO)
       iv0(:,1) = (/ 0, 0, 0 /)
       iv0(:,2) = (/ 1, 0, 0 /)
       iv0(:,3) = (/ 0, 1, 0 /)
       iv0(:,4) = (/ 0, 0, 1 /)

       ! testing gradient_mode, parallelization does not enter these.
       if (gradient_mode == -1) then
          ndiff = 0
          ngrd1 = 0
          ngrd2 = 0
          write (uout,'(A)') "* Comparing grid point assign. to attractors"
          write (uout,'(A)') "* trm1 = actual ode, full gradient | trm2 = DP45 embedded, full gradient "
          write (uout,'(A2,"         ",A,"                     ",A,"              ",A4," ",A4)') &
             "tt", "--- x_convex ---", "--- x_crys ---", "trm1", "trm2"
       else if (gradient_mode == -2) then
          ndiff = 0
          ngrd1 = 0
          ngrd2 = 0
          write (uout,'(A)') "* Comparing grid point assign. to attractors"
          write (uout,'(A)') "* trm1 = actual ode, color | trm2 = DP45 embedded, full gradient "
          write (uout,'(A2,"         ",A,"                     ",A,"              ",A4," ",A4)') &
             "tt", "--- x_convex ---", "--- x_crys ---", "trm1", "trm2"
       else if (gradient_mode == -3) then
          ndiff = 0
          ngrd1 = 0
          ngrd2 = 0
          write (uout,'(A)') "* Comparing grid point assign. to attractors"
          write (uout,'(A)') "* trm1 = actual ode, qtree | trm2 = DP45 embedded, full gradient "
          write (uout,'(A2,"         ",A,"                     ",A,"              ",A4," ",A4)') &
             "tt", "--- x_convex ---", "--- x_crys ---", "trm1", "trm2"
       end if

       ! recurisve subdivision
       call tetrah_subdivide(tt,iv0,0,acum_atprop,trm,fgr,lapgr,vgr)

       ! deferred corner integration -- sum
       if (intcorner_deferred) then
          !$omp critical (IO)
          write (uout,'(I3,"|",A)')  tt, " Deferred CORNER integration, accumulating properties."
          !$omp end critical (IO)
          call integ_corner_sum(tt,trm,vgr,acum_atprop)
       end if

       ! output gradient differences in debug mode, multi-threads do not enter this.
       if (gradient_mode < 0) then
          write (uout,'(A,I7,A,I7)') "* Number of term. differences : ", ndiff, " of ", siz
          write (uout,'(A,I7)') "* ngrd of method 1 : ", ngrd1
          write (uout,'(A,I7)') "* ngrd of method 2 : ", ngrd2
          write (uout,*)
       end if

       ! Copy trm, fgr and lapgr across tetrahedron faces. Only single-thread.
       if (docontacts .and. color_allocate == 1) then
          write (uout,'(I3,"|",A)')  tt, " Copying contacts."
          do ff = 1, 4
             if (tcontact(tt,ff) == tcontact_void) cycle
             al1 = 24 * leqv
             al2 = 6 * leqv
             al3 = leqv
             al4 = leqv
             idum = tcontact(tt,ff)
             invop = (idum < 0)
             idum = abs(idum) - 1
             tto = idum / al1 + 1
             idum = idum - (tto-1)*al1
             ffo = idum / al2 + 1
             idum = idum - (ffo-1)*al2
             pp = idum / al3 + 1
             idum = idum - (pp-1)*al3
             cid = idum / al4 + 1
             idum = idum - (cid-1)*al4
             opid = idum + 1
             if (tt == tto .and. ff == ffo) cycle

             ! Determine origin face
             if (ff == 1) then
                xface_orig(:,1) = 0d0
                xface_orig(:,2) = xface_orig(:,1) + bvec(:,1,tt) * l2
                xface_orig(:,3) = xface_orig(:,1) + bvec(:,2,tt) * l2
             else if (ff == 2) then
                xface_orig(:,1) = 0d0
                xface_orig(:,2) = xface_orig(:,1) + bvec(:,1,tt) * l2
                xface_orig(:,3) = xface_orig(:,1) + bvec(:,3,tt) * l2
             else if (ff == 3) then
                xface_orig(:,1) = 0d0
                xface_orig(:,2) = xface_orig(:,1) + bvec(:,2,tt) * l2
                xface_orig(:,3) = xface_orig(:,1) + bvec(:,3,tt) * l2
             else
                xp = 0d0
                xface_orig(:,1) = xp + bvec(:,1,tt) * l2
                xface_orig(:,2) = xp + bvec(:,2,tt) * l2
                xface_orig(:,3) = xp + bvec(:,3,tt) * l2
             end if
             xface_orig = matmul(xface_orig,perm3(:,:,pp)) ! permute before calculating convex coordinates
             rface_orig = matmul(cmat(:,:,tt),xface_orig)
             do i = 1, 3
                xface_orig(:,i) = sy%c%c2x(xface_orig(:,i))
                xface_orig(:,i) = xface_orig(:,i) + torig(:,tt)
                if (.not.invop) then
                   xface_orig(:,i) = matmul(lrotm(1:3,1:3,opid),xface_orig(:,i))
                   xface_orig(:,i) = xface_orig(:,i) + sy%c%cen(:,cid)
                else
                   xface_orig(:,i) = xface_orig(:,i) + sy%c%cen(:,cid)
                   xface_orig(:,i) = matmul(lrotm(1:3,1:3,opid),xface_orig(:,i))
                end if
             end do

             ! Determine end face
             if (ffo == 1) then
                xface_end(:,1) = 0d0
                xface_end(:,2) = xface_end(:,1) + bvec(:,1,tto) * l2
                xface_end(:,3) = xface_end(:,1) + bvec(:,2,tto) * l2
             else if (ffo == 2) then
                xface_end(:,1) = 0d0
                xface_end(:,2) = xface_end(:,1) + bvec(:,1,tto) * l2
                xface_end(:,3) = xface_end(:,1) + bvec(:,3,tto) * l2
             else if (ffo == 3) then
                xface_end(:,1) = 0d0
                xface_end(:,2) = xface_end(:,1) + bvec(:,2,tto) * l2
                xface_end(:,3) = xface_end(:,1) + bvec(:,3,tto) * l2
             else
                xp = 0d0
                xface_end(:,1) = xp + bvec(:,1,tto) * l2
                xface_end(:,2) = xp + bvec(:,2,tto) * l2
                xface_end(:,3) = xp + bvec(:,3,tto) * l2
             end if
             rface_end = matmul(cmat(:,:,tto),xface_end)

             do i = 1, 3
                xface_end(:,i) = sy%c%c2x(xface_end(:,i))
                xface_end(:,i) = xface_end(:,i) + torig(:,tto)
             end do

             ok = .true.
             do i = 1, 3
                face_diff = xface_orig(:,i) - xface_end(:,i)
                face_diff = face_diff - nint(face_diff)
                ok = ok .and. all(abs(face_diff) < eps_tetrah_contact)
             end do
             ! write (uout,'("+ Copying from tetrah = ",I3," face = ",I2)') tt, ff
             ! write (uout,'("            to tetrah = ",I3," face = ",I2)') tto, ffo
             ! write (uout,'("  Permutation = ",I2)') pp
             ! write (uout,'("  Rotation matrix = ",I2)') opid
             ! write (uout,'("  Centering vector = ",I2)') cid
             ! write (uout,'("  Invert rotation/centering? = ",L2)') invop
             ! write (uout,'("  Face consistency check : ",L2)') ok
             ! write (uout,*)
             if (.not. ok) then
                write (uout,'(" xface_orig : ")')
                do i = 1, 3
                   write (uout,'("    ",1p,3(E20.12,"  "))') xface_orig(:,i)
                end do
                write (uout,'(" xface_end : ")')
                do i = 1, 3
                   write (uout,'("    ",1p,3(E20.12,"  "))') xface_end(:,i)
                end do
                call ferror('qtree','faces in contact have non-consistent coordinates',warning)
                cycle
             end if

             do i = 0, l2
                do j = 0, l2-i
                   vin = nint((l2-i-j)*rface_orig(:,1) + i*rface_orig(:,2) + j*rface_orig(:,3))
                   vino = nint((l2-i-j)*rface_end(:,1) + i*rface_end(:,2) + j*rface_end(:,3))
                   idx = cindex(vin,maxl)
                   idxo = cindex(vino,maxl)
                   if (idx > size(trm,1) .or. idx < 1) call ferror('qtree_integration','trm exceeded',faterr)
                   if (idxo > size(trm,1) .or. idx < 1) call ferror('qtree_integration','trm exceeded',faterr)
                   trm(idxo,tto) = trm(idx,tt)
                   if (savefgr) fgr(idxo,tto) = fgr(idx,tt)
                   if (savelapgr) lapgr(idxo,tto) = lapgr(idx,tt)
                end do
             end do
          end do
       end if

       ! copy values from the accumulator
       ! nnuc+1 are on-CP grid points
       ! nnuc+2 are points belonging to inactive basins
       ! nnuc+3 are gradient path errors
       !$omp critical (IO)
       write (uout,'(I3,"|",A)')  tt, " Summing atomic contributions. "
       !$omp end critical (IO)
       !$omp critical (atprop1)
       do i = 1, nnuc+3
          atprop(i,:) = atprop(i,:) + acum_atprop(i,:)
       end do
       !$omp end critical (atprop1)
    end do ! tetrahedra color assign
    !$omp end parallel do

    write (uout,'("-------------------------------------------------------"/)')

    ! close the differences tess file
    if (gradient_mode < 0) then
       write (roottess,'(A,A,I2.2)') trim(fileroot),"_diffterm",maxl
       call close_difftess(roottess)
       write (uout,'(A,A)') "* DIFF. TESS file generated : ", trim(roottess) // ".tess"
       write (uout,*)
    end if
    ! close the sticks files
    if (plot_mode > 0 .and. plotsticks) then
       do i = 0, maxl
          write (roottess,'(A,A,I2.2,A,I2.2)') trim(fileroot),"_level",maxl,".",i
          call fclose(lustick(i))
          write (uout,'(A,A)') "* STICK file generated : ", trim(roottess) // ".stick"
       end do
       write (uout,*)
    end if

    if (plot_mode > 0) then
       if (color_allocate == 1) then
          ! whole cell
          write (roottess,'(A,A,I2.2)') trim(fileroot),"_level",maxl
          call small_writetess(roottess,0,trm)
          write (uout,'(A,A)') "* TESS file generated : ", trim(roottess) // ".tess"
          write (uout,*)
          ! per-atom
          if (plot_mode == 4 .or. plot_mode == 5) then
             do i = 1, nnuc+1
                write (roottess,'(A,A,I2.2,A,I2.2)') trim(fileroot),"_level",maxl,"_",i
                call small_writetess(roottess,i,trm)
                write (uout,'(A,A)') "* TESS file generated : ", trim(roottess) // ".tess"
                write (uout,*)
             end do
          end if
       else
       end if
    end if

    ! scale integrals and sum spheres
    do i = 1, nnuc
       atprop(i,:) = atprop(i,:) * leqvf / sy%f(sy%iref)%cp(i)%mult
       atprop(i,2:sy%npropi) = atprop(i,2:sy%npropi) + sphereprop(2:sy%npropi,i)
    end do

    ! output the results
    bas%nattr = nnuc
    allocate(bas%icp(nnuc),bas%xattr(3,nnuc),bas%docelatom(0:sy%f(sy%iref)%ncpcel))
    bas%docelatom = .true.
    k = 0
    do i = 1, nnuc
       k = k + 1
       do j = 1, sy%npropi
          res(j)%psum(i) = atprop(i,j)
       end do
       do j = 1, sy%f(sy%iref)%ncpcel
          if (sy%f(sy%iref)%cpcel(j)%idx == i) then
             bas%icp(k) = j
             bas%xattr(:,k) = sy%f(sy%iref)%cpcel(j)%x
             exit
          end if
       end do
    end do

    call int_output_header(bas,res,.true.,.true.)
    call int_output_fields(bas,res,.true.,.true.)
    deallocate(bas%icp,bas%xattr,bas%docelatom,res)

    ! clean up
    if (.not.quiet) call tictac("End QTREE")
    call qtree_cleanup()
    if (allocated(trm)) deallocate(trm)
    if (allocated(fgr)) deallocate(fgr)
    if (allocated(lapgr)) deallocate(lapgr)
    if (allocated(vgr)) deallocate(vgr)
    if (allocated(acum_atprop)) deallocate(acum_atprop)

  end subroutine qtree_integration

  !xx! private procedures

  !> Set sphere sizes according to user's input or, alternatively, calculate them
  !> by analyzing the system at a smaller level (lvl).
  subroutine qtree_setsph1(lvl,verbose)
    use systemmod, only: sy
    use qtree_tetrawork, only: term_rec
    use qtree_basic, only: qtreeidx, qtreei, qtreer, torig, tvec, maxlen, nnuc,&
       nt_orig, r_betagp, tvol, cindex, r_betaint, find_beta_rodriguez, get_tlengths,&
       qtree_initialize, qtree_checksymmetry, qtree_cleanup, sphfactor, sphintfactor,&
       color_allocate, vcutoff
    use fieldmod, only: type_elk, type_wien
    use global, only: rbetadef
    use tools_io, only: ferror, faterr, uout, string

    integer, intent(in) :: lvl
    logical, intent(in) :: verbose

    integer :: l2, i, j, k, l, tt, tto, idum
    integer :: nid, nid0, vin(3)
    integer(qtreeidx) :: idx
    real*8 :: xx(3), dist, mdist, rmin, rmax, rdist, rmt, rnn2
    real*8, allocatable :: rref(:)
    logical, allocatable :: nucmask(:), nfrozen(:)
    integer :: niter, icolor
    logical :: docalc
    integer :: save_cl
    ! data arrays
    integer(qtreei), allocatable :: trm(:,:)
    real(qtreer), allocatable :: fgr(:,:), lapgr(:,:), vgr(:)
    real*8, allocatable :: acum_atprop(:,:)

    real*8, parameter :: initial_f = 1.0d0
    real*8, parameter :: shrink = 0.95d0
    integer, parameter :: miter = 400

    ! save the color_allocate variable
    icolor = color_allocate
    color_allocate = 1

    ! get the nnuc
    nnuc = 0
    do i = 1, sy%f(sy%iref)%ncp
       if (sy%f(sy%iref)%cp(i)%typ == sy%f(sy%iref)%typnuc) nnuc = nnuc + 1
    end do

    ! allocate arrays
    allocate(nucmask(nnuc),nfrozen(nnuc),rref(nnuc))
    l2 = 2**lvl

    ! Calculate sphere sizes (r_betagp and r_betaint)
    nfrozen = .false.
    do i = 1, nnuc
       rnn2 = sy%c%get_rnn2(i)
       if (i<=sy%c%nneq .and. sy%f(sy%iref)%type == type_elk) then
          rmt = sy%f(sy%iref)%elk%rmt(sy%c%at(i)%is)
       elseif (i<=sy%c%nneq .and. sy%f(sy%iref)%type == type_wien) then
          rmt = sy%f(sy%iref)%wien%rmt_atom(sy%c%at(i)%x)
       else
          rmt = rnn2
       end if

       if (i <= sy%c%nneq) then
          rref(i) = rnn2
       else
          xx = sy%f(sy%iref)%cp(i)%x
          call sy%f(sy%iref)%nearest_cp(xx,idum,rref(i),type=sy%f(sy%iref)%typnuc,nozero=.true.)
          rref(i) = rref(i) / 2d0
       end if

       if (sphfactor(i) > 1d-12) then
          r_betagp(i) = min(sphfactor(i) * rref(i),rmt)
          nfrozen(i) = .true.
       else if (sphfactor(i) < -1d-12) then
          r_betagp(i) = min(abs(sphfactor(i)) * rref(i),rmt)
          call find_beta_rodriguez(i,r_betagp(i))
          nfrozen(i) = .true.
       else
          ! calculate using the tetrahedral grid
          r_betagp(i) = min(initial_f * rref(i),rmt)
          nfrozen(i) = .false.
       end if
    end do
    sphfactor(i) = r_betagp(i) / rref(i)

    docalc = any(.not.nfrozen)
    if (docalc) then
       ! initialize, avoid allocation of only 1 iwst
       save_cl = color_allocate
       call qtree_initialize(lvl,0,acum_atprop,trm,fgr,lapgr,vgr,.false.)
       call qtree_checksymmetry()
       color_allocate = save_cl

       ! Minimum sphere size
       mdist = Rbetadef

       ! minlen and maxlen
       call get_tlengths(rmin,rmax)
       write (uout,'("+ Calculated min. length: ",F12.6)') rmin
       write (uout,'("+ Calculated max. length: ",F12.6)') rmax
       write (uout,'("+ Estimated max. length: ",F12.6)') maxlen/l2
       write (uout,'("+ Minimum beta-sphere radius: ",F12.6)') mdist

       ! use the geometric mean of rmax and maxlen/l2
       !rdist = sqrt(rmax * maxlen / l2)
       rdist = rmax
       write (uout,'("+ Using safe-shell distance: ",F12.6)') rdist
       write (uout,'("  This distance does NOT ensure good spheres, but works most of the times.")')
       write (uout,'("  If a small sphere error is obtained, either reduce the spheres manually")')
       write (uout,'("  or use a higher qtree level.")')
       write (uout,*)

       ! Calculate un-frozen radii
       if (verbose) &
          write (uout,'("* Calculating BETA-SPHERE radii")')

       niter = 0
       do while (any(.not.nfrozen))
          niter = niter + 1
          if (niter > miter) call ferror('qtree_setsph1','too many iterations',faterr)
          if (verbose) &
             write (uout,'("+ Pass ",I3)') niter

          nucmask = .not.nfrozen
          do tt = 1, nt_orig
             if (tvol(tt) < vcutoff) cycle
             if (color_allocate == 0) then
                tto = 1
             else
                tto = tt
             end if

             do i = 0, l2
                do j = 0, l2-i
                   do k = 0, l2-i-j
                      vin = (/i,j,k/)
                      xx = torig(:,tt)
                      do l = 1, 3
                         xx = xx + tvec(:,l,tt) * real(vin(l),8) / real(l2,8)
                      end do

                      do nid = 1, nnuc
                         if (.not.nucmask(nid)) cycle
                         ! check if this point is inside a shell around the nucleus
                         call sy%f(sy%iref)%nearest_cp(xx,nid0,dist,nid0=nid)
                         if (dist > r_betagp(nid)-mdist .and. dist < r_betagp(nid)+min(rdist,r_betagp(nid))) then
                            idx = cindex(vin,lvl)
                            trm(idx,tto) = int(term_rec(tt,vin,lvl,trm,fgr,lapgr),1)
                            if (abs(trm(idx,tto)) /= nid) then
                               r_betagp(nid) = r_betagp(nid) * shrink
                               sphfactor(nid) = sphfactor(nid) * shrink
                               if (verbose) write (uout,'("+ Shrinking sphere of nuc ",I4," to ",F14.6)') nid, r_betagp(nid)
                               if (r_betagp(nid) < mdist) then
                                  r_betagp(nid) = mdist
                                  sphfactor(nid) = r_betagp(nid) / rref(nid)
                                  if (verbose) write (uout,&
                                     '("+ Freezing nuc ",I4," and setting small beta-sphere -- r= ",F14.6)')&
                                     nid, r_betagp(nid)
                                  nfrozen(nid) = .true.
                               end if
                               nucmask(nid) = .false.
                               cycle
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
          nfrozen = nfrozen .or. nucmask
       end do
       write (uout,*)
    end if

    if (verbose) then
       write (uout,'("* BETA-SPHERE sizes calculation, final radii")')
    end if
    do i = 1, nnuc
       r_betaint(i) = sphintfactor(i) * r_betagp(i)
       if (verbose) then
          rnn2 = sy%c%get_rnn2(i)
          write (uout,'("+ Sphfactor/rbeta/rnn2 of nuc ",A," (",A,") :",3("  ",A))') &
             string(i), string(sy%c%at(i)%name), string(sphfactor(i),'f',decimal=6), &
             string(r_betagp(i),'f',decimal=6), string(rnn2,'f',decimal=6)
       end if
    end do
    if (verbose) then
       write (uout,'("+ Sphfactor values to be used in input")')
       do i = 1, nnuc
          write (uout,'("sphfactor ",A,"  ",A)') string(i), string(sphfactor(i),'f',decimal=6)
       end do
       write (uout,*)
    end if

    ! recover the color allocate
    color_allocate = icolor

    ! clean up
    deallocate(nucmask,nfrozen,rref)

    if (docalc) call qtree_cleanup()
    if (allocated(trm)) deallocate(trm)
    if (allocated(fgr)) deallocate(fgr)
    if (allocated(lapgr)) deallocate(lapgr)
    if (allocated(vgr)) deallocate(vgr)
    if (allocated(acum_atprop)) deallocate(acum_atprop)

  end subroutine qtree_setsph1

  !> Set sphere sizes according to user's input or, alternatively, calculate them
  !> by analyzing the system at a smaller level (lvl).
  subroutine qtree_setsph2(verbose)
    use systemmod, only: sy
    use surface, only: minisurf
    use qtree_basic, only: nnuc, r_betagp, r_betaint, find_beta_rodriguez, sphfactor, sphintfactor
    use fieldmod, only: type_elk, type_wien
    use tools_io, only: uout, string
    logical, intent(in) :: verbose

    type(minisurf) :: srf
    integer :: i, j, idum, ii
    real*8 :: xx(3), dist, plen, rnn2
    real*8 :: rref, x0(3), unit(3), rmt
    integer :: ier
    integer :: ndo, nstep
    integer, allocatable :: ido(:)
    logical :: doagain

    real*8, parameter :: initial_f = 1.0d0
    real*8, parameter :: shrink = 0.95d0
    integer, parameter :: nleb = 38
    real*8, parameter :: dthres = 0.5d0

    ! get the nnuc
    nnuc = 0
    do i = 1, sy%f(sy%iref)%ncp
       if (sy%f(sy%iref)%cp(i)%typ == sy%f(sy%iref)%typnuc) nnuc = nnuc + 1
    end do

    ! allocate arrays
    allocate(ido(nnuc))

    ! Calculate sphere sizes (r_betagp and r_betaint)
    ndo = 0
    do i = 1, nnuc
       rnn2 = sy%c%get_rnn2(i)
       if (i<=sy%c%nneq .and. sy%f(sy%iref)%type == type_elk) then
          rmt = sy%f(sy%iref)%elk%rmt(sy%c%at(i)%is)
       elseif (i<=sy%c%nneq .and. sy%f(sy%iref)%type == type_wien) then
          rmt = sy%f(sy%iref)%wien%rmt_atom(sy%c%at(i)%x)
       elseif (i<=sy%c%nneq) then
          rmt = rnn2
       else
          rmt = 1d30
       end if

       if (i <= sy%c%nneq) then
          rref = rnn2
       else
          xx = sy%f(sy%iref)%cp(i)%x
          call sy%f(sy%iref)%nearest_cp(xx,idum,rref,type=sy%f(sy%iref)%typnuc,nozero=.true.)
          rref = rref / 2d0
       end if

       if (sphfactor(i) > 1d-12) then
          r_betagp(i) = min(sphfactor(i) * rref,rmt)
       else if (sphfactor(i) < -1d-12) then
          r_betagp(i) = min(abs(sphfactor(i)) * rref,rmt)
          call find_beta_rodriguez(i,r_betagp(i))
       else
          ! calculate using the tetrahedral grid
          r_betagp(i) = min(initial_f * rref,rmt)
          ndo = ndo + 1
          ido(ndo) = i
       end if
       sphfactor(i) = r_betagp(i) / rref
    end do

    call srf%init(nleb,0)
    call srf%clean()
    call srf%lebedev_nodes(nleb)
    !$omp parallel do &
    !$omp private(i,x0,rref,idum,doagain,j,unit,xx,nstep,ier,dist,plen)
    do ii = 1, ndo
       i = ido(ii)
       if (i <= sy%c%nneq) then
          x0 = sy%c%at(i)%r
          rref = sy%c%get_rnn2(i)
       else
          x0 = sy%f(sy%iref)%cp(i)%r
          call sy%f(sy%iref)%nearest_cp(xx,idum,rref,type=sy%f(sy%iref)%typnuc,nozero=.true.)
          rref = rref / 2d0
       end if

       doagain = .true.
       do while(doagain)
          doagain = .false.
          do j = 1, srf%nv
             unit = (/ sin(srf%th(j)) * cos(srf%ph(j)),&
                sin(srf%th(j)) * sin(srf%ph(j)),&
                cos(srf%th(j)) /)
             xx = x0 + r_betagp(i) * unit
             call sy%f(sy%iref)%gradient(xx,+1,nstep,ier,.true.,plen)
             dist = norm2(xx - x0)
             if (dist > dthres) then
                !$omp critical (write)
                r_betagp(i) = r_betagp(i) * shrink
                sphfactor(i) = r_betagp(i) / rref
                !$omp end critical (write)
                doagain = .true.
                exit
             end if
          end do
       end do
    end do
    !$omp end parallel do
    call srf%end()

    if (verbose) then
       write (uout,'("* BETA-SPHERE sizes calculation, final radii")')
    end if
    do i = 1, nnuc
       r_betaint(i) = sphintfactor(i) * r_betagp(i)
       if (verbose) then
          if (i <= sy%c%nneq) then
             rnn2 = sy%c%get_rnn2(i)
             write (uout,'("+ Sphfactor/rbeta/rnn2 of nuc ",A," (",A,") :",3(A,"  "))') &
                string(i), string(sy%c%at(i)%name), string(sphfactor(i),'f',decimal=6), &
                string(r_betagp(i),'f',decimal=6), string(rnn2,'f',decimal=6)
          else
             write (uout,'("+ Sphfactor/rbeta/rnn2 of nuc ",A," (",A,") :",3(A,"  "))') &
                string(i), "nnm", string(sphfactor(i),'f',decimal=6), &
                string(r_betagp(i),'f',decimal=6), string(0d0,'f',decimal=6)
          end if
       end if
    end do
    if (verbose) then
       write (uout,'("+ Sphfactor values to be used in input")')
       do i = 1, nnuc
          write (uout,'("sphfactor ",A,"  ",A)') string(i), string(sphfactor(i),'f',decimal=6)
       end do
       write (uout,'("+ End input")')
       write (uout,*)
    end if

    ! clean up
    deallocate(ido)

  end subroutine qtree_setsph2

end submodule proc
