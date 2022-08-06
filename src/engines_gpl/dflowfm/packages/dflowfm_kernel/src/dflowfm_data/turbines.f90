
module m_rdturbine
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2013.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  Se e the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: turbines.f90 66081 2020-02-27 17:01:24Z ccchart.x $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/SANDIA/fm_tidal_v3/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/turbines.f90 $
!!--declarations----------------------------------------------------------------
    implicit none
    
contains

subroutine rdturbine(filtrb, lundia, turbines, error)
!!--description-----------------------------------------------------------------
!
! Read turbine data from input file
!
!!--declarations----------------------------------------------------------------
    use precision
    use properties ! includes tree_structures
    !use grid_dimens_module, only: griddimtype
    use m_structures, only: structure_turbines, allocate_turbines
    use message_module, only: write_error, write_warning, FILE_NOT_FOUND, FILE_READ_ERROR, PREMATURE_EOF
    use table_handles, only: readtable, gettable, GETTABLE_NAME
    use mathconsts, only: pi
    use m_sferic, only: dg2rd
    !
    implicit none
!
! Call variables
!
!    type (griddimtype)          , pointer     :: griddim
    character(*)                , intent(in)  :: filtrb
    integer                     , intent(in)  :: lundia
    type(structure_turbines)    , intent(inout) :: turbines
    logical                     , intent(out) :: error
!
! Local variables
!
    integer                     , external    :: newunit
    !
    integer                                   :: istat
    integer                                   :: j
    integer                                   :: itrb
    integer                                   :: nval
    real(fp), dimension(2)                    :: xyloc
    character(10)                             :: tempstring
    character(20)                             :: parname
    character(256)                            :: curvename
    character(256)                            :: trbthrustfil
    character(1024)                           :: message
    type (tree_data)            , pointer     :: aturbine
    type (tree_data)            , pointer     :: turbine_tree
!
!! executable statements -------------------------------------------------------
!
    error = .false.
    call tree_create  ( "Turbine input", turbine_tree )
    !
    if (filtrb == ' ') then
       call write_error('Empty turbine file name specified.',unit=lundia)
       error = .true.
       return
    endif     
    !! TODO FIX HACK
    !! CCC DEBUG
    open(unit=5979,file='Turbine_Power.dat',form='formatted')
    !
    ! Read turbine-file
    !
    call prop_file('ini', trim(filtrb), turbine_tree, istat)
    if (istat /= 0) then
       select case (istat)
       case(1)
          call write_error(FILE_NOT_FOUND//trim(filtrb), unit=lundia)
       case(3)
          call write_error(PREMATURE_EOF//trim(filtrb), unit=lundia)
       case default
          call write_error(FILE_READ_ERROR//trim(filtrb), unit=lundia)
       endselect
       error = .true.
       return
    endif
    !
    call prop_get_string(turbine_tree, 'TurbineFileInformation', 'FileVersion', tempstring)
    if (trim(tempstring) /= '01.00') then
       call write_error('FileVersion should match 01.00 for turbine file.',unit=lundia)
       error = .true.
       return
    endif
    !
    trbthrustfil = ' '
    call prop_get_string(turbine_tree,'General','CurvesFil',trbthrustfil)
    if (trbthrustfil == ' ') then
       call write_error('Missing Turbine Loss Coefficients file.',unit=lundia)
       error = .true.
       return
    else
       message = ' '
       call readtable(turbines%curves, trbthrustfil, 0, message)
       if (message /= ' ') then
          call write_error(message,unit=lundia)
          error = .true.
          return
       endif
    endif
    !
    itrb = 0
    do j = 1, size(turbine_tree%child_nodes)
       !
       ! Does turbine_tree contain any child with name 'turbine' (converted to lower case)?
       !
       aturbine => turbine_tree%child_nodes(j)%node_ptr
       parname = tree_get_name( aturbine )
       if (parname == 'turbine') then
          itrb = itrb+1
       endif
    enddo
    ! 
    call allocate_turbines(turbines,itrb,lundia,error)
    if (error) return
    !
    itrb = 0
    do j = 1, size(turbine_tree%child_nodes)
       !
       ! Does turbine_tree contain any child with name 'turbine' (converted to lower case)?
       !
       aturbine => turbine_tree%child_nodes(j)%node_ptr
       parname = tree_get_name( aturbine )
       if (parname == 'turbine') then
          itrb = itrb+1
          call prop_get_string(aturbine, '*', 'Name', turbines%nr(itrb)%name)
          if (turbines%nr(itrb)%name == ' ') write(turbines%nr(itrb)%name,'(A,I0)') 'Turbine ',itrb
          !
          call prop_get(aturbine, '*', 'TurbType', turbines%nr(itrb)%turbtype) ! BJ 20150326
          call prop_get(aturbine, '*', 'Diameter', turbines%nr(itrb)%diam)     ! BJ 20150326 Need to keep this for NDiaDist4Vel
          if (turbines%nr(itrb)%turbtype == 0) then ! BJ 20150326
           turbines%nr(itrb)%turbarea = 0.25_fp*pi*turbines%nr(itrb)%diam**2
          elseif (turbines%nr(itrb)%turbtype == 1) then
             call prop_get(aturbine, '*', 'Width', turbines%nr(itrb)%width)
             call prop_get(aturbine, '*', 'Height', turbines%nr(itrb)%height)
             turbines%nr(itrb)%turbarea = turbines%nr(itrb)%width*turbines%nr(itrb)%height
          endif
          call prop_get(aturbine, '*', 'XYLoc', xyloc, 2)
          turbines%nr(itrb)%xyz(1:2) = xyloc
          call prop_get(aturbine, '*', 'Orientation', turbines%nr(itrb)%angle)
          turbines%nr(itrb)%csturb = cos(turbines%nr(itrb)%angle*dg2rd)
          turbines%nr(itrb)%snturb = sin(turbines%nr(itrb)%angle*dg2rd) 

          call prop_get(aturbine, '*', 'NDiaDist4Vel', turbines%nr(itrb)%ndiamu)
          call prop_get(aturbine, '*', 'TurbineModel', turbines%nr(itrb)%turbinemodel)
          call prop_get(aturbine, '*', 'TurbulenceModel', turbines%nr(itrb)%turbulencemodel)
          call prop_get(aturbine, '*', 'Beta_p', turbines%nr(itrb)%beta_p)
          call prop_get(aturbine, '*', 'Beta_d', turbines%nr(itrb)%beta_d)
          call prop_get(aturbine, '*', 'Cep4', turbines%nr(itrb)%cep4)
          call prop_get(aturbine, '*', 'Cep5', turbines%nr(itrb)%cep5)
          !
          tempstring = ' '
          call prop_get_string(aturbine, '*', 'VertPos', tempstring)
          call small(tempstring,10)
          if (tempstring=='fixed' .or. tempstring==' ') then
              turbines%nr(itrb)%vertpos = 0
              call prop_get(aturbine, '*', 'AxisLevel', turbines%nr(itrb)%xyz(3))
          elseif (tempstring=='floating') then
              turbines%nr(itrb)%vertpos = 1
              call prop_get(aturbine, '*', 'AxisDepth', turbines%nr(itrb)%xyz(3))
          else
              write(message,'(4A)') 'Invalid vertical position "',trim(tempstring),'" for ',turbines%nr(itrb)%name
              call write_error(message,unit=lundia)
              error = .true.
              return
          endif
          !
          call prop_get_string(aturbine, '*', 'ThrustCurve', turbines%nr(itrb)%thrustcrvname)
          call gettable(turbines%curves, turbines%nr(itrb)%thrustcrvname, 'thrust coefficient', &
                & turbines%nr(itrb)%thrustcrvnr(1), turbines%nr(itrb)%thrustcrvnr(2), &
                & nval, 1, message, GETTABLE_NAME)
          if (nval/=1) then
              write(message,'(3A)') 'Unable to find table for thrust curve "',trim(turbines%nr(itrb)%thrustcrvname),'"'
              call write_error(message,unit=lundia)
              error = .true.
              return
          endif
          !
          call prop_get_string(aturbine, '*', 'PowerCurve', turbines%nr(itrb)%powercrvname)
          if (turbines%nr(itrb)%powercrvname /= ' ') then
              call gettable(turbines%curves, turbines%nr(itrb)%powercrvname, 'power coefficient', &
                    & turbines%nr(itrb)%powercrvnr(1), turbines%nr(itrb)%powercrvnr(2), &
                    & nval, 1, message, GETTABLE_NAME)
              if (nval/=1) then
                  write(message,'(3A)') 'Unable to find table for power curve "',trim(turbines%nr(itrb)%powercrvname),'"'
                  call write_error(message,unit=lundia)
                  error = .true.
                  return
              endif
          endif
       endif
    enddo
    !
    call tree_destroy(turbine_tree)
end subroutine rdturbine


subroutine echoturbine(turbines, lundia)
!!--description-----------------------------------------------------------------
!
! Display turbine data
!
!!--declarations----------------------------------------------------------------
    use precision
    use m_structures, only: structure_turbines
    !
    implicit none
!
! Call variables
!
    integer                     , intent(in)  :: lundia
    type(structure_turbines)    , intent(in)  :: turbines
!
! Local variables
!
    integer                                   :: j
    character(30)                             :: txtput1
    character(50)                             :: txtput2
!
!! executable statements -------------------------------------------------------
!
    if (associated(turbines%nr)) then
       write (lundia, '(a)')   '*** Start  of turbine input'
       do j = 1, size(turbines%nr)
          txtput1 = 'Turbine'
          write (lundia, '(3a)') txtput1, ': ', trim(turbines%nr(j)%name)
          if (turbines%nr(j)%turbtype==0) then ! BJ 20150326
             txtput1 = '  Diameter'
             write (lundia, '(2a,f12.3)') txtput1, ': ', turbines%nr(j)%diam
          else
             txtput1 = '  Width'
             write (lundia, '(2a,f12.3)') txtput1, ': ', turbines%nr(j)%width
             txtput1 = '  Height'
             write (lundia, '(2a,f12.3)') txtput1, ': ', turbines%nr(j)%height
          endif
          !
          txtput1 = '  X coordinate'
          write (lundia, '(2a,e12.4)') txtput1, ': ', turbines%nr(j)%xyz(1)
          txtput1 = '  Y coordinate'
          write (lundia, '(2a,e12.4)') txtput1, ': ', turbines%nr(j)%xyz(2)
          txtput1 = '  Orientation of turbine'
          write (lundia, '(2a,f12.3)') txtput1, ': ', turbines%nr(j)%angle
          !
          txtput1 = '  Z position of turbine axis'
          select case(turbines%nr(j)%vertpos)
             case (0)
                write(txtput2,'(a,e12.4)') 'fixed at level =',turbines%nr(j)%xyz(3)
             case (1)
                write(txtput2,'(a,e12.4)') 'floating at depth =',turbines%nr(j)%xyz(3)
          end select
          write (lundia, '(3a)') txtput1, ': ', trim(txtput2)
          !
          txtput1 = '  Loss Curve Name'
          write (lundia, '(3a)') txtput1, ': ', trim(turbines%nr(j)%thrustcrvname)
          !
          txtput1 = '  Power Curve Name'
          if (turbines%nr(j)%powercrvname /= ' ') then
              write (lundia, '(3a)') txtput1, ': ', trim(turbines%nr(j)%powercrvname)
          else
              write (lundia, '(3a)') txtput1, ': ', 'N/A'
          endif
       enddo
       write (lundia, '(a)')   '*** End    of turbine input'
    endif
end subroutine echoturbine

subroutine mapturbine(turbines,error)
!subroutine mapturbine(turbines, xcor, ycor, guu, gvv, kcu, kcv, error, gdp)
!!--description-----------------------------------------------------------------
!
! Determine locations of turbines on grid
!
!!--declarations----------------------------------------------------------------
    use precision
    use m_structures, only: structure_turbines, structure_turbine
    use m_GlobalParameters, only: INDTP_1D, INDTP_2D, INDTP_ALL
    !use m_kdtree2
    use kdtree2Factory
    use mathconsts
    use m_flowgeom, only: Lnx,dxi,lncn
    use m_flow, only: kmx
    use network_data, only: numL, xk, yk, kn
    use geometry_module, only: dbdistance
    use m_missing, only: dmiss
    use m_sferic, only: jsferic, jasfer3D, dg2rd
    use m_alloc
    !
    implicit none
!
! Call variables
!
    type(structure_turbines)                                                        , intent(inout) :: turbines
    !real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: xcor
    !real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: ycor
    !real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: guu
    !real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: gvv
    !integer      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: kcu
    !integer      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(in)    :: kcv
    logical                                                                         , intent(out)   :: error
!
! Local variables
!
    logical                                            :: inside
    integer                                            :: dn
    integer                                            :: dm
    integer                                            :: i
    integer                                            :: istat
    integer                                            :: j
    integer                                            :: n
    integer                                            :: nedges
    integer                                            :: nf
    integer                                            :: nt
    integer                                            :: m
    integer                                            :: mf
    integer                                            :: mt
    real(fp)                                           :: ang_nline
    real(fp)                                           :: ang_mline
    real(fp)                                           :: ang_ax
    real(fp)                                           :: df
    real(fp)                                           :: dt
    real(fp)                                           :: radius
    real(fp)                                           :: rn
    real(fp)                                           :: rm
    real(fp)                                           :: xu
    real(fp)                                           :: yu
    character(1)                                       :: normalflowdir
    character(10)                                      :: angn_str
    character(10)                                      :: angm_str
    character(10)                                      :: angx_str
    type(structure_turbine)                  , pointer :: turbine
    !
    integer, pointer                                   :: kmax
    integer, pointer                                   :: lundia
    logical, pointer                                   :: spheric

    integer                                            :: num, nlinks, ierror
    integer,          dimension(:),   allocatable      :: iLink, ipol
    double precision, dimension(:),   allocatable      :: xx, yy
    double precision, dimension(:),   allocatable      :: dSL
    
    character(len=40), dimension(1)                    :: namturb      !< names of points
    double precision,  dimension(1)                    :: xturb, yturb
    integer,           dimension(1)                    :: kturb
    double precision                                   :: xi1, xi2
    integer                                            :: jakdtree, jaoutside
    
    double precision                                   :: xh, yh, xe, ye, R, dL
    integer                                            :: iL, k1, k2, LL, L
    
    !double precision, external                         :: dbdistance, dprodin
    double precision, external                         :: dprodin

!
!! executable statements -------------------------------------------------------
!
    
!   get maximum mesh width
    dL = 0d0
    do L=1,numL
       k1 = kn(1,L)
       k2 = kn(2,L)
       dL = max(dL,dbdistance(xk(k1),yk(k1),xk(k2),yk(k2),jsferic, jasfer3D, dmiss))
    end do
    
    error = .false.
    if (associated(turbines%nr)) then
        !kmax    => gdp%d%kmax
        !lundia  => gdp%gdinout%lundia
        !spheric => gdp%gdtricom%sferic
        do j = 1, size(turbines%nr)
            turbine => turbines%nr(j)

            !! TODO locate velocity reference points using ndiamup


            !! calculate intersected flow links
            num = 2  ! just draw a line segment
            allocate(xx(num), yy(num))

            allocate(iLink(Lnx))
            iLink = 0
            allocate(ipol(Lnx))
            ipol = 0
            allocate(dSL(Lnx))
            dSL = 0d0
            
            if (turbine%turbtype == 1) then
               R = 0.5d0*turbine%width
            else
               R = 0.5d0*turbine%diam
            endif

            xx(1) = turbine%xyz(1) - sin(turbine%angle*degrad) * (R+dL)
            yy(1) = turbine%xyz(2) + cos(turbine%angle*degrad) * (R+dL)
            xx(2) = turbine%xyz(1) + sin(turbine%angle*degrad) * (R+dL)
            yy(2) = turbine%xyz(2) - cos(turbine%angle*degrad) * (R+dL)

            ! netlinks or flowlinks???   Flowlinks
            call find_crossed_links_kdtree2(treeglob,num,xx,yy,2,Lnx,0,nlinks, iLink, ipol, dSL, ierror)
            !subroutine find_crossed_links_kdtree2(treeinst,NPL,xpl,ypl,itype,nLinks,jaboundarylinks,numcrossedLinks, iLink, iPol, dSL, ierror)
            !   type(kdtree_instance),               intent(inout) :: treeinst
            !   integer,                             intent(in)    :: NPL                !< polyline length
            !   double precision, dimension(NPL),    intent(in)    :: xpl, ypl           !< polyline node coordinates
            !   integer,                             intent(in)    :: itype              !< netlinks (1) or flowlinks(2)
            !   integer,                             intent(in)    :: nLinks             !< number of links ( Lnx for flowlinks, numL for netlinks)
            !   integer,                             intent(in)    :: jaboundarylinks    !< include boundary links (1) or not (0), flowlinks only
            !   integer,                             intent(out)   :: numcrossedLinks    !< number of crossed flowlinks
            !   integer,          dimension(nLinks), intent(inout) :: iLink              !< crossed flowlinks
            !   integer,          dimension(nLinks), intent(inout) :: iPol               !< polygon section
            !   double precision, dimension(nLinks), intent(inout) :: dSL                !< polygon section cross location
            !   integer,                             intent(out)   :: ierror             !< ierror (1) or not (0)

            deallocate(xx, yy)
            write(*,*) 'nlinks:',nlinks
            write(*,*) 'ilinks:',ilink(1:nlinks)
            allocate(turbine%edgelist(nlinks)       , stat=istat)
!            if (istat==0) allocate(turbine%reldist  (nlinks+1)     , stat=istat)
!            if (istat==0) allocate(turbine%zlevel   (nlinks,0:kmax), stat=istat)
!            if (istat==0) allocate(turbine%area     (nlinks,kmax)  , stat=istat)
            if (istat==0) allocate(turbine%blockfrac(nlinks,kmx)  , stat=istat)
            if (istat/=0) then
               ! error: memory allocation
               !call write_error('Memory allocation error in MAPTURBINE',unit=lundia)
               error = .true.
               return
            endif
            
!           get projected coordinates
            if ( allocated(turbine%xi1) ) deallocate(turbine%xi1)
            if ( allocated(turbine%xi2) ) deallocate(turbine%xi2)
            allocate(turbine%xi1(nlinks), stat=istat)
            allocate(turbine%xi2(nlinks), stat=istat)
            
!           get hub coordinates
            xh = turbine%xyz(1)
            yh = turbine%xyz(2)
            
            num = 0

            do iL = 1,nlinks
               LL = turbine%edgelist(iL)
               LL = ilink(iL)
               
!              compute corner coordinates in rotor plane
               k1 = lncn(1,LL)
               k2 = lncn(2,LL)
               
               if (turbine%turbtype == 1) then
                  R =  0.5d0*turbine%width                  
               else
                  R =  0.5d0*turbine%diam
               endif
               xe = xh - R*sin(turbine%angle*dg2rd)
               ye = yh + R*cos(turbine%angle*dg2rd)
               
               xi1 = dprodin(xh,yh,xe,ye,xh,yh,xk(k1),yk(k1))/R
               xi2 = dprodin(xh,yh,xe,ye,xh,yh,xk(k2),yk(k2))/R
               
               if ( abs(xi1).le.R .or. abs(xi2).le.R ) then
                  num = num+1
               
                  turbine%xi1(num) =  xi1
                  turbine%xi2(num) =  xi2

            !     store the list of crossed flow links
                  turbine%edgelist(num) = LL
               end if
            end do
            
!           de/reallocate            
            deallocate(ilink,ipol,dSL)
!            call realloc(turbine%edgelist,num,keepExisting=.true.)
            call realloc(turbine%xi1,num,keepExisting=.true.)
            call realloc(turbine%xi2,num,keepExisting=.true.)
            
            turbine%numedges = num
            
            
            ! store center flownode
            namturb = ''
            xturb(1) = turbine%xyz(1)  ! note: hub can be in two cells neighboring the rotor, one of them is selected
            yturb(1) = turbine%xyz(2)
            jakdtree = 1
            jaoutside = 0
            call find_flownode(1, xturb, yturb, namturb, kturb, jakdtree, jaoutside, INDTP_ALL)
            turbine%cellnr = kturb(1)

            ! store center upstream reference node
            xturb(1) = turbine%xyz(1) + turbine%ndiamu*turbine%diam*turbine%csturb
            yturb(1) = turbine%xyz(2) + turbine%ndiamu*turbine%diam*turbine%snturb
            call find_flownode(1, xturb, yturb, namturb, kturb, jakdtree, jaoutside, INDTP_ALL)
            turbine%cellu(1) = kturb(1)

            ! store center downstream reference node
            xturb(1) = turbine%xyz(1) - turbine%ndiamu*turbine%diam*turbine%csturb
            yturb(1) = turbine%xyz(2) - turbine%ndiamu*turbine%diam*turbine%snturb
            call find_flownode(1, xturb, yturb, namturb, kturb, jakdtree, jaoutside, INDTP_ALL)
            turbine%cellu(2) = kturb(1)
            
        enddo
    endif
end subroutine mapturbine


!subroutine updturbine(turbines, dzu, dzv, dpu, dpv, hu, hv, zw, thick, &
!                    & u, v, alfas, hdt, nmaxddb)
subroutine updturbine(turbines)
!!--description-----------------------------------------------------------------
!
! Apply turbines
!
!!--declarations----------------------------------------------------------------
    use precision
    use mathconsts
    use m_structures, only: structure_turbines, structure_turbine, turbinecurve, intersect_turbine !, turbines
    !use m_kdtree2
    use kdtree2Factory
    use m_flowtimes, only: dts,time0
    use m_flowgeom, only: dxi
    use m_flow, only: adve, advi, u0, u1, Ltop, Lbot, hu, kmx, s1, zws, ucx, ucy
    use m_flowgeom, only: bob
    use unstruc_messages
    use m_sferic, only: dg2rd
    !
    implicit none
!
! Call variables
!
    type(structure_turbines)                                   , intent(inout) :: turbines
    !real(fp)                                                   , intent(in)    :: hdt
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: alfas
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub, gdp%d%kmax), intent(in)    :: dzu
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub, gdp%d%kmax), intent(in)    :: dzv
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: dpu
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: dpv
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: hu
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: hv
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)            , intent(in)    :: zw
    !real(fp)     , dimension(gdp%d%kmax)                       , intent(in)    :: thick
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub, gdp%d%kmax), intent(in)    :: u
    !real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub, gdp%d%kmax), intent(in)    :: v
    !integer                                                    , intent(in)    :: nmaxddb
!
! Local variables
!
    integer                                            :: j, il, l
    !integer                                            :: k
    !integer                                            :: n
    !integer                                            :: nm
    real(fp)                                           :: rhow
    !real(fp)                                           :: unm1
    !real(fp)                                           :: unm2
    !real(fp)                                           :: z0
    !character(256)                                     :: errorstring
    type(structure_turbine)                  , pointer :: turbine
    !!
    !integer, pointer                                   :: kmax
    !integer, pointer                                   :: lundia
    !logical, pointer                                   :: zmodel
    logical                                            :: error
    integer                                            :: numcrossedlinks
    double precision                                   :: Ct, Cd, aaa
    
!    double precision, dimension(2,4)                   :: pt   ! corner coordinates in rotor plane
    
    double precision                                   :: totalForce, uavg, dxavg, favg, pavg, ud   ! temp debug variables
    
    double precision                                   :: xi1, xi2
    double precision                                   :: xh, yh, zh, zb
    integer                                            :: k1, k2, LL, Lb, Lt, k, kk, kb, kt
    
    double precision, dimension(2)                     :: reldist
    double precision, dimension(2)                     :: zlevel
    double precision                                   :: blockfrac
    double precision                                   :: area
    double precision                                   :: uref, udx, udy
    character(256)                                     :: errorstring

    double precision, external                         :: dprodin
!
!! executable statements -------------------------------------------------------
!
    double precision :: arat

    rhow = 1000.0_fp;
    if (associated(turbines%nr)) then
        !kmax   => gdp%d%kmax
        !lundia => gdp%gdinout%lundia
        !zmodel => gdp%gdprocs%zmodel
        do j = 1, size(turbines%nr)
            turbine => turbines%nr(j)

            numcrossedlinks = turbine%numedges  ! size(turbine%edgelist)
            
!           get hub coordinates
            xh = turbine%xyz(1)
            yh = turbine%xyz(2)
            if ( turbine%vertpos.eq.1 ) then
               zh = s1(turbine%cellnr) - turbine%xyz(3)
            else
               zh = turbine%xyz(3)
            end if

          if (turbine%cellnr>0) then

            !
            ! Local Velocity Thrust calculation
            !
            ! F = 0.5*A*rho*Ct*Uinf^2
            ! Ud = Uinf*(1-a)
            ! F = 0.5*A*rho*Ud^2/(1-a)^2
            !   Ct = 4.*a*(1.-a)  => a = (1.0 - sqrt(1.0-Ct))/2.0
            !   Ct/(1-a)^2 = 4.*a/(1.-a)
            !   Cd = 4.0*Ct/(1-a)^2 = 4.0*(1.0 - sqrt(1.0-Ct))/(1.0 + sqrt(1.0-Ct))
            ! F = 0.5*A*rho*Cd*Ud^2
            
            
            ! disk velocity
            kk = turbine%cellnr                 ! turbine location
            call getkbotktop(kk,kb,kt)          ! get 3D cell bed- and surface layer cell numbers
            do k=kb,kt                          ! loop over the layers
                if ( zws(k-1)<= zh .and. zh < zws(k) ) then                     ! use zws, the cell/layer-interface vertical coordinate, layer k is between zws(k-1) and zws(k)
                    udx = ucx(k)
                    udy = ucy(k)
                    exit
                end if
            end do
            
            if ( udx*turbine%csturb + udy*turbine%snturb < 0.0 ) then   ! dot disk velocity with turbine orientation to get upwind side
                kk = turbine%cellu(1)
            else
                kk = turbine%cellu(2)
            endif
            call getkbotktop(kk,kb,kt)          ! get 3D cell bed- and surface layer cell numbers
            do k=kb,kt                          ! loop over the layers
                if ( zws(k-1)<= zh .and. zh < zws(k) ) then                     ! use zws, the cell/layer-interface vertical coordinate, layer k is between zws(k-1) and zws(k)
                    uref = abs(ucx(k)*turbine%csturb + ucy(k)*turbine%snturb)   ! use projected reconstructed velocity vector (ucx,ucy) in turbine normal direction
                    exit
                end if
            end do

            turbine%thrustcoef = turbinecurve(turbines%curves, turbine%thrustcrvnr, uref, errorstring)
            
            if (turbine%powercrvnr(1)>0 .and. errorstring == '  ') then
                turbine%powercoef      = turbinecurve(turbines%curves, turbine%powercrvnr, uref, errorstring)
            else
                turbine%powercoef      = 0.0_fp
            endif
            
            !TODO put in stop conditions for ill defined coefficients, Ct > 1.0, Ct < 0.0
            !if (errorstring /= ' ') then
            !    write(lundia,'(A)') trim(errorstring)
            !    call d3stop(1, gdp)
            !endif

            Ct = turbine%thrustcoef
            aaa = (1.0 - sqrt(1.0-Ct))/2.0
            Cd = 4.0*(1.0 - sqrt(1.0-Ct))/(1.0 + sqrt(1.0-Ct))
            
            turbine%current_power = 0.0
            do il = 1,numcrossedlinks
                LL = turbine%edgelist(il)

                if ( kmx.eq.0 ) then   ! 2D
                  advi(LL) = advi(LL) + 0.5d0*Cd*dxi(LL)*abs(u1(LL))
                else
                  ! get horizontal corner coordinates in rotor plane
                  reldist(1) = turbine%xi1(iL)
                  reldist(2) = turbine%xi2(iL)
                  
                  ! get bedlevel
                  zb = 0.5d0*(bob(1,LL)+bob(2,LL))
                  call getlbotltop(LL,Lb,Lt)

                  do L=Lb,Lt
                     k = L - Lb + 1
                     ! get vertical corner coordinates
                     zlevel(1) = zb + hu(L-1)
                     zlevel(2) = zb + hu(L)
                    
                     call intersect_turbine(reldist,zlevel,zh,turbine%diam,turbine%width,turbine%height,turbine%turbtype,area,blockfrac)
                     turbine%blockfrac(il,k) = blockfrac

                     if (turbine%turbinemodel == 1) then
                        advi(L) = advi(L) + 0.5d0*Cd*dxi(LL)*abs(u0(L))*blockfrac
                     else
                        adve(L) = adve(L) + 0.5d0*Ct*dxi(LL)*uref**2*blockfrac
                     endif

                     turbine%current_power  = turbine%current_power + 0.5_fp * turbine%powercoef * rhow * area * abs(uref**3)*blockfrac
                  end do
                end if
            enddo
            turbine%cumul_power = turbine%cumul_power + dts*turbine%current_power
            !! CCC DEBUG
            !! TODO fix this ugly hack
            write(5979,'(a,f,a,i,a,e,a,e)') 'Time: ',time0,'  turbine # ',j,'  Power ',turbine%current_power, &
                                        '  Total Power  ',turbine%cumul_power
          endif
        enddo

        !call add_loss_due_to_turbines(turbines, v, u, 1, nmaxddb, 4, gdp)
        !call add_loss_due_to_turbines(turbines, u, v, nmaxddb, 1, 4, gdp)
    endif
end subroutine updturbine

function wrturbine_cnst(turbines, fds, grpnam, grpind) result (ierror)
!!--description-----------------------------------------------------------------
!
! Write NEFIS elements of constants group
!
!!--declarations----------------------------------------------------------------
    use precision
    use m_structures, only: structure_turbines
    !
    implicit none
!
! Call variables
!
    type(structure_turbines)                                   , intent(in)    :: turbines
    integer                                                    , intent(in)    :: fds
    character(*)                                               , intent(in)    :: grpnam
    integer, dimension(3,5)                                    , intent(in)    :: grpind
    integer                                                                    :: ierror
!
! Local variables
!
    integer                                           :: i
    integer                                           :: nturb
    integer                            , external     :: putelt
    real(sp), dimension(:), allocatable               :: sbuff1
    real(sp), dimension(:,:), allocatable             :: sbuff2
    character(256), dimension(:), allocatable         :: cbuff
!
!! executable statements -------------------------------------------------------
!
    nturb = size(turbines%nr)
    if (nturb==0) return
    !
    allocate(cbuff(nturb),sbuff1(nturb),sbuff2(nturb,2))
    !
    !TODO FIXME update write routines to use NETCDF not NEFIS
    !do i = 1,nturb
    !    cbuff(i)    = turbines%nr(i)%name
    !    sbuff1(1)   = real(turbines%nr(i)%angle,sp)
    !    sbuff2(i,1) = real(turbines%nr(i)%xyz(1),sp)
    !    sbuff2(i,2) = real(turbines%nr(i)%xyz(2),sp)
    !enddo
    !!
    !ierror = putelt(fds, grpnam, 'NTURBINES', grpind, 1, nturb)
    !ierror = putelt(fds, grpnam, 'NAMTURBINES', grpind, 1, cbuff)
    !ierror = putelt(fds, grpnam, 'ANGTURBINES', grpind, 1, sbuff1)
    !ierror = putelt(fds, grpnam, 'XYTURBINES', grpind, 1, sbuff2)
    !
    deallocate(cbuff, sbuff1, sbuff2)
    
    ierror = 0
end function wrturbine_cnst


function wrturbine_time(turbines, fds, grpnam, grpind) result (ierror)
!!--description-----------------------------------------------------------------
!
! Write NEFIS elements of time-varying group
!
!!--declarations----------------------------------------------------------------
    use precision
    use m_structures, only: structure_turbines
    !
    implicit none
!
! Call variables
!
    type(structure_turbines)                                   , intent(in)    :: turbines
    integer                                                    , intent(in)    :: fds
    character(*)                                               , intent(in)    :: grpnam
    integer, dimension(3,5)                                    , intent(in)    :: grpind
    integer                                                                    :: ierror
!
! Local variables
!
    integer                                           :: i
    integer                                           :: nturb
    integer                            , external     :: putelt
    real(sp), dimension(:), allocatable               :: sbuff1
!
!! executable statements -------------------------------------------------------
!
    nturb = size(turbines%nr)
    ierror = 0
    if (nturb==0) return
    !
    allocate(sbuff1(nturb))
    !
    do i = 1,nturb
        sbuff1(i) = real(turbines%nr(i)%current_zlevel,sp)
    enddo    

    !TODO FIXME update write routines to use NETCDF not NEFIS
    !ierror = putelt(fds, grpnam, 'ZTURBINES', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%current_uref,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'UTURBINES', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%powercoef,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'POWERCOEF', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%thrustcoef,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'THRUSTCOEF', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%friccoef,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'LOSSCOEF', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%current_sim_thrust,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'INST_SIMTHRUST', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%cumul_sim_thrust,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'CUM_SIMTHRUST', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%current_thrust,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'INST_THRUST', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%cumul_thrust,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'CUM_THRUST', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%current_power,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'INST_POWER', grpind, 1, sbuff1)
    !!
    !do i = 1,nturb
    !    sbuff1(i) = real(turbines%nr(i)%cumul_power,sp)
    !enddo    
    !ierror = putelt(fds, grpnam, 'CUM_POWER', grpind, 1, sbuff1)
    !
    deallocate(sbuff1)
end function wrturbine_time

end module m_rdturbine
