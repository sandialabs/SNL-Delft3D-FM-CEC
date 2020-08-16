module m_structures

!----- AGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2017-2020.!
!  This file is part of Delft3D (D-Flow Flexible Mesh component).
!
!  Delft3D is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Affero General Public License as
!  published by the Free Software Foundation version 3.
!
!  Delft3D  is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Affero General Public License for more details.
!
!  You should have received a copy of the GNU Affero General Public License
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.
!
!  contact: delft3d.support@deltares.nl
!  Stichting Deltares
!  P.O. Box 177
!  2600 MH Delft, The Netherlands
!
!  All indications and logos of, and references to, "Delft3D",
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!
!-------------------------------------------------------------------------------
   
! $Id: unstruc_structures.f90 65986 2020-02-16 13:04:44Z jagers $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/SANDIA/fm_tidal_v3/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/unstruc_structures.f90 $

use properties
use m_GlobalParameters
use unstruc_channel_flow, only: network
use precision
use table_handles, only: handletype, MAXTABLECLENGTH, cleartable
implicit none

type(tree_data), pointer, public :: strs_ptr !< A property list with all input structure specifications of the current model. Not the actual structure set.
integer :: jaoldstr !< tmp backwards comp: we cannot mix structures from EXT and from structure-input files. Use one or the other.

 ! Structure Parameters
 double precision, dimension(:,:), allocatable :: valpump     !< Array for pump;      (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) pump discharge w.r.t. structure orientation
                                                              !<                      (3,:) pump structure water level up
                                                              !<                      (4,:) pump structure water level down
                                                              !<                      (5,:) pump structure head
                                                              !<                      (6,:) pump capacity
                                                              !<                      (7,:) actual pump stage
                                                              !<                      (8,:) pump head
                                                              !<                      (9,:) pump reduction factor
                                                              !<                      (10,:) pump water level at delivery side
                                                              !<                      (11,:) pump water level at suction side
                                                              !<                      (12,:) pump discharge w.r.t. pumping orientation (same sign as capacity)

 double precision, dimension(:,:), allocatable :: valgate     !< Array for gate;      (1,:) discharge through gate
 double precision, dimension(:,:), allocatable :: valcdam     !< Array for cdam;      (1,:) discharge through controlable dam
                                                              !<                      (2,:) Upstream average water levels
                                                              !<                      (3,:) downstream average water level
                                                              !<                      (4,0) width of dam
 double precision, dimension(:,:), allocatable :: valgategen  !< Array for gate(new), (1,:) discharge through gate
                                                              !<                      (2,:) Upstream average water level
                                                              !<                      (3,:) gate width
 double precision, dimension(:,:), allocatable :: valweirgen  !< Array for weir;      (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through weir
                                                              !<                      (3,:) weir structure water level up
                                                              !<                      (4,:) weir structure water level down
                                                              !<                      (5,:) weir structure head
                                                              !<                      (6,:) weir flow area
                                                              !<                      (7,:) weir velocity
                                                              !<                      (8,:) water level on crest
                                                              !<                      (9,:) weir crest level
                                                              !<                      (10,:) weir crest width
                                                              !<                      (11,:) weir state (0: closed, 1: free weir, 2: drowned/submerged weir)
                                                              !<                      (12,:) weir force difference per unit width
                                                              !<                      (13,:) weir counters of partitions for parallel
 double precision, dimension(:,:), allocatable :: valcgen     !< Array for general structure (old ext), (1,:) discharge
 double precision, dimension(:,:), allocatable :: valgenstru  !< Array for general structure (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through general structure
                                                              !<                      (3,:) general structure water level up
                                                              !<                      (4,:) general structure water level down
                                                              !<                      (5,:) general structure head
                                                              !<                      (6,:) general structure flow area
                                                              !<                      (7,:) general structure velocity
                                                              !<                      (8,:) general structure water level on crest
                                                              !<                      (9,:) general structure crest level
                                                              !<                      (10,:) general structure crest width
                                                              !<                      (11,:) general structure state (0: closed, 1: free weir, 2: drowned/submerged weir)
                                                              !<                      (12,:) general structure force difference per unit width
                                                              !<                      (13,:) general structure gate opening width
                                                              !<                      (14,:) general structure gate lower edge level
                                                              !<                      (15,:) general structure gate opening height
                                                              !<                      (16,:) general structure gate upper edge level
                                                              !<                      (17,:) general structure discharge through gate opening
                                                              !<                      (18,:) general structure discharge over gate
                                                              !<                      (19,:) general structure discharge under gate
                                                              !<                      (20,:) general structure flow area in gate opening
                                                              !<                      (21,:) general structure flow area over gate
                                                              !<                      (22,:) general structure flow area under gate
                                                              !<                      (23,:) general structure velocity through gate opening
                                                              !<                      (24,:) general structure velocity over gate
                                                              !<                      (25,:) general structure velocity under gate
                                                              !<                      (26,:) general structure counters of partitions for parallel
 double precision, dimension(:,:), allocatable, target :: valdambreak !< Array for dambreak, (1,:)  flow link width
                                                              !<                      (2,:) instantanuous discharge
                                                              !<                      (3,:) dambreak water level up
                                                              !<                      (4,:) dambreak water level down
                                                              !<                      (5,:) dambreak structure head
                                                              !<                      (6,:) dambreak flow area
                                                              !<                      (7,:) dambreak normal velocity
                                                              !<                      (8,:) dambreak crest level
                                                              !<                      (9,:) dambreak crest width
                                                              !<                      (10,:) dambreak water level jump
                                                              !<                      (11,:) dambreak breach width time derivative
                                                              !<                      (12,:) cumulative discharge
 double precision, dimension(:,:), allocatable :: valorifgen  !< Array for orifice (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through orifice
                                                              !<                      (3,:) orifice water level up
                                                              !<                      (4,:) orifice water level down
                                                              !<                      (5,:) orifice head
                                                              !<                      (6,:) orifice flow area
                                                              !<                      (7,:) orifice velocity
                                                              !<                      (8,:) orifice water level on crest
                                                              !<                      (9,:) orifice crest level
                                                              !<                      (10,:) orifice crest width
                                                              !<                      (11,:) orifice state (0: closed, 1: free weir, 2: drowned/submerged weir)
                                                              !<                      (12,:) orifice force difference per unit width
                                                              !<                      (13,:) orifice gate opening width (not applicable)
                                                              !<                      (14,:) orifice gate lower edge level
                                                              !<                      (15,:) orifice gate opening height
                                                              !<                      (16,:) orifice counters of partitions for parallel
 double precision, dimension(:,:), allocatable :: valbridge   !< Array for bridge;    (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through bridge
                                                              !<                      (3,:) bridge water level up
                                                              !<                      (4,:) bridge water level down
                                                              !<                      (5,:) bridge head
                                                              !<                      (6,:) bridge flow area
                                                              !<                      (7,:) bridge velocity
                                                              !<                      (8,:) bridge bed level up
                                                              !<                      (9,:) bridge bed level down
                                                              !<                      (10,:) bridge actual bed level (crest)
 double precision, dimension(:,:), allocatable :: valculvert  !< Array for culvert;   (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through culvert
                                                              !<                      (3,:) culvert water level up
                                                              !<                      (4,:) culvert water level down
                                                              !<                      (5,:) culvert structure head
                                                              !<                      (6,:) culvert flow area
                                                              !<                      (7,:) culvert velocity
                                                              !<                      (8,:) culvert crest level
                                                              !<                      (9,:) culvert state (0: closed, 1: free weir, 2: drowned/submerged weir)
                                                              !<                      (10,:) culvert gate lower edge level
                                                              !<                      (11,:) culvert gate opening height
 double precision, dimension(:,:), allocatable :: valuniweir  !< Array for universal weir; (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through universal weir
                                                              !<                      (3,:) universal weir water level up
                                                              !<                      (4,:) universal weir water level down
                                                              !<                      (5,:) universal weir head
                                                              !<                      (6,:) universal weir flow area
                                                              !<                      (7,:) universal weir velocity
                                                              !<                      (8,:) universal weir crest level
 double precision, dimension(:,:), allocatable :: valcmpstru  !< Array for compound structure; (1,:) flow link width, used for averaging.
                                                              !<                      (2,:) discharge through compound structure
                                                              !<                      (3,:) compound structure water level up
                                                              !<                      (4,:) compound structure water level down
                                                              !<                      (5,:) compound structure head
                                                              !<                      (6,:) compound structure flow area
                                                              !<                      (7,:) compound structure velocity
 integer                           :: NUMVALS_PUMP = 12       !< Number of variables for pump
 integer                           :: NUMVALS_GATE = 5        !< Number of variables for gate
 integer                           :: NUMVALS_CDAM = 4        !< Number of variables for controble dam
 integer                           :: NUMVALS_CGEN = 4        !< Number of variables for general structure (old ext file)
 integer                           :: NUMVALS_GATEGEN = 9     !< Number of variables for gate (new)
 integer                           :: NUMVALS_WEIRGEN = 13    !< Number of variables for weir
 integer                           :: NUMVALS_GENSTRU = 26    !< Number of variables for general structure( new exe file)
 integer                           :: NUMVALS_DAMBREAK = 12   !< Number of variables for dambreak
 integer                           :: NUMVALS_ORIFGEN = 16    !< Number of variables for orific
 integer                           :: NUMVALS_BRIDGE  = 10    !< Number of variables for bridge
 integer                           :: NUMVALS_CULVERT = 11    !< Number of variables for culvert
 integer                           :: NUMVALS_UNIWEIR = 8     !< Number of variables for univeral weir
 integer                           :: NUMVALS_CMPSTRU = 7     !< Number of variables for compound structure
 
 integer                           :: jahiscgen               !< Write structure parameters to his file, 0: n0, 1: yes
 integer                           :: jahispump               !< Write pump      parameters to his file, 0: n0, 1: yes
 integer                           :: jahisgate               !< Write gate      parameters to his file, 0: n0, 1: yes
 integer                           :: jahiscdam               !< Write dam       parameters to his file, 0: n0, 1: yes
 integer                           :: jahisweir               !< Write weir      parameters to his file, 0: n0, 1: yes
 integer                           :: jahisdambreak           !< Write dambreak  parameters to his file, 0: n0, 1: yes
 integer                           :: jahisorif               !< Write orifice   parameters to his file, 0: no, 1: yes
 integer                           :: jahisbridge             !< Write bridge    parameters to his file, 0: no, 1: yes
 integer                           :: jahisculv               !< Write culvert   parameters to his file, 0: no, 1: yes
 integer                           :: jahisuniweir            !< Write univeral weir parameters to his file, 0: no, 1: yes
 integer                           :: jahiscmpstru            !< Write compound structure parameters to his file, 0: no, 1: yes
 
 integer, parameter :: IOPENDIR_FROMLEFT  = -1 !< Gate door opens/closes from left side.
 integer, parameter :: IOPENDIR_FROMRIGHT =  1 !< Gate door opens/closes from right side.
 integer, parameter :: IOPENDIR_SYMMETRIC =  0 !< Gate door opens/closes symmetrically (from center).

 type tgate                                          !< Gate structure type, before it gets evaluated as a general structure.
    !double precision :: sill_level       !< Not used: stored in zcgen(1,igen)
    !double precision :: lower_edge_level !< Not used: stored in zcgen(2,igen)
    !double precision :: opening_width    !< Not used: stored in zcgen(3,igen)
    double precision :: door_height       !< Height of the door, used for 'double-barrier' overflow. Time-INDEPENDENT.
    double precision :: sill_width        !< Width of the sill, may be larger than the opening width, such that in open part we have weir flow and in closed part we have gate flow. Time-INDEPENDENT.
    integer          :: opening_direction !< Direction from which the gate opens/closes, IOPENDIR_FROMLEFT|FROMRIGHT|SYMMETRIC.
 end type tgate

 type(tgate), allocatable :: gates(:)

 ! TIDAL TURBINES: Insert allocatable of type structure_turbines here

 type structure_turbine                                  !< turbine structure type containing data pertaining to tidal current turbines
     integer, dimension(2)             :: thrustcrvnr    !< index of turbine thurst curve
     integer, dimension(2)             :: powercrvnr     !< index of turbine power curve
     integer                           :: vertpos        !< vertical positioning type of turbine: 0 = fixed z level, 1 = turbine axis at fixed depth relative to water surface
     integer                           :: turbtype       !< BJ 20150326: turbine type: 0 = circular, 1 = rectangular
     real(fp)                          :: angle          !< orientation of turbine indicated by angle between turbine axis and positive x/lon direction
     real(fp)                          :: csturb         !< cosine of orientation angle
     real(fp)                          :: snturb         !< sine of orientation angle
     real(fp)                          :: diam           !< diameter of turbine (turbtype=0)
     real(fp)                          :: width          !< BJ 20150326: width of turbine (turbtype=1)
     real(fp)                          :: height         !< BJ 20150326: height of turbine (turbtype=1)
     real(fp)                          :: turbarea       !< turbine area (0.25*pi*D**2)
     real(fp), dimension(3)            :: xyz            !< (relative) position of turbine
     character(256)                    :: name           !< name of turbine
     character(MAXTABLECLENGTH)        :: thrustcrvname  !< name of turbine thrust curve
     character(MAXTABLECLENGTH)        :: powercrvname   !< name of turbine power curve
     !
     integer                           :: turbinemodel   !< turbine model used 0: default deltares model, 1: SNL model
     integer                           :: turbulencemodel!< turbulence model used at turbine 0: none 1: Rethore implementation 2: Rados implementation
     real(fp)                          :: beta_p         !< Beta_p coefficient
     real(fp)                          :: beta_d         !< Beta_d coefficient
     real(fp)                          :: cep4           !< Cepsilon4 coefficient
     real(fp)                          :: cep5           !< Cepsilon5 coefficient
     real(fp)                          :: friccoef       !< latest friction coefficient
     real(fp)                          :: powercoef      !< latest power coefficient
     real(fp)                          :: thrustcoef     !< latest thrust coefficient
     !
     real(fp)                          :: current_uref       !< latest reference velocity
     real(fp)                          :: current_zlevel     !< latest vertical position
     real(fp)                          :: current_sim_thrust !< latest simulated thrust
     real(fp)                          :: cumul_sim_thrust   !< integrated simulated thrust
     real(fp)                          :: current_thrust     !< latest analytical thrust
     real(fp)                          :: cumul_thrust       !< integrated analytical thrust
     real(fp)                          :: current_power      !< latest analytical power
     real(fp)                          :: cumul_power        !< integrated analytical power
     !
     integer                           :: cellnr         !< index of cell in which turbine is located
     real(fp), dimension(2)            :: relhpos        !< relative location within cell
     !
     real(fp)                          :: ndiamu         !< distance of velocity measurement point "up- & downstream"
     integer, dimension(2)             :: cellu          !< indices of "up- & downstream" cells for velocity
     !
     character(1)                      :: edgetype       !< type of edge/link on which structure is located
     integer                           :: numedges       !< number of flowlinks in rotorplane
     integer , dimension(:)  , pointer :: edgelist       !< indices of edges on which structure is located
     real(fp), dimension(:)  , pointer :: reldist        !< distance along edges measured relative to location of turbine
     real(fp), dimension(:,:), pointer :: zlevel         !< work array of right size to store levels of layer interfaces per cell
     real(fp), dimension(:,:), pointer :: area           !< work array to store full    exchange area          per layer and cell
     real(fp), dimension(:,:), pointer :: blockfrac      !< work array to store blocked exchange area fraction per layer and cell
     
     double precision, dimension(:), allocatable    :: xi1, xi2
 end type structure_turbine

!< type for storing array of turbine structure types and performance curves
 type structure_turbines
     type(structure_turbine)   , dimension(:)    , pointer :: nr         !< pointer to individual turbines
     type(handletype)                                      :: curves     !< thrust and power curve tables
 end type structure_turbines

 type (structure_turbines) :: turbines

   contains

!> @brief Initialize turbine data records
!!
!! @param[out]      turbines    a turbine data structure
subroutine init_turbines(turbines)

    implicit none
!
! Call variables
!
    type(structure_turbines)     , intent(out) :: turbines

!
! executable statements
!
    nullify(turbines%nr)

end subroutine init_turbines

!> @brief Allocate turbine data records
!!
!! @param[in]      lundia       diagnostic file unit number
!! @param[in]      nturbines    number of turbines required
!! @param[out]     %error       returns .false. if memory could not be allocated
!! @param[inout]   turbines     structure_turbines type
subroutine allocate_turbines(turbines,nturbines,lundia,error)
!
! Declarations
!
    use message_module, only: write_error
    implicit none
!
! Call variables
!
    integer                     , intent(in)    :: lundia
    integer                     , intent(in)    :: nturbines
    logical                     , intent(out)   :: error
    type(structure_turbines)    , intent(inout) :: turbines
!
! Local variables
!
    integer :: i
    integer :: istat
!
! Executable statements
!
    error = .false.
    allocate(turbines%nr  (nturbines)  , stat=istat)
    if (istat/=0) then
       call write_error('Memory allocation error in ALLOCATE_TURBINES.',unit=lundia)
       error = .true.
       return
    endif
    !
    do i = 1, nturbines
       turbines%nr(i)%thrustcrvnr        = -999
       turbines%nr(i)%powercrvnr         = -999
       turbines%nr(i)%vertpos            = 0
       turbines%nr(i)%turbtype           = 0       ! BJ 20150326
       turbines%nr(i)%angle              = 0.0_fp
       turbines%nr(i)%csturb             = 0.0_fp
       turbines%nr(i)%snturb             = 0.0_fp
       turbines%nr(i)%width              = 0.0_fp  ! BJ 20150326
       turbines%nr(i)%height             = 0.0_fp  ! BJ 20150326
       turbines%nr(i)%diam               = 0.0_fp
       turbines%nr(i)%turbarea           = 0.0_fp
       turbines%nr(i)%xyz                = 0.0_fp
       turbines%nr(i)%name               = ' '
       turbines%nr(i)%thrustcrvname      = ' '
       turbines%nr(i)%powercrvname       = ' '
       !
       turbines%nr(i)%cellnr             = -999
       turbines%nr(i)%relhpos            = 0.0_fp
       !
       turbines%nr(i)%ndiamu             = 5.0_fp
       turbines%nr(i)%cellu              = 0
       !
       turbines%nr(i)%turbinemodel       = 0
       turbines%nr(i)%turbulencemodel    = 0
       turbines%nr(i)%beta_p             = -999.0_fp
       turbines%nr(i)%beta_d             = -999.0_fp
       turbines%nr(i)%cep4               = -999.0_fp
       turbines%nr(i)%cep5               = -999.0_fp
       turbines%nr(i)%friccoef           = -999.0_fp
       turbines%nr(i)%powercoef          = -999.0_fp
       turbines%nr(i)%thrustcoef         = -999.0_fp
       !
       turbines%nr(i)%current_uref       = -999.0_fp
       turbines%nr(i)%current_sim_thrust = -999.0_fp
       turbines%nr(i)%cumul_sim_thrust   =    0.0_fp
       turbines%nr(i)%current_thrust     = -999.0_fp
       turbines%nr(i)%cumul_thrust       =    0.0_fp
       turbines%nr(i)%current_power      = -999.0_fp
       turbines%nr(i)%cumul_power        =    0.0_fp
       !
       turbines%nr(i)%edgetype           = ' '
       nullify(turbines%nr(i)%edgelist)
       nullify(turbines%nr(i)%reldist)
       nullify(turbines%nr(i)%zlevel)
       nullify(turbines%nr(i)%area)
       nullify(turbines%nr(i)%blockfrac)
    enddo
end subroutine allocate_turbines

!> @brief Deallocate turbine data records
!!
!! @param[inout]   turbines     structure_turbines type
subroutine deallocate_turbines(turbines)

!
! declarations
!
    implicit none
!
! Call variables
!
    type(structure_turbines)    , intent(inout) :: turbines
    
!
! Local variables
!
    integer :: i
    integer :: istat
    
!
! executable statements
!
    if (associated(turbines%nr)) then
       do i = 1,size(turbines%nr)
          if (associated(turbines%nr(i)%edgelist))  deallocate(turbines%nr(i)%edgelist, stat=istat)
          if (associated(turbines%nr(i)%reldist))   deallocate(turbines%nr(i)%reldist, stat=istat)
          if (associated(turbines%nr(i)%zlevel))    deallocate(turbines%nr(i)%zlevel, stat=istat)
          if (associated(turbines%nr(i)%area))      deallocate(turbines%nr(i)%area, stat=istat)
          if (associated(turbines%nr(i)%blockfrac)) deallocate(turbines%nr(i)%blockfrac, stat=istat)
          if (allocated(turbines%nr(i)%xi1)) deallocate(turbines%nr(i)%xi1)
          if (allocated(turbines%nr(i)%xi2)) deallocate(turbines%nr(i)%xi2)
       enddo
       deallocate(turbines%nr  , stat=istat)
    endif
    !
    call cleartable(turbines%curves)
end subroutine deallocate_turbines

!> @brief Calculate proportion of exchange area that is blocked by the turbine
subroutine intersect_turbine(reldist,zlevel,z0,diam,width,height,turbtype,area,blockfrac)
    use precision
    implicit none
!
! Call variables
!
    real(fp), dimension(2)   , intent(in)  :: reldist   !< horizontal distance relative to turbine
    real(fp), dimension(2)   , intent(in)  :: zlevel    !< vertical position of layer interfaces
    real(fp)                 , intent(in)  :: z0        !< z coordinate of turbine axis
    real(fp)                 , intent(in)  :: diam      !< turbine diameter
    real(fp)                 , intent(in)  :: width     !< turbine width
    real(fp)                 , intent(in)  :: height    !< turbine height
    integer                  , intent(in)  :: turbtype  !< turbine type
    real(fp)                 , intent(out) :: area      !< full    exchange area
    real(fp)                 , intent(out) :: blockfrac !< blocked exchange area fraction
!
! Local variables
!
    real(fp) :: rad
!
! executable statements
!
    rad = diam/2.0_fp
    !
    area     = abs(zlevel(1) - zlevel(2)) * abs(reldist(1) - reldist(2))
    if (turbtype == 0) then   ! BJ 20150326, existing case of circular turbine
       blockfrac = intersect_circle_rect(reldist(1),reldist(2),zlevel(1),zlevel(2),rad,0.0_fp,z0) / area
    elseif (turbtype == 1) then   ! BJ 20150326, new case of rectangular turbine
       blockfrac = abs(min(max(reldist(1),-width/2.0_fp),width/2.0_fp) - min(max(reldist(2),-width/2.0_fp),width/2.0_fp)) * &
                 & abs(min(max(zlevel(2)-z0,-height/2.0_fp),height/2.0_fp) - min(max(zlevel(1)-z0,-height/2.0_fp),height/2.0_fp)) / &
                 & area
    endif   ! BJ 20150326
end subroutine intersect_turbine

!> @brief Calculate area of the intersection of rectangle with circle
function intersect_circle_rect(x1,x2,z1,z2,rad,x0,z0) result (area)
    use precision
    implicit none
!
! Call variables
!
    real(fp), intent(in)  :: x1    !< x coordinate of corner 1
    real(fp), intent(in)  :: x2    !< x coordinate of corner 2
    real(fp), intent(in)  :: z1    !< z coordinate of corner 1
    real(fp), intent(in)  :: z2    !< z coordinate of corner 2
    real(fp), intent(in)  :: rad   !< radius of circle
    real(fp), intent(in)  :: x0    !< x coordinate of circle centre
    real(fp), intent(in)  :: z0    !< z coordinate of circle centre
    real(fp)              :: area  !< intersection area
!
! Local variables
!
    real(fp)              :: zero = 0.0_fp ! as the name says: zero
    !
    real(fp)              :: dxmin ! minimum of x coordinates
    real(fp)              :: dxmax ! maximum of x coordinates
    real(fp)              :: dzmin ! minimum of z coordinates
    real(fp)              :: dzmax ! maximum of z coordinates
    real(fp)              :: area1 ! area of intersection with circle part 1 (  0 -  90 deg)
    real(fp)              :: area2 ! area of intersection with circle part 1 ( 90 - 180 deg)
    real(fp)              :: area3 ! area of intersection with circle part 1 (180 - 270 deg)
    real(fp)              :: area4 ! area of intersection with circle part 1 (270 - 360 deg)
!
! executable statements
!
    !
    ! convert to coordinates relative to centre of circle
    !
    dxmin = min(x1,x2) - x0
    dxmax = max(x1,x2) - x0
    dzmin = min(z1,z2) - z0
    dzmax = max(z1,z2) - z0
    !
    ! The total area of the intersection of the rectangle with the circle
    ! equals the sum of the areas of the intersections of each quarter of the
    ! circle with the rectangle.
    !
    ! In each case we only have to consider part of the rectangle in the same
    ! quadrant (i.e. the associated sign of x and z). We implement only the
    ! first quarter and map the other quarters to the first one.
    !
    area1 = intersect_circle_rect4( max(zero,dxmin), max(zero,dxmax), max(zero,dzmin), max(zero,dzmax),rad) ! upper right (Q1)
    area2 = intersect_circle_rect4(-min(zero,dxmax),-min(zero,dxmin), max(zero,dzmin), max(zero,dzmax),rad) ! upper left  (Q2) - flip x
    area3 = intersect_circle_rect4(-min(zero,dxmax),-min(zero,dxmin),-min(zero,dzmax),-min(zero,dzmin),rad) ! lower left  (Q3) - flip x - flip z
    area4 = intersect_circle_rect4( max(zero,dxmin), max(zero,dxmax),-min(zero,dzmax),-min(zero,dzmin),rad) ! lower right (Q4) - flip z
    area  = area1 + area2 + area3 + area4
end function intersect_circle_rect

!> @brief Calculate area of the intersection of rectangle with circle quadrant
function intersect_circle_rect4(dxmin,dxmax,dzmin,dzmax,rad) result (area)
    use precision
    implicit none
!
! Call variables
!
    real(fp), intent(in)  :: dxmin !< minimum relative x coordinate (x=0 at circle centre)
    real(fp), intent(in)  :: dxmax !< maximum relative x coordinate (x=0 at circle centre)
    real(fp), intent(in)  :: dzmin !< minimum relative z coordinate (z=0 at circle centre)
    real(fp), intent(in)  :: dzmax !< maximum relative z coordinate (z=0 at circle centre)
    real(fp), intent(in)  :: rad   !< radius of circle
    real(fp)              :: area  !< intersection area
!
! Local variables
!
    real(fp)              :: zero = 0.0_fp ! as the name says: zero
    !
    real(fp)              :: rad2  ! squared radius
    real(fp)              :: a1    ! angle 1 (radian) associated with first intersection of circle and box
    real(fp)              :: a2    ! angle 2 (radian) associated with second intersection of circle and box
    real(fp)              :: dx1   ! x coordinate of intersection point 1
    real(fp)              :: dx2   ! x coordinate of intersection point 2
    real(fp)              :: dz1   ! z coordinate of intersection point 1
    real(fp)              :: dz2   ! z coordinate of intersection point 2
!
! executable statements
!
    rad2 = rad*rad
    !
    if ((dxmin*dxmin + dzmin*dzmin) > rad2) then
        ! rectangle completely outside circle
        area = zero
    elseif ((dxmax*dxmax + dzmax*dzmax) < rad2) then
        ! rectangle completely inside circle
        area = (dxmax-dxmin)*(dzmax-dzmin)
    else
        ! rectangle partly inside circle
        if ((dxmax*dxmax + dzmin*dzmin) > rad2) then
            ! edge (dxmin,dzmin) - (dxmax,dzmin) intersects the circle
            a1 = asin(dzmin/rad)
            dx1 = sqrt(rad2-dzmin*dzmin)
            dz1 = dzmin
        else
            ! edge (dxmax,dzmin) - (dxmax,dzmax) intersects the circle
            a1 = acos(dxmax/rad)
            dx1 = dxmax
            dz1 = sqrt(rad2-dxmax*dxmax)
        endif
        if ((dxmin*dxmin + dzmax*dzmax) > rad2) then
            ! edge (dxmin,dzmin) - (dxmin,dzmax) intersects the circle
            a2 = acos(dxmin/rad)
            dx2 = dxmin
            dz2 = sqrt(rad2-dxmin*dxmin)
        else
            ! edge (dxmin,dzmax) - (dxmax,dzmax) intersects the circle
            a2 = asin(dzmax/rad)
            dx2 = sqrt(rad2-dzmax*dzmax)
            dz2 = dzmax
        endif
        area = 0.5_fp*(a2-a1)*rad2 &
             & + 0.5_fp*dx1*dz1 &
             & + 0.5_fp*dx2*dz2 &
             & - dxmin*dz2 &
             & - (dx1-dxmin)*dzmin
    endif
end function intersect_circle_rect4

!> @brief Retrieve turbine thrust coefficient for a given velocity
function turbinecurve(curves, curvenr, u, errorstring) result(thrustcoef)
    use precision
    use table_handles, only: gettabledata  
    implicit none
!
! Call variables
!
    integer, dimension(2), intent(in) :: curvenr     !< index of turbine thrust curve
    type(handletype)     , intent(in) :: curves      !< thrust and power curve tables
    real(fp)             , intent(in) :: u           !< velocity u
    real(fp)                          :: thrustcoef  !< resulting thrust coefficient
    character(256)       ,intent(out) :: errorstring !< error string
!
! Local variables
!
    real(fp), dimension(1)            :: thrustvar   ! resulting thrust coefficient
!
!! executable statements -------------------------------------------------------
!
    call gettabledata(curves, curvenr(1), curvenr(2), 1, 1, &
                    & thrustvar, u, 0, errorstring)
    thrustcoef = thrustvar(1) 
end function turbinecurve

   !> Allocates and initializes all "valstruct"(:,:) arrays.
   !! Used for history output and/or restart file output for hydraulic structures.
   subroutine init_structure_hisvalues()
      use m_flowexternalforcings , only: npumpsg, ncgensg, ngatesg, ncdamsg, ngategen, ngenstru, nweirgen, ndambreaksg
      !use m_structures, only: NUMVALS_PUMP, NUMVALS_GATE, NUMVALS_CDAM, NUMVALS_CGEN, &
      !                        NUMVALS_GATEGEN, NUMVALS_WEIRGEN, NUMVALS_GENSTRU
      use m_alloc
      use m_flowtimes, only: ti_rst
      implicit none

      if( jahispump > 0 .and. npumpsg > 0) then
         if( allocated( valpump ) ) deallocate( valpump )
         allocate( valpump(NUMVALS_PUMP,npumpsg) ) ; valpump = 0d0
      endif
      if( jahiscgen > 0 ) then
         if( ncgensg > 0 ) then
            if( allocated( valcgen ) ) deallocate( valcgen )
            allocate( valcgen(NUMVALS_CGEN,ncgensg) ) ; valcgen = 0d0
         endif
         
         if (ngenstru == 0) then ! If it is new general structure, then it is stored in the network type
            ngenstru = network%sts%numGeneralStructures
         end if
         if( ngenstru > 0 ) then
            if( allocated( valgenstru ) ) deallocate( valgenstru )
            allocate( valgenstru(NUMVALS_GENSTRU,ngenstru) ) ; valgenstru  = 0d0
         endif
      endif
      if( jahisgate > 0 ) then
         if( ngatesg > 0 ) then
            if( allocated( valgate ) ) deallocate( valgate )
            allocate( valgate(NUMVALS_CGEN,ngatesg) ) ; valgate = 0d0
         endif
         if( ngategen > 0 ) then
            if( allocated( valgategen ) ) deallocate( valgategen )
            allocate( valgategen(NUMVALS_GATEGEN,ngategen) ) ; valgategen = 0d0
         endif
      endif
      if( jahiscdam > 0 .and. ncdamsg > 0) then
         if( allocated( valcdam) ) deallocate( valcdam )
         allocate( valcdam(NUMVALS_CDAM,ncdamsg) ) ; valcdam = 0d0
      endif
      if (nweirgen == 0) then ! If it is new 1D weir, the weir is stored in the network type
         nweirgen = network%sts%numWeirs
      end if
      
      if( jahisweir > 0 .and. nweirgen > 0) then
         if( allocated( valweirgen) ) deallocate( valweirgen )
         allocate( valweirgen(NUMVALS_WEIRGEN,nweirgen) ) ; valweirgen = 0d0
      endif
      if( jahisdambreak > 0 .and. ndambreaksg > 0) then
         if( allocated( valdambreak ) ) deallocate( valdambreak )
         allocate( valdambreak(NUMVALS_DAMBREAK,ndambreaksg) ) ; valdambreak = 0d0
      endif
      if( jahisorif > 0 .and. network%sts%numOrifices > 0) then
         if( allocated( valorifgen) ) deallocate( valorifgen )
         allocate( valorifgen(NUMVALS_ORIFGEN,network%sts%numOrifices) ) ; valorifgen = 0d0
      endif
      if( jahisbridge > 0 .and. network%sts%numBridges > 0) then
         if( allocated( valbridge) ) deallocate( valbridge )
         allocate( valbridge(NUMVALS_BRIDGE,network%sts%numBridges) ) ; valbridge = 0d0
      endif
      if( (ti_rst > 0 .or. jahisculv > 0) .and. network%sts%numCulverts > 0) then
         if( allocated( valculvert) ) deallocate( valculvert )
         allocate( valculvert(NUMVALS_CULVERT,network%sts%numCulverts) ) ; valculvert = 0d0
      endif
      if( jahisuniweir > 0 .and. network%sts%numUniWeirs > 0) then
         if( allocated( valuniweir ) ) deallocate( valuniweir )
         allocate( valuniweir(NUMVALS_UNIWEIR,network%sts%numUniWeirs) ) ; valuniweir = 0d0
      endif
      if( jahiscmpstru > 0 .and. network%cmps%count > 0) then
         if( allocated( valcmpstru ) ) deallocate( valcmpstru )
         allocate( valcmpstru(NUMVALS_CMPSTRU,network%cmps%count) ) ; valcmpstru = 0d0
      endif

! TIDAL TURBINES: Insert init_turbines here

   end subroutine init_structure_hisvalues

!> Sets ALL (scalar) variables in this module to their default values.
!! For a reinit prior to flow computation, only call reset_structures() instead.
subroutine default_structures()

call tree_destroy(strs_ptr)

call reset_structures()

! TIDAL TURBINES: Insert calls to deallocate_turbines and init_turbines here

   ! default settings for structure output to history file
   jahiscgen = 1
   jahispump = 1
   jahisgate = 1
   jahiscdam = 1
   jahisweir = 1
   jahisorif = 1
   jahisculv = 1
   jahisbridge   = 1
   jahisdambreak = 1
   jahisuniweir = 1
   jahiscmpstru = 1

end subroutine default_structures


!> Resets only structures variables intended for a restart of an existing flow simulation (same MDU).
!! Upon loading of new model/MDU, call default_structures() instead.
subroutine reset_structures()
   if (allocated(gates)) deallocate(gates)
end subroutine reset_structures

!> Fills the valstruct array for one given structure on a given link LL.
!! All values are filled, both the generic ones, as well as the type-specific ones.
!! Note: old-style structures may call this with istrtypein = ST_UNSET.
subroutine fill_valstruct_perlink(valstruct, L, dir, istrtypein, istru, L0)
   use m_missing, only: dmiss
   use m_flow, only: q1, s1, au, hu
   use m_flowgeom, only: wu, ln, teta, bl
   use m_1d_structures, only: get_discharge_under_compound_struc
   use m_General_Structure
   implicit none
   double precision, dimension(:), intent(inout) :: valstruct   !< Output values on structure (e.g. valweirgen(:)):
                                                                !< (1) total width
                                                                !< (2) structure discharge
                                                                !< (3) structure water level up
                                                                !< (4) structure water level down
                                                                !< (5) structure head
                                                                !< (6) flow area
                                                                !< (7) velocity
                                                                !< (8) water level on crest
                                                                !< (9) crest level
                                                                !< (10) crest width
                                                                !< (11) state
                                                                !< (12) force difference per unit width
                                                                !< (13) gate opening width
                                                                !< (14) gate lower edge level
                                                                !< (15) gate opening height
                                                                !< (16) gate upper edge level
                                                                !< (17) discharge through gate opening
                                                                !< (18) discharge over gate
                                                                !< (19) discharge under gate
                                                                !< (20) flow area in gate opening
                                                                !< (21) flow area over gate
                                                                !< (22) flow area under gate
                                                                !< (23) velocity through gate opening
                                                                !< (24) velocity over gate
                                                                !< (25) velocity under gate
   integer,                        intent(in   ) :: L           !< Flow link number.
   double precision,               intent(in   ) :: dir         !< Direction of flow link w.r.t. structure orientation (1.0 for same direction, -1.0 for opposite).
   integer,                        intent(in   ) :: istrtypein  !< The type of the structure. May differ from the struct%type, for example:
                                                                !< an orifice should be called with istrtypein = ST_ORIFICE, whereas its struct(istru)%type = ST_GENERAL_ST.
   integer,                        intent(in   ) :: istru       !< Structure index in network%sts set.
   integer,                        intent(in   ) :: L0          !< Local flow link index in the struct%linknumbers array.

   integer :: ku, kd, k1, k2
   type(t_GeneralStructure), pointer :: genstr
   double precision :: qcmp

   if (dir > 0) then
      ku = ln(1,L)
      kd = ln(2,L)
   else
      ku = ln(2,L)
      kd = ln(1,L)
   end if

   ! 1. Generic values that apply to all structure types
   valstruct(1) = valstruct(1) + wu(L)
   
   if (istru > 0) then ! When it is not old weir and not old general structure and not a compound structure
      if (network%sts%struct(istru)%compound > 0) then ! for a structure that belongs to a compound structure
         k1 = ln(1,L)
         k2 = ln(2,L)
         if (hu(L) > 0) then
            qcmp = get_discharge_under_compound_struc(network%sts%struct(istru), L0, s1(k1), s1(k2), teta(L))
         else
            qcmp = 0d0
         end if
         valstruct(2) = valstruct(2) + qcmp*dir
      else
         valstruct(2) = valstruct(2) + q1(L)*dir
      end if
   else
      valstruct(2) = valstruct(2) + q1(L)*dir
   end if
   
   valstruct(3) = valstruct(3) + s1(ku)*wu(L)
   valstruct(4) = valstruct(4) + s1(kd)*wu(L)
   valstruct(5) = valstruct(5) + (s1(ku) - s1(kd))*wu(L)

   if (istrtypein /= ST_PUMP) then ! Compute flow area for structures except for pump
      if (istru > 0) then ! When it is not old weir and not old general structure and not a compound structure
         if (network%sts%struct(istru)%compound > 0) then ! for a structure that belongs to a compound structure
            valstruct(6) = valstruct(6) + network%sts%struct(istru)%au(L0)
         else
            valstruct(6) = valstruct(6) + au(L)
         end if
      else
         valstruct(6) = valstruct(6) + au(L)
      end if
   end if

   ! 2. More specific valus that apply to certain structure types only

   ! General structure-based structures with a crest.
   if (any(istrtypein == (/ ST_GENERAL_ST, ST_WEIR, ST_ORIFICE /))) then ! TODO: ST_GATE
      valstruct(8)  = valstruct(8) + network%sts%struct(istru)%generalst%sOnCrest(L0)*wu(L)
      valstruct(12) = valstruct(12) + get_force_difference(istru, L0)*wu(L)
   end if
   
   ! General structure-based structures with a (gate) door.
   if (any(istrtypein == (/ ST_GENERAL_ST /))) then ! TODO: ST_GATE
      k1 = ln(1,L)
      k2 = ln(2,L)

      genstr => network%sts%struct(istru)%generalst
      if (hu(L) > 0) then
         valstruct(17) = valstruct(17) + get_discharge_through_gate_opening(genstr, L0, s1(k1), s1(k2))*dir
         valstruct(18) = valstruct(18) + get_discharge_over_gate(genstr, L0, s1(k1), s1(k2))*dir
         valstruct(19) = valstruct(19) + get_discharge_under_gate(genstr, L0, s1(k1), s1(k2))*dir
      
         valstruct(20) = valstruct(20) + genstr%au(3,L0) ! flow area through gate opening
         valstruct(21) = valstruct(21) + genstr%au(2,L0) ! flow area over gate
         valstruct(22) = valstruct(22) + genstr%au(1,L0) ! flow area under gate
      end if
   end if
   
   ! 3. More specific values that apply to bridge
   if (istrtypein == ST_BRIDGE) then
      valstruct(8)  = valstruct(8) + bl(ku)*wu(L)
      valstruct(9)  = valstruct(9) + bl(kd)*wu(L)
      valstruct(10) = valstruct(10) + network%sts%struct(istru)%bridge%bedLevel_actual*wu(L)
   end if

end subroutine fill_valstruct_perlink


!> Averages the values on one structure across all links,
!! where needed taking care of partition models.
!! Note 1: fill_valstructs_perlink must have been called in
!! a loop prior to calling this averaging routine.
!! Note 2: if it is a general structure (jagenst == 1), then (6)-(12) are computed as well.
subroutine average_valstruct(valstruct, istrtypein, istru, nlinks, icount)
   use m_missing, only: dmiss
   use m_partitioninfo, only: jampi
   use m_1d_structures
   use m_General_Structure, only: t_GeneralStructure
   implicit none
   double precision, dimension(:), intent(inout) :: valstruct   !< Output values on structure (e.g. valpump(:)):
                                                                !< (1) total width (unchanged)
                                                                !< (2) structure discharge (unchanged)
                                                                !< (3) structure water level up (averaged)
                                                                !< (4) structure water level down (averaged)
                                                                !< (5) structure head (averaged)
                                                                !< (6) flow area (unchanged)
                                                                !< (7) velocity (computed)
                                                                !< (8) water level on crest (averaged)
                                                                !< (9) crest level (computed)
                                                                !< (10) crest width (computed)
                                                                !< (11) state (if all links have the same state, then write it. Otherwise it is missing value)
                                                                !< (12) force difference per unit width (averaged)
                                                                !< (13) gate opening width
                                                                !< (14) gate lower edge level
                                                                !< (15) gate opening height
                                                                !< (16) gate upper edge level
                                                                !< (17) discharge through gate opening
                                                                !< (18) discharge over gate
                                                                !< (19) discharge under gate
                                                                !< (20) flow area in gate opening
                                                                !< (21) flow area over gate
                                                                !< (22) flow area under gate
                                                                !< (23) velocity through gate opening
                                                                !< (24) velocity over gate
                                                                !< (25) velocity under gate
                                                                !< (icount) counters of partitions for parallel
   integer,                        intent(in   ) :: istrtypein  !< The type of the structure. May differ from the struct%type, for example:
                                                                !< an orifice should be called with istrtypein = ST_ORIFICE, whereas its struct(istru)%type = ST_GENERAL_ST.
   integer,                        intent(in   ) :: istru       !< Structure index in network%sts set.
   integer,                        intent(in   ) :: nlinks      !< Number of flow links for this structure (on the current partition)
   integer,                        intent(in   ) :: icount      !< Index of the counter element in valstruct array,
                                                                !! it is the last element of the array.
   
   integer:: i, tmp, jadif
   type(t_structure), pointer :: pstru
   type(t_GeneralStructure), pointer :: genstr
   
   ! 1. Generic values that apply to all structure types
   if (jampi == 0) then
      if (valstruct(1) == 0d0 ) then ! zero width
         valstruct(2) = dmiss  ! discharge
         valstruct(3) = dmiss  ! s1up
         valstruct(4) = dmiss  ! s1down
         valstruct(5) = dmiss  ! head
         if (istrtypein /= ST_PUMP) then
            valstruct(6) = dmiss ! flow area
            valstruct(7) = dmiss ! velocity
         end if

         if (any(istrtypein == (/ ST_GENERAL_ST, ST_WEIR, ST_ORIFICE /))) then ! TODO: ST_GATE
            valstruct(8) = dmiss ! water level on crest
            valstruct(9) = dmiss ! crest level
            valstruct(10)= dmiss ! crest width
            valstruct(11)= dmiss ! state
            valstruct(12)= dmiss ! force difference per unit width
         end if
      else
         ! valstruct(2): keep discharge at the summed value
         ! Average the remaining values:
         valstruct(3) = valstruct(3) / valstruct(1)        ! s1up
         valstruct(4) = valstruct(4) / valstruct(1)        ! s1down
         valstruct(5) = valstruct(5) / valstruct(1)        ! head

         if (istrtypein /= ST_PUMP) then
            if (valstruct(6) > 0d0) then ! non-zero flow area
               valstruct(7) = valstruct(2) / valstruct(6)  ! velocity
            else
               valstruct(7) = 0d0
            end if
         end if

         if (any(istrtypein == (/ ST_GENERAL_ST, ST_WEIR, ST_ORIFICE /))) then ! TODO: ST_GATE
            pstru => network%sts%struct(istru)
            valstruct(8) = valstruct(8) / valstruct(1)     ! water level on crest
            valstruct(12)= valstruct(12)/ valstruct(1)     ! force difference per unit width
            
         end if
      endif
   endif

   ! 2. More specific valus that apply to certain structure types only

   ! General structure-based structures with a crest.
   if (any(istrtypein == (/ ST_GENERAL_ST, ST_WEIR, ST_ORIFICE /)) & ! TODO: ST_GATE
       .and. nlinks > 0) then ! If it is a new general structure, and there are links
      valstruct(icount) = 1                     ! count the current partition
      valstruct(9) = get_crest_level(pstru)     ! crest level
      valstruct(10)= get_width(pstru)           ! crest width
      ! determine state
      tmp = maxval(pstru%generalst%state(1:3,1))
      jadif = 0
      do i = 2, nlinks
         if (tmp /= maxval(pstru%generalst%state(1:3,i))) then
            jadif = 1
            exit
         end if
      end do
      if (jadif == 0) then
         valstruct(11) = dble(tmp)
      else
         valstruct(11) = dmiss
      end if
   end if

   ! General structure-based structures with a (gate) door.
   if (any(istrtypein == (/ ST_GENERAL_ST, ST_ORIFICE /))) then ! TODO: ST_GATE
      if (nlinks > 0) then ! If it is a new general structure, and there are links
         genstr => network%sts%struct(istru)%generalst
         valstruct(13) = genstr%gateopeningwidth_actual           ! gate opening width
         valstruct(14) = get_gle(pstru)                           ! gate lower edge level
         valstruct(15) = get_opening_height(pstru)                ! gate opening height
         valstruct(16) = valstruct(14) + genstr%gatedoorheight    ! gate upper edge level
         valstruct(icount) = 1
      end if

      if (jampi == 0 ) then
         if (valstruct(1) == 0d0) then ! zero width
            valstruct(13:) = dmiss
         else
            ! only for general structure
            if (istrtypein == ST_GENERAL_ST) then 
               if (valstruct(20) > 0) then ! flow area in gate opening
                  valstruct(23) = valstruct(17) / valstruct(20) ! velocity through gate opening
               end if
               if (valstruct(21) > 0) then ! flow area over gate
                  valstruct(24) = valstruct(18) / valstruct(21) ! velocity over gate
               end if
               if (valstruct(22) > 0) then ! flow area under gate
                  valstruct(25) = valstruct(19) / valstruct(22) ! velocity under gate
               end if
            end if
         end if
      end if 
   end if
   
   ! 3. More specific values that apply to bridge
   if (istrtypein == ST_BRIDGE) then
      if (jampi == 0) then
         if (valstruct(1) == 0d0 ) then ! zero width
            valstruct(8) = dmiss
            valstruct(9) = dmiss
            valstruct(10)= dmiss
         else
            valstruct(8)  = valstruct(8) / valstruct(1)
            valstruct(9)  = valstruct(9) / valstruct(1)
            valstruct(10) = valstruct(10)/ valstruct(1)
         end if
      end if
   end if

end subroutine average_valstruct


!!> Gets force difference per unit width over structure (weir, gate, general structure) per link.
double precision function get_force_difference(istru, L)
   use m_missing
   use m_flowgeom, only: ln
   use m_flow, only: s1
   use m_1d_structures, only: get_crest_level
   implicit none   
   integer, intent(in   )   :: istru !< structure index
   integer, intent(in   )   :: L     !< current link L
   
   double precision  :: s1up   !< water level up
   double precision  :: s1dn   !< water level down
   double precision  :: crestl
   integer           :: k1, k2
   double precision  :: rholeft, rhoright
   
   crestl = get_crest_level(network%sts%struct(istru))
  
   k1 = ln(1,L)
   k2 = ln(2,L)
   s1up = max(s1(k1), s1(k2))
   s1dn = min(s1(k1), s1(k2))
   if (crestl > dmiss + 0.1d0) then
      rholeft  = 1000.0d0
      rhoright = 1000.0d0
      
      get_force_difference =  max((s1up - crestl), 0.0d0)**2 * rholeft  * gravity / 2.0d0 -  &
                            max((s1dn - crestl), 0.0d0)**2 * rhoright * gravity / 2.0d0
   else
      get_force_difference = dmiss
   end if

end function get_force_difference


!> Gets discharge through gate opening per link.
double precision function get_discharge_through_gate_opening(genstr, L0, s1m1, s1m2)
   use m_missing
   use m_General_Structure
   implicit none   
   type(t_GeneralStructure), pointer, intent(in   ) :: genstr !< Derived type containing general structure information.
   integer,                           intent(in   ) :: L0     !< Local link index in genstr%..(:) link-based arrays.
   double precision,                  intent(in   ) :: s1m1   !< (geometrical) upstream water level.
   double precision,                  intent(in   ) :: s1m2   !< (geometrical) downstream water level.
   double precision  :: u1L, dsL, gatefraction
   
   dsL = s1m2 - s1m1 
   gatefraction = genstr%gateclosedfractiononlink(L0)
   
   if (gatefraction > gatefrac_eps) then
      u1L = genstr%ru(3,L0) - genstr%fu(3,L0)*dsL
      get_discharge_through_gate_opening = genstr%au(3,L0) * u1L
   else
      get_discharge_through_gate_opening = 0d0
   end if

end function get_discharge_through_gate_opening

!> Gets discharge over gate opening per link.
double precision function get_discharge_over_gate(genstr, L0, s1m1, s1m2)
   use m_missing
   use m_General_Structure
   implicit none   
   type(t_GeneralStructure), pointer, intent(in   ) :: genstr !< Derived type containing general structure information
   integer,                           intent(in   ) :: L0     !< Local link index in genstr%..(:) link-based arrays.
   double precision,                  intent(in   ) :: s1m1   !< (geometrical) upstream water level.
   double precision,                  intent(in   ) :: s1m2   !< (geometrical) downstream water level.
   double precision  :: u1L, dsL, gatefraction
   
   dsL = s1m2 - s1m1
   gatefraction = genstr%gateclosedfractiononlink(L0)
   
   if (gatefraction > gatefrac_eps) then
      u1L = genstr%ru(2,L0) - genstr%fu(2,L0)*dsL
      get_discharge_over_gate = genstr%au(2,L0) * u1L
   else
      get_discharge_over_gate = 0d0
   end if

end function get_discharge_over_gate

!> Gets discharge under gate per link.
double precision function get_discharge_under_gate(genstr, L0, s1m1, s1m2)
   use m_missing
   use m_General_Structure
   implicit none   
   type(t_GeneralStructure), pointer, intent(in   ) :: genstr !< Derived type containing general structure information
   integer,                           intent(in   ) :: L0     !< Local link index in genstr%..(:) link-based arrays.
   double precision,                  intent(in   ) :: s1m1   !< (geometrical) upstream water level.
   double precision,                  intent(in   ) :: s1m2   !< (geometrical) downstream water level.
   double precision  :: u1L, dsL, gatefraction
   
   dsL = s1m2 - s1m1
   gatefraction = genstr%gateclosedfractiononlink(L0)
   
   if (gatefraction > gatefrac_eps) then
      u1L = genstr%ru(1,L0) - genstr%fu(1,L0)*dsL
      get_discharge_under_gate = genstr%au(1,L0) * u1L
   else
      get_discharge_under_gate = 0d0
   end if

end function get_discharge_under_gate

!> Updates structure parameters for the output to restart file.
!! Only computes the needed values, and
!! only when they are not computed for history output.
!! Values are stored in the val*(:,:) arrays, shared with history output.
subroutine structure_parameters_rst()
   use m_1d_structures, only: get_opening_height
   implicit none
   integer :: n, istru

   do n = 1, network%sts%numCulverts
      istru = network%sts%culvertIndices(n)
      valculvert(11,n) = get_opening_height(network%sts%struct(istru))
   end do

end subroutine structure_parameters_rst
end module m_structures