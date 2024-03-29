module m_dredge_data
      !
      ! Contains parameters and types related to dredging and dumping facilities
      !
      use precision
      use handles
      
      implicit none
      
      public dredtype
      public dumptype
      
      integer, parameter :: DEPTHDEF_REFPLANE        = 1
      integer, parameter :: DEPTHDEF_WATERLVL        = 2
      integer, parameter :: DEPTHDEF_MAXREFWL        = 3
      integer, parameter :: DEPTHDEF_MINREFWL        = 4
      integer, parameter :: DEPTHDEF_MAX             = 4
   !
   ! Dredge type enumeration
   !
      integer, parameter :: DREDGETYPE_DREDGING      = 1
      integer, parameter :: DREDGETYPE_SANDMINING    = 2
      integer, parameter :: DREDGETYPE_NOURISHMENT   = 3
   !
   ! Dredge trigger type enumeration
   !
      integer, parameter :: DREDGETRIG_POINTBYPOINT  = 1
      integer, parameter :: DREDGETRIG_ALLBYONE      = 2
      integer, parameter :: DREDGETRIG_ALLBYAVG      = 3
   !
   ! Dredge distribution enumeration
   !
      integer, parameter :: DREDGEDISTR_UNIFORM      = 1
      integer, parameter :: DREDGEDISTR_HIGHEST      = 2
      integer, parameter :: DREDGEDISTR_PROPORTIONAL = 3
      integer, parameter :: DREDGEDISTR_HIGHFIRST    = 4
      integer, parameter :: DREDGEDISTR_LOWFIRST     = 5
      integer, parameter :: DREDGEDISTR_SHALLOWEST   = 6
      integer, parameter :: DREDGEDISTR_SHALLOWFIRST = 7
      integer, parameter :: DREDGEDISTR_DEEPFIRST    = 8
      integer, parameter :: DREDGEDISTR_MAX          = 8
   !
   ! Dredge to dump distribution enumeration
   !
      integer, parameter :: DR2DUDISTR_PERCENTAGE    = 1
      integer, parameter :: DR2DUDISTR_SEQUENTIAL    = 2
      integer, parameter :: DR2DUDISTR_PROPORTIONAL  = 3
      integer, parameter :: DR2DUDISTR_MAX           = 3
   !
   ! Dump distribution enumeration
   !
      integer, parameter :: DUMPDISTR_UNIFORM        = 1
      integer, parameter :: DUMPDISTR_LOWEST         = 2
      integer, parameter :: DUMPDISTR_DEEPEST        = 3
      integer, parameter :: DUMPDISTR_PROPORTIONAL   = 4
      integer, parameter :: DUMPDISTR_MAX            = 4 
   !
   ! Location check whether cell is included in dredge/dump activity
   !
      integer, parameter :: CHKLOC_ALLCORNER         = 1
      integer, parameter :: CHKLOC_CENTRE            = 2
      integer, parameter :: CHKLOC_ANYCORNER         = 3
   !
   ! collection of dredging area related parameters
   !
   type dredtype
      real(fp)                                :: alpha_dh        ! multiplication factor between 0.0 and 0.5 
                                                                 ! when including duneheights in the trigger height
      real(fp)                                :: dredge_depth    ! dredge depth of dredge area
      real(fp)                                :: maxvolrate      ! maximum dredging volume rate
      real(fp)                                :: clearance       ! extra depth when dredging
      real(fp)                                :: plough_effic    ! efficiency of dune ploughing in reducing dune height
      real(fp)      , dimension(2)            :: dredgeloc       ! dredgelocation (x,y)
      real(fp)                                :: totalvolsupl    ! total volume to be nourished.
      real(fp)      , dimension(:)  , pointer :: area            ! (npnt)   grid cell area at nm point (gsqs)
      real(fp)      , dimension(:)  , pointer :: hdune           ! (npnt)   dune height at nm point
      real(fp)      , dimension(:)  , pointer :: dz_dredge       ! (npnt)   thickness to be removed at nm point
                                                                 !          from the sediment layer due to dredging
      real(fp)      , dimension(:)  , pointer :: reflevel        ! (npnt)   current reference level at nm point
      real(fp)      , dimension(:)  , pointer :: dunetoplevel    ! (npnt)   dune top level at nm point
      real(fp)      , dimension(:)  , pointer :: triggerlevel    ! (npnt)   trigger level at nm point for dredging
      real(fp)      , dimension(:)  , pointer :: bedlevel        ! (npnt)   bed level at nm point
      real(fp)      , dimension(:)  , pointer :: troughlevel     ! (npnt)   dune trough level at nm point
      real(fp)      , dimension(:)  , pointer :: sedimentdepth   ! (npnt)   sediment depth at nm point
      real(fp)      , dimension(:)  , pointer :: sortvar         ! (npnt)   variable used for sorting
      !
      integer       , dimension(4)            :: paractive       ! parameters for getting on/off status from time-series; active(1)==0 means always on
      integer                                 :: depthdef        ! 1 = depth relative to reference plane
                                                                 ! 2 = depth relative to water level
                                                                 ! 3 = depth relative to maximum of reference plane and water level
                                                                 ! 4 = depth relative to minimum of reference plane and water level
      integer                                 :: triggertype     ! 1 = trigger point by point (previous trigger_all = .false.)
                                                                 ! 2 = one point triggers all (previous trigger_all = .true. )
                                                                 ! 3 = average value triggers dredging (of all points)
      integer                                 :: dredgedistr     ! 1 = dredge uniformly
                                                                 ! 2 = dredge high locations first
                                                                 ! 3 = dredge proportionally (default)
                                                                 ! 4 = dredge highest first
                                                                 ! 5 = dredge lowest first
                                                                 ! 6 = dredge shallow locations first
                                                                 ! 7 = dredge shallowest first
                                                                 ! 8 = dredge deepest first
      integer                                 :: ichkloc         ! 1 = all cell corners within polygon
                                                                 ! 2 = cell center within polygon
                                                                 ! 3 = any one of the cell corners within polygon
      integer                                 :: idx_type        ! dredging number (only for itype=1 or 2)
                                                                 ! suppletion number (only for itype=3)
      integer                                 :: itype           ! 1 = critical level
                                                                 ! 2 = sand mining
                                                                 ! 3 = sediment nourishment
      integer                                 :: dumpdistr       ! distribution of sediment over dump areas
                                                                 ! 1 = percentage (only option if not dumplimited)
                                                                 ! 2 = sequential (default for dumplimited)
                                                                 ! 3 = proportional to capacity
      integer                                 :: outletlink      ! link number of outlet (sediment removed from model)
      integer       , dimension(:)  , pointer :: inm             ! (npnt)   index 1 to npnt into nm array
      integer       , dimension(:)  , pointer :: nm              ! (npnt)   nm points
      integer                                 :: npnt            ! number of points in dredging area
      !
      logical                                 :: active          ! T: dredge this time step
      logical       , dimension(:)  , pointer :: triggered       ! (npnt)   flag to indicate whether dredging was triggered at nm point
      logical                                 :: stilldredging   ! T: continue dredging at next time step
      logical                                 :: dredgewhendry   ! T: dredge also at dry points (kfsed /= 1)
      logical                                 :: dumplimited     ! T: dredging limited to the amount that can be dumped
      logical                                 :: if_morfac_0     ! T: allow instanteneous dredging while morfac equals 0
      logical                                 :: in1domain       ! T: dredge area located in only 1 domain
      logical                                 :: obey_cmp        ! T: stop dredging if bed composition is empty
      logical                                 :: use_dunes       ! T: include effect of dune heights on dredging and vice versa
      !
      character( 80)                          :: name            ! name identifying dredge area
   end type dredtype
   
   !
   ! collection of dump area related parameters
   !
   type dumptype
      real(fp), dimension(2)                  :: dumploc         ! dredgelocation (x,y)
      real(fp)                                :: mindumpdepth    ! Minimum Dump Depth
      real(fp)      , dimension(:)  , pointer :: area            ! (npnt)   grid cell area at nm point (gsqs)
      real(fp)      , dimension(:)  , pointer :: hdune           ! (npnt)   dune height at nm point
      real(fp)      , dimension(:)  , pointer :: reflevel        ! (npnt)   current reference level at nm point
      real(fp)      , dimension(:)  , pointer :: bedlevel        ! (npnt)   bed level at nm point
      real(fp)      , dimension(:)  , pointer :: dz_dump         ! (npnt)   thickness of dumped sediment at nm point
      real(fp)      , dimension(:)  , pointer :: sortvar         ! (npnt)   variable used for sorting
      !
      integer                                 :: depthdef        ! 1 = depth relative to reference plane
                                                                 ! 2 = depth relative to water level
                                                                 ! 3 = depth relative to maximum of reference plane and water level
                                                                 ! 4 = depth relative to minimum of reference plane and water level
      integer       , dimension(:)  , pointer :: inm             ! (npnt)   index 1 to npnt into nm array
      integer       , dimension(:)  , pointer :: nm              ! nm points
      integer                                 :: npnt            ! number of nm points
      integer                                 :: dumpdistr       ! 1 = dump sediment uniformly (default)
                                                                 ! 2 = dump in lowest locations first
                                                                 ! 3 = dump in deepest locations first
      integer                                 :: ichkloc         ! 1 = all cell corners within polygon
                                                                 ! 2 = cell center within polygon
                                                                 ! 3 = any one of the cell corners within polygon
      !
      logical                                 :: dumpwhendry     ! T: dump also at dry points (kfsed /= 1)
      logical                                 :: dumpcapaflag    ! T: dump area is capacity limited (mindumpdepth /= -999)
      logical                                 :: in1domain       ! T: dump area located in only 1 domain
      logical                                 :: use_dunes       ! T: include effect of dune heights on dumping and vice versa
      !
      character( 80)                          :: name            ! name identifying dredge area
   end type dumptype

   type sv_dredge
      type (handletype)                       :: tseriesfile     ! table  containing dredging time-series
      !
      real(fp)      , dimension(:,:), pointer :: link_percentage ! (nalink,lsedtot) distribution of dredged material
                                                                 !          from dredge to dump areas
      real(fp)      , dimension(:)  , pointer :: link_distance   ! (nalink) distance from dredge to dump area
      real(fp)      , dimension(:,:), pointer :: link_sum        ! (nalink,lsedtot) cumulative dredged sediment
                                                                 !          transported via this link
      real(fp)      , dimension(:)  , pointer :: dzdred          ! (nmmax)  thickness to be removed from the sediment layer
                                                                 !          due to dredging/sandmining
                                                                 !          dumping is handled via array dbodsd
      real(fp)      , dimension(:)  , pointer :: refplane        ! (nmmax) reference plane for the depth (default: 0)
      real(fp)      , dimension(:,:), pointer :: voldred         ! (nadred,lsedtot+1) volume dredged sediment at a timestep
                                                                 ! lsedtot+1 used to store dredged "non-modelled subsoil sediment"
      real(fp)      , dimension(:)  , pointer :: voldune         ! (nmmax) volume to be dredged when considering dunes
      real(fp)      , dimension(:)  , pointer :: totvoldred      ! (nadred) total volume dredged
      real(fp)      , dimension(:,:), pointer :: totvoldredfrac  ! (nadred+nasupl,lsedtot+1) total volume dredged per fraction
      real(fp)      , dimension(:)  , pointer :: globalareadred  ! (nadred) global area for dredging (over all domains)
      real(fp)      , dimension(:,:), pointer :: voldump         ! (nadump,lsedtot) volume dredged sediment at a timestep
      real(fp)      , dimension(:,:), pointer :: percsupl        ! (nasupl,lsedtot) percentage of sediment fraction for nourishment 
      real(fp)      , dimension(:)  , pointer :: totvoldump      ! (nadump) total volume dumped
      real(fp)      , dimension(:,:), pointer :: totvoldumpfrac  ! (nadump,lsedtot+1) total volume dumped
      real(fp)      , dimension(:)  , pointer :: localareadump   ! (nadump) local area for dumping (only this domain)
      real(fp)      , dimension(:)  , pointer :: globalareadump  ! (nadump) global area for dumping (over all domains)
      real(fp)      , dimension(:)  , pointer :: globaldumpcap   ! (nadump) global dump capacity (over all domains)
      !
      integer                                 :: dredge_domainnr ! domain number of dredge
      integer                                 :: dredge_ndomains ! number of domains that use dredge&dump functionality
      integer                                 :: nadred          ! number of dredge areas
      integer                                 :: nadump          ! number of dump areas
      integer                                 :: nasupl          ! number of nourishment areas
      integer                                 :: nalink          ! number of links
      real(fp)                                 :: ntimaccum       ! number of (half) timesteps over which ndredge and nplough has been accumulated
      !
      integer       , dimension(:,:), pointer :: link_def        ! (nalink,2) actual transports from dredge (1st column)
                                                                 !            to dump (2nd column) areas
      real(fp)       , dimension(:)  , pointer :: ndredged        ! (nadred) number of (half) timesteps dredged
      real(fp)       , dimension(:)  , pointer :: nploughed       ! (nadred) number of (half) timesteps ploughed
      !
      logical                                 :: tsmortime       ! T: time scale of time-series is morphological time
      logical                                 :: firstdredge     ! T: first time in computational dredge routine
      !
      character(256)                          :: dredgefile      ! name of dredge input file
      character( 80), dimension(:)  , pointer :: dredge_areas    ! (nadred) names identifying dredge areas
      character( 80), dimension(:)  , pointer :: dump_areas      ! (nadump) names identifying dump   areas
      type (dredtype), dimension(:) , pointer :: dredge_prop     ! (nadred) dredge area properties
      type (dumptype), dimension(:) , pointer :: dump_prop       ! (nadump) dump area properties
       
   end type sv_dredge

end module m_dredge_data
