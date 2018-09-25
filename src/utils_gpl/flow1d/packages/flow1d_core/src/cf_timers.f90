module cf_timers
!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2018.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  This program is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
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
!  $Id: cf_timers.f90 8044 2018-01-24 15:35:11Z mourits $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/utils_gpl/flow1d/packages/flow1d_core/src/cf_timers.f90 $
!-------------------------------------------------------------------------------

      implicit none

!     Function       : Handles everything to conveniently time a FORTRAN program
!
!         contains                 meaning
!         timini   ( )           : Should be called once first to initialize. Logical 'timon'
!                                  must be switched to .true. afterwards by the caller !
!         timinc   ( )           : Used internally if the arrays run out of space
!         timstrt(subrou,ihandl) : Starts timing for this subroutine or program part.
!                                  'subrou' is a max. 40 character ID-string.
!                                  'ihandl' must be saved by the caller
!         timstop(ihandl)        : stops timing for this handle and accumulates the result
!         timdump(filename)      : writes the results to the report file 'filename'
!
!     Example:
!
!     *** in the highest level calling program: ***
!
!     use timers                     ! to make the logical 'timon' available  (before implicit none)
!         .......
!     integer(4) ithndl              ! handle to time this subroutine
!     data       ithndl / 0 /
!     call timini ( )
!     timon = .true.
!     call timstrt( "your main program name", ithndl )
!         .......
!     if ( timon ) then
!        call timstop ( ithndl )
!        call timdump ( ifnam(1:index(ifnam,".",.true.)-1)//'-timers.out' )
!     endif                          ! the file name is constructed here from the model file name
!
!     *** in each subroutine: ***
!
!     use timers                     ! to make the logical 'timon' available  (before implicit none)
!         .......
!     integer(4) ithndl              ! handle to time this subroutine
!     integer(4) ithnd2              ! handle to time a critical section of this subroutine
!     data       ithndl / 0 /        ! initialize it to zero for the first call
!     data       ithnd2 / 0 /        ! initialize it to zero for the first call
!     if ( timon ) call timstrt( "your subroutine name", ithndl )
!         .......
!        if ( timon ) call timstrt( "the name of your critical section", ithnd2 )
!         .......
!           timon = .false.
!                 OMP parallel section with subroutine calls that may contain timer calls
!           timon = .true.
!         .......
!        if ( timon ) call timstop ( ithnd2 )    ! here your critical section stops
!         .......
!     if ( timon ) call timstop ( ithndl )
!     return
!     end
!
!     NOTE1: You must take care that you do not leave your subroutine before you close its timer
!
!     NOTE2: There is no maximum to the number of subroutines or program parts that can be timed.
!            The maximum number of times that a timer for the same subroutine / program part
!            can appear in a different call tree (called context) is 99.
!
!     NOTE3: Since the call-tree is introduced, the computation time of a subroutine or program
!            part may be spilt up into parts corresponding to their different location in the call
!            tree. This may result in multiple lines for the same subroutine or program part in the
!            report file. Each line gives the computation time of the subroutine at that location
!            in the call tree. The total time in the subroutine or program part is then not printed.
!
!     NOTE4: You may start a timer with your own explanatory indentification string at any point
!            in your program, you only must take care that always a corresponding stop of the
!            timer occurs. So you are not limited to the timing of whole subroutines only.
!            The ease of use facilitates to easily split up a routine into several parts that are
!            each timed with an own timer to identify bottle necks within a routine.
!
!     NOTE5: Since the use of subroutines is only advantageous if a subroutine is called multiple
!            times, it is more efficient to time your subroutines from inside like in the example
!            above, than from outside in the calling program. An additional advantage of this
!            policy is that the location of your timing commands is limited to the top and the
!            bottom of your subroutine and not somewhere else throughout the code.
!
!     NOTE6: These timers are NOT treadsafe (yet), so always set 'timon' to .false. before you enter an
!            OMP parallel section. You switch timon to .true. after you left the parallel section.
!            (The instructions to switch 'timon' off and on cannot be nested)
!            Furthermore it is assumed that your handles are saved by the program parts, so you
!            must insert the 'save' instruction for the handle if that is not done automatically.

      logical                             :: timon          ! is the timer switched on or off
!      
      integer  ( 4), private              :: nohmax         ! current maximum size of the timer arrays
      integer  ( 4)                       :: nohandl        ! current highest timer handle
      integer  ( 4), private              :: noshndl        ! current highest subroutine handle
      integer  ( 4), private              :: prevhnd        ! previous timer handle
      integer  ( 4), private              :: dlevel         ! current level in the call tree
      integer  ( 4), private              :: maxlvl         ! maximum level of call trees
      integer  ( 4),          allocatable :: ntimcal(:)     ! call frequency
      integer  ( 4), private, allocatable :: level  (:)     ! call level
      integer  ( 4), private, allocatable :: context(:,:,:) ! call context
      integer  ( 4), private, allocatable :: ncontxt(:)     ! number of contexts per subroutine
      integer  ( 8), private              :: count          ! system clock count
      integer  ( 8), private              :: rate           ! ticks per second
      real     ( 8), private, allocatable :: cpstart(:)     ! to save cpu startimes
      real     ( 8),          allocatable :: cptime (:)     ! to accumulate cpu times
      real     ( 8), private, allocatable :: wcstart(:)     ! to save wall clock startimes
      real     ( 8),          allocatable :: wctime (:)     ! to accumulate wall clock times
      character(40),          allocatable :: tmsubnm(:)     ! name of the subroutine

      integer  ( 4),          allocatable :: ntimcal_prev(:) ! call frequency previous timer his output
      real     ( 8),          allocatable :: cptime_prev(:)  ! accumulated cpu times previous timer his output
      real     ( 8),          allocatable :: wctime_prev(:)  ! accumulated wall clock times previous timer his output

contains

!***************

subroutine timini  ( )

      nohmax  = 500
      nohandl =   0
      noshndl =   0
      prevhnd =   0
      timon   = .false.
      dlevel  =   0
      maxlvl  =   0
      if ( .not. allocated ( ntimcal ) ) allocate ( ntimcal(     nohmax) )
      if ( .not. allocated ( level   ) ) allocate ( level  (     nohmax) )
      if ( .not. allocated ( cpstart ) ) allocate ( cpstart(     nohmax) )
      if ( .not. allocated ( cptime  ) ) allocate ( cptime (     nohmax) )
      if ( .not. allocated ( wcstart ) ) allocate ( wcstart(     nohmax) )
      if ( .not. allocated ( wctime  ) ) allocate ( wctime (     nohmax) )
      if ( .not. allocated ( tmsubnm ) ) allocate ( tmsubnm(     nohmax) )
      if ( .not. allocated ( ncontxt ) ) allocate ( ncontxt(     nohmax) )
      if ( .not. allocated ( context ) ) allocate ( context(99,2,nohmax) )

      if (.not. allocated (ntimcal_prev)) allocate (ntimcal_prev(nohmax))
      if (.not. allocated (cptime_prev)) allocate (cptime_prev(nohmax))
      if (.not. allocated (wctime_prev)) allocate (wctime_prev(nohmax))

      ncontxt =   0
end subroutine timini

!***************

subroutine timinc  ( )

      integer  ( 4), allocatable :: ntimtmp(:)     !  call frequency
      integer  ( 4), allocatable :: levltmp(:)     !  call level
      integer  ( 4), allocatable :: contemp(:,:,:) !  context
      integer  ( 4), allocatable :: ncontmp(:)     !  nr of contexts per subroutine
      real     ( 8), allocatable :: cpstemp(:)     !  to save cp startimes
      real     ( 8), allocatable :: cpttemp(:)     !  to save cp times
      real     ( 8), allocatable :: wcstemp(:)     !  to save wc startimes
      real     ( 8), allocatable :: wcttemp(:)     !  to save wc times
      character(20), allocatable :: tmsutmp(:)     !  name of the subroutine

      integer  ( 4), allocatable :: ntimtmp_prev(:) !  to save previous call frequency
      real     ( 8), allocatable :: cpttemp_prev(:) !  to save previous cp times
      real     ( 8), allocatable :: wcttemp_prev(:) !  to save previous wc times

      if ( allocated(ntimtmp) ) deallocate( ntimtmp )
      if ( allocated(levltmp) ) deallocate( levltmp )
      if ( allocated(contemp) ) deallocate( contemp )
      if ( allocated(cpstemp) ) deallocate( cpstemp )
      if ( allocated(cpttemp) ) deallocate( cpttemp )
      if ( allocated(wcstemp) ) deallocate( wcstemp )
      if ( allocated(wcttemp) ) deallocate( wcttemp )
      if ( allocated(tmsutmp) ) deallocate( tmsutmp )

      if (allocated(ntimtmp_prev)) deallocate(ntimtmp_prev)
      if (allocated(cpttemp_prev)) deallocate(cpttemp_prev)
      if (allocated(wcttemp_prev)) deallocate(wcttemp_prev)

      allocate ( ntimtmp(     nohmax) )
      allocate ( levltmp(     nohmax) )
      allocate ( contemp(99,2,nohmax) )
      allocate ( ncontmp(     nohmax) )
      allocate ( cpstemp(     nohmax) )
      allocate ( cpttemp(     nohmax) )
      allocate ( wcstemp(     nohmax) )
      allocate ( wcttemp(     nohmax) )
      allocate ( tmsutmp(     nohmax) )

      allocate (ntimtmp_prev(nohmax))
      allocate (cpttemp_prev(nohmax))
      allocate (wcttemp_prev(nohmax))

      ntimtmp = ntimcal
      levltmp = level
      contemp = context
      ncontmp = ncontxt
      cpstemp = cpstart
      cpttemp = cptime
      wcstemp = wcstart
      wcttemp = wctime
      tmsutmp = tmsubnm

      ntimtmp_prev = ntimcal_prev
      cpttemp_prev = cptime_prev
      wcttemp_prev = wctime_prev

      deallocate ( ntimcal )
      deallocate ( level   )
      deallocate ( context )
      deallocate ( ncontxt )
      deallocate ( cpstart )
      deallocate ( cptime  )
      deallocate ( wcstart )
      deallocate ( wctime  )
      deallocate ( tmsubnm )

      deallocate (ntimcal_prev)
      deallocate (cptime_prev)
      deallocate (wctime_prev)

      nohmax  = nohmax + 100
      allocate ( ntimcal(     nohmax) )
      allocate ( level  (     nohmax) )
      allocate ( context(99,2,nohmax) )
      allocate ( ncontxt(     nohmax) )
      allocate ( cpstart(     nohmax) )
      allocate ( cptime (     nohmax) )
      allocate ( wcstart(     nohmax) )
      allocate ( wctime (     nohmax) )
      allocate ( tmsubnm(     nohmax) )

      allocate (ntimcal_prev(nohmax))
      allocate (cptime_prev(nohmax))
      allocate (wctime_prev(nohmax))

      ntimcal(    1:nohandl) = ntimtmp
      level  (    1:nohandl) = levltmp
      context(:,:,1:nohandl) = contemp
      ncontxt(    1:nohandl) = ncontmp
      cpstart(    1:nohandl) = cpstemp
      cptime (    1:nohandl) = cpttemp
      wcstart(    1:nohandl) = wcstemp
      wctime (    1:nohandl) = wcttemp
      tmsubnm(    1:nohandl) = tmsutmp

      ntimcal_prev(1:nohandl) = ntimtmp_prev
      cptime_prev(1:nohandl) = cpttemp_prev
      wctime_prev(1:nohandl) = wcttemp_prev

      deallocate ( ntimtmp )
      deallocate ( levltmp )
      deallocate ( contemp )
      deallocate ( ncontmp )
      deallocate ( cpstemp )
      deallocate ( cpttemp )
      deallocate ( wcstemp )
      deallocate ( wcttemp )
      deallocate ( tmsutmp )

      deallocate (ntimtmp_prev)
      deallocate (cpttemp_prev)
      deallocate (wcttemp_prev)
end subroutine timinc

!***************

subroutine timstrt ( subrou, ihandl )

      character*(*), intent(in   ) :: subrou    !  name of (part of) subroutine to monitor
      integer(4)   , intent(inout) :: ihandl    !  handle of the section
      integer(4)                      handle    !  handle of the timer
      integer(4)                      i         !  loop counter
      integer(4)                      ival(8)
      real   (4)                      time

      if (.not. timon) then
         return
      endif
      
      handle = 0
      if ( ihandl == 0 ) then                              !  first time that the timer is called for
         noshndl = noshndl + 1                               !  this handle
         ihandl  = noshndl
      else                                                   !  find its occurence in the call trees
         do i = 1, ncontxt(ihandl)
            if ( context(i,1,ihandl) == prevhnd ) then
               handle = context(i,2,ihandl)
               exit
            endif
         enddo
      endif
      if ( handle == 0 ) then                              !  this is new call tree entry
         i               = ncontxt(ihandl) + 1               !  increase context counter for this ihandl
         ncontxt(ihandl) = i
         context(i,1,ihandl) = prevhnd                       !  save unique timer handle of the caller,
         if ( nohandl == nohmax ) call timinc ( )            !  to allow to find the calling context
         nohandl = nohandl + 1
         handle  = nohandl                                   !  make a new timer handle
         context(i,2,ihandl) = handle                        !  save this handle for this context
         tmsubnm(handle) = subrou                            !  save the ID of this timer
         ntimcal(handle) = 0                                 !  zero the accumulators
         cptime (handle) = 0.0d00
         wctime (handle) = 0.0d00
      endif
      dlevel  = dlevel + 1                                   !  level is 1 deeper than previous level
      level  (handle) = dlevel                               !  levels are only used to indent the
      maxlvl = max (maxlvl,dlevel)                           !  reported output
      prevhnd = handle                                       !  now this timer may become the caller

      call system_clock  ( count, rate )
      call cpu_time      (          time )                   !  this is straight forward timing
      cpstart(handle) = time
      call date_and_time ( values = ival )
      wcstart(handle) = real( count, 8 ) / real( rate, 8 )
      ntimcal(handle) = ntimcal(handle) + 1
end subroutine timstrt

!***************

subroutine timstop ( ihandl )

      integer(4)   , intent(in   ) :: ihandl    !  handle of the subroutine
      integer(4)                      handle    !  handle of the timer
      integer(4)                      i         !  loop counter
      real   (8)                      stopt
      real   (4)                      time

      if (.not. timon) then
         return
      endif
      
      dlevel = dlevel-1                                      !  we return, decrease the level
      handle = -1
      do i = 1, ncontxt(ihandl)                              !  find the context of the timer handle
         if ( context(i,2,ihandl) == prevhnd ) then          !  (prevhnd) that we close now
            handle  = prevhnd                                !  this is the handle that we close
            prevhnd = context(i,1,ihandl)                    !  we retrieve the calling handle from
            exit                                             !  context
         endif
      enddo

      if ( handle == -1 ) then
         write( *, * ) 'Programming error: unbalanced calls to timstart/timstop'
         write( *, * ) 'Found in the context of handle ', ihandl, ' - no reliable routine name avaiable'
         return
      endif

      call cpu_time      (        time )                   !  this is straight forward timing
      call system_clock  ( count, rate )

      cptime(handle) = cptime(handle) + time  - cpstart(handle)
      stopt = real( count, 8 ) / real( rate, 8 )
      wctime(handle) = wctime(handle) + stopt - wcstart(handle)
end subroutine timstop

!***************

subroutine timdump ( fileunit )

      integer       fileunit
      integer(4)    i, j                       !   loop accross timer handles
      integer(4)    ihandl                     !   loop accross subroutine handles
      integer(4)    k                          !   loop accross contexts
      real   (8)    cpfact, wcfact
      character(62) forchr
      data          forchr / '(i5,i11,x,ES13.6,x,f7.2,x,ES13.6,'  /

      if (.not. timon) then
         return
      endif
      
      write (fileunit, '(a,98(a ))' ) ' nr.     times     cpu time      cpu    wall clock      wc', &
                                    '  level',('     ',i=3,maxlvl),' routine name'
      write (fileunit, '(a,98(i5))' ) '        called    in seconds      %     in seconds       %', &
                                 ( i, i=2,maxlvl)

      cpfact = 100.0d00/cptime(1)
      wcfact = 100.0d00/wctime(1)
      do i = 1, nohandl
         if ( level(i) == -1 ) cycle
         write (forchr(34:), '(i5,''x,f8.2,'',i6,''x,a40)'')' ) (level(i)-1)*5+2,(maxlvl-level(i))*5+1
         write (fileunit, forchr ) i,ntimcal(i),cptime(i),cptime(i)*cpfact,wctime(i),wctime(i)*wcfact, tmsubnm(i)
         if ( i == nohandl ) cycle
         if ( level(i+1) .lt. level(i) ) then            !  we are going up
            do ihandl = 1, noshndl                       !  find the context of this timer
               do k = 1, ncontxt(ihandl)
                  j = context(k,2,ihandl)
                  if ( j == i ) then
                     prevhnd = context(k,1,ihandl)       !  it is actually the calling timer we go for
                     goto 10                             !  because the calling timer need to be closed
                  endif                                  !  since we are going one level up
               enddo
            enddo
   10       do ihandl = k+1, noshndl                     !  now find subsequent timers with the same caller
               do k = 1, ncontxt(ihandl)                 !  and print them
                  j = context(k,2,ihandl)
                  if ( context(k,1,ihandl) == prevhnd .and. level(j) .ne. -1 .and. j .ne. i ) then
                     write (forchr(32:), '(i5,''x,f8.2,'',i6,''x,a40)'')' ) (level(j)-1)*5+2,(maxlvl-level(j))*5+1
                     write (fileunit, forchr ) j,ntimcal(j),cptime(j),cptime(j)*cpfact, wctime(j),wctime(j)*wcfact,tmsubnm(j)
                     level(j) = -1
                  endif
               enddo
            enddo
         endif
         level(i) = -1
      enddo
end subroutine timdump

!***************

subroutine timfinalize  ( )

      implicit none

      if ( allocated ( ntimcal ) ) deallocate ( ntimcal )
      if ( allocated ( level   ) ) deallocate ( level   )
      if ( allocated ( cpstart ) ) deallocate ( cpstart )
      if ( allocated ( cptime  ) ) deallocate ( cptime  )
      if ( allocated ( wcstart ) ) deallocate ( wcstart )
      if ( allocated ( wctime  ) ) deallocate ( wctime  )
      if ( allocated ( tmsubnm ) ) deallocate ( tmsubnm )
      if ( allocated ( ncontxt ) ) deallocate ( ncontxt )
      if ( allocated ( context ) ) deallocate ( context )

      timon = .false.
end subroutine timfinalize

!***************

end module cf_timers
