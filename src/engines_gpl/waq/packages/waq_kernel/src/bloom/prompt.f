!!  Copyright (C)  Stichting Deltares, 2012-2018.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

!    Date:       13 Dec 1989
!    Time:       08:29
!    Program:    PROMPT.FOR
!    Version:    1.1
!    Programmer: ??????
!    Previous version(s):
!    1.0 -- 13 Dec 1989 -- 07:49 -- Operating System: DOS
!    0.0 -- 12 Dec 1989 -- 10:19 -- Operating System: DOS
!
!  *********************************************************************
!  *SUBROUTINE PROMPT TO REQUEST INPUT FROM THE USER IN AN INTERACTIVE *
!  *         RUN, INDICATING IN WHICH PART OF THE PROGRAM HE IS.       *
!  *********************************************************************
!
! PC program version for Microsoft fortran.
! Use escape sequences to position the cursor after the last prompting
! message.
!
      SUBROUTINE BLPROMPT (INDEX, ARG)
      CHARACTER*80 OUTSTR
      CHARACTER*1 OUTST2 (80), ESCAPE, BLANK
      INTEGER STOS, IESCAP, ARG
      INCLUDE 'ioblck.inc'
      DATA BLANK  /' '/
      DATA IESCAP /27/
!
      ESCAPE = CHAR (IESCAP)
      CALL PROMES (INDEX, ARG, OUTSTR)
      IF (INDEX .LT. 1000) RETURN
      LENOUT = LENST (OUTSTR, 80)
      IRC = STOS (OUTSTR, 1, LENOUT, OUTST2, LENOU2)
      CALL VIDEO (7)
      WRITE (OUUNI,10) (OUTST2(I),I=1,LENOUT),BLANK
10    FORMAT (80A1)
      CALL VIDEO (0)
!
! Position the cursor 2 lines backwards and afterwards shift it right
! LENOUT positions.
!
      WRITE (OUUNI,20) ESCAPE
 20   FORMAT (1X,A1,'[2A')
      IF (LENOUT .LT. 10) THEN
         WRITE (OUUNI,30) ESCAPE, LENOUT
 30      FORMAT (1X,A1,'[',I1,'C')
      ELSE
         WRITE (OUUNI,40) ESCAPE, LENOUT
 40      FORMAT (1X,A1,'[',I2,'C')
      END IF
      RETURN
      END
