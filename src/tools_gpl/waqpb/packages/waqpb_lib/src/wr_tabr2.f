!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2022.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
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
!  $Id: wr_tabr2.f 140618 2022-01-12 13:12:04Z klapwijk $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/tools_gpl/waqpb/packages/waqpb_lib/src/wr_tabr2.f $

      SUBROUTINE WR_TABR2 ( DEFFDS      ,
     +                      NO_CNSB     , R2_CID      ,
     +                      R2_SID      , LUNREP      ,
     +                      IERROR      )
C
C     Deltares
C
C     CREATED            :  june 1999 by Jan van Beek
C
C     FUNCTION           :  Write TABLE_R2 group to NEFIS file
C
C     FILES              :  NEFIS file assumed opened
C
C     SUBROUTINES CALLED :
C
C     ARGUMENTS
C
C     NAME    TYPE     LENGTH     FUNCT.  DESCRIPTION
C     ----    -----    ------     ------- -----------
C     DEFFDS       INT    2993    I/O     Definition file descriptor
C     DATFDS       INT    999     I/O     Data file descriptor
C     NO_CNSB      INT            I       number of rows in table r2
C     R1_CID       CHA*10 NO_CNSB I       configuration identification
C     R1_SID       CHA*10 NO_CNSB I       substance identification
C     LUNREP       INT    1       I       Unit number report file
C     IERROR       INT    1       O       Error
C
C     IMPLICIT NONE for extra compiler checks
C     SAVE to keep the group definition intact
C
      IMPLICIT NONE
      SAVE
C
C     declaration of arguments
C
      INTEGER       NO_CNSB     , LUNREP     ,
     +              IERROR
      INTEGER       DEFFDS
      CHARACTER*10  R2_CID      (NO_CNSB)
      CHARACTER*10  R2_SID      (NO_CNSB)
C
C     Local variables
C
C     GRPNAM  CHAR*16     1       LOCAL   group name (table)
C     NELEMS  INTEGER     1       LOCAL   number of elements in group (=cell)
C     ELMNMS  CHAR*16  NELEMS     LOCAL   name of elements on file
C     ELMTPS  CHAR*16  NELEMS     LOCAL   type of elements
C     ELMDMS  INTEGER  6,NELEMS   LOCAL   dimension of elements
C     NBYTSG  INTEGER  NELEMS     LOCAL   length of elements (bytes)
C
      INTEGER       NELEMS
      PARAMETER   ( NELEMS = 3 )
C
      INTEGER       I               , IELM
      INTEGER       ELMDMS(2,NELEMS), NBYTSG(NELEMS),
     +              UINDEX(3)
      CHARACTER*16  GRPNAM
      CHARACTER*16  ELMNMS(NELEMS)  , ELMTPS(NELEMS)
      CHARACTER*64  ELMDES(NELEMS)
C
C     External NEFIS Functions
C
      INTEGER   CREDAT
     +         ,DEFCEL
     +         ,DEFELM
     +         ,DEFGRP
     +         ,FLSDAT
     +         ,FLSDEF
     +         ,PUTELS
     +         ,PUTELT
      EXTERNAL  CREDAT
     +         ,DEFCEL
     +         ,DEFELM
     +         ,DEFGRP
     +         ,FLSDAT
     +         ,FLSDEF
     +         ,PUTELS
     +         ,PUTELT
C
C     element names
C
      DATA  GRPNAM  /'TABLE_R2'/
      DATA
     + (ELMNMS(I),ELMTPS(I),NBYTSG(I),ELMDMS(1,I),ELMDMS(2,I),ELMDES(I),
     +  I = 1 , NELEMS)
     +/'NO_CNSB'  ,'INTEGER'  , 4,1,1,'number of rows  in table R2'    ,
     + 'R2_CID'   ,'CHARACTER',10,1,0,'configuration identification'   ,
     + 'R2_SID'   ,'CHARACTER',10,1,0,'substance identification'       /
C
C     Set dimension of table
C
      DO IELM = 2 , NELEMS
         ELMDMS(2,IELM) = NO_CNSB
      ENDDO
C
C     Define elements
C
      WRITE(LUNREP,*) ' WRITING GROUP:',GRPNAM
      DO IELM = 1 , NELEMS
         IERROR = DEFELM (DEFFDS        , ELMNMS(IELM)  ,
     +                    ELMTPS(IELM)  , NBYTSG(IELM)  ,
     +                    ' '           , ' '           ,
     +                    ELMDES(IELM)  , ELMDMS(1,IELM),
     +                    ELMDMS(2,IELM))
         IF ( IERROR .NE. 0 ) THEN
            WRITE(LUNREP,*) 'ERROR defining element:',ELMNMS(IELM)
            WRITE(LUNREP,*) 'ERROR number:',IERROR
            GOTO 900
         ENDIF
      ENDDO
C
C     Define group
C
      IERROR = DEFCEL (DEFFDS, GRPNAM, NELEMS, ELMNMS)
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR defining cell for group',GRPNAM
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      IERROR = DEFGRP (DEFFDS, GRPNAM, GRPNAM, 1, 1, 1)
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR defining group',GRPNAM
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      IERROR = CREDAT (DEFFDS, GRPNAM, GRPNAM)
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR creating data',GRPNAM
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      IERROR = FLSDEF(DEFFDS)
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR flushing definition file'
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      IERROR = FLSDAT(DEFFDS)
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR flushing data file'
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
C
C     Nu het schrijven
C
      UINDEX(1) = 1
      UINDEX(2) = 1
      UINDEX(3) = 1
      WRITE(LUNREP,*) ' WRITING ELEMENT:',ELMNMS(1)
      IERROR = PUTELT (DEFFDS ,
     +                 GRPNAM , ELMNMS(1),
     +                 UINDEX , 1        ,
     +                 NO_CNSB           )
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR writing element',ELMNMS(1)
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      WRITE(LUNREP,*) ' WRITING ELEMENT:',ELMNMS(2)
      IERROR = PUTELS (DEFFDS ,
     +                 GRPNAM , ELMNMS(2),
     +                 UINDEX , 1        ,
     +                 R2_CID            )
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR writing element',ELMNMS(2)
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
      WRITE(LUNREP,*) ' WRITING ELEMENT:',ELMNMS(3)
      IERROR = PUTELS (DEFFDS ,
     +                 GRPNAM , ELMNMS(3),
     +                 UINDEX , 1        ,
     +                 R2_SID            )
      IF ( IERROR .NE. 0 ) THEN
         WRITE(LUNREP,*) 'ERROR writing element',ELMNMS(3)
         WRITE(LUNREP,*) 'ERROR number:',IERROR
         GOTO 900
      ENDIF
C
  900 CONTINUE
      RETURN
C
      END
