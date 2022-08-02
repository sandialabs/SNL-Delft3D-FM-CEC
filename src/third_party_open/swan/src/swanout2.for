!
!     SWAN/OUTPUT       file 2 of 2
!
!  Contents of this file:
!     SWBLOK
!     SBLKPT
!     SWBLKP
!     SRAWPT
!     SWTABP
!     SUHEAD
!     SWSPEC
!     SWCMSP
!     SWRMAT
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWBLOK ( RTYPE, OQI , OQR , IVTYP, FAC, PSNAME,          41.40
     &                    MXK  , MYK , IRQ , VOQR , VOQ        )          40.51 40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3, ONLY: NSTATM                                           41.62
      USE SWCOMM4, ONLY: KSPHER                                           41.62
      USE OUTP_DATA                                                       40.13
      USE swn_outnc                                                       41.40
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.81: Annette Kieftenburg
!     34.01: Jeroen Adema
!     40.03: Nico Booij
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     41.62: Andre van der Westhuysen
!
!  1. UPDATE
!
!     30.81, Jan. 99: Replaced variable FROM by FROM_ (because FROM is
!                     a reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Nov. 99: NVAR in write statement replaced by OREQ(18)
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.30, May  03: extension to write block output to Matlab files
!     40.31, Jul. 03: small correction w.r.t. length of OVSNAM in
!                     call SWRMAT
!     40.31, Dec. 03: removing POOL construction
!     40.41, Jun. 04: some improvements with respect to MATLAB
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.62, Nov. 15: included wave partitioning output (raw partition file)
!
!  2. PURPOSE
!
!       Preparing output in the form of a block that is printed by
!         subroutine SBLKPT
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       RTYPE   ch*4   input    type of output request:
!                               'BLKP' for output on paper,
!                               'BLKD' and 'BLKL' for output to datafile
!       PSNAME  ch*8   input    name of outpu frame
!       MXK     int    input    number of grid points in x-direction
!       MYK     int    input    number of grid points in y-direction
!       VOQR
!       VOQ
!
!  5. SUBROUTINES CALLING
!
!       SWOUTP (SWAN/OUTP)
!
!  6. SUBROUTINES USED
!
!       SBLKPT, SCUNIT, SFLFUN (all SWAN/OUTP), TABHED,
!       MSGERR, COPYCH, FOR
!       SWRMAT, TXPBLA
!
      LOGICAL STPNOW                                                      34.01
!
!  7. ERROR MESSAGES
!
!       If the point set is not of the type frame, an error message
!       is printed and control returns to subroutine OUTPUT
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       If output is on paper then
!           Call TABHED to print heading
!       Else
!           Call FOR to open file
!       ----------------------------------------------------------------
!       For each required variable do
!           Determine type of variable and factor of multiplication
!           Call SBLKPT to write block output to printer or datafile
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!                                                                             30.81
      CHARACTER (LEN=8) :: PSNAME       ! name of output locations        40.13
      CHARACTER (LEN=4) :: RTYPE        ! output type                     40.13
      INTEGER   VOQR(*), IPD
      INTEGER   OQI(4), IVTYP(OQI(3))                                     40.31
      REAL*8    OQR(2)                                                    41.40 40.31
      REAL      VOQ(MXK*MYK,*), FAC(OQI(3))                               41.40 40.31
      INTEGER IF, IL                                                      40.41 40.30
      INTEGER, SAVE :: IREC(MAX_OUTP_REQ)=0                               40.51 40.41
      LOGICAL, SAVE :: MATLAB=.FALSE.                                     40.41 40.30
      LOGICAL, SAVE :: NCF   =.FALSE.                                     41.40
      LOGICAL       :: EXIST = .FALSE.                                    41.43
      LOGICAL, SAVE :: RAWPRT=.FALSE.                                     41.62
      CHARACTER (LEN=20) :: CTIM                                          40.41
      CHARACTER (LEN=30) :: NAMVAR                                        40.41
      CHARACTER (LEN=80) :: HTXT(3)                                       41.62

      INTEGER, SAVE :: IENT=0                                             40.13
      IF (LTRACE) CALL STRACE (IENT,'SWBLOK')
!
!     **** obtain destination and number of variables from array OUTR ***
      NREF = OQI(1)                                                       40.31
      IF (RTYPE .EQ. 'BLKP') THEN
!       printer type output with header
        IPD   = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWAN', PRINTF)                  30.20
      ELSE IF (RTYPE .EQ. 'BLKD') THEN
!       output to datafile without header
        IPD = 2
      ELSE
        IPD = 3
      ENDIF
      IF (ITEST.GE.90) WRITE (PRTEST, 21)  RTYPE,NREF, OQI(3)             40.31 40.03
  21  FORMAT (' Test SWBLOK: RTYPE NREF NVAR ',A4,2(1X,I6))
      FILENM = OUTP_FILES(OQI(2))                                         40.31 40.13
      MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                          40.41 40.30
     &         INDEX (FILENM, '.mat' ).NE.0                               40.41 40.30
      NCF    = INDEX( FILENM, '.NC'  ).NE.0 .OR.                          41.40
     &         INDEX (FILENM, '.nc'  ).NE.0                               41.40
      RAWPRT = INDEX( FILENM, '.RAW' ).NE.0 .OR.                          41.62
     &         INDEX (FILENM, '.raw' ).NE.0                               41.62
!NNCF      IF (NREF.EQ.0) THEN
      IF (.NOT.NCF .AND. NREF.EQ.0) THEN                                  41.40
        IOSTAT = -1                                                       20.75
        CALL FOR (NREF, FILENM, 'UF', IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
        OQI(1) = NREF                                                     40.31 30.00
        OUTP_FILES(OQI(2)) = FILENM                                       40.41
        IF (MATLAB) THEN                                                  40.30
           CLOSE(NREF)                                                    40.30
           OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',               40.30
     &          STATUS='REPLACE',
     &          ACCESS='DIRECT', RECL=1)                                  40.30
           IREC(IRQ) = 1                                                  40.51
        END IF
        IF (RAWPRT.AND.IPD.EQ.1) THEN                                     41.62
           IF (NSTATM.EQ.1) THEN
              WRITE (HTXT(1),'(a)') ' yyyymmdd hhmmss'
           ELSE
              WRITE (HTXT(1),'(a)') ''
           ENDIF
           IF (KSPHER.EQ.0) THEN
              WRITE (HTXT(2),'(a)') '         x             y'
           ELSE
              WRITE (HTXT(2),'(a)') '       lat           lon'
           ENDIF
           WRITE (HTXT(3),'(a)')
     &              '       name       nprt depth uabs  udir cabs  cdir'
!
           WRITE(NREF,'(A26)') 'SWAN PARTITIONED DATA FILE'
           WRITE(NREF,'(A16,A24,A50)')  TRIM(HTXT(1)), TRIM(HTXT(2)),
     &                                  TRIM(HTXT(3))
           WRITE(NREF,'(A31,A20)') '        hs     tp     lp       ',
     &                             'theta     sp      wf'
        END IF
      ELSE IF (NCF .AND. NCOFFSET(IRQ).EQ.0) THEN                         41.40
        ! reserve free unit number
        IOSTAT = -1
        INQUIRE(FILE=FILENM, EXIST=EXIST)                                 41.43
        CALL FOR (NREF, FILENM, 'UU', IOSTAT)
        IF (STPNOW()) RETURN
        IF (.NOT.EXIST) CLOSE(NREF, STATUS='DELETE')                      41.43
        OQI(1) = NREF
        CALL swn_outnc_openblockfile(FILENM, MYK, MXK,
     &                               OVLNAM, OQI, OQR, IVTYP, IRQ,
     &                               VOQ(:,VOQR(1)),VOQ(:,VOQR(2)))
      ENDIF
      IDLA = OQI(4)                                                       40.31 30.00
      NVAR = OQI(3)                                                       40.31 30.00
!
      IF (ITEST.GE.90) WRITE (PRTEST, 22)  NREF, FILENM
  22  FORMAT (' Test SWBLOK: NREF FILENM  ', I6, A40)
!
      CTIM = CHTIME                                                       40.41
      CALL TXPBLA(CTIM,IF,IL)                                             40.41
      CTIM(9:9)='_'                                                       40.41
!
      IF (RAWPRT) THEN                                                    41.62
!        generate a dump of the raw partition data
         CALL SRAWPT ( NREF, VOQR, VOQ, MXK, MYK )
         GOTO 900
      END IF
!
      DO JVAR = 1, NVAR
        IVTYPE = IVTYP(JVAR)                                              40.31 30.00
        DFAC   = FAC(JVAR)                                                40.31 30.00
!
        IF (IPD.EQ.1) THEN
          IF (DFAC.LE.0.) THEN
!           determine default factor for print output
            IF (OVHEXP(IVTYPE) .LT. 0.5E10) THEN
              IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13                 30.20
            ELSE
              IF (OVSVTY(IVTYPE).EQ.1) THEN
                FMAX = 1.E-8
                DO 10 IP = 1, MXK*MYK
                  FTIP = ABS(VOQ(IP,VOQR(IVTYPE)))
                  FMAX = MAX (FMAX, FTIP)
  10            CONTINUE
              ELSE IF (OVSVTY(IVTYPE).EQ.2) THEN
                FMAX = 1000.
              ELSE IF (OVSVTY(IVTYPE).EQ.3) THEN
                FMAX = 1.E-8
                DO 11 IP = 1, MXK*MYK
                  FTIP1 = ABS(VOQ(IP,VOQR(IVTYPE)))
                  FTIP2 = ABS(VOQ(IP,VOQR(IVTYPE)+1))
                  FMAX  = MAX (FMAX, FTIP1, FTIP2)
  11            CONTINUE
              ENDIF
              IFAC = INT (10.+LOG10(FMAX)) - 13
            ENDIF
            DFAC = 10.**IFAC
          ENDIF
        ELSE
          IF (DFAC.LE.0.) DFAC = 1.
        ENDIF
!
        IF (ITEST .GE. 80) WRITE(PRTEST, 6020) JVAR, IVTYPE, DFAC,
     &    COSCQ, SINCQ
 6020   FORMAT(' Test SWBLOK: jvar, ivtype, dfac, coscq, sincq',
     &          2I10,3E12.5)
!
        IF (OVSVTY(IVTYPE) .LT. 3) THEN
!                      scalar quantities
          IF (MATLAB) THEN                                                40.30
             IF (IL.EQ.1 .OR. IVTYPE.LT.3 .OR. IVTYPE.EQ.52) THEN         40.94 40.41
                NAMVAR = OVSNAM(IVTYPE)                                   40.41
             ELSE                                                         40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_'//CTIM                                        40.41
             END IF                                                       40.41
             CALL SWRMAT( MYK, MXK, NAMVAR,                               40.41
     &                    VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),           40.51
     &                    IDLA, OVEXCV(IVTYPE) )
          ELSE IF (NCF) THEN                                              41.40
             IF (IVTYPE.GT.2.AND.IVTYPE.NE.40) THEN
                CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
     &                                     IRQ, VOQ(1,VOQR(IVTYPE)),
     &                                     OVEXCV(IVTYPE), 1)
             END IF
          ELSE
             CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &       MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,VOQR(IVTYPE)))
          END IF
        ELSE
!                     vectorial quantities
          IF (MATLAB) THEN                                                40.30
             IF (IL.EQ.1) THEN                                            40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_x'                                             40.41
             ELSE                                                         40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_x_'//CTIM                                      40.41
             END IF                                                       40.41
             CALL SWRMAT( MYK, MXK, NAMVAR,                               40.41
     &                 VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),              40.51
     &                 IDLA, OVEXCV(IVTYPE))
             IF (IL.EQ.1) THEN                                            40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_y'                                             40.41
             ELSE                                                         40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_y_'//CTIM                                      40.41
             END IF                                                       40.41
             CALL SWRMAT( MYK, MXK, NAMVAR,                               40.41
     &                 VOQ(1,VOQR(IVTYPE)+1), NREF, IREC(IRQ),            40.51
     &                 IDLA, OVEXCV(IVTYPE))
          ELSE IF (NCF) THEN                                              41.40
             IF ( IVTYPE.GT.3 ) THEN
                CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
     &                                    IRQ, VOQ(1,VOQR(IVTYPE)),
     &                                    OVEXCV(IVTYPE), 1)
                CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
     &                                    IRQ, VOQ(1,VOQR(IVTYPE)+1),
     &                                    OVEXCV(IVTYPE), 2)
             END IF
          ELSE
             CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &       MXK, MYK, IDLA, OVLNAM(IVTYPE)//'X-comp',
     &       VOQ(1,VOQR(IVTYPE)))
             CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &       MXK, MYK, IDLA, OVLNAM(IVTYPE)//'Y-comp',
     &       VOQ(1,VOQR(IVTYPE)+1))
          END IF
        ENDIF
!
      END DO
  900 CONTINUE
      IF ( NCF ) CALL swn_outnc_close_on_end(OQI(1), IRQ)                 41.40
      IF (IPD.EQ.1 .AND. NREF.EQ.PRINTF) WRITE (PRINTF, 6030)
 6030 FORMAT (///)
!
      RETURN
! * end of subroutine SWBLOK *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SBLKPT (IPD, NREF, DFAC, PSNAME, QUNIT,
     &                    MXK, MYK, IDLA, STRING, OQVALS)
!                                                                      *
!***********************************************************************

      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE OUTP_DATA                                                       40.13
      USE TIMECOMM                                                        40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     00.00, Mar. 87: subroutine heading added, some variable names
!                     line numbers changed, layout modified
!     00.04, Feb. 90: lay-out of output changed according to IDLA=1
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.82, Nov. 98: Corrected syntax format statement
!     40.13, July 01: variable formats introduced, using module OUTP_DATA
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Writing the block output either on paper or to datafile
!
!  3. Method
!
!     ---
!
!  4. PARAMETERLIST
!
!     IPD     INT    input    switch for printing on paper (IPD=1)
!                             or writing to datafile (IPD = 2 or 3)
!     NREF    INT    input    unit reference number of output file
!     DFAC    REAL   input    multiplication factor of block output
!     IVTYPE  INT    input    type of the output quantity
!                             Note: IVTYPE=0 for Y-component of a
!                             vectorial quantity
!     PSNAME  CH*8   input    name of output point set (frame)
!     QUNIT   CH*6   input    physical unit (dimension) of variable
!     MXK     int    input    number of points in x-direction of frame
!     MYK     int    input    number of points in y-direction of frame
!     IDLA    INT    input    controls lay-out of output (see user manual)
!     STRING  CH*(*) input    description of output variable
!
!  8. Subroutines used
!
!       ---
!
!  9. Subroutines calling
!
!       SWBLOK (SWAN/OUTP)
!
! 10. Error messages
!
!       ---
!
! 11. Remarks
!
!       ---
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       If IPD = 1 (output on paper) then
!           If DFAC < 0 (DFAC not given by the user) then
!               Compute maximum value of output variable
!               Compute multiplication factor DFAC
!           ------------------------------------------------------------
!           Print block heading
!           For each IX of the output frame do
!               Print IX and for every IY the value of the outputvariable
!           ------------------------------------------------------------
!       Else
!           If DFAC < 0 then DFAC = 1.
!           Write output variable line by line to datafile
!       ----------------------------------------------------------------
!
! 13. Source text
!

      CHARACTER (LEN=20) :: WFORM1 = '(A1, 2X, 151(I6))'                  40.13
      CHARACTER (LEN=21) :: WFORM2 = '(1X,I4,1X, 151(F6.0))'              41.41 40.13
      CHARACTER (LEN=20) :: WFORM3 = '(5X, 151(F6.0))'                    40.13

      CHARACTER PSNAME*8, STRING*(*), QUNIT*(*)                           40.00
      REAL      DFAC, OQVALS(*)
      INTEGER   NREF, MXK, MYK, IPD
      LOGICAL   BPRN
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL STRACE (IENT,'SBLKPT')
!
      IF (ITEST.GE.150) WRITE (PRTEST, 10) NREF,IPD,MXK,MYK
  10  FORMAT (' SBLKPT', 4(I6))                                           30.82
!
!
!     divide all output values by the given factor (DFAC)
!
      IF (ABS(DFAC-1.) .GT. 0.001) THEN
        RPDFAC=1./DFAC
        DO 15 IP = 1, MXK*MYK
           OQVALS(IP) = OQVALS(IP)*RPDFAC
   15   CONTINUE
      ENDIF
!
!
!      IFF = VOQR(IVTYPE)
!
      IF (IPD.EQ.1) THEN
!
!       ***** output on paper *****
!
        WRITE (NREF, 20) OUT_COMMENT                                      40.13
        WRITE (NREF, 20) OUT_COMMENT                                      40.13
  20    FORMAT (A)
        WRITE (NREF, 22) OUT_COMMENT, PROJNR, PSNAME, STRING,             40.13
     &                   DFAC, QUNIT                                      40.13
  22    FORMAT (A,' Run:', A4, '  Frame:  ',A8,' **  ',A,', Unit:',       40.13
     &          E12.4, 1X, A)                                             40.13
        IF (NSTATM .GT. 0) THEN
          WRITE (NREF, 24) OUT_COMMENT, CHTIME                            40.13
  24      FORMAT (A,' Time:', A)                                          40.13
        ELSE
          WRITE (NREF, 20) OUT_COMMENT                                    40.13
        ENDIF                                                             40.13
        WRITE (NREF, 20) OUT_COMMENT                                      40.13

        ISP = 151                                                         30.21
        DO  31  IXP1 = 1, MXK, ISP                                        30.72
          IXP2 = IXP1+ISP-1
          IF (IXP2.GT.MXK) IXP2=MXK

          WRITE (WFORM1(15:15), '(I1)') DEC_BLOCK                         40.13
          WRITE (WFORM2(17:17), '(I1)') DEC_BLOCK                         40.13
          WRITE (WFORM3(11:11), '(I1)') DEC_BLOCK                         40.13
          IF (ITEST.GE.80) WRITE (PRTEST, 25) WFORM1, WFORM2, WFORM3      40.13
  25      FORMAT (' SBLKPT Formats: ', A, /, 17X, A, /, 17X, A)           40.13

          WRITE (NREF, 26) OUT_COMMENT                                    40.13
  26      FORMAT (A1,'         X --->')                                   40.13
          WRITE (NREF, 20) OUT_COMMENT                                    40.13
          WRITE (NREF, WFORM1) OUT_COMMENT, (II-1,II=IXP1,IXP2)           40.13
          WRITE (NREF, 99030) OUT_COMMENT                                 40.13
99030     FORMAT (A1, 'Y')                                                40.13

          BPRN = .TRUE.
          DO 30 IYK = MYK, 1, -1
            IP = (IYK-1)*MXK
            IF (BPRN) THEN
              WRITE (NREF, WFORM2) IYK-1,
     &        (OQVALS(IP+IXK), IXK=IXP1,IXP2)
            ELSE
              WRITE (NREF, WFORM3)
     &        (OQVALS(IP+IXK), IXK=IXP1,IXP2)
            ENDIF
!!!            BPRN = .NOT. BPRN
   30     CONTINUE                                                        30.72
   31   CONTINUE                                                          30.72
      ELSE
!
!       ***** output to datafile *****
!
        ISP=6
        IF (IDLA.EQ.4) THEN
          WRITE (NREF, FLT_BLOCK) (OQVALS(IP), IP=1, MXK*MYK)             40.13
        ELSE
          DO 50 IYK = 1, MYK                                              13/FEB
            IF (IDLA.EQ.3) THEN
              IP = (IYK-1)*MXK
            ELSE
              IP = (MYK-IYK)*MXK
            ENDIF
            WRITE (NREF, FLT_BLOCK) (OQVALS(IP+IXK), IXK=1,MXK)           40.13
   50     CONTINUE
        ENDIF
      ENDIF
!
      RETURN
! * end of subroutine SBLKPT *
      END
!****************************************************************
!
      SUBROUTINE SWBLKP ( OQI   , IVTYP, MXK  , MYK, VOQR, VOQ,
     &                    IONOD )                                         40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.51
!PUN      USE SIZES, ONLY: MYPROC
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: further optimization
!
!  2. Purpose
!
!     Write block data to process output file
!
!  3. Method
!
!     Write data without header using the format FLT_BLKP
!
!  4. Argument variables
!
!     IONOD       array indicating in which subdomain output points       40.51
!                 are located                                             40.51
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         array containing output request data
!     VOQ         output variables
!     VOQR        array containing information for output
!
      INTEGER MXK, MYK, OQI(4), IVTYP(OQI(3)), VOQR(*)
      INTEGER IONOD(*)                                                    40.51
      REAL    VOQ(MXK*MYK,*)
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IOSTAT:     status of input/output
!     IP    :     pointer
!     IPROC :     processor number
!     IVTYPE:     type number output variable
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     NREF  :     unit reference number
!     NVAR  :     number of variables
!
      INTEGER IENT, IOSTAT, IP, IPROC,IVTYPE, IXK, IYK, JVAR, NREF, NVAR
!
!  8. Subroutines used
!
!     FOR              General open file routine
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLKP')

      NREF = OQI(1)
      IF (NREF.EQ.0) THEN
         FILENM = OUTP_FILES(OQI(2))
         IOSTAT = -1
         CALL FOR (NREF, FILENM, 'UU', IOSTAT)
         IF (STPNOW()) RETURN
         OQI(1) = NREF
         OUTP_FILES(OQI(2)) = FILENM                                      40.41
      END IF
      NVAR = OQI(3)

      IPROC = INODE
!PUN      IPROC = MYPROC

      DO JVAR = 1, NVAR
         IVTYPE = IVTYP(JVAR)
         DO IYK = 1, MYK                                                  40.51
            IP = (IYK-1)*MXK                                              40.51
            DO IXK = 1, MXK                                               40.51
               IF ( IONOD(IP+IXK).EQ.IPROC )                              40.51
     &            WRITE (NREF) VOQ(IP+IXK,VOQR(IVTYPE))                   40.51
            END DO                                                        40.51
         END DO                                                           40.51
         IF ( OVSVTY(IVTYPE).GE.3 ) THEN
            DO IYK = 1, MYK                                               40.51
               IP = (IYK-1)*MXK                                           40.51
               DO IXK = 1, MXK                                            40.51
                  IF ( IONOD(IP+IXK).EQ.IPROC )                           40.51
     &               WRITE (NREF) VOQ(IP+IXK,VOQR(IVTYPE)+1)              40.51
               END DO                                                     40.51
            END DO                                                        40.51
         END IF
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SRAWPT ( NREF, VOQR, VOQ, MXK, MYK )
!
!****************************************************************
!
      USE SWCOMM1
      USE SWCOMM3, ONLY: BNAUT, NSTATM, PI
      USE SWCOMM4, ONLY: KSPHER
      USE OCPCOMM4
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     41.62: Andre van der Westhuysen
!
!  1. Updates
!
!     41.62, Nov. 15: New subroutine
!
!  2. Purpose
!
!     Generates a dump of the raw partition data at all grid points
!
!  4. Argument variables
!
!     MXK     int    input    number of points in x-direction of frame
!     MYK     int    input    number of points in y-direction of frame
!     NREF    int    input    unit reference number of output file
!     VOQR
!     VOQ
!
      INTEGER NREF, MXK, MYK
      INTEGER VOQR(*)
      REAL    VOQ(MXK*MYK,*)
!
!  6. Local variables
!
!     CABS    :     magnitude of current
!     CDIR    :     direction of current
!     IP      :     pointer
!     IXK     :     counter in x-direction
!     IYK     :     counter in y-direction
!     NPT     :     actual number of partitions
!     UABS    :     magnitude of wind
!     UDIR    :     direction of wind
!
      INTEGER IENT, II, IP, IXK, IYK
      INTEGER NPT
      REAL    UABS, UDIR, CABS, CDIR
      REAL    HS, TP, DIR, XP, YP, DEP, DSPR, WL
!
      REAL    HSPT(10), TPPT(10), WLPT(10), DIRPT(10),
     &        DSPT(10), WFPT(10), STPT(10)
!
      REAL    DEGCNV
!
!  9. Subroutines calling
!
!     SWBLOK
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SRAWPT')
!
      DO IYK = MYK, 1, -1
         DO IXK = 1, MXK
            IP = (IYK-1)*MXK+IXK
!
            HS  = VOQ(IP,VOQR(10))
            TP  = VOQ(IP,VOQR(12))
            DIR = VOQ(IP,VOQR(13))
!
!           --- if not an exception value of wave height, peak period
!               and wave direction, write values to file
!
            IF ( (HS  .NE. OVEXCV(10)) .AND.
     &           (TP  .NE. OVEXCV(12)) .AND.
     &           (DIR .NE. OVEXCV(13))      ) THEN
!
               XP   = VOQ(IP,VOQR(1))
               YP   = VOQ(IP,VOQR(2))
               DEP  = VOQ(IP,VOQR(4))
               DSPR = VOQ(IP,VOQR(16))
               WL   = VOQ(IP,VOQR(17))
!
               NPT = INT(VOQ(IP,VOQR(171)))
!
!              --- compute magnitude and direction of wind and ambient current
!
               UABS = SQRT(VOQ(IP,VOQR(26))**2+VOQ(IP,VOQR(26)+1)**2)
               UDIR = ATAN2(VOQ(IP,VOQR(26)+1),VOQ(IP,VOQR(26)))*180./PI
               IF (.NOT.BNAUT) UDIR = UDIR + ALCQ * 180./PI
               IF (UDIR.LT.0.) UDIR = UDIR + 360.
               UDIR = DEGCNV(UDIR)
!
               CABS = SQRT(VOQ(IP,VOQR(5))**2+VOQ(IP,VOQR(5)+1)**2)
               CDIR = ATAN2(VOQ(IP,VOQR(5)+1),VOQ(IP,VOQR(5)))*180./PI
               IF (CDIR.GT.360.) CDIR = CDIR - 360.
               IF (CDIR.LT.0.  ) CDIR = CDIR + 360.
!
!              store partition parameters per grid point
!
               DO II = 0, 9
                  HSPT (II+1) = VOQ(IP,VOQR(100+II))
                  TPPT (II+1) = VOQ(IP,VOQR(110+II))
                  WLPT (II+1) = VOQ(IP,VOQR(120+II))
                  DIRPT(II+1) = VOQ(IP,VOQR(130+II))
                  DSPT (II+1) = VOQ(IP,VOQR(140+II))
                  WFPT (II+1) = VOQ(IP,VOQR(150+II))
!                  STPT (II+1) = VOQ(IP,VOQR(160+II))
               END DO
!
!              write the partition data
!
               IF (NSTATM.EQ.1) THEN
                  IF (KSPHER.EQ.0) THEN
                     WRITE(NREF,'(A9,X,A6,2F14.4,A14,I3,F7.1,
     &                 F5.1,F6.1,F5.1,F6.1)')
     &                 CHTIME(1:8),CHTIME(10:16),YP,XP,'''grid_point''',
     &                 NPT, DEP, UABS, UDIR, CABS, CDIR
                  ELSE
                     WRITE(NREF,'(A9,X,A6,2F12.6,A14,I3,F7.1,
     &                 F5.1,F6.1,F5.1,F6.1)')
     &                 CHTIME(1:8),CHTIME(10:16),YP,XP,'''grid_point''',
     &                 NPT, DEP, UABS, UDIR, CABS, CDIR
                  ENDIF
               ELSE
                  IF (KSPHER.EQ.0) THEN
                     WRITE(NREF,'(16X,2F14.4,A14,I3,F7.1,
     &                 F5.1,F6.1,F5.1,F6.1)')
     &                 YP,XP,'''grid_point''',
     &                 NPT, DEP, UABS, UDIR, CABS, CDIR
                  ELSE
                     WRITE(NREF,'(16X,2F12.6,A14,I3,F7.1,
     &                 F5.1,F6.1,F5.1,F6.1)')
     &                 YP,XP,'''grid_point''',
     &                 NPT, DEP, UABS, UDIR, CABS, CDIR
                  ENDIF
               ENDIF
!
               WRITE(NREF,'(I3,F8.2,F8.2,F8.2,F9.2,F9.2,F7.2)')
     &                 0, HS, TP, WL, DIR, DSPR, 999.99
!
               DO II = 1, NPT
                  WRITE(NREF,'(I3,F8.2,F8.2,F8.2,F9.2,F9.2,F7.2)')
     &                 II       , HSPT(II), TPPT(II), WLPT(II),
     &                 DIRPT(II), DSPT(II), WFPT(II)
               END DO
!
            END IF
!
         END DO
      END DO
!
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWTABP (RTYPE , OQI  , OQR , IVTYP, PSNAME, MIP, VOQR,
!NNCF      SUBROUTINE SWTABP (RTYPE , OQI  , IVTYP, PSNAME, MIP, VOQR,         40.31
     &                   VOQ, IONOD)                                      40.51
!                                                                      *
!***********************************************************************

      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.13
      USE TIMECOMM                                                        40.41
      USE M_PARALL                                                        40.51
      USE swn_outnc, only: swn_outnc_openblockfile,
     &                     swn_outnc_appendblock,
     &                     swn_outnc_close_on_end
!PUN      USE SIZES, ONLY: MYPROC, MNPROC
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.62: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.80: Nico Booij
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     32.01: Roeland Ris & Cor van der Schelde
!     34.01: Jeroen Adema
!     40.00: Nico Booij (Non-stationary boundary conditions)
!     40.03, 40.13: Nico Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!
!  1. Updates
!
!     30.50, Sep. 96: option TABI (indexed file) added
!     30.62, Jul. 97: corrected initialisation of table output
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     32.01, Jan. 98: Extended initialisation of NUMDEC for SETUP
!     30.80, Apr. 98: number of decimals for setup from 2 to 3
!     40.00, June 98: severely revised
!     30.82, Oct. 98: Header information is now also printed in PRINT file
!     30.81, Jan. 99: Replaced variable FROM by FROM_ (because FROM is
!                     a reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Mar. 00: number of decimals (NUMDEC) is made larger
!     40.13, Jan. 01: program version now written into table heading
!            Mar. 01: XOFFS and YOFFS were incorrectly added to coordinates
!                     (they are already included in VOQ values)
!     40.13, July 01: variable formats introduced, using module OUTP_DATA
!                     comment sign in front of heading lines
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.31, Dec. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: further optimization
!
!  2. Purpose
!
!     Printing of output in the form of a table for any type of
!     output point set
!
!  3. Method
!
!     A table is made in which for each point the required output
!     variables are printed in the order given by the user. If more
!     variables are required than one line can contain, writing is
!     continued on the next line before output for the next point is
!     started.
!
!  4. Argument variables
!
!     PSNAME
! i   RTYPE : Type of output request
!             ='TABD'; Output to datafile (no header information)
!             ='TABI'; Indexed output for table in ArcView format
!             ='TABP'; Output to paper (with header information)
!             ='TABS';
!             ='TABT';
!             ='TABC'; NETCDF output
!
      CHARACTER RTYPE*4, PSNAME*8
!
!     MIP
!     VOQR
!
      INTEGER   MIP, VOQR(*), OQI(4), IVTYP(OQI(3))
      INTEGER   IONOD(*)                                                  40.51
!
!     VOQ
!
      REAL      VOQ(MIP,*)
      REAL*8    OQR(2)
!
!  5. Parameter variables
!
!     MXOUTL
!
      INTEGER MXOUTL
!
      PARAMETER (MXOUTL=720)
!
!  6. Local variables
!
!     NUMDEC
!
      INTEGER NUMDEC
      LOGICAL EXIST
!
!  8. Subroutines used
!
!     SUHEAD (all SWAN/OUTP)
!     FOR
!     TABHED (all Ocean Pack)
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     OUTPUT (SWAN/OUTP)
!
! 10. ERROR MESSAGES
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       If unit ref. number = 0
!       Then Read filename from array IOUTR
!            Call FOR to open datafile
!            If Rtype = 'TABP' or 'TABI'
!            Then Print heading for required table
!       ----------------------------------------------------------------
!       If Rtype = 'TABS'
!       Then write time into file
!       ----------------------------------------------------------------
!       Make Output line blank
!       Make Linkar = 1
!       If Rtype = 'TABI'
!       Then write index into output line
!            update Linkar
!       ----------------------------------------------------------------
!       For every output point do
!           For all output quantities do
!               Get value for array VOQ
!               If output quantity is TIME
!               Then make Format='(A18)'
!                    write time into output line
!                    make lfield = 18
!               Else if Rtype = 'TABD'
!                    Then Format = '(E12.4)'
!                         make lfield = 12
!                    Else Format = '(F13.X)'
!                         make lfield = 13
!                         determine number of decimals and write into
!                         Format
!                    ---------------------------------------------------
!                    Write value into output line according to Format
!               --------------------------------------------------------
!               Make Linkar = Linkar + lfield + 1
!           ------------------------------------------------------------
!           Write Output line to file
!       ----------------------------------------------------------------
!
! 13. Source text
!
      CHARACTER (LEN=15)     :: FSTR                                      40.22 40.13
      CHARACTER (LEN=MXOUTL) :: OUTLIN                                    40.13
      CHARACTER (LEN=8)      :: CRFORM = '(2F14.4)'                       40.03

      INTEGER, SAVE :: IENT=0                                             40.13
      IF (LTRACE) CALL STRACE(IENT,'SWTABP')
!
      NREF = OQI(1)                                                       40.31 30.00
      NVAR = OQI(3)                                                       40.31 30.00
!
!     Header information is printed once for each data file and for
!     each entry in this routine when the table is written to the PRINT file
!
      IF (NREF .EQ. 0 .OR. NREF.EQ.PRINTF) THEN                           30.82
        IF (NREF.EQ.0) THEN                                               30.82
          FILENM = OUTP_FILES(OQI(2))                                     40.31 40.13
          IOSTAT = -1                                                     20.75
          INQUIRE(FILE=FILENM, EXIST=EXIST)
          CALL FOR (NREF, FILENM, 'UF', IOSTAT)
          IF (STPNOW()) RETURN                                            34.01
          OQI(1) = NREF                                                   40.31 30.00
          OUTP_FILES(OQI(2)) = FILENM                                     40.41
          IF ( RTYPE .EQ. 'TABC' ) THEN
            IF (.NOT.EXIST) CLOSE(NREF, STATUS='DELETE')
            CALL swn_outnc_openblockfile(FILENM, 1, MIP,
     &                                   OVLNAM, OQI, OQR, IVTYP,OQI(2),
     &                                   VOQ(:,VOQR(1)),VOQ(:,VOQR(2)))
          ENDIF
        END IF                                                            30.82
        IF (RTYPE .NE. 'TABD') THEN
          OUTLIN = '    '
!
!         write heading into file
!
          IF (RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI') THEN
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
  20        FORMAT (A)
            WRITE (NREF, 43) OUT_COMMENT, PROJNR, PSNAME, TRIM(VERTXT)    40.13
  43        FORMAT (A1, ' Run:', A4,'  Table:',A8, 10X,                   40.13
     &      'SWAN version:', A)                                           40.13
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
!           write (short) names of output quantities
            IF (RTYPE.EQ.'TABI') THEN
              OUTLIN(1:12) = OUT_COMMENT // '           '                 40.13
              LINKAR = 12
            ELSE
              OUTLIN(1:3) = OUT_COMMENT                                   40.13
              LINKAR = 4                                                  40.13
            ENDIF
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (IVTYPE.EQ.40) THEN
                LFIELD = 18
              ELSE
                LFIELD = 13
              ENDIF
              IF (OVSVTY(IVTYPE).LE.2) THEN
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     '//OVSNAM(IVTYPE)//'              '
              ELSE
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '   X-'//OVSNAM(IVTYPE)//'              '
                LINKAR = LINKAR+LFIELD+1
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '   Y-'//OVSNAM(IVTYPE)//'              '
              ENDIF
              LINKAR = LINKAR+LFIELD+1
            ENDDO
            WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
!           write units of output quantities
            OUTLIN = '    '
            IF (RTYPE.EQ.'TABI') THEN
              OUTLIN(1:12) = OUT_COMMENT // '           '                 40.13
              LINKAR = 12
            ELSE
              OUTLIN(1:3) = OUT_COMMENT                                   40.13
              LINKAR = 4                                                  40.13
            ENDIF
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (IVTYPE.EQ.40) THEN
                LFIELD = 18
              ELSE
                LFIELD = 13
              ENDIF
              DO ISTR = LEN(OVUNIT(IVTYPE)), 1, -1
                IF (OVUNIT(IVTYPE)(ISTR:ISTR) .NE. ' ') THEN
                  LSTR = ISTR
                  GOTO 51
                ENDIF
              ENDDO
              LSTR = 1
  51          OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     ['//OVUNIT(IVTYPE)(1:LSTR)//']            '
              IF (OVSVTY(IVTYPE).GT.2) THEN
                LINKAR = LINKAR+LFIELD+1
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     ['//OVUNIT(IVTYPE)(1:LSTR)//']            '
              ENDIF
              LINKAR = LINKAR+LFIELD+1
            ENDDO
            WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
          ELSE IF (RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS') THEN
            WRITE (NREF, 101) 1
 101        FORMAT ('SWAN', I4, T41, 'Swan standard file, version')
            WRITE (NREF, 111) OUT_COMMENT, VERTXT                         40.13
 111        FORMAT (A1, '   Data produced by SWAN version ', A)           40.13
            WRITE (NREF, 113) OUT_COMMENT, PROJID, PROJNR                 40.13
 113        FORMAT (A1, '   Project: ', A, ';  run number: ', A)          40.13
 102        FORMAT (A, T41, A)
 103        FORMAT (I6, T41, A)
            IF (RTYPE.EQ.'TABT') THEN
              WRITE (NREF,102) 'TABLE'
            ELSE
              IF (NSTATM.EQ.1) THEN
                WRITE (NREF, 102) 'TIME', 'time-dependent data'
                WRITE (NREF, 103) ITMOPT, 'time coding option'
              ENDIF
              IF (KSPHER.EQ.0) THEN                                       33.09
                WRITE (NREF, 102) 'LOCATIONS', 'locations in x-y-space'
                CRFORM = '(2F14.4)'                                       40.03
              ELSE                                                        33.09
                WRITE (NREF, 102) 'LONLAT',
     &                    'locations in spherical coordinates'            33.09
                CRFORM = '(2F12.6)'                                       40.13
              ENDIF                                                       33.09
              WRITE (NREF, 103) MIP, 'number of locations'
              DO 110 IP = 1, MIP
                WRITE (NREF, FMT=CRFORM)                                  40.03
     &                DBLE(VOQ(IP,VOQR(1))), DBLE(VOQ(IP,VOQR(2)))        40.13
 110          CONTINUE
            ENDIF
            WRITE (NREF, 102) 'QUANT', 'description of quantities'
            NKOLS = NVAR
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (OVSVTY(IVTYPE).GT.2) NKOLS = NKOLS + 1
            ENDDO
            WRITE (NREF, 103) NKOLS, 'number of quantities in table'
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (OVSVTY(IVTYPE).LE.2) THEN
                WRITE (NREF, 102) OVSNAM(IVTYPE), OVLNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
 104            FORMAT (E14.4, T41, A)
              ELSE
                WRITE (NREF, 102) 'X-'//OVSNAM(IVTYPE), OVLNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
                WRITE (NREF, 102) 'Y-'//OVSNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!
!     ***** printing of the table *****
!
      IF (RTYPE.EQ.'TABC') THEN
        DO JVAR = 1, NVAR
          IVTYPE = IVTYP(JVAR)                                            40.31
          IF (IVTYPE.GT.3.AND.IVTYPE.NE.40) THEN
            CALL swn_outnc_appendblock(1, MIP, IVTYPE, OQI(1),
     &                                 OQI(2), VOQ(1,VOQR(IVTYPE)),
     &                                 OVEXCV(IVTYPE), 1)
            IF ( OVSVTY(IVTYPE).EQ.3 ) THEN
              CALL swn_outnc_appendblock(1, MIP, IVTYPE, OQI(1),
     &                                   OQI(2), VOQ(1,VOQR(IVTYPE)+1),
     &                                   OVEXCV(IVTYPE), 2)
            ENDIF
          ENDIF
        ENDDO
        IF (NVAR > 0) CALL swn_outnc_close_on_end(OQI(1), OQI(2))
        RETURN
      ENDIF
      IF (RTYPE.EQ.'TABS') THEN
        IF (NSTATM.EQ.1) WRITE (NREF, 102) CHTIME, 'date and time'
      ENDIF
      DO 70 IP = 1, MIP
       IF ( .NOT.PARLL .OR. IONOD(IP).EQ.INODE ) THEN                     40.51
!PUN       IF ( MNPROC.EQ.1 .OR. IONOD(IP).EQ.MYPROC ) THEN                   41.36
        LINKAR = 1
        OUTLIN = '    '
        IF (RTYPE.EQ.'TABI') THEN
!         write point sequence number as first column
          WRITE (OUTLIN(1:8), '(I8)') IP
          LINKAR = 9
          OUTLIN(LINKAR:LINKAR) = ' '
        ENDIF
        DO 60 JVAR = 1, NVAR
          IVTYPE = IVTYP(JVAR)                                            40.31
          IF (IVTYPE.EQ.40) THEN                                          40.00
!           For time 18 characters are needed
            FSTR = '(A18)'
            LFIELD = 18
            OUTLIN(LINKAR:LINKAR+LFIELD-1) = CHTIME
          ELSE
            IF (RTYPE.EQ.'TABD') THEN
              FSTR = FLT_TABLE                                            40.13
              LFIELD = FLD_TABLE                                          40.13
            ELSE
              FSTR = '(F13.X)'
              LFIELD = 13
!             NUMDEC is number of decimals in the table for each output quantity
              NUMDEC = MAX (0, 6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))        40.03
              IF (NUMDEC.GT.9) NUMDEC = 9                                 40.00
              WRITE (FSTR(6:6), '(I1)') NUMDEC                            40.00
            ENDIF
!           write value into OUTLIN
            CALL WRITENUM (OUTLIN(LINKAR:LINKAR+LFIELD-1), FSTR,
     &             VOQ(IP,VOQR(IVTYPE)))
            IF (OVSVTY(IVTYPE).EQ.3) THEN
              LINKAR = LINKAR + LFIELD + 1
!             write second component of a vectorial quantity
              CALL WRITENUM (OUTLIN(LINKAR:LINKAR+LFIELD-1), FSTR,
     &               VOQ(IP,VOQR(IVTYPE)+1))
            ENDIF
          ENDIF
          LINKAR = LINKAR + LFIELD + 1
          OUTLIN(LINKAR-1:LINKAR) = '  '
  60    CONTINUE
        WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
       END IF                                                             40.51
  70  CONTINUE
!
      RETURN
! * end of subroutine SWTABP *
      CONTAINS
      SUBROUTINE WRITENUM(STRING, FSTR, NUMBER)
         ! helper function to convert a number, and avoid printing '****'
         CHARACTER(LEN=*), INTENT(OUT) :: STRING
         REAL            , INTENT(IN)  :: NUMBER
         CHARACTER(LEN=*), INTENT(IN)  :: FSTR

         CHARACTER(LEN=7) :: NWFMT
         CHARACTER(LEN=*), PARAMETER :: FMT1 = "('(E',I1,'.',I1,')')"
         CHARACTER(LEN=*), PARAMETER :: FMT2 = "('(E',I2,'.',I1,')')"
         INTEGER          :: LENSTR, NUMDIG

         WRITE(STRING, FSTR) NUMBER
         IF (INDEX(STRING, '*') > 0) THEN
            LENSTR = MIN(99, LEN(STRING))
            NUMDIG = MAX(1, MIN(9, LENSTR - 7))
            IF (LENSTR > 9) THEN
               WRITE(NWFMT, FMT2) LENSTR, NUMDIG
            ELSE
               WRITE(NWFMT, FMT1) LENSTR, NUMDIG
            ENDIF
            WRITE(STRING, NWFMT) NUMBER
         ENDIF
      END SUBROUTINE WRITENUM
      END
!***********************************************************************
!                                                                      *
      CHARACTER *8 FUNCTION SUHEAD (QUNIT)
!                                                                      *
!***********************************************************************
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  1. UPDATE
!
!        3 JAN 1990 : 0.0, first draft
!
!  2. PURPOSE
!
!       Preparation of unit for the table print output
!       in the form:   [unit]
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       QUNIT   CH*(*) input    unit of the variable to be printed
!                               in the table headings
!
!  5. SUBROUTINES CALLING
!
!       SWTABP (SWAN/OUTP)
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ---
!
! 10. SOURCE TEXT
!
      CHARACTER QUNIT *(*), TEXT1 *6, TEXT2 *8
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SUHEAD')
!
      TEXT2 = '        '
!      L = LEN (QUNIT)
!      IF (L.EQ.0 .OR. L.GT.6) THEN
!        IF (ITEST.GE.10) WRITE(PRINTF,9910) L
! 9910   FORMAT(' ** Error SUHEAD, length of unit in heading =', i2,
!     &        ' out of range')
!      ENDIF
      TEXT1 = QUNIT(1:6)
!     determine the position of the last non-blank character
      DO 10 I = 6,1,-1
        IF (TEXT1(I:I).NE.' ') GOTO 20
 10   CONTINUE
 20   IEND = I
      IF (IEND.EQ.0) THEN
        TEXT1 = '-'
        IEND = 1
      ENDIF
!     shift the unit-string one position to the right
      DO 30 I=1,IEND
        J = I+1
        TEXT2(J:J) = TEXT1(I:I)
 30   CONTINUE
!     enclose the unit by brackets
      TEXT2(1:1) = '['
      TEXT2(IEND+2:IEND+2) = ']'
      SUHEAD = TEXT2
!
 100  RETURN
!     end of subroutine SUHEAD
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWSPEC (RTYPE, OQI, OQR, MIP, VOQR, VOQ, AC2, ACLOC,     41.40 40.31 20.28
     &                   SPCSIG, SPCDIR, DEP2, KGRPNT, CROSS, IONOD)      40.90 40.31 30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.80
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.13
      USE M_PARALL                                                        40.31
      USE swn_outnc, only: swn_outnc_spec                                 41.40
!PUN      USE SIZES, ONLY: MYPROC, MNPROC
!PUN      use SwanGriddata, only: ivertg
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     32.01: Roeland Ris & Cor van der Schelde
!     34.01: Jeroen Adema
!     40.00, 40.03, 40.13: Nico Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.90: Nico Booij
!
!  1. Updates
!
!     20.28         : completely new version
!     20.43         : arguments ECOS and ESIN replaced by SPCDIR
!     32.01, Jan. 98: Introduced nautical convention (project h3268)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Jan. 99: Replaced variable FROM by FROM_ (because FROM is
!                     a reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.00, Aug. 99: new file structure introduced
!     40.03, May  00: correct time coding option written to heading of file
!            Oct. 00: write 'LOCATION' in upper case
!     40.13, Mar. 01: format for writing coordinates different for Cartesian
!                     and spherical coordinates
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.31, Dec. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: further optimization
!     40.90, June 08: argument CROSS added to enable subroutine SWCMSP to
!                     take into account obstacles
!
!  2. Purpose
!
!     Printing of action density spectrum in the form of a table
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IONOD : array indicating in which subdomain output points           40.51
!             are located                                                 40.51
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      LOGICAL, INTENT(IN) :: CROSS(1:4,1:MIP) ! true if obstacle is       40.90
      ! between output point and computational grid point                 40.90
!
!     RTYPE   ch*4   input    type of output request: 'SPEC' for 2-D spectral
!                             output, 'SPE1' for 1-D freq. spectrum
!     MIP     int    input    number of output points in set PSNAME
!     ACLOC   real   local    case SPEC: 2-D spectrum at one output location
!                             case SPE1: 1-D spectra at output locations
!     AK      real   input    wavenumber array at output location
!     UX, UY  real   input    current velocities at output location
!
!  8. Subroutines used
!
!     DEGCNV: Transforms dir. from nautical to cartesian or vice versa    32.01
!     ANGDEG: Transforms degrees to radians                               32.01
!     SWCMSP
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!       ---
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       get NREF from array OREQ  (output requests)
!       if NREF = 0
!       then open output file
!            write heading into the file                                  40.00
!       ----------------------------------------------------------------
!
! 13. Source text
!
      CHARACTER (LEN=*) :: RTYPE                                          40.13 30.81
      CHARACTER (LEN=8) :: CRFORM = '(2F14.4)'                            40.13
      INTEGER       :: VOQR(*), OQI(4), OTYPE, KGRPNT(MXC,MYC)            40.31 40.13
      INTEGER       :: IONOD(*)                                           40.31
      INTEGER, SAVE :: IVERF = 1                                          40.13
      REAL*8    OQR(2)
      REAL      VOQ(MIP,*), AC2(MDC,MSC,MCGRD),
     &          ACLOC(*), DEP2(MCGRD)                                     40.00
      LOGICAL   EQREAL
      LOGICAL, SAVE :: NCF =.FALSE.                                       41.40

      INTEGER, SAVE :: IENT=0                                             40.13
      IF (LTRACE) CALL STRACE(IENT,'SWSPEC')
!
      FILENM = OUTP_FILES(OQI(2))                                         41.40
      NCF = INDEX( FILENM, '.NC' ).NE.0 .OR.                              41.40
     &      INDEX (FILENM, '.nc' ).NE.0                                   41.40
!
      IF ( NCF ) THEN                                                     41.40
         ! When PARALLEL, write intermediate binary files
         CALL swn_outnc_spec ( RTYPE, OQI, OQR, MIP, VOQR,
     &                         VOQ, AC2, SPCSIG, SPCDIR,
     &                         DEP2, KGRPNT, CROSS, IONOD )
         RETURN
      ENDIF
!
      NREF = OQI(1)                                                       40.31 30.00
      IF (INRHOG.EQ.1) THEN
        OFAC = RHO * GRAV                                                 30.20
      ELSE
        OFAC = 1.
      ENDIF
!
      IF (NREF .EQ. 0) THEN
        IERR = -1
        FILENM = OUTP_FILES(OQI(2))                                       40.31 40.13
        IOSTAT = -1                                                       20.75
        CALL FOR (NREF, FILENM, 'UF', IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
        OQI(1) = NREF                                                     40.31 30.00
        OUTP_FILES(OQI(2)) = FILENM                                       40.41
!       IF (IOSTAT.NE.0) WRITE (PRINTF, 6020) FILENM, IOSTAT
!6020   FORMAT (' Open error: ', A36, I6)
!       write heading into the file                                       40.00
!       write keyword SWAN and version number
        WRITE (NREF, 101) IVERF
 101    FORMAT ('SWAN', I4, T41, 'Swan standard spectral file, version')  40.00
        WRITE (NREF, 111) VERTXT                                          40.03
 111    FORMAT ('$   Data produced by SWAN version ', A)                  40.03
        WRITE (NREF, 113) PROJID, PROJNR
 113    FORMAT ('$   Project: ', A, ';  run number: ', A)
        IF (NSTATM.EQ.1) THEN
          WRITE (NREF, 102) 'TIME', 'time-dependent data'
 102      FORMAT (A, T41, A)                                              40.00
          WRITE (NREF, 103) ITMOPT, 'time coding option'                  40.03
 103      FORMAT (I6, T41, A)                                             40.00
        ENDIF
        IF (KSPHER.EQ.0) THEN                                             33.09 NB!
          WRITE (NREF, 102) 'LOCATIONS', 'locations in x-y-space'
          CRFORM = '(2F14.4)'                                             40.13
        ELSE                                                              33.09 NB!
          WRITE (NREF, 102) 'LONLAT',
     &                    'locations in spherical coordinates'            33.09 nb!
          CRFORM = '(2F12.6)'                                             40.13
        ENDIF                                                             33.09 NB!
!PUN        IF ( MNPROC.EQ.1 .OR. .NOT.LCOMPGRD ) THEN
        WRITE (NREF, 103) MIP, 'number of locations'
        DO 110 IP = 1, MIP
          WRITE (NREF, FMT=CRFORM) DBLE(VOQ(IP,VOQR(1))),                 40.13
     &                             DBLE(VOQ(IP,VOQR(2)))                  40.13
 110    CONTINUE
!PUN        ENDIF
        IF (RTYPE(3:3) .EQ. 'R') THEN
          WRITE (NREF, 102) 'RFREQ', 'relative frequencies in Hz'         40.00
        ELSE
          WRITE (NREF, 102) 'AFREQ', 'absolute frequencies in Hz'         40.00
        ENDIF
        WRITE (NREF, 103) MSC, 'number of frequencies'                    40.00
        DO 120 IS = 1, MSC
          WRITE (NREF, 114) SPCSIG(IS)/PI2
 114      FORMAT (F10.4)
 120    CONTINUE
        IF (RTYPE(4:4).EQ.'C') THEN
!         full 2-D spectrum
          IF (BNAUT) THEN
            WRITE (NREF, 102) 'NDIR',
     &                        'spectral nautical directions in degr'      40.00
          ELSE
            WRITE (NREF, 102) 'CDIR',
     &                        'spectral Cartesian directions in degr'     40.00
          ENDIF
          WRITE (NREF, 103) MDC, 'number of directions'
          DO 130 ID = 1, MDC
            IF (BNAUT) THEN
              WRITE (NREF, 124) 180. + DNORTH - SPCDIR(ID,1)*180./PI
            ELSE
              WRITE (NREF, 124) SPCDIR(ID,1)*180./PI
            ENDIF
 124        FORMAT (F10.4)
 130      CONTINUE
          WRITE (NREF, 132) 1
 132      FORMAT ('QUANT', /, I6, T41, 'number of quantities in table')   40.00
          IF (INRHOG.EQ.1) THEN
            WRITE (NREF, 102) 'EnDens',
     &                        'energy densities in J/m2/Hz/degr'          40.00
            WRITE (NREF, 102) 'J/m2/Hz/degr', 'unit'
            WRITE (NREF, 104) OVEXCV(22), 'exception value'               40.00
 104        FORMAT (E14.4, T41, A)
          ELSE
            WRITE (NREF, 102) 'VaDens',
     &                        'variance densities in m2/Hz/degr'          40.00
            WRITE (NREF, 102) 'm2/Hz/degr', 'unit'
            WRITE (NREF, 104) OVEXCV(22), 'exception value'               40.00
          ENDIF
        ELSE
!         1-D spectrum
          WRITE (NREF, 132) 3
          IF (INRHOG.EQ.1) THEN
            WRITE (NREF, 102) 'EnDens',  'energy densities in J/m2/Hz'
            WRITE (NREF, 102) 'J/m2/Hz', 'unit'
            WRITE (NREF, 104) OVEXCV(22), 'exception value'               40.00
          ELSE
            WRITE (NREF, 102) 'VaDens', 'variance densities in m2/Hz'
            WRITE (NREF, 102) 'm2/Hz',  'unit'
            WRITE (NREF, 104) OVEXCV(22), 'exception value'               40.00
          ENDIF
          IF (BNAUT) THEN
            WRITE (NREF, 102) 'NDIR',
     &                        'average nautical direction in degr'
          ELSE
            WRITE (NREF, 102) 'CDIR',
     &                        'average Cartesian direction in degr'
          ENDIF
          WRITE (NREF, 102) OVUNIT(13), 'unit'                            40.00
          WRITE (NREF, 104) OVEXCV(13), 'exception value'                 40.00
          WRITE (NREF, 102) 'DSPRDEGR', OVLNAM(16)                        40.00
          WRITE (NREF, 102) OVUNIT(16), 'unit'                            40.00
          WRITE (NREF, 104) OVEXCV(16), 'exception value'                 40.00
        ENDIF
      ENDIF
!
!     writing of heading is completed, write time if nonstationary
!
      IF (NSTATM.EQ.1) THEN
        WRITE (NREF, 202) CHTIME                                          40.00
 202    FORMAT (A18, T41, 'date and time')
      ENDIF
!
      IF (RTYPE(4:4).EQ.'C') THEN
        IF (RTYPE.EQ.'SPEC') THEN
          OTYPE = -2
        ELSE
          OTYPE = 2
        ENDIF
      ELSE
        IF (RTYPE.EQ.'SPE1') THEN
          OTYPE = -1
        ELSE
          OTYPE = 1
        ENDIF
      ENDIF
!
      DO 290 IP = 1, MIP
        IF (OPTG.NE.5) THEN                                               40.80
           XC = VOQ(IP,VOQR(24))
           YC = VOQ(IP,VOQR(25))
        ELSE                                                              40.80
           IF (.NOT.EQREAL(VOQ(IP,1),OVEXCV(1))) XC = VOQ(IP,1) - XOFFS   40.80
           IF (.NOT.EQREAL(VOQ(IP,2),OVEXCV(2))) YC = VOQ(IP,2) - YOFFS   40.80
        ENDIF                                                             40.80
        DEP = VOQ(IP,VOQR(4))
        IF (DEP.LE.0. .OR. EQREAL(DEP,OVEXCV(4))) THEN
          IF ( .NOT.PARLL .OR. IONOD(IP).EQ.INODE )                       40.51
!PUN          IF ( MNPROC.EQ.1 .OR. IONOD(IP).EQ.MYPROC )                     41.36
     &       WRITE (NREF, 220) 'NODATA'                                   40.00
 220      FORMAT (A6)
          GOTO 290
        ENDIF
        IF ( PARLL .AND. IONOD(IP).NE.INODE ) GOTO 290                    40.51
!PUN        IF ( MNPROC.GT.1 .AND. IONOD(IP).NE.MYPROC ) GOTO 290             41.36
        IF (ICUR.GT.0) THEN
          UX = VOQ(IP,VOQR(5))
          UY = VOQ(IP,VOQR(5)+1)
        ELSE
          UX = 0.
          UY = 0.
        ENDIF
!
        CALL SWCMSP (OTYPE       ,XC         ,YC          ,               40.00
     &               AC2         ,ACLOC      ,SPCSIG      ,
     &               DEP         ,DEP2       ,UX          ,               40.00
     &               UY          ,SPCDIR(1,2),SPCDIR(1,3) ,
     &               OFAC        ,KGRPNT     ,CROSS(1,IP) ,IERR        )  40.90 40.00
!
        IF (IERR.GT.0) THEN
          WRITE (NREF, 220) 'NODATA'
        ELSE
          IF (ABS(OTYPE).EQ.2) THEN
!           write 2d spectrum
            CALL WRSPEC (NREF, ACLOC)
          ELSE
!           write 1d spectrum
            WRITE (NREF, 115) IP
!PUN            IF ( MNPROC.GT.1 .AND. LCOMPGRD ) THEN
!PUN               WRITE (NREF, 115) ivertg(IP)
!PUN            ELSE
!PUN               WRITE (NREF, 115) IP
!PUN            ENDIF
 115        FORMAT ('LOCATION', I6)                                       40.03
            DO IFR = 1, MSC
!             write frequency spectra to file
              WRITE (NREF, 116) (ACLOC(IFR+JJ*MSC), JJ=0,2)               40.00
 116          FORMAT (E12.4, 2F7.1)
            ENDDO
          ENDIF
        ENDIF
 290  CONTINUE
!
      RETURN
! * end of subroutine SWSPEC *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWCMSP (OTYPE     ,XC        ,YC        ,                40.00
     &                   AC2       ,ACLOC     ,SPCSIG    ,
     &                   DEP       ,DEP2      ,UX        ,                40.00
     &                   UY        ,ECOS      ,ESIN      ,
     &                   OFAC      ,KGRPNT    ,CROSS     ,IERR         )  40.90 30.21
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.80
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     40.00: Nico Booij
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.90: Nico Booij
!
!  1. Updates
!
!     20.xx         : New subroutine
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.81, Nov. 98: Adjustment for 1-D case of new boundary conditions
!     30.81, Dec. 98: Argument list KSCIP1 adjusted
!     40.00, Jan. 98: number of output points is always 1
!                     subr produces 3 parameters in case 1D
!                     in 2D cases, loops over ID and ISIGM swapped
!                     interpolation changed: if a corner of the mesh is dry,
!                     exception values are written
!                     argument DEP2 added
!     30.82, Apr. 99: Conversion from m^2/rad/s to m^2/Hz correctly implemented
!     30.82, July 99: Corrected argumentlist KSCIP1
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!     40.90, June 08: interpolation near obstacles improved
!
!  2. Purpose
!
!     Computation of energy density spectrum 1-D or 2-D
!
!  3. Method
!
!     Energy is assumed to be distributed evenly over the interval
!     from Sigma/Frinth to Sigma*Frinth
!     This energy is tranferred to the Omega axis; it is determined
!     how much energy is to be assigned to each interval
!
!     Bilinear interpolation within a mesh of the computational grid
!     to obtain action density in an output location.
!     To transform from relative to absolute frequency, an interval in
!     sigma-space is partitioned, the energy in a submesh is determined
!     and transferred to omega-space where it is added to the energy for
!     a grid step; to obtain energy density this value is divided
!     by the length of the interval.
!     To obtain frequency spectra (1-D) energy density is integrated
!     over theta (spectral direction)
!
!  4. Argument variables
!
!     SPCSIG: input  Relative frequencies in computational domain in      30.72
!                    sigma-space                                          30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!       OTYPE   int    input    type of spectrum wanted: 2 or -2 for 2-D
!                               spectrum, 1 or -1 for 1-D freq. spectrum
!                               positive: relative freq, negative: abs. fr.
!       MPP     int    input    number of output points in set PSNAME
!       XC, YC  real   input    coordinates of output location(s)
!       ACLOC   real   local    |OTYPE|=2: 2-D spectrum at one output location
!                               |OTYPE|=1: 1-D spectra at output locations
!       DEP     real   input    depths at output location
!       UX, UY  real   input    current velocities at output location
!       ECOS  real   input    cosines of spectral directions
!       ESIN  real   input    sines of spectral directions
!       OFAC    real   input    output factor (if INRHOG=1, equal to Rho*Grav)
!
      LOGICAL, INTENT(IN) :: CROSS(1:4) ! true if obstacle is between     40.90
              ! output point and computational grid point                 40.90
!
!  5. SUBROUTINES CALLING
!
!       SWSPEC (SWAN/OUTP)
!
!  6. SUBROUTINES USED
!
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!
      INTEGER   OTYPE           ,KGRPNT(MXC,MYC)                          30.21
      REAL      XC, YC, UX, UY, DEP, DEP2(MCGRD), AC2(MDC,MSC,MCGRD),     40.00
     &          ACLOC(*)        ,ECOS(MDC)         ,ESIN(MDC)             40.00
      REAL      CG(1), K1(1), K2(1), N(1), ND(1), SIG1(1), SIG2(1)        30.82
      LOGICAL EXCPT                                                       40.80
      LOGICAL EQREAL                                                      41.54
      REAL, ALLOCATABLE :: ACL(:,:)                                       40.80
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE(IENT,'SWCMSP')
!
!     make initial value of energy density 0
!
      IF (ABS(OTYPE).EQ.1) THEN
!       1-D spectra
        DO II = 1, 3*MSC                                                  40.00
          ACLOC(II) = 0.
        ENDDO
      ELSE
!       2-D spectra
        DO II = 1, MDC*MSC
          ACLOC(II) = 0.
        ENDDO
      ENDIF
!
!     ***** determine energy densities *****
!
      IF (DEP.LE.0.) GOTO 800
!
!     the action density spectrum is interpolated
!
      ALLOCATE(ACL(MDC,MSC))                                              40.90
      IF (OPTG.EQ.5) THEN                                                 40.90
        CALL SwanInterpolateAc ( ACL, XC, YC, AC2, EXCPT )                40.80
      ELSE                                                                40.90
        IF (KREPTX.EQ.0) THEN                                             40.90
!         non-repeating grid                                              40.90
          IF (XC .LT. -0.01)            GOTO 800                          40.90
          IF (XC .GT. REAL(MXC-1)+0.01) GOTO 800                          40.90
        ENDIF                                                             40.90
        IF (YC .LT. -0.01)            GOTO 800                            40.90
        IF (YC .GT. REAL(MYC-1)+0.01) GOTO 800                            40.90
        CALL SWOINA (XC, YC, AC2, ACL, KGRPNT, DEP2, CROSS(1), EXCPT)     40.90
      ENDIF                                                               40.90
      IF (EXCPT) GOTO 800                                                 40.90
!
      DO 110 ID = 1, MDC
        IF (ICUR.GT.0 .AND. OTYPE.LT.0) THEN
          UDIR = UX * ECOS(ID) + UY * ESIN(ID)
        ENDIF
        DO 100 ISIGM = 1, MSC
!
          ACLL = ACL(ID,ISIGM)                                            40.80
!
!         energy density interpolated in space:
          ECLL = OFAC * ACLL * SPCSIG(ISIGM)                              40.00
!
          IF (ICUR.EQ.0 .OR. OTYPE.GT.0
     &                                    ) THEN
!
!           spectrum as function of relative frequency (SPCSIG)
!
            IF (ABS(OTYPE).EQ.2) THEN
               ACLOC(ID+(ISIGM-1)*MDC) = ECLL                             40.00
            ELSE
!             1-D spectrum of rel. frequency
              ECLL = ECLL * DDIR                                          40.00
              ACLOC(ISIGM) = ACLOC(ISIGM) + ECLL
              ACLOC(ISIGM+  MSC) = ACLOC(ISIGM+  MSC) + ECLL * ECOS(ID)
              ACLOC(ISIGM+2*MSC) = ACLOC(ISIGM+2*MSC) + ECLL * ESIN(ID)
              IF (ITEST.GE.250 .OR. IOUTES .GE. 40) WRITE (PRTEST, 83)
     &           ISIGM, ECLL, (ACLOC(ISIGM+JJ*MSC),JJ=0,2)
 83           FORMAT (' Test SWCMSP ', I6, 4(1X,E12.4))
            ENDIF
          ELSE
!
!           spectrum as function of absolute frequency (OMEGA)
!           WK is wavenumber
!           energy density is assumed constant over the interval from
!           SIG1 to SIG2
!
            SIG1(1)  = SPCSIG(ISIGM) / FRINTH                             30.82
            CALL KSCIP1 (1, SIG1, DEP, K1, CG, N, ND)                     30.82
            OMEG1 = SIG1(1) + K1(1) * UDIR                                30.82
            SIG2(1)  = SPCSIG(ISIGM) * FRINTH                             30.82
            CALL KSCIP1 (1, SIG2, DEP, K2, CG, N, ND)                     30.82
            OMEG2 = SIG2(1) + K2(1) * UDIR                                30.82
            DSIG  = FRINTF * SPCSIG(ISIGM)
!
!           EE is energy density in Omega:
!
            IF ( .NOT.EQREAL(OMEG1,OMEG2) ) THEN                          41.54
               EE = ECLL * DSIG / ABS(OMEG2-OMEG1)
            ELSE                                                          41.54
               EE = ECLL * DSIG                                           41.54
            ENDIF                                                         41.54
            IF (ITEST.GE.250 .OR. IOUTES .GE. 40) WRITE (PRTEST, 86)
     &          ID, ISIGM, SIG1(1), K1(1), OMEG1, SIG2(1), K2(1), OMEG2   30.82
  86        FORMAT (' Test SWCMSP/86 ', 2I6, 8(1X,E12.4))
!
!           assign the energy to omega interval
!
            IF (OMEG1.GT.OMEG2) THEN
!             swap the two values
              RR    = OMEG2
              OMEG2 = OMEG1
              OMEG1 = RR
            ENDIF
            DO 90 IOM = 1, MSC
              OMEGA = SPCSIG(IOM) / FRINTH
              OMEGB = SPCSIG(IOM) * FRINTH
              IF (OMEG1.LT.OMEGB) THEN
                RLOW = MAX (OMEG1,OMEGA)
              ELSE
                GOTO 90
              ENDIF
              IF (OMEG2.GT.OMEGA) THEN
                RUPP = MIN (OMEG2,OMEGB)
              ELSE
                GOTO 90
              ENDIF
              IF (RUPP.LT.RLOW) THEN
                WRITE (PRINTF, 88) ISIGM, IOM, OMEG1, OMEG2,
     &                             OMEGA, OMEGB, RUPP, RLOW
  88            FORMAT (' error SWCMSP:', 2I4, 8(1X,E12.4))
              ELSE
                IF ( .NOT.EQREAL(OMEG1,OMEG2) ) THEN                      41.54
                   DOMEG = RUPP - RLOW                                    41.54
                ELSE                                                      41.54
                   DOMEG = 1.                                             41.54
                ENDIF                                                     41.54
                IF (OTYPE.EQ.-2) THEN
                  ACLOC(ID+(IOM-1)*MDC) =
     &               ACLOC(ID+(IOM-1)*MDC) + EE * DOMEG                   41.54 40.00
                ELSE
                  EADD = EE * DDIR * DOMEG                                41.54 40.00
                  ACLOC(IOM) = ACLOC(IOM) + EADD
                  ACLOC(IOM+  MSC) = ACLOC(IOM+  MSC) + EADD * ECOS(ID)
                  ACLOC(IOM+2*MSC) = ACLOC(IOM+2*MSC) + EADD * ESIN(ID)
                ENDIF
              ENDIF
  90        CONTINUE
          ENDIF
 100    CONTINUE
        IF (OTYPE.EQ.-2 .AND. ICUR.GT.0) THEN                             40.00
          DO 150 IOM = 1, MSC
            DOMEG = FRINTF * SPCSIG(IOM)
            ACLOC(ID+(IOM-1)*MDC) =
     &             ACLOC(ID+(IOM-1)*MDC) / DOMEG                          40.00
            IF (ITEST.GE.250 .OR. IOUTES .GE. 40) WRITE (PRTEST, 135)
     &        ID, IOM, ACLOC(ID+(IOM-1)*MDC), DOMEG                       40.00
 135        FORMAT (' Test SWCMSP ', 2I6, 4(1X,E12.4))
 150      CONTINUE
        ENDIF
 110  CONTINUE
      IF (ABS(OTYPE).EQ.1) THEN
!       1-D spectrum
        IF (ICUR.GT.0 .AND. OTYPE.EQ.-1) THEN
          DO 250 IOM = 1, MSC
            DOMEG = FRINTF * SPCSIG(IOM)
            DO 249 JJ = 0, 2
              ACLOC(IOM+JJ*MSC) = ACLOC(IOM+JJ*MSC) / DOMEG
 249        CONTINUE
 250      CONTINUE
        ENDIF
        DO 270 IOM = 1, MSC
          IF (ACLOC(IOM).GT.1.E-12) THEN                                  40.41
            EX = ACLOC(IOM+MSC) / ACLOC(IOM)
            EY = ACLOC(IOM+2*MSC) / ACLOC(IOM)
            ACLOC(IOM+MSC) = DEGCNV (ATAN2(EY,EX) * 180./PI)
            FF = MIN (1.,SQRT(EX**2+EY**2))
            ACLOC(IOM+2*MSC) = SQRT(2.-2.*FF)*180./PI
!
!           To convert ACLOC from m^2/rad/s to m^2/Hz (1D spectrum)
!
            ACLOC(IOM) = ACLOC(IOM) * PI2                                 40.00
!
          ELSE
!           exception values for VaDens, DIR and DSPR
            ACLOC(IOM)       = OVEXCV(22)                                 40.41
            ACLOC(IOM+MSC)   = OVEXCV(13)
            ACLOC(IOM+2*MSC) = OVEXCV(16)
          ENDIF
 270    CONTINUE
      ENDIF
!
 190  IERR = 0
      GOTO 900
!
!     point is outside grid
 800  IERR = 1
!
 900  CONTINUE
      DEALLOCATE(ACL)                                                     40.90 40.80
      RETURN
! * end of subroutine SWCMSP *
      END
!MatL4!****************************************************************
!MatL4!
!MatL4      SUBROUTINE SWRMAT ( MROWS , NCOLS, MATNAM, RDATA,
!MatL4     &                    IOUTMA, IREC , IDLA  , DUMVAL )
!MatL4!
!MatL4!****************************************************************
!MatL4!
!MatL4      USE OCPCOMM4                                                        40.41
!MatL4!
!MatL4      IMPLICIT NONE
!MatL4!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!MatL4!
!MatL4!  0. Authors
!MatL4!
!MatL4!     40.30: Marcel Zijlema
!MatL4!     40.41: Marcel Zijlema
!MatL4!
!MatL4!  1. Updates
!MatL4!
!MatL4!     40.30, May 03: New subroutine
!MatL4!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!MatL4!
!MatL4!  2. Purpose
!MatL4!
!MatL4!     Writes block output to a binary file in the MAT-format
!MatL4!     to be used in MATLAB
!MatL4!
!MatL4!  3. Method
!MatL4!
!MatL4!     1) The binary file BINFIL must be opened with the following
!MatL4!        statement:
!MatL4!
!MatL4!        OPEN(UNIT=IOUTMA, FILE=BINFIL, FORM='UNFORMATTED',
!MatL4!             ACCESS='DIRECT', RECL=1)
!MatL4!
!MatL4!        Furthermore, initialize record counter to IREC = 1
!MatL4!
!MatL4!     2) Be sure to close the binary file when there are no more
!MatL4!        matrices to be saved
!MatL4!
!MatL4!     3) The matrix may contain signed infinity and/or Not a Numbers.
!MatL4!        According to the IEEE standard, on a 32-bit machine, the real
!MatL4!        format has an 8-bit biased exponent (=actual exponent increased
!MatL4!        by bias=127) and a 23-bit fraction or mantissa. The leftmost
!MatL4!        bit is the sign bit. Let a fraction, biased exponent and sign
!MatL4!        bit be denoted as F, E and S, respectively. The following
!MatL4!        formats adhere to IEEE standard:
!MatL4!
!MatL4!          S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!MatL4!          S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!MatL4!          S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!MatL4!
!MatL4!        Hence, the representation of +Inf equals 2**31 - 2**23. A
!MatL4!        representation of a NaN equals the representation of +Inf
!MatL4!        plus 1.
!MatL4!
!MatL4!        The Cray machine C916 at SARA Information Centre does not
!MatL4!        support the IEEE standard.
!MatL4!
!MatL4!     4) The NaN's or Inf's are indicated by a dummy value as given
!MatL4!        by dumval
!MatL4!
!MatL4!     For more information consult "Appendix - MAT-File Structure"
!MatL4!     of the MATLAB External Data Reference guide (Version 4.2)
!MatL4!
!MatL4!  4. Argument variables
!MatL4!
!MatL4!     DUMVAL      a dummy value meant for indicating NaN
!MatL4!     IDLA        controls lay-out of output (see user manual)
!MatL4!     IOUTMA      unit number of binary MAT-file
!MatL4!     IREC        direct access file record counter
!MatL4!     MATNAM      character array holding the matrix name
!MatL4!     MROWS       a 4-byte integer representing the number of
!MatL4!                 rows in matrix
!MatL4!     NCOLS       a 4-byte integer representing the number of
!MatL4!                 columns in matrix
!MatL4!     RDATA       real array consists of MROWS * NCOLS real
!MatL4!                 elements stored column wise
!MatL4!
!MatL4      INTEGER       MROWS, NCOLS, IDLA, IOUTMA, IREC
!MatL4      REAL          RDATA(*), DUMVAL
!MatL4      CHARACTER*(*) MATNAM
!MatL4!
!MatL4!  5. Parameter variables
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4!  6. Local variables
!MatL4!
!MatL4!     BVAL  :     a byte value
!MatL4!     CHARS :     array to pass character info to MSGERR
!MatL4!     I     :     loop variable
!MatL4!     IENT  :     number of entries
!MatL4!     IF    :     first non-character in string
!MatL4!     IL    :     last non-character in string
!MatL4!     IMAGF :     a 4-byte imaginary flag. Possible values are:
!MatL4!                 0: there is only real data
!MatL4!                 1: the data has also an imaginary part
!MatL4!     IOS   :     auxiliary integer with iostat-number
!MatL4!     ITYPE :     the type flag containing a 4-byte integer whose
!MatL4!                 decimal digits encode storage information.
!MatL4!                 If the integer is represented as ABCD then:
!MatL4!                 "A" indicates the format to write the binary
!MatL4!                 data to a file on the machine. Possible values are:
!MatL4!                   0: Intel based machines (PC 386/486, Pentium)
!MatL4!                   1: Motorola 68000 based machines (Macintosh,
!MatL4!                      HP 9000, SPARC, Apollo, SGI)
!MatL4!                   2: VAX-D format
!MatL4!                   3: VAX-G format
!MatL4!                   4: Cray
!MatL4!                 "B" is always zero
!MatL4!                 "C" indicates which format the data is stored.
!MatL4!                  Possible values are:
!MatL4!                   0: double precision (64 bit) floating point numbers
!MatL4!                   1: single precision (32 bit) floating point numbers
!MatL4!                   2: 32-bit signed integers
!MatL4!                   3: 16-bit signed integers
!MatL4!                   4: 16-bit unsigned integers
!MatL4!                   5: 8-bit unsigned integers
!MatL4!                 "D" indicates the type of data (matrix).
!MatL4!                  Possible values:
!MatL4!                   0: numeric matrix
!MatL4!                   1: textual matrix
!MatL4!                   2: sparse  matrix
!MatL4!     J     :     index
!MatL4!     M     :     loop variable
!MatL4!     MSGSTR:     string to pass message to call MSGERR
!MatL4!     N     :     loop variable
!MatL4!     NAMLEN:     a 4-byte integer representing the number of
!MatL4!                 characters in matrix name plus 1
!MatL4!     NANVAL:     an integer representing Not a Number
!MatL4!
!MatL4      INTEGER I, J, IENT, IF, IL, IOS, M, N
!MatL4      INTEGER BVAL(4), IMAGF, ITYPE, NAMLEN, NANVAL
!MatL4      CHARACTER*20 INTSTR, CHARS
!MatL4      CHARACTER*80 MSGSTR
!MatL4!
!MatL4!  8. Subroutines used
!MatL4!
!MatL4!     INTSTR           Converts integer to string
!MatL4!     MSGERR           Writes error message
!MatL4!     TXPBLA           Removes leading and trailing blanks in string
!MatL4!     SWI2B            Calculates 32-bit representation of an
!MatL4!                      integer number
!MatL4!     SWR2B            Calculates 32-bit representation of a
!MatL4!                      floating-point number
!MatL4!
!MatL4!  9. Subroutines calling
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 10. Error messages
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 11. Remarks
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 12. Structure
!MatL4!
!MatL4!     set Not a Number
!MatL4!
!MatL4!     set some flags
!MatL4!
!MatL4!     write header consisting of ITYPE, MROWS, NCOLS, IMAGF, NAMLEN and
!MatL4!     name of matrix MATNAM
!MatL4!
!MatL4!     write matrix
!MatL4!
!MatL4!     if necessary, give message that error occurred while writing file
!MatL4!
!MatL4! 13. Source text
!MatL4!
!MatL4      SAVE IENT
!MatL4      DATA IENT/0/
!MatL4      IF (LTRACE) CALL STRACE (IENT,'SWRMAT')
!MatL4
!MatL4!     --- set Not a Number
!MatL4
!MatL4      NANVAL = 255 * 2**23 + 1
!MatL4
!MatL4!     --- set some flags
!MatL4
!MatL4      ITYPE = 1010
!MatL4      IMAGF = 0
!MatL4      IOS   = 0
!MatL4
!MatL4!     --- write header consisting of ITYPE, MROWS, NCOLS, IMAGF,
!MatL4!         NAMLEN and name of matrix MATNAM
!MatL4!         the name should be ended by zero-byte terminator
!MatL4
!MatL4      CALL SWI2B ( ITYPE, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( MROWS, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( NCOLS, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( IMAGF, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL TXPBLA(MATNAM,IF,IL)
!MatL4      NAMLEN = IL - IF + 2
!MatL4      CALL SWI2B ( NAMLEN, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      DO I = IF, IL
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MATNAM(I:I)
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4      IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(0)
!MatL4      IREC = IREC + 1
!MatL4
!MatL4!     --- write matrix
!MatL4
!MatL4      DO M = 1, NCOLS
!MatL4         DO N = 1, MROWS
!MatL4            IF ( IDLA.EQ.1 ) THEN
!MatL4               J = (MROWS-N)*NCOLS + M
!MatL4            ELSE
!MatL4               J = (N-1)*NCOLS + M
!MatL4            END IF
!MatL4            IF (RDATA(J).NE.DUMVAL ) THEN
!MatL4               CALL SWR2B ( RDATA(J), BVAL )
!MatL4            ELSE IF (.NOT. DUMVAL.NE.0. ) THEN
!MatL4               CALL SWR2B ( RDATA(J), BVAL )
!MatL4            ELSE
!MatL4               CALL SWI2B ( NANVAL, BVAL )
!MatL4            END IF
!MatL4            DO I = 1, 4
!MatL4              IF (IOS.EQ.0)
!MatL4     &                   WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4              IREC = IREC + 1
!MatL4            END DO
!MatL4         END DO
!MatL4      END DO
!MatL4
!MatL4!     --- if necessary, give message that error occurred while writing file
!MatL4
!MatL4      IF ( IOS.NE.0 ) THEN
!MatL4         CHARS = INTSTR(IOS)
!MatL4         CALL TXPBLA(CHARS,IF,IL)
!MatL4         MSGSTR = 'Error while writing binary MAT-file - '//
!MatL4     &            'IOSTAT number is '//CHARS(IF:IL)
!MatL4         CALL MSGERR ( 4, MSGSTR )
!MatL4         RETURN
!MatL4      END IF
!MatL4
!MatL4      RETURN
!MatL4      END
!****************************************************************
!
      SUBROUTINE SWRMAT ( MROWS , NCOLS, MATNAM, RDATA,
     &                    IOUTMA, IREC , IDLA  , DUMVAL )
!
!****************************************************************
!
      USE OCPCOMM2                                                        41.08
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2020  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     41.08: Pieter Smit
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.08, Aug. 09: adapted to write binary file in the Level 5 MAT-file format
!
!  2. Purpose
!
!     Writes block output to a binary file in the MAT-format
!     to be used in MATLAB
!
!  3. Method
!
!     1) The binary file BINFIL must be opened with the following
!        statement:
!
!        OPEN(UNIT=IOUTMA, FILE=BINFIL, FORM='UNFORMATTED',
!             ACCESS='DIRECT', RECL=4)
!
!        Furthermore, initialize record counter to IREC = 1
!
!     2) Be sure to close the binary file when there are no more
!        matrices to be saved
!
!     3) The matrix may contain signed infinity and/or Not a Numbers.
!        According to the IEEE 754 standard, on a 32-bit machine, the real
!        format has an 8-bit biased exponent (=actual exponent increased
!        by bias=127) and a 23-bit fraction or mantissa. The leftmost
!        bit is the sign bit. Let a fraction, biased exponent and sign
!        bit be denoted as F, E and S, respectively. The following
!        formats adhere to IEEE standard:
!
!          S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!          S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!          S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!
!        Hence, the representation of +Inf equals 2**31 - 2**23. A
!        representation of a NaN equals the representation of +Inf
!        plus 1.
!
!     4) The NaN's or Inf's are indicated by a dummy value as given
!        by dumval
!
!     For more information on the Level 5 MAT-file format consult
!     document "MAT-File Format" of MathWorks
!
!  4. Argument variables
!
!     DUMVAL      a dummy value meant for indicating NaN
!     IDLA        controls lay-out of output (see user manual)
!     IOUTMA      unit number of binary MAT-file
!     IREC        direct access file record counter
!     MATNAM      character array holding the matrix name
!     MROWS       a 4-byte integer representing the number of
!                 rows in matrix
!     NCOLS       a 4-byte integer representing the number of
!                 columns in matrix
!     RDATA       real array consists of MROWS * NCOLS real
!                 elements stored column wise
!
      INTEGER       MROWS, NCOLS, IDLA, IOUTMA, IREC
      REAL          RDATA(*), DUMVAL
      CHARACTER*(*) MATNAM
!
!  5. Parameter variables
!
!     BlockSize   size of matlab data segment
!     DataSize    number of bytes written per write statement
!     HeaderSize  size of the header in bytes
!     mChar       character data
!     mInt32      signed   integer*4
!     mUInt32     unsigned integer*8
!     mSingle     real
!
!     --- standard sizes
!
      INTEGER, PARAMETER :: DataSize   = 4
      INTEGER, PARAMETER :: HeaderSize = 128
      INTEGER, PARAMETER :: BlockSize  = 8
!
!     --- Matlab data types
!
      INTEGER, PARAMETER :: mChar      = 1
      INTEGER, PARAMETER :: mInt32     = 5
      INTEGER, PARAMETER :: mUInt32    = 6
      INTEGER, PARAMETER :: mSingle    = 7
!
!  6. Local variables
!
!     CTMP  :     a temporary character array
!     HEADER:     header of binary MAT-file
!     I     :     loop variable
!     IENT  :     number of entries
!     IOS   :     auxiliary integer with iostat-number
!     IRECS :     size of array including tags and flags
!     J     :     index
!     M     :     loop variable
!     MSGSTR:     string to pass message to call MSGERR
!     N     :     loop variable
!     NAMLEN:     a 4-byte integer representing the number of
!                 characters in matrix name
!     NANVAL:     an integer representing Not a Number
!     NTOT  :     size of data array
!
      INTEGER I, J, IENT, IOS, M, N, NTOT
      INTEGER NAMLEN, NANVAL
      INTEGER, SAVE :: IRECS
      CHARACTER*80 MSGSTR
      CHARACTER(LEN=HeaderSize) HEADER
      CHARACTER(LEN=BlockSize) CTMP

!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     set Not a Number
!     length of name matrix
!     size of data array
!     write header once
!     array name
!     write matrix
!     write the size of the array
!     if necessary, give message that error occurred while writing file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRMAT')

      IOS = 0

!     --- set Not a Number

      NANVAL = 255 * 2**23 + 1

!     --- length of name matrix

      NAMLEN = LEN_TRIM(MATNAM)

!     --- size of data array

      NTOT = MROWS * NCOLS

!     --- descriptive header

      WRITE (HEADER, '(6A)') 'Data produced by SWAN version ',
     &                       TRIM(VERTXT),'; project: ',TRIM(PROJID),
     &                       '; run number: ',PROJNR

!     --- data offset

      HEADER(117:124) = CHAR(ICHAR(' '))

!     --- version

      HEADER(125:126) = CHAR(0) // CHAR(1)

!     --- endian indicator

      WRITE(HEADER(127:128),'(A)') INT(19785,KIND=2)

!     --- write header once

      IF ( IREC.EQ.1 ) THEN
         DO I = 1, HeaderSize/DataSize
            J = DataSize*(I-1) + 1
            IF ( IOS.EQ.0 )
     &         WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) HEADER(J:J+DataSize-1)
            IREC = IREC + 1
         END DO
      END IF

!     --- array tag

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 14                ! Matlab array
      IREC = IREC + 1
      IRECS = IREC
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
      IREC = IREC + 1

!     --- array flags

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 7                 ! single precision array
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
      IREC = IREC + 1

      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
         IREC = IREC + 1
      END IF

!     --- dimensions array

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MROWS
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NCOLS
      IREC = IREC + 1

      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
         IREC = IREC + 1
      END IF

!     --- array name

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mChar
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NAMLEN
      IREC = IREC + 1

      I = 1
      DO
        CTMP(1:8) = CHAR(ICHAR(' '))
        IF ( I.GT.NAMLEN ) THEN
           EXIT
        ELSE IF ( I+BlockSize.LE.NAMLEN ) THEN
           CTMP(1:8) = MATNAM(I:I+BlockSize-1)
        ELSE
           CTMP(1:NAMLEN-I+1) = MATNAM(I:NAMLEN)
        END IF

        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(1:4)
        IREC = IREC + 1
        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(5:8)
        I = I + BlockSize
        IREC = IREC + 1
      END DO

!     --- write matrix

      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) mSingle
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) NTOT*DataSize
      IREC = IREC + 1

      DO M = 1, NCOLS
         DO N = 1, MROWS
            IF ( IDLA.EQ.1 ) THEN
               J = (MROWS-N)*NCOLS + M
            ELSE
               J = (N-1)*NCOLS + M
            END IF
            IF (RDATA(J).NE.DUMVAL ) THEN
               WRITE (IOUTMA,REC=IREC) RDATA(J)
            ELSE IF (.NOT. DUMVAL.NE.0. ) THEN
               WRITE (IOUTMA,REC=IREC) RDATA(J)
            ELSE
               WRITE (IOUTMA,REC=IREC) NANVAL
            END IF
            IREC = IREC + 1
         END DO
      END DO

      IF ( MOD(NTOT,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) 0.
         IREC = IREC + 1
      END IF

!     --- write the size of the array

      IF ( IOS.EQ.0 )
     &   WRITE (IOUTMA,REC=IRECS,IOSTAT=IOS) (IREC-IRECS-1)*DataSize

!     --- if necessary, give message that error occurred while writing file

      IF ( IOS.NE.0 ) THEN
         WRITE (MSGSTR, '(A,I5)')
     &            'Error while writing binary MAT-file - '//
     &            'IOSTAT number is ', IOS
         CALL MSGERR( 4, TRIM(MSGSTR) )
         RETURN
      END IF

      RETURN
      END
