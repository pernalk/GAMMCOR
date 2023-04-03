!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
      SUBROUTINE DUMP_RHO2(CREF,XNDXCI,WORK,KFREE,LFREE)
!
!     Dump RHO2 to GammCor program
!
      IMPLICIT NONE

      INTEGER KFREE,LFREE
      REAL*8  XNDXCI(*), WORK(*)
      REAL*8  CREF(NCONF)
C
C KP 23/02/18
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*), DOAPSG
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFVAR : NCONF
#include "dummy.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "priunit.h"

      INTEGER RHOTYP
      INTEGER ISPIN1,ISPIN2
      INTEGER KRHO1, KRHO2, NNRHO2
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
      LOGICAL TRIP_SAV

      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
      TDM    = .FALSE.
      NORHO2 = .FALSE.

      ISPIN1 = 0
      ISPIN2 = 0
      CSFEXP = .NOT.FLAG(27)

      RHOTYP = 2
      KRHO1  = 1
      KFREE  = KRHO1  + N2ASHX
      KRHO2  = KFREE
      NNRHO2 = NASHT*NASHT*NASHT*NASHT
      KFREE  = KFREE  + NNRHO2
      CALL SETVEC(WORK(KRHO1),0.D0,N2ASHX)
      CALL SETVEC(WORK(KRHO2),0.D0,NNRHO2)
      CALL DENSID(LSYM,LSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),
     & WORK(KRHO2),
     & RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     & XNDXCI,WORK,KFREE,LFREE)

      CALL DISK_RHO2(WORK(KRHO2),NASHT,ISPIN1,ISPIN2)

      TRIP_SAV = .TRUE.
      IF (TRIP_SAV) THEN

      ISPIN1 = 1
      ISPIN2 = 1

      CALL SETVEC(WORK(KRHO1),0.D0,N2ASHX)
      CALL SETVEC(WORK(KRHO2),0.D0,NNRHO2)
      CALL DENSID(LSYM,LSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),
     & WORK(KRHO2),
     & RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     & XNDXCI,WORK,KFREE,LFREE)

      CALL DISK_RHO2(WORK(KRHO2),NASHT,ISPIN1,ISPIN2)

C     for Ms = 0 RHO2(1,0) and RHO2(0,1) should be always 0
      ISPIN1 = 0
      ISPIN2 = 1

      CALL SETVEC(WORK(KRHO1),0.D0,N2ASHX)
      CALL SETVEC(WORK(KRHO2),0.D0,NNRHO2)
      CALL DENSID(LSYM,LSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),
     & WORK(KRHO2),
     & RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     & XNDXCI,WORK,KFREE,LFREE)

      CALL DISK_RHO2(WORK(KRHO2),NASHT,ISPIN1,ISPIN2)

      ISPIN1 = 1
      ISPIN2 = 0

      CALL SETVEC(WORK(KRHO1),0.D0,N2ASHX)
      CALL SETVEC(WORK(KRHO2),0.D0,NNRHO2)
      CALL DENSID(LSYM,LSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),
     & WORK(KRHO2),
     & RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     & XNDXCI,WORK,KFREE,LFREE)

      CALL DISK_RHO2(WORK(KRHO2),NASHT,ISPIN1,ISPIN2)

C     1-electron spin-density
      RHOTYP = 1
      ISPIN1 = 1
      ISPIN2 = 0
      NORHO2 = .TRUE.

      CALL SETVEC(WORK(KRHO1),0.D0,N2ASHX)

      CALL DENSID(LSYM,LSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),
     & DUMMY,
     & RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     & XNDXCI,WORK,KFREE,LFREE)

      CALL DISK_RHO1(WORK(KRHO1),NASHT)
c     CALL DISK_RHOS12(WORK(KRHO1),WORK(KRHO2),NASHT)

      END IF

      END

      SUBROUTINE DISK_RHO2(RHO2,NASHT,ISPIN1,ISPIN2)
#include "priunit.h"
      INTEGER NASHT
      INTEGER ISPIN1,ISPIN2
      INTEGER I,J,K,L
      REAL*8  RHO2(NASHT,NASHT,NASHT,NASHT)
      CHARACTER*12 RDMFILE

      IF(ISPIN1.EQ.0 .AND. ISPIN2.EQ.0) RDMFILE = 'rdm2.dat'
      IF(ISPIN1.EQ.1 .AND. ISPIN2.EQ.1) RDMFILE = 'rdms2.dat'
      IF(ISPIN1.EQ.0 .AND. ISPIN2.EQ.1) RDMFILE = 'rdms201.dat'
      IF(ISPIN1.EQ.1 .AND. ISPIN2.EQ.0) RDMFILE = 'rdms210.dat'

      write(lupri,*) 'rdmfile',RDMFILE
C
C KP 23/02/18
C
      OPEN(938,FILE=RDMFILE)
      DO 111 I=1,NASHT
      DO 111 J=1,NASHT
      DO 111 K=1,NASHT
      DO 111 L=1,NASHT
  111 IF(ABS(RHO2(I,J,K,L)).GT.1.D-8) WRITE(938,'(4I4,F19.12)')
     $ I,J,K,L,RHO2(I,J,K,L)
      CLOSE(938)

      END

      SUBROUTINE DISK_RHO1(RHO1,NASHT)

#include "priunit.h"
      INTEGER NASHT
      INTEGER I,J
      REAL*8  RHO1(NASHT,NASHT)
      CHARACTER*9 RDMFILE

      RDMFILE = 'rdms1.dat'

      write(lupri,*) 'nasht',NASHT
      write(lupri,*) 'rho1',RHO1

      OPEN(938,FILE=RDMFILE)
      DO J=1,NASHT
      DO I=1,NASHT
         WRITE(938,'(2I4,F19.12)') I,J,RHO1(I,J)
      END DO
      END DO
      CLOSE(938)

      END

      SUBROUTINE DISK_RHOS12(RHO1,RHO2,NASHT)
C     dump spin RHO1 and RHO2
C     to a binary file
#include "priunit.h"
      INTEGER NASHT
      REAL*8  RHO1(NASHT,NASHT),RHO2(NASHT,NASHT,NASHT,NASHT)
      INTEGER I,J,IUNIT

      OPEN(NEWUNIT=IUNIT,FILE='2RDM10',FORM='UNFORMATTED')

      WRITE(IUNIT) NASHT
      WRITE(IUNIT) RHO1(1:NASHT,1:NASHT)
      WRITE(IUNIT) RHO2(1:NASHT,1:NASHT,1:NASHT,1:NASHT)

      CLOSE(IUNIT)

      END

      SUBROUTINE DUMP_Jsr(Jsr)
#include "inforb.h"
#include "priunit.h"
      double precision :: Jsr(N2BASX)
      character(8)     :: LAB123(3)
      data LAB123/'********','********','********'/

      OPEN(newunit=iunit,FILE='dftSRfile.dat',
     $     form='unformatted')
      REWIND(IUNIT)
      WRITE(IUNIT) LAB123,'JSR     '
      write(iunit) INT(N2BASX,KIND=4)
      write(iunit) Jsr
      close(iunit)
C
c     CALL QUIT('DUMP_Jsr not implemented yet')

      END

      SUBROUTINE DUMP_Vsrxc(Vsr_xc,NSRAO)
#include "inforb.h"
#include "dfterg.h"
#include "priunit.h"
      integer :: NSRAO
      double precision :: Vsr_xc(N2BASX)
      character(8)     :: LAB123(3)
      data LAB123/'********','********','********'/

      OPEN(newunit=iunit,FILE='dftSRfile.dat',
     $     position='append',action='write',form='unformatted')
      WRITE(IUNIT) LAB123,'VSRXC   '
      WRITE(IUNIT) INT(N2BASX,KIND=4)
      WRITE(IUNIT) Vsr_xc(1:N2BASX)
      WRITE(IUNIT) ESRDFT
      CLOSE(IUNIT)

c     CALL QUIT('DUMP_Vsrxc not implemented yet')

      END

      SUBROUTINE DUMP_DFTGRID_INIT(NTOT_GRID_POINTS, DOGGA_in)
#include "priunit.h"
!      implicit none

      integer, intent(in) :: NTOT_GRID_POINTS
      logical, intent(in) :: DOGGA_in

      logical :: DOGGA
      common /cb_dump_dftgrid/ DOGGA

      DOGGA = DOGGA_in

      OPEN(938,FILE='dftgrid.dat',status='new',form='unformatted')
      WRITE(938) 'DFTGRID '
      WRITE(938) NTOT_GRID_POINTS
      WRITE(938) DOGGA
c      write(lupri,*) 'TEST NGRID', NTOT_GRID_POINTS,DOGGA
c     CALL QUIT('DUMP_DFTGRID_INIT not implemented yet')

      ! open new file for dumping grid info. and save NTOT_GRID_POINTS
      ! save NTOT_GRID_POINTS and DOGGA to this file
      ! save DOGGA in a common b

      END

      SUBROUTINE DUMP_DFTGRID_POINT(WEIGHT,COR,GSO)

!     Input:
!        WEIGHT   : weight of current grid point
!        COR(1:3) : coordinates of current grid point
!        GSO(:)   : amplitude in grid point of basis functions as symmetry orbitals
!           GSO(        i) contains amplitude of xhi(i)
!           GSO(  NBAST+i) contains value of d/dx xhi(i)
!           GSO(2*NBAST+i) contains value of d/dy xhi(i)
!           GSO(3*NBAST+i) contains value of d/dz xhi(i)
!

      implicit none
      real*8, intent(in) :: WEIGHT, COR(3), GSO(*)

! inforb.h : NBAST
#include "inforb.h"

      logical :: DOGGA
      common /cb_dump_dftgrid/ DOGGA

      IF (DOGGA) THEN
         WRITE(938) WEIGHT,COR,GSO(1:4*NBAST)
      ELSE ! LDA
         WRITE(938) WEIGHT,COR,GSO(1:NBAST)
      END IF

      END

      SUBROUTINE DUMP_DFTGRID_END()
#include "priunit.h"

      write(lupri,*) 'Closing grid dump: DUMP_DFTGRID_END()'
      CLOSE(938)
c     CALL QUIT('DUMP_DFTGRID_END not implemented yet')

      ! close file used for dumping grid info. and save NTOT_GRID_POINTS

      END
