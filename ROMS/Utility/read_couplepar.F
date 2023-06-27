#include "cppdefs.h"
#if defined MODEL_COUPLING && !defined MCT_LIB
      SUBROUTINE read_CouplePar (model)
!
!git $Id$
!svn $Id: read_couplepar.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports multiple model coupling input        !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_iounits
      USE mod_scalars
!
      USE inp_decode_mod
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: Lwrite
      logical :: Lvalue(1)
!
      integer :: Npts, Nval, i, ic, io_err, j, inp, ng, out, status
      integer :: Ivalue(1)
!
      real(r8), dimension(nRval) :: Rval

      real(r8), allocatable :: MyRval(:)
!
      character (len=40 ) :: KeyWord
      character (len=256) :: io_errmsg, line
      character (len=256) :: Cname
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Determine coupling standard input file name.  In distributed-memory,
!  this name is assigned at the executtion command line and processed
!  with the Unix routine GETARG.  The ROMS/TOMS input parameter script
!  name is specified in this coupling script.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      Lwrite=Master
      inp=1
      out=stdout
!
      IF (MyRank.eq.0) CALL my_getarg (1,Cname)
      CALL mp_bcasts (1, model, Cname)
      IF (MyRank.eq.0) THEN
        WRITE(stdout,*) 'Coupled Input File name = ', TRIM(Cname)
      END IF
      OPEN (inp, FILE=TRIM(Cname), FORM='formatted', STATUS='old',      &
     &      IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        IF (MyRank.eq.0) WRITE (stdout,10) TRIM(io_errmsg)
        exit_flag=5
        RETURN
 10     FORMAT (/,' READ_COUPLEPAR - Unable to open coupling input',    &
     &          ' script.',/18x,'ERROR: ',a,/,                          &
     &          /,18x,'In distributed-memory applications, the input',  &
     &          /,18x,'script file is processed in parallel. The Unix', &
     &          /,18x,'routine GETARG is used to get script file name.',&
     &          /,18x,'For example, in MPI applications make sure that',&
     &          /,18x,'command line is something like:',/,              &
     &          /,18x,'mpirun -np 4 romsM coupling.in',/,               &
     &          /,18x,'and not',/,                                      &
     &          /,18x,'mpirun -np 4 romsM < coupling.in',/)
      END IF
# else
      Lwrite=Master
      inp=stdinp
      out=stdout
# endif
!
!-----------------------------------------------------------------------
!  Read in multiple models coupling parameters. Then, load input
!  data into module. Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Nmodels')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nmodels=Ivalue(1)
              IF (.not.allocated(MyRval) ) THEN
                allocate ( MyRval(Nmodels) )
              END IF
              IF (.not.allocated(OrderLabel) ) THEN
                allocate ( OrderLabel(Nmodels) )
              END IF
              IF (.not.allocated(Nthreads) ) THEN
                allocate ( Nthreads(Nmodels) )
                Nthreads=0
              END IF
              IF (.not.allocated(TimeInterval) ) THEN
                allocate ( TimeInterval(Nmodels,Nmodels) )
                TimeInterval=0.0_r8
              END IF
              IF (.not.allocated(INPname) ) THEN
                allocate ( INPname(Nmodels) )
              END IF
              IF (.not.allocated(Nexport) ) THEN
                allocate ( Nexport(Nmodels) )
                Nexport=0
              END IF
              IF (.not.allocated(Nimport) ) THEN
                allocate ( Nimport(Nmodels) )
                Nimport=0
              END IF
            CASE ('Lreport')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Lreport=Lvalue(1)
            CASE ('OrderLabel')
              DO i=1,Nmodels
                IF (i.eq.Nval) THEN
                  OrderLabel(i)=TRIM(ADJUSTL(Cval(Nval)))
                  IF (INDEX(TRIM(OrderLabel(i)),'ocean').ne.0) THEN
                    Iocean=i
                  ELSE IF (INDEX(TRIM(OrderLabel(i)),'waves').ne.0) THEN
                    Iwaves=i
                  ELSE IF (INDEX(TRIM(OrderLabel(i)),'atmos').ne.0) THEN
                    Iatmos=i
                  END IF
                END IF
              END DO
            CASE ('Nthreads(ocean)')
              IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nthreads(Iocean)=Ivalue(1)
              END IF
            CASE ('Nthreads(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nthreads(Iwaves)=Ivalue(1)
              END IF
            CASE ('Nthreads(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nthreads(Iatmos)=Ivalue(1)
              END IF
            CASE ('TimeInterval')
              Npts=load_r(Nval, Rval, Nmodels, MyRval)
              ic=0
              DO j=1,Nmodels
                DO i=1,Nmodels
                  IF (i.gt.j) THEN
                    ic=ic+1
                    TimeInterval(i,j)=MyRval(ic)
                    TimeInterval(j,i)=MyRval(ic)
                  END IF
                END DO
              END DO
            CASE ('INPname(ocean)')
              IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
                INPname(Iocean)=TRIM(ADJUSTL(Cval(Nval)))
                Iname=TRIM(INPname(Iocean))
              END IF
            CASE ('INPname(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                INPname(Iwaves)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('INPname(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                INPname(Iatmos)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('CPLname')
              CPLname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('Nexport(ocean)')
              IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nexport(Iocean)=Ivalue(1)
              END IF
            CASE ('Nexport(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nexport(Iwaves)=Ivalue(1)
              END IF
            CASE ('Nexport(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nexport(Iatmos)=Ivalue(1)
              END IF
            CASE ('Export(ocean)')
              IF (.not.allocated(Export)) THEN
                allocate ( Export(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Export(i)%code(MAX(1,Nexport(i))) )
                  Export(i)%code=' '
                END DO
              END IF
              IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
                IF (Nval.le.Nexport(Iocean)) THEN
                  Export(Iocean)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Export(waves)')
              IF (.not.allocated(Export)) THEN
                allocate ( Export(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Export(i)%code(MAX(1,Nexport(i))) )
                  Export(i)%code=' '
                END DO
              END IF
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                IF (Nval.le.Nexport(Iwaves)) THEN
                  Export(Iwaves)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Export(atmos)')
              IF (.not.allocated(Export)) THEN
                allocate ( Export(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Export(i)%code(MAX(1,Nexport(i))) )
                  Export(i)%code=' '
                END DO
              END IF
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                IF (Nval.le.Nexport(Iatmos)) THEN
                  Export(Iatmos)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Nimport(ocean)')
              IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nimport(Iocean)=Ivalue(1)
              END IF
            CASE ('Nimport(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nimport(Iwaves)=Ivalue(1)
              END IF
            CASE ('Nimport(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Ivalue)
                Nimport(Iatmos)=Ivalue(1)
              END IF
            CASE ('Import(ocean)')
              IF (.not.allocated(Import)) THEN
                allocate ( Import(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Import(i)%code(MAX(1,Nimport(i))) )
                  Import(i)%code=' '
                END DO
              END IF
            IF ((0.lt.Iocean).and.(Iocean.le.Nmodels)) THEN
              IF (Nval.le.Nimport(Iocean)) THEN
                Import(Iocean)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            END IF
            CASE ('Import(waves)')
              IF (.not.allocated(Import)) THEN
                allocate ( Import(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Import(i)%code(MAX(1,Nimport(i))) )
                  Import(i)%code=' '
                END DO
              END IF
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                IF (Nval.le.Nimport(Iwaves)) THEN
                  Import(Iwaves)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Import(atmos)')
              IF (.not.allocated(Import)) THEN
                allocate ( Import(Nmodels) )
                DO i=1,Nmodels
                  allocate ( Import(i)%code(MAX(1,Nimport(i))) )
                  Import(i)%code=' '
                END DO
              END IF
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                IF (Nval.le.Nimport(Iatmos)) THEN
                  Import(Iatmos)%code(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
          END SELECT
        END IF
      END DO
  20  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  30  CLOSE (inp)
!
  40  FORMAT (/,' READ_CouplePar - Error while processing line: ',/,a)
!
      RETURN
      END SUBROUTINE read_CouplePar
#else
      SUBROUTINE read_CouplePar
      END SUBROUTINE read_CouplePar
#endif
