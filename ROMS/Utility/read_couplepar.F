#include "cppdefs.h"
#ifdef MODEL_COUPLING
      SUBROUTINE read_CouplePar (model)
!
!svn $Id: read_couplepar.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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

      integer :: Npts, Nval, i, ic, j, inp, ng, out, status
      integer :: Ivalue(1)

      integer :: decode_line, load_i, load_l, load_r

      real(r8), dimension(200) :: Rval

      real(r8), allocatable :: MyRval(:)

      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256) :: Cname
      character (len=256), dimension(200) :: Cval
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
     &      ERR=10)
      GO TO 30
 10   WRITE (stdout,20)
      IF (Master) WRITE(stdout,*) 'MyRank = ', MyRank, TRIM(Cname)
      exit_flag=4
      RETURN
 20   FORMAT (/,' INP_PAR - Unable to open coupling input script.',     &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 masterM coupling.in',/,               &
     &        /,11x,'and not',/,                                        &
     &        /,11x,'mpirun -np 4 masterM < coupling.in',/)
 30   CONTINUE
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
        READ (inp,'(a)',ERR=40,END=50) line
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
                Npts=load_i(Nval, Rval, 1, Nthreads(Iocean))
              END IF
            CASE ('Nthreads(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nthreads(Iwaves))
              END IF
            CASE ('Nthreads(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nthreads(Iatmos))
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
                Npts=load_i(Nval, Rval, 1, Nexport(Iocean))
              END IF
            CASE ('Nexport(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nexport(Iwaves))
              END IF
            CASE ('Nexport(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nexport(Iatmos))
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
                Npts=load_i(Nval, Rval, 1, Nimport(Iocean))
              END IF
            CASE ('Nimport(waves)')
              IF ((0.lt.Iwaves).and.(Iwaves.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nimport(Iwaves))
              END IF
            CASE ('Nimport(atmos)')
              IF ((0.lt.Iatmos).and.(Iatmos.le.Nmodels)) THEN
                Npts=load_i(Nval, Rval, 1, Nimport(Iatmos))
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
  40  IF (Master) WRITE (out,60) line
      exit_flag=4
      RETURN
  50  CLOSE (inp)

  60  FORMAT (/,' READ_CouplePar - Error while processing line: ',/,a)

      RETURN
      END SUBROUTINE read_CouplePar
#else
      SUBROUTINE read_CouplePar
      END SUBROUTINE read_CouplePar
#endif
