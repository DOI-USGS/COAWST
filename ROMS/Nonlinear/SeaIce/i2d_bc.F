#include "cppdefs.h"
      MODULE i2d_bc_mod
#ifdef ICE_MODEL
! 
!***********************************************************************
!  Compute lateral boundary conditions for any 2D ice variable.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC i2d_bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE i2d_bc_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, linew,                             &
     &                        ai_west, ai_east, ai_north, ai_south,     &
     &                        ui, vi, ai, S)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
      TYPE(T_LBC), intent(in) :: S(4)

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(in)    :: ai_west(LBj:)
      real(r8), intent(in)    :: ai_east(LBj:)
      real(r8), intent(in)    :: ai_north(LBi:)
      real(r8), intent(in)    :: ai_south(LBi:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: ai_west(LBj:UBj)
      real(r8), intent(in)    :: ai_east(LBj:UBj)
      real(r8), intent(in)    :: ai_north(LBi:UBi)
      real(r8), intent(in)    :: ai_south(LBi:UBi)
# endif

!
!  Local variable declarations.
!
      integer :: i, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau 

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (S(iwest)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=ai(Istr-1,j,know)-ai(Istr-1,j-1,know)
# ifdef MASKING
            grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
# endif
            grad(Istr,j)=ai(Istr,j,know)-ai(Istr,j-1,know)
# ifdef MASKING
            grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
# endif
          END DO
          DO j=Jstr,Jend
            dTdt=ai(Istr,j,know)-ai(Istr  ,j,linew)
            dTdx=ai(Istr,j,linew)-ai(Istr+1,j,linew)
            IF (S(iwest)%nudging) THEN
              tau=Tobc_out(1,ng,iwest)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,iwest)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
              dTde=grad(Istr,j  )
            ELSE
              dTde=grad(Istr,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
# ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
# else
            Ce=0.0_r8
# endif
            ai(Istr-1,j,linew)=(cff*ai(Istr-1,j,know)+                  &
     &                          Cx *ai(Istr  ,j,linew)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Istr-1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/       &
     &                              (cff+Cx)
            IF (S(iwest)%nudging) THEN
              ai(Istr-1,j,linew)=ai(Istr-1,j,linew)+                    &
     &                       tau*(ai_west(j)-ai(Istr-1,j,know))
            END IF
# ifdef MASKING
            ai(Istr-1,j,linew)=ai(Istr-1,j,linew)*                      &
     &                              GRID(ng)%rmask(Istr-1,j)
# endif
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (S(iwest)%clamped) THEN
          DO j=Jstr,Jend
            ai(0,j,linew)=ai_west(j)
# ifdef MASKING
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                 GRID(ng)%rmask(0,j)
# endif
# ifdef WET_DRY
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                 GRID(ng)%rmask_wet(0,j)
# endif
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(iwest)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(1,j,linew).ge.0._r8) THEN
              ai(0,j,linew)=ai_west(j)
# ifdef MASKING
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask(0,j)
# endif
# ifdef WET_DRY
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask_wet(0,j)
# endif
            ELSE
              ai(0,j,linew)=ai(1,j,liold)
# ifdef MASKING
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask(0,j)
# endif
# ifdef WET_DRY
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask_wet(0,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (S(iwest)%closed) THEN
          DO j=Jstr,Jend
            ai(0,j,linew)=ai(1,j,linew)
# ifdef MASKING
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
# endif
# ifdef WET_DRY
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
# endif
          END DO
        END IF

# if defined CHUKCHI && defined OUTFLOW_MASK
        DO j=Jstr,Jend
          ai(1,j,linew)=ai_west(j)
#  ifdef MASKING
          ai(1,j,linew)=ai(1,j,linew)*                                  &
     &                 GRID(ng)%rmask(1,j)
#  endif
#  ifdef WET_DRY
          ai(1,j,linew)=ai(1,j,linew)*                                  &
     &                 GRID(ng)%rmask_wet(1,j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (S(ieast)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend,j)=ai(Iend,j,know)-ai(Iend,j-1,know)
# ifdef MASKING
            grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
# endif
            grad(Iend+1,j)=ai(Iend+1,j,know)-ai(Iend+1,j-1,know)
# ifdef MASKING
            grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
# endif
          END DO
          DO j=Jstr,Jend
            dTdt=ai(Iend,j,know)-ai(Iend  ,j,linew)
            dTdx=ai(Iend,j,linew)-ai(Iend-1,j,linew)
            IF (S(ieast)%nudging) THEN
              tau=Tobc_out(1,ng,ieast)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,ieast)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dTde=grad(Iend,j  )
            ELSE
              dTde=grad(Iend,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
# ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
# else
            Ce=0.0_r8
# endif
            ai(Iend+1,j,linew)=(cff*ai(Iend+1,j,know)+                  &
     &                          Cx *ai(Iend  ,j,linew)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                              (cff+Cx)
            IF (S(ieast)%nudging) THEN
              ai(Iend+1,j,linew)=ai(Iend+1,j,linew)+                    &
     &             tau*(ai_east(j)-ai(Iend+1,j,know))
            END IF
# ifdef MASKING
            ai(Iend+1,j,linew)=ai(Iend+1,j,linew)*                      &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (S(ieast)%clamped) THEN
          DO j=Jstr,Jend
            ai(Lm(ng)+1,j,linew)=ai_east(j)
# ifdef MASKING
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                        GRID(ng)%rmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                        GRID(ng)%rmask_wet(Lm(ng)+1,j)
# endif
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(iwest)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
              ai(Lm(ng)+1,j,linew)=ai_east(j)
# ifdef MASKING
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
# endif
            ELSE
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,liold)
# ifdef MASKING
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (S(ieast)%closed) THEN
          DO j=Jstr,Jend
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
# ifdef MASKING
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
# endif
          END DO
        END IF

# if defined CHUKCHI && defined OUTFLOW_MASK
        DO j=Jstr,Jend
          ai(Lm(ng),j,linew)=ai_east(j)
#  ifdef MASKING
          ai(Lm(ng),j,linew)=ai(Lm(ng),j,linew)*                        &
     &                        GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng),j,linew)=ai(Lm(ng),j,linew)*                        &
     &                        GRID(ng)%rmask_wet(Lm(ng),j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (S(isouth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr)=ai(i,Jstr,know)-ai(i-1,Jstr,know)
# ifdef MASKING
            grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
# endif
            grad(i,Jstr-1)=ai(i,Jstr-1,know)-ai(i-1,Jstr-1,know)
# ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
# endif
          END DO
          DO i=Istr,Iend
            dTdt=ai(i,Jstr,know)-ai(i,Jstr  ,linew)
            dTde=ai(i,Jstr,linew)-ai(i,Jstr+1,linew)
            IF (S(isouth)%nudging) THEN
              tau=Tobc_out(1,ng,isouth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,isouth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jstr)
            ELSE
              dTdx=grad(i+1,Jstr)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
# ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
# else
            Cx=0.0_r8
# endif
            Ce=dTdt*dTde
            ai(i,Jstr-1,linew)=(cff*ai(i,Jstr-1,know)+                  &
     &                          Ce *ai(i,Jstr  ,linew)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/       &
     &                              (cff+Ce)
            IF (S(isouth)%nudging) THEN
              ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)+                    &
     &           tau*(ai_south(i)-ai(i,Jstr-1,know))
            END IF
# ifdef MASKING
            ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (S(isouth)%clamped) THEN
          DO i=Istr,Iend
            ai(i,0,linew)=ai_south(i)
# ifdef MASKING
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
# endif
# ifdef WET_DRY
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
# endif
          END DO
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(isouth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0._r8) THEN
              ai(i,0,linew)=ai_south(i)
# ifdef MASKING
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask(i,0)
# endif
# ifdef WET_DRY
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask_wet(i,0)
# endif
            ELSE
              ai(i,0,linew)=ai(i,1,liold)
# ifdef MASKING
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask(i,0)
# endif
# ifdef WET_DRY
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask_wet(i,0)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (S(isouth)%closed) THEN
          DO i=Istr,Iend
            ai(i,0,linew)=ai(i,1,linew)
# ifdef MASKING
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
# endif
# ifdef WET_DRY
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
# endif
          END DO
        END IF

# if defined CHUKCHI && defined OUTFLOW_MASK
        DO i=Istr,Iend
          ai(i,1,linew)=ai_south(i)
#  ifdef MASKING
          ai(i,1,linew)=ai(i,1,linew)*                                  &
     &                   GRID(ng)%rmask(i,1)
#  endif
#  ifdef WET_DRY
          ai(i,1,linew)=ai(i,1,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,1)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (S(inorth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend)=ai(i,Jend,know)-ai(i-1,Jend,know)
# ifdef MASKING
            grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
# endif
            grad(i,Jend+1)=ai(i,Jend+1,know)-ai(i-1,Jend+1,know)
# ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
# endif
          END DO
          DO i=Istr,Iend
            dTdt=ai(i,Jend,know)-ai(i,Jend  ,linew)
            dTde=ai(i,Jend,linew)-ai(i,Jend-1,linew)
            IF (S(inorth)%nudging) THEN
              tau=Tobc_out(1,ng,inorth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,inorth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jend)
            ELSE
              dTdx=grad(i+1,Jend)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
# ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
# else
            Cx=0.0_r8
# endif
            Ce=dTdt*dTde
            ai(i,Jend+1,linew)=(cff*ai(i,Jend+1,know)+                  &
     &                          Ce *ai(i,Jend  ,linew)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/       &
     &                              (cff+Ce)
            IF (S(inorth)%nudging) THEN
              ai(i,Jend+1,linew)=ai(i,Jend+1,linew)+                    &
     &              tau*(ai_north(i)-ai(i,Jend+1,know))
            END IF
# ifdef MASKING
            ai(i,Jend+1,linew)=ai(i,Jend+1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jend+1)
# endif
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (S(inorth)%clamped) THEN
          DO i=Istr,Iend
            ai(i,Mm(ng)+1,linew)=ai_north(i)
# ifdef MASKING
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
# endif
# ifdef WET_DRY
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
# endif
          END DO
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(inorth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
              ai(i,Mm(ng)+1,linew)=ai_north(i)
# ifdef MASKING
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
# endif
# ifdef WET_DRY
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
# endif
            ELSE
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),liold)
# ifdef MASKING
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
# endif
# ifdef WET_DRY
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (S(inorth)%closed) THEN
          DO i=Istr,Iend
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
# ifdef MASKING
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
# endif
# ifdef WET_DRY
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
# endif
          END DO
        END IF

# if defined CHUKCHI && defined OUTFLOW_MASK
        DO i=Istr,Iend
          ai(i,Mm(ng),linew)=ai_north(i)
#  ifdef MASKING
          ai(i,Mm(ng),linew)=ai(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng),linew)=ai(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng))
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          ai(0,0,linew)=0.5_r8*(ai(1,0,linew)+                          &
     &                         ai(0,1,linew))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ai(Lm(ng)+1,0,linew)=0.5_r8*(ai(Lm(ng)+1,1,linew)+            &
     &                                ai(Lm(ng)  ,0,linew))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ai(0,Mm(ng)+1,linew)=0.5_r8*(ai(0,Mm(ng)  ,linew)+            &
     &                                ai(1,Mm(ng)+1,linew))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ai(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                           &
     &             (ai(Lm(ng)+1,Mm(ng)  ,linew)+                        &
     &              ai(Lm(ng)  ,Mm(ng)+1,linew))
        END IF
      END IF
      RETURN
      END SUBROUTINE i2d_bc_tile
#endif

      END MODULE i2d_bc_mod
