!
! Basal respiration of phytoplankton 
! Basal Respiration is a function of phytoplankton biomass and temperature
!
! Parameters and computation
!
                br20   = 0.025_r8               ! Basal respiration at 20 degress
                brthta = 1.047_r8               ! Basal respiration temperature coefficient
                cff1=br20*brthta**(itemp)       ! Respiration coefficient at temperature (/d) 
                cff2=dtdays*cff1
!
                N_Flux_BaseResp=cff2*max(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
!
! Update other model systems (In future, consider if basal resp goes to iSDeN and iSDeC)
!
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_BaseResp
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_BaseResp

#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       rOxNH4*(N_Flux_BaseResp)
#endif

#ifdef CARBON
!
!  Total inorganic carbon (CO2) released during phytoplankton basal respiration.
!
                cff1=PhyCN(ng)*(N_Flux_BaseResp)
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-cff1
#endif

