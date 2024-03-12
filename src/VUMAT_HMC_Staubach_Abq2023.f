!=======================================================================================================
!This file is part of VUMAT_HMC_Staubach_Abq2023.
!
!VUMAT_HMC_Staubach_Abq2023 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License !as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!VUMAT_HMC_Staubach_Abq2023 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with VUMAT_HMC_Staubach_Abq2023. If not, see <https://www.gnu.org/licenses/>. 
!=======================================================================================================
!
! SUBROUTINE: VUMAT_HMC_Staubach
!
!> @author Patrick Staubach, patrick.staubach@yahoo.de
!          Bauhaus University Weimar, Ruhr-University Bochum
!
! DESCRIPTION:
!> @brief Contains the hydro-mechanically coupled explicit vumat interface
!> @brief Implementation according to "Vibratory pile driving in water-saturated sand:
!          Back-analysis of model tests using a hydro-mechanically coupled CEL method"
!  https://doi.org/10.1016/j.sandf.2020.11.005
!
! REVISION HISTORY
!> @date 02.02.2021 - Initial version
!=======================================================================================================
      subroutine vumat(
c Read only (unmodifiable) variables -
     &  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     &  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     &  props, density, strainInc, relSpinInc,
     &  tempOld, stretchOld, defgradOld, fieldOld,
     &  stressOld, stateOld, enerInternOld, enerInelasOld,
     &  tempNew, stretchNew, defgradNew, fieldNew,
c Write only (modifiable) variables -
     &  stressNew, stateNew, enerInternNew, enerInelasNew)
c 
      include 'vaba_param.inc'
      !implicit none
c 
c Variable declaration
c 
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     &  charLength(nblock),dtArray(2*(nblock)+1), strainInc(nblock, ndir+nshr),
     &  relSpinInc(nblock,nshr), tempOld(nblock),
     &  stretchOld(nblock,ndir+nshr),
     &  defGradOld(nblock, ndir+nshr+nshr),
     &  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     &  stateOld(nblock,nstatev), enerInternOld(nblock),
     &  enerInelasOld(nblock), tempNew(nblock),
     &  stretchNew(nblock,ndir+nshr),
     &  defgradNew(nblock,ndir+nshr+nshr),
     &  fieldNew(nblock,nfieldv),
     &  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     &  enerInternNew(nblock), enerInelasNew(nblock)
 
      character*80 cmname
      
      !double precision stepTime, totalTime

      integer ntens, kinc, kstep
      ! ndir nstatev, nprops, lanneal
      integer i, iblock
      !, nblock, nfieldv, nshr

      double precision time(2), predef(1), drot(3,3), dpred(1)
      double precision dfgrd0(3,3), dfgrd1(3,3), div_acc, water_table
      double precision tr_dstran,density2,KPerm,gravity,KO,dir_grav
      double precision dot_pw,viscosity,pore_pressure,gamma_w,cavitation
      double precision stress(ndir+nshr), statev(nstatev)
      double precision ddsdde(ndir+nshr,ndir+nshr), ddsddt(ndir+nshr)
      double precision drplde(ndir+nshr), dstran(ndir+nshr)
      double precision stran(ndir+nshr), coords(3)

      double precision drpldt
      double precision dtime, dt
      double precision mat_props(nprops)

      integer isfirstinc, ndi, nshrmat, nmat_props, nstatev_mat
      
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       !Variables for pore water pressure
       KPerm       = 1.00d-10 ! permeability
       viscosity   = 1.00d-6  ! dynamic viscosity
       density2    = 1.8649d0 ! density
       gamma_w     = 10.0d0   ! dead weight water
       gravity     = 10.0d0   ! gravity
       KO          = 0.4d0    ! lateral stress coefficient
       water_table = -1.0d0   ! water table height
       dir_grav    = 3        ! direction of gravity
       cavitation  = -100.0d0 ! Total pore water pressure at which cavitation occurs
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !! Output variables
      !! statev(23) : tr_dstran
      !! statev(24) : pore water pressure
      !! statev(25) : excess pore water pressure
      !! statev(26) : divergence of acceleration
      !! statev(27:26+ntens) : effective stress
      !! statev(33) : 2: below water table, 1 above
!*************************************************************************

      ntens       = ndir + nshr
      mat_props   = props
      time(1)     = stepTime
      time(2)     = totalTime
      dtime       = dtArray(1)
      dt          = dtArray(1)
      ndi         = ndir
      nshrmat     = nshr
      nmat_props  = nprops
      nstatev_mat = nstatev
        
      if ((time(1).le.0.0d0).and.(time(2).le.0.10d0)) then
        isfirstinc=1
        KSTEP=1
        KINC=1
      else
        KSTEP=10
        KINC=10
      endif
      
!*************************************************************************
               
      do iblock = 1, nblock     

        statev(:) = stateOld(iblock,:)
        
        do i = 1, ndi
          stress(i)     = stressOld(iblock,i)
          dstran(i)     = strainInc(iblock,i)
          coords(i)     = coordMp(iblock,i)
        enddo
        
        stress(4) = stressOld(iblock,4)
        dstran(4) = 2.0d0 * strainInc(iblock,4)
        
        if (nshr > 1) then      
          stress(6)  = stressOld(iblock,5)
          stress(5) = stressOld(iblock,6)
          dstran(5) = 2.0d0 * strainInc(iblock,6)
          dstran(6) = 2.0d0 * strainInc(iblock,5)
        endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      tr_dstran = 0.0d0
      do i = 1, ndi
        tr_dstran = tr_dstran + dstran(i)
      enddo
     
      div_acc = 0.0d0
      dot_pw = 0.0d0
      !! Assign new div acceleration
      if(dt>1d-11) then
        div_acc = (tr_dstran-statev(23)) / dt /2.0d0
     &               + statev(26)/2.0d0
        dot_pw =  (tempNew(iblock) - tempOld(iblock))/dt
      else
        div_acc = 0.0d0
        dot_pw  = 0.0d0
      endif
     
      !! Consider hydrostatic pressure from ground water 
      statev(24) = 0.0d0
      statev(24) = (coords(int(dir_grav)) + water_table)*gamma_w
      if(statev(24)<0.0d0) statev(24) = 0.0d0
     
       !! Effective stress
      if(int(dir_grav) == 3) then
       do i = 1, ndi
         if (i == 1) stress(i) = stress(i) + statev(25)
     &                         + statev(24)*KO
         if (i == 2) stress(i) = stress(i) + statev(25)
     &                         + statev(24)*KO
         if (i == 3) stress(i) = stress(i) + statev(25) 
     &                         + statev(24)
       enddo
       else
       do i = 1, ndi
         if (i == 1) stress(i) = stress(i) + statev(25)
     &                         + statev(24)*KO
         if (i == 2) stress(i) = stress(i) + statev(25)
     &                         + statev(24)
         if (i == 3) stress(i) = stress(i) + statev(25) 
     &                         + statev(24)*KO
       enddo
       endif
       !! Assign excess pore water pressure change
       statev(25) = statev(25) + dot_pw*dt

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        if (isfirstinc == 1) then
           
           mat_props(1) = 8000.0d0
           mat_props(2) = 0.3d0
            
           call UMATElastic(stress,statev,ddsdde,0,0,0,
     &       0,ddsddt,drplde,drpldt,
     &       stran,dstran,time,dtime,0,0,predef,dpred,cmname,
     &       ndi,nshrmat,ntens,nstatev_mat,mat_props,nmat_props,coords,
     &       drot,1,0,dfgrd0,dfgrd1,1,0,0,0,kstep,kinc)
 
        else
           

          call umat(stress,statev,ddsdde,0,0,0,
     &      0,ddsddt,drplde,drpldt,
     &      stran,dstran,time,dtime,0,0,predef,dpred,cmname,
     &      ndi,nshrmat,ntens,nstatev_mat,mat_props,nmat_props,coords,
     &      drot,1,0,dfgrd0,dfgrd1,1,0,0,0,kstep,kinc)
           
        endif
     
    
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        statev(23) = tr_dstran 
        !! Assign new div acc
        statev(26) = div_acc
        !! Assign inelastic energy
        if((statev(25)+ statev(24))<cavitation) then
          if((-tr_dstran/density2
     &       + div_acc*KPerm/viscosity/density2)<0.0d0) then
            statev(34) = 1.0d0
          else
            enerInelasNew(iblock) = enerInelasOld(iblock)
     &                        - tr_dstran/density2
     &                        + div_acc*KPerm/viscosity/density2
          endif
        else
          enerInelasNew(iblock) = enerInelasOld(iblock)
     &                      - tr_dstran/density2
     &                      + div_acc*KPerm/viscosity/density2
        endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !! Assign effective stress to statev
        statev(27:26+ntens) = stress

       if(int(dir_grav) == 3) then
       do i = 1, ndi
         if (i == 1) stressNew(iblock,i) = stress(i) - statev(25)
     &                         - statev(24)*KO
         if (i == 2) stressNew(iblock,i) = stress(i) - statev(25)
     &                         - statev(24)*KO
         if (i == 3) stressNew(iblock,i) = stress(i) - statev(25) 
     &                         - statev(24)
       enddo
       else
       do i = 1, ndi
         if (i == 1) stressNew(iblock,i) = stress(i) - statev(25)
     &                         - statev(24)*KO
         if (i == 2) stressNew(iblock,i) = stress(i) - statev(25)
     &                         - statev(24)
         if (i == 3) stressNew(iblock,i) = stress(i) - statev(25) 
     &                         - statev(24)*KO
       enddo
       endif
       
        stressNew(iblock,4) = stress(4)
      
        if( nshr > 1 ) then
          stressNew(iblock,6)  = stress(5)
          stressNew(iblock,5)  = stress(6)
        endif

        stateNew(iblock,:) = statev(:)
   
      enddo

      end
    
      SUBROUTINE UMATElastic(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPSU,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      real*8 PROPSU(nprops),STRESS(NTENS),DSTRAN(NTENS)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
      E=PROPSU(1)
      ANU=PROPSU(2)
      ALAMDA=ANU*E/ (ONE+ANU)/(ONE-TWO*ANU)
      AMU=E/TWO/(ONE+ANU)
  
      DO I=1,NTENS
       DO J=1,NTENS
        DDSDDE(I,J)=0.0D0
       ENDDO
      ENDDO
      DDSDDE(1,1)=ALAMDA+TWO*AMU
      DDSDDE(2,2)=DDSDDE(1,1)
      DDSDDE(3,3)=DDSDDE(1,1)
      DDSDDE(4,4)=AMU
      if (ntens==6) then
      DDSDDE(5,5)=AMU
      DDSDDE(6,6)=AMU
      endif
      DDSDDE(1,2)=ALAMDA
      DDSDDE(1,3)=ALAMDA
      DDSDDE(2,3)=ALAMDA
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,3)
      DDSDDE(3,2)=DDSDDE(2,3)
C
      DO I=1,NTENS
        DO J=1,NTENS
        STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
        ENDDO
      ENDDO

      RETURN
      END 
