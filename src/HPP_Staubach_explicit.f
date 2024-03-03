!=======================================================================================================
!This file is part of VUMAT_HMC_Staubach_Abq2020.
!
!VUMAT_HMC_Staubach_Abq2020 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License !as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!VUMAT_HMC_Staubach_Abq2020 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with VUMAT_HMC_Staubach_Abq2020. If not, see <https://www.gnu.org/licenses/>. 
!=======================================================================================================
!
! SUBROUTINE: HPP_Staubach
!
!> @author Patrick Staubach, patrick.staubach@yahoo.de     
!          Bauhaus University Weimar, formerly Karlsruhe (IBF)
!
! DESCRIPTION:
!> @brief Contains the hypoplastic constitutive model with the extension of the intergranular strain
!> @brief Uses an adaptive Newton scheme to evaluate the stress increment and the jacobian
! REVISION HISTORY
!> @date 20.05.2018 - Initial version   
!=======================================================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      use tools 
      
      implicit none

      CHARACTER*80 CMNAME
      
      integer nprops,ntens,kinc,
     1 ndi,nshr,nstatv,noel,npt,layer,kspt,KSTEP
      
      real(8) STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      double precision celent
      real(8) temp,sse,spd,scd,rpl,dtemp,drpldt,pnewdt,dtime

      
      !Constitutive
      ! ------------------------------------
      real(8) :: dstranU(ntens)                             !subincrement of strain (geotechnical sign)
      real(8) :: doth(3,3)                                  !increment of intergranular strain
      real(8) :: T(3,3),D(3,3)                              !subincrement of strain (geotechnical sign)
      real(8) :: Mstatev(NSTATV)                            !copy of statev vector, updated in iteration
      real(8) :: Mstiff0(ntens,ntens),Mstiff1(ntens,ntens)  !Stiffness matrices in first iteration respectively second
      real(8) :: Tcorrected(3,3)                            !Stress after projection on the bounding surface

      !Variables needed for adaptive newton scheme
      ! ------------------------------------
      real(8) :: mDTime,mTime,Error,Dstress_0(ntens),Dstress_1(ntens)
      real(8) :: Delta_stress_0(ntens),Delta_stress_1(ntens)
      real(8) :: stress0(ntens),stressCorrector(ntens),stress1(ntens)
      
      !Helpful variables
      ! ------------------------------------
      integer :: i_proj

      !Parameters
      ! ------------------------------------
      real(8), parameter :: errorTol      = 0.01d0 ! relative error tolerance for stress (should be about 0,1-0,5 %)
      real(8), parameter :: zero =0.0d0
      
      !Viscosity
      ! ------------------------------------
      logical            :: useVisco      = .true.
      integer, parameter :: npropsVisco   = 6
      real(8)            :: propsVisco(6)
      real(8)            :: stress_visc(ntens),CMATRIX(ntens,ntens)

      !Bounding projection
      ! ------------------------------------
      logical            :: useprojection = .true.
      real(8), parameter :: pmin          = -0.01d0  !minimum allowed mean pressure
      
      !Material properties (see Subroutine hypoStiff for detailed list)
      ! ------------------------------------

      !Geotechnical sign convention
      ! ------------------------------------
      stress          = -stress
      dstran          = -dstran

      !Substepping
      ! ------------------------------------
      Dstress_0       = zero                                            !current substep increment
      Dstress_1       = zero                                            !from previous substep increment
      stress0         = stress                                          
      stress1         = stress                                          
      Delta_stress_0  = zero                                            !first substep inc
      Delta_stress_1  = zero                                            !second substep inc
      stressCorrector = zero
      dstranU         = zero
      mTime           = zero

	  
      !propsVisco
      ! ------------------------------------
      propsVisco      = 0.0d0
      propsVisco(1)   = 0.5d0             !minimum of lower stress point (about 0.5 kPa, positive!)
      propsVisco(2)   = 1.0d0              !maximum of lower stress point (about 2.0 kPa, positive!)
      propsVisco(3)   = 2.0d0              !maximum of lambda
      propsVisco(4)   = 0.8d0              !minimum of lambda
      propsVisco(5)   = 2.0d0              !maximum of mou
      propsVisco(6)   = 0.8d0              !minimum of mou

      !Initialize
      ! ------------------------------------
      DDSDDE          = zero
      stress_visc     = zero
      i_proj          = 0
      T               = zero
      D               = zero
      
      !State variables (see Subroutine hypoStiff for detailed list)
      ! ------------------------------------  
      Mstatev         = zero
      Mstatev         = statev !copy of the statevs vector
      
      !Initial time step
      ! ------------------------------------
      mDTime          = 1.0d0  !Try to achieve convergence in just one step (mDtime=1.0d0)
      
      !Add viscous stress stress_visco
      ! ------------------------------------
      if(useVisco) then   
          
          stress = stress + Mstatev(2+ntens:1+ntens+ntens) !viscosity (pressure is negative in the statev vector!)
          
      endif

      
      !Get full tensor notations
      ! ------------------------------------
      T = map2T(stress,ntens)

      !Projection of vanishing stresses on the bounding surface
      ! ------------------------------------
      if (useprojection) then
          
          !Compute corrected stress state Tcorrected
          ! ------------------------------------
          call proj_bounding(-T,Tcorrected,props(1),pmin,i_proj) !uses mechanical sign convention
          
          T       = -Tcorrected !back to geotechnical sign convention
          stress  = map2stress(T,ntens)
          stress0 = stress
          stress1 = stress
      endif
	  
      if(statev(1)>1.05*props(6)) statev(1)=1.05*props(6)
      if(statev(1)<0.95*props(5)) statev(1)=0.95*props(5)
      
      !Start Subincrement
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      do while (mTime < 1.0d0) 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Substep increments of strain
      ! ------------------------------------
      dstranU = dstran*mDTime                            
      
      !Get full tensor notations
      ! ------------------------------------
      T = map2T(stress,ntens)
      D = map2D(dstranU,ntens)

      !Get hypoplastic stiffness (Mstiff0) using T, hypoStiff does not update Mstatev or T! It only gives Mstiff and doth (rate of intergranular strain)
      ! ------------------------------------
      call hypoStiff(T,ntens,props,nprops,D,Mstatev,nstatv,
     &               Mstiff0,doth,noel)
     
      !Now calculate stress increment
      !Substepping stress increments
      ! ------------------------------------
      Delta_stress_0 = matmul(Mstiff0,dstranU)
      
      ! First approximation of stress increment
      ! ------------------------------------
      Dstress_1 = Dstress_0 + Delta_stress_0
      
      !Updated stress from  first approximation 
      ! ------------------------------------
      stress1  = stress0 + Dstress_1       
      
      !===============================================================================================
      !Repeat everything
      !Stiffness will be calculated again using the updated stress state (stress1)
      !===============================================================================================

      !Get full tensor notations
      ! ------------------------------------
      T = map2T(stress1,ntens)
      
      !Projection of vanishing stresses on the bounding surface
      ! ------------------------------------
      if (useprojection) then
          
          !Compute corrected stress state Tcorrected
          ! ------------------------------------
          call proj_bounding(-T,Tcorrected,props(1),pmin,i_proj) !uses mechanical sign convention
          
          T = -Tcorrected !back to geotechnical sign convention
      endif
      
      !Get hypoplastic stiffness (Mstiff) using T
      ! ------------------------------------
      call hypoStiff(T,ntens,props,nprops,D,Mstatev,nstatv,
     &               Mstiff1,doth,noel)
      
      !Now calculate stress increment
      !Substepping stress increments
      ! ------------------------------------
      Delta_stress_1 = matmul(Mstiff1,dstranU)
      
      ! Second improved approximation
      ! ------------------------------------
      stressCorrector = Dstress_0+0.5d0*(Delta_stress_0+Delta_stress_1) !Mid-point-rule
      
      ! Compute the relative error of stress
      ! ------------------------------------
      if(norm2(map2T(stressCorrector,ntens))>0.00000000001d0) then
        Error  = norm2(map2T(stressCorrector-Dstress_1,ntens))
     &         / norm2(map2T(stressCorrector,ntens))
      else
        Error = 1d-14
      endif

      if (Error > errorTol) then
      ! Error is too large, compute new smaller time increment and try again
      ! ------------------------------------
          mDTime              = mDtime*min(0.9d0*(errorTol/Error)
     &                          **(0.5d0),0.1d0)
          
          if(mDTime < 0.0000000000001d0) then
               write(*,*) '-----------------------------------'
               write(*,*) '--------$HPP_UMAT_Staubach$--------'
               write(*,*) '-----------------------------------'
               write(*,*) 'Substep didnt converge'
               write(*,*) '-> Needed Time increment is too small'
               write(*,*) 'Element=',NOEL
               write(*,*) 'The calculation has been terminated'
               write(*,*) 'Current State (2D case):'
               write(*,*) '-----------------------------------'
               write(*,*) 'STRESS11:',T(1,1),'STRESS22:',T(2,2)
               write(*,*) 'STRESS33:',T(3,3),'STRESS12:',T(1,2)
               write(*,*) 'Trace Stress:',tr(T)
               write(*,*) 'Epor:',Mstatev(1)
               write(*,*) 'Strain (11,22,33,12):'
               write(*,*)  dstran(1:ntens) 
               write(*,*) 'Intergranular strain (11,22,33,12):'
               write(*,*)  Mstatev(1+1:1+ntens) 
               write(*,*) 'Projection on bounding surface:'
               write(*,*)  i_proj,'0: no projection made'
               write(*,*) '-----------------------------------'
               write(*,*) 'Geotechnical sign convention!'
               
              !Abaqus should try again with 1/4 of the previous step size
              ! ------------------------------------ 
              pnewdt = 0.25d0
              return
          endif
      
      else
      !Step is accepted, the time increment (mDtime) is added to the total time 
      !and state variables, stress and time increment are updated
      ! ------------------------------------ 
          
          !Update integranular strain
          ! ------------------------------------ 
          Mstatev(2:1+ntens)   = Mstatev(2:1+ntens)
     &                         + map2stran(doth,ntens)
          
          
          !Update void ratio
          ! ------------------------------------ 
          Mstatev(1)           = Mstatev(1)-(dstranU(1)+dstranU(2)
     &                         + dstranU(3))*(1.0d0+statev(1))

          !Cumulative change in stress 
          ! ------------------------------------ 
          Dstress_0            = stressCorrector
          
          !Updated stress from previous substep
          ! ------------------------------------
          stress               = stress0 + Dstress_0

          !Projection of vanishing stresses on the bounding surface
          ! ------------------------------------
          if (useprojection) then

              !Get full tensor notations
              ! ------------------------------------
              T = map2T(stress,ntens)
          
              !Compute corrected stress state Tcorrected
              ! ------------------------------------
              call proj_bounding(-T,Tcorrected,props(1),pmin,i_proj) !uses mechanical sign convention
              
              T      = -Tcorrected !back to geotechnical sign convention
              
              stress = map2stress(T,ntens)
            
          endif

          !Update total time
          ! ------------------------------------ 
          mTime                = mTime  + mDtime
          
          !New time increment (increase computational speed by increasing the second factor below)
          ! ------------------------------------ 
          mDTime               = mDtime*min(0.9d0*(errorTol/Error)
     &                           **(0.5d0),3.0d0)
          
          !To make sure that not too much strain is considered check if Dtime may be to large
          ! ------------------------------------ 
          if(mTime+mDTime>1.0d0) then
              mDTime           = 1.0d0-mTime
          endif

      endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      enddo !Newton has finished
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Recover statevs
      ! ------------------------------------ 
      statev(1:1+ntens) = Mstatev(1:1+ntens) !statev(1)=epor, statev(2:1+ntens)=intergranular strain(1:ntens)
      
      !Add viscous stress stress_visc
      ! ------------------------------------
      if(useVisco) then
          
        !Get full tensor notations
        ! ------------------------------------
        T           = map2T(stress,ntens)

        !Compute stress_visc and CMATRIX, GEOTECHNICAL sign convention (dstran is negative for tension)
        ! ------------------------------------
        call PhantomVisco(T,ntens,dstran,propsVisco,npropsVisco,
     &                          stress_visc,CMATRIX)
        
        if (dtime<1.0d-8) dtime=1.0d0
        
        stress_visc = stress_visc/dtime
        CMATRIX     = CMATRIX/dtime
        
        stress      = stress + stress_visc
        DDSDDE      = DDSDDE + CMATRIX
          
        !Save viscous stress for beginning of next increment
        ! ------------------------------------
        statev(2+ntens:1+ntens+ntens) = -stress_visc(1:ntens)      !viscosity, now MECHANICAL sign convention to match stress sign
      endif

      ! Mean pressure
      ! ------------------------------------
      statev(3+ntens+ntens) = (stress(1)+stress(2)+stress(3))/3.0d0 
      
      statev(4+ntens+ntens) = i_proj
      
      statev(5+ntens+ntens:5+3*ntens) = stress
      
      ! Back to MECHANICAL sign convention
      ! ------------------------------------
      stress = - stress
	  
      if(statev(1)>1.05*props(6)) statev(1)=1.05*props(6)
      if(statev(1)<0.95*props(5)) statev(1)=0.95*props(5)
      return
      end subroutine
      
      !===============================================================================================  
      ! openGEOFEA
      !===============================================================================================  
      !
      ! SUBROUTINE: hypoStiff                                                               
      !---------------                                                                    
      !> @author Patrick Staubach, patrick.staubach@yahoo.de                                       
      !                                                                                   
      ! DESCRIPTION:                                                                      
      !> @brief Returns the hypoplastic stiffness Mstiff
      !                                                                         
      ! REVISION HISTORY                                                        
      !> @date 20.05.2018 - Initial version               
      !---------------                                                          
      !                                                                         
      !> @param[in] T         : stress in GEOTECHNICAL sign convention and full tensor notation
      !> @param[in] D         : strain increment in GEOTECHNICAL sign convention and full tensor notation
      !> @param[in] props     : material parameters as defined by the user (see below for list of props)
      !> @param[in] statev    : state variables (see below for list of statev)
      !> @param[out] DDSDDE   : Hypoplastic stiffness in voigt-notation
      !> @param[out] doth     : increment of intergranular strain in full tensor notation
      !> @param[out] Mstiff   : Hypoplastic stiffness in voigt-notation in full tensor notation
      !===============================================================================================  
      subroutine hypoStiff(T,ntens,props,nprops,D,statev,
     &                     nstatv,DDSDDE,doth,noel)

      use tools 
      implicit none

      integer,intent(in)    :: ntens,nprops,nstatv,noel
      real(8),intent(in)    :: props(nprops),D(3,3)
      real(8),intent(in)    :: statev(nstatv)
      real(8),intent(in)    :: T(3,3)
      
      real(8),intent(out)   :: DDSDDE(ntens,ntens),doth(3,3)

      real(8),parameter     :: zero=0.0d0,onethird=0.3333333333333d0
      
      !Intern, needed for vHP hypoplastic model
      ! ------------------------------------
      real(8) :: hatT(3,3),devT(3,3),NN(3,3),Stressp,trdevT2,epor
      real(8) :: LLhat(3,3,3,3),LLhat2(3,3,3,3),Mstiff(3,3,3,3)
      real(8) :: fb,fe,fd,fs,aphi,edi,eci,eii,bigF,psi,cos3theta
      real(8) :: mtd,med,fdq,eporII,stiff_e,determinate,elast
      real(8) :: Tcorrected(3,3)
      ! ------------------------------------
      
      !Intern, needed for intergranular strain concept
      ! ------------------------------------
      real(8) :: h_inter_strain(3,3)
      real(8) :: h_inter_strainD(3,3),rho,dirh
	  
      
      !Material properties as defined in the input order
      ! ------------------------------------
      type MatProps
      real(8) ::
     & phi,
     & nu, 
     & bauerHs,
     & bauerN,
     & ed0,
     & ec0,
     & ei0,
     & alpha,
     & beta,
     & m_T,
     & m_R,
     & Rmax,
     & betaR,
     & Chi, 
     & Kw 
      end type MatProps
      
      type(MatProps) :: MatP
      
      !Material properties as defined in the input  file (same order as Niemunis UMAT)
      ! ------------------------------------
      MatP%phi      = props(1)        ! friction angle (degree)
      MatP%nu       = props(2)        ! poisson's ratio for increased shear stiffness
      MatP%bauerHs  = props(3)        ! granular hardness
      MatP%bauerN   = props(4)        ! exponent in Bauer's compression rule
      MatP%ed0      = props(5)        ! min. void ratio at p=0
      MatP%ec0      = props(6)        ! crit. void ratio at p=0
      MatP%ei0      = props(7)        ! max void ratio at p=0 (isotropic compresison)
      MatP%alpha    = props(8)        ! baro-pykno exponent for  f_d
      MatP%beta     = props(9)        ! baro-pykno exponent for f_e
      MatP%m_T      = props(10)       ! stiffness increase 90 deg of the previous strain path
      MatP%m_R      = props(11)       ! stiffness increase 180 deg of the previous strain path
      MatP%Rmax     = props(12)       ! Maximum radius of the elastic integranular strain
      MatP%betaR    = props(13)       ! Exponent
      MatP%Chi      = props(14)       ! Exponent
      MatP%Kw       = props(15)       ! Bulk modulus of pore water for ideal undrained simulations
      
      epor          = statev(1)       ! Void ratio
	  
      !intergranular strain [statev(2:1+ntens)]
      ! ------------------------------------
      h_inter_strain = zero
      h_inter_strain = map2D(statev(2:1+ntens),ntens)
      
      
      !Mean stress
      ! ------------------------------------
      Stressp = tr(T)/3.0d0
      Tcorrected = T
      if (Stressp < 3.0d0) then
        Tcorrected(1,1)   = T(1,1)+3.0d0/3.0d0
        Tcorrected(2,2)   = T(2,2)+3.0d0/3.0d0
        Tcorrected(3,3)   = T(3,3)+3.0d0/3.0d0
        Stressp           = 3.0d0
      elseif(Stressp > 1d3) then
        Stressp           = 1d3
      endif

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !This box contains som e material specific definitions as defined in the habilitation of Niemunis 
      ! ------------------------------------

      edi = MatP%ed0*exp(-(3.0d0*Stressp/MatP%bauerHs)**MatP%bauerN)
      eci = MatP%ec0*exp(-(3.0d0*Stressp/MatP%bauerHs)**MatP%bauerN)
      eii = MatP%ei0*exp(-(3.0d0*Stressp/MatP%bauerHs)**MatP%bauerN)
      
      !In case of epor>eii take special treatment into account according to Niemunis UMAT
      ! ------------------------------------
      eporII     = 1.5d0*eii 
      stiff_e    = 1.0d0
      
      if((epor .gt. eii).and. (epor .lt. eporII))  then        
        stiff_e  = (1.0d0-((epor-eii)/(eporII-eii)))
        epor     = eii
      endif
      
      if(epor .ge. eporII) then      
          
        stiff_e =0.0d0
        epor    = eii
        
      endif
      

      aphi = sq3*(3.0d0-dsin(MatP%phi))/(sq2*2.0d0*dsin(MatP%phi))    ! (d-) operator means that the values given and returned are of type real(8)
      
      fb   = MatP%bauerHs/(MatP%bauerN)/(MatP%ec0/MatP%ei0)**MatP%beta 
     &     * (1.0d0+eii)/eii*(3.0d0*Stressp/MatP%bauerHs)
     &     **(1.0d0-MatP%bauerN)/(3.0d0+aphi**2.0d0-aphi*sq3*((MatP%ei0
     &     - MatP%ed0)/(MatP%ec0-MatP%ed0))**MatP%alpha)
      

      fe   = (eci/epor)**MatP%beta 
      
      if (epor .lt. edi) then
          fd    = -((edi-epor)/(eci-edi))**MatP%alpha
      else      
          fd    = ((epor-edi)/(eci-edi))**MatP%alpha
      endif
      

      if(fd.lt.1.0d0) then
        mtd     = -edi/MatP%bauerHs*MatP%bauerN
     &          * (3.0d0*Stressp/MatP%bauerHs)**(MatP%bauerN-1.0d0)
                
        med     = 1.0d0
        fdq     = -(med*sq3*(1.0d0+epor)
     &          + mtd*fb*fe*sq3*(3.0d0+aphi**2))/(mtd*3.0d0*fb*fe*aphi)
        fd      = fd + (1-fd)**4*fdq  
      endif
      
      hatT      = zero
      hatT      = hated(Tcorrected)

      fs        = 1.0d0/(hatT.xx.hatT) ! fs has originally been introduced by Niemunis and is not used in the commonly hypoplasticity
      
      devT      = zero
      devT      = hatT-onethird*delta
       
      trdevT2   = zero
      trdevT2   = tr(matmul(devT,devT))
      
      cos3theta = 1.0d0
      if (trdevT2.gt.1.d-10) then
          cos3theta = -sq6*tr(matmul(matmul(devT,devT),devT))
     &              / (trdevT2)**1.5d0
      endif
      
      if (cos3theta.gt.1.0d0)  cos3theta =  1.0d0
      if (cos3theta.lt.-1.0d0) cos3theta = -1.0d0
      
      psi   = norm2(dev(hatT))*sq3
      
      If(abs(2.0d0+sq2*psi*cos3theta)<1d-8) then
        bigF  = dsqrt(abs(0.125d0*psi**2.0d0+(2.0d0-psi**2)     !The abs() was added as otherwise NAN occurs in some specific stress state
     &        / (1d-8)))-psi/(sq2*2.0d0)     !This procedure was applied by Niemunis as well (function GET_FM, computation FM)
      else
        bigF  = dsqrt(abs(0.125d0*psi**2.0d0+(2.0d0-psi**2)     !The abs() was added as otherwise NAN occurs in some specific stress state
     &        / (2.0d0+sq2*psi*cos3theta)))-psi/(sq2*2.0d0)     !This procedure was applied by Niemunis as well (function GET_FM, computation FM)
      endif
      
      LLhat = zero
      LLhat = stiff_e*fs*fb*fe*
     &        (bigF**2*Jdelta+aphi**2.0d0*(hatT.out.hatT))    
      
      NN    = zero
      NN    = stiff_e*fs*fb*fe*fd*bigF*aphi*(hatT+devT)
      

   !   !increased shear stiffness (comment next 7 lines out if it is not wished, code works fine without)
   !   ! ------------------------------------
   !  LLhat2 = zero
   !  LLhat2 = LLhat                                                    
   ! &       + ((1.0d0 + aphi**2/3.0d0 + aphi/sq3)*(1.0d0-2.0d0*MatP%nu)
   ! &       / (1.0d0+MatP%nu)-1.0d0)
   ! &       * (Jdelta - (delta.out.delta)/3.0d0)*fs*fb*fe
   !  
   !  NN     = LLhat2 .xx. (Inv(LLhat,ok) .xx. NN)                      
   !  LLhat  = LLhat2
      
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !Get intergranular strain (explicit); implicit-> call intStiffimpl (same args)
      ! ------------------------------------
      call intStiffimpl(h_inter_strain,D,MatP%Rmax,MatP%betaR,
     &               doth,dirh,rho,h_inter_strainD)
      !Compute Mstiff in dependence of the direction of integranular strain
      ! ------------------------------------
      Mstiff = zero

      if (dirh>zero) then
          Mstiff=(rho**MatP%Chi*MatP%m_T+(1.0d0-rho**MatP%Chi)       
     &          *MatP%m_R)* LLhat + rho**MatP%Chi*(1.0d0-MatP%m_T)
     &          *(LLhat.xx.(h_inter_strainD.out.h_inter_strainD))  
     &          + rho**MatP%Chi*(NN.out.h_inter_strainD)  
      else
          
          Mstiff=(rho**MatP%Chi*MatP%m_T+(1.0d0-rho**MatP%Chi)
     &          *MatP%m_R)* LLhat + rho**MatP%Chi*(MatP%m_R-MatP%m_T)
     &          *(LLhat.xx.(h_inter_strainD .out.h_inter_strainD))  
      endif
      
      
      !Obtain DDSDDE by mapping Mstiff into voigt notation
      ! ------------------------------------
      DDSDDE = zero
      DDSDDE = map2ddsdde(Mstiff,ntens)
	  
       !determinate=get_det(ddsdde,ntens)
       !if (determinate.le.0.0d0) then
       !    elast=(ddsdde(1,1)+ddsdde(2,2)+ddsdde(3,3))/20.0d0
       !    do while (determinate.le.0.0d0)
       !    do i=1,ntens 
       !    ddsdde(i,i)=ddsdde(i,i)+elast
       !    enddo
       !    determinate=get_det(ddsdde,ntens)
       !    enddo
       !endif
      
      if(any(isnan(DDSDDE))) then
           write(*,*) '-----------------------------------'
           write(*,*) '--------$HPP_UMAT_Staubach$--------'
           write(*,*) '-----------------------------------'
           write(*,*) 'NaN in DDSDDE'
               write(*,*) '-----------------------------------'
               write(*,*) '--------$HPP_UMAT_Staubach$--------'
               write(*,*) '-----------------------------------'
               write(*,*) 'Substep didnt converge'
               write(*,*) '-> Needed Time increment is too small'
               write(*,*) 'Element=',NOEL
               write(*,*) 'The calculation has been terminated'
               write(*,*) 'Current State (2D case):'
               write(*,*) '-----------------------------------'
               write(*,*) 'STRESS11:',T(1,1),'STRESS22:',T(2,2)
               write(*,*) 'STRESS33:',T(3,3),'STRESS12:',T(1,2)
               write(*,*) 'Trace Stress:',tr(T)
               write(*,*) 'Epor:',statev(1)
               write(*,*) 'Strain (11,22,33,12):'
               write(*,*)  D(1:3,1:3) 
               write(*,*) 'Intergranular strain (11,22,33,12):'
               write(*,*)  statev(1+1:1+ntens) 
               write(*,*) 'Projection on bounding surface:'
               write(*,*) '-----------------------------------'
               write(*,*) 'Geotechnical sign convention!'
      endif
      
      return
      end subroutine

      !===============================================================================================  
      ! openGEOFEA
      !===============================================================================================  
      !
      ! SUBROUTINE: intStiff                                                               
      !---------------                                                                    
      !> @author Patrick Staubach, patrick.staubach@yahoo.de                                       
      !                                                                                   
      ! DESCRIPTION:                                                                      
      !> @brief Returns the updated intergranular state variables using an EXPLICIT procedure
      !                                                                         
      ! REVISION HISTORY                                                        
      !> @date 20.05.2018 - Initial version               
      !---------------                                                          
      !                                                                         
      !> @param[in] h_inter_strain   : intergranular strain tensor 
      !> @param[in] D                : strain sub increment in GEOTECHNICAL sign convention and full tensor notation
      !> @param[in] Rmax             : Max radius of intergranular strain
      !> @param[in] betaR            : Intergranular interpolator
      !> @param[out] doth            : increment of intergranular strain in full tensor notation
      !> @param[out] dirh            : Direction of intergranular strain
      !> @param[out] rho             : Mobilization of intergranular strain
      !> @param[out] h_inter_strainD : Normalized intergranular strain
      !===============================================================================================  
      subroutine intStiff(h_inter_strain,D,Rmax,betaR,doth,
     &                    dirh,rho,h_inter_strainD)
      use tools 
      implicit none

      real(8),intent(in) :: h_inter_strain(3,3),D(3,3)
      real(8),intent(in) :: betaR,Rmax
      real(8),intent(out):: doth(3,3),dirh,rho,h_inter_strainD(3,3)
      
      real(8) :: h_inter_strainNorm,aux(3,3,3,3)
      
      real(8), parameter:: zero=0.0d0

      h_inter_strainD    = zero
      h_inter_strainD    = normalized33(h_inter_strain)

      ! Inside this routine the MECHANICAL sign convention is used for D
      ! ------------------------------------
      dirh = zero
      dirh = (h_inter_strainD.xx.-D)
      
      rho = zero
      rho = MAX(h_inter_strainNorm/Rmax,1.0d-10)         
      
      if (dirh>zero) then
          
        h_inter_strainNorm = zero
        h_inter_strainNorm = norm2(h_inter_strain)
      
        
        aux  = zero
        aux  = Jdelta -(h_inter_strainD.out.h_inter_strainD)
     &           * rho**betaR
            
        doth = zero
        doth = (-D.xx.aux)
        
      else
      
        doth = zero
        doth = -D
        
      endif
      
      return
      end subroutine

      !===============================================================================================  
      ! openGEOFEA
      !===============================================================================================  
      !
      ! SUBROUTINE: intStiffimpl                                                             
      !---------------                                                                    
      !> @author Patrick Staubach, patrick.staubach@yahoo.de                                       
      !                                                                                   
      ! DESCRIPTION:                                                                      
      !> @brief Returns the updated intergranular state variables using an IMPLICIT procedure
      !> @brief Based on the subroutine of Niemunis
      !                                                                         
      ! REVISION HISTORY                                                        
      !> @date 20.05.2018 - Initial version               
      !---------------                                                          
      !                                                                         
      !> @param[in] h_inter_strain   : intergranular strain tensor 
      !> @param[in] D                : strain sub increment in GEOTECHNICAL sign convention and full tensor notation
      !> @param[in] Rmax             : Max radius of intergranular strain
      !> @param[in] betaR            : Intergranular interpolator
      !> @param[out] doth            : increment of intergranular strain in full tensor notation
      !> @param[out] dirh            : Dircetion of intergranular strain
      !> @param[out] rho             : Mobilization of intergranular strain
      !> @param[out] h_inter_strainD : Normalized intergranular strain
      !===============================================================================================  
      subroutine intStiffimpl(h_inter_strain,D,Rmax,betaR,doth,
     &                    dirh,rho,h_inter_strainD)
      use tools  
      implicit none

      real(8),intent(in) :: h_inter_strain(3,3),D(3,3)
      real(8),intent(in) :: betaR,Rmax
      real(8),intent(out):: doth(3,3),dirh,rho,h_inter_strainD(3,3)
      
      real(8) :: h_inter_strainNorm
      real(8) :: rhobetax
      real(8) :: c1
      real(8) :: h_inter_strain2(3,3)
      
      real(8), dimension (1:3,1:3,1:3,1:3) :: U,h,hi,hiU
      real(8), dimension (1:3,1:3,1:3,1:3,1:3,1:3) ::   dUijkldamn
      
      real(8), parameter:: zero=0.0d0
      integer :: i,j,k,l,m,n
      
      h_inter_strainD   = zero
      h_inter_strainD   = normalized33(h_inter_strain)
      
      dirh = zero
      dirh = (h_inter_strainD.xx.-D)
      

      if (dirh>zero) then
          
        h_inter_strainNorm = zero
        h_inter_strainNorm = norm2(h_inter_strain)
      
        rho=h_inter_strainNorm/Rmax 
        
        rhobetax          = zero
        c1                = zero

        
        if (h_inter_strainNorm > zero) then
            
           rhobetax           = rho**betaR
           c1                 = 1.0d0/h_inter_strainNorm*rhobetax
           
        endif
        
      U = zero  
      U = Jdelta  - rhobetax * (h_inter_strainD .out. h_inter_strainD)
      
      dUijkldamn = zero
      
      do 20 i=1,3
      do 20 j=1,3
      do 20 k=1,3
      do 20 l=1,3
      do 20 m=1,3
      do 20 n=1,3
          
  20  dUijkldamn(i,j,k,l,m,n)=-c1*((betaR-2.0d0)
     &       *h_inter_strain(i,j)*h_inter_strain(k,l)
     &       *h_inter_strain(m,n)+h_inter_strainNorm*h_inter_strain(k,l)
     &       *delta(i,m)*delta(j,n)+h_inter_strainNorm
     &       *h_inter_strain(i,j)*delta(k,m)*delta(l,n))
      
      h = zero
      h = Jdelta - ( dUijkldamn .xx. -D)*0.5d0
      hi = zero
      hi = inv3333(h,ok)

      hiU  = zero
      hiU  = hi .xx. U
      
      doth = zero
      doth = hiU .xx. -D 
      
      else
          
          
        h_inter_strainNorm = zero
        h_inter_strainNorm = norm2(h_inter_strain)
        
        rho                = h_inter_strainNorm/Rmax 
        
        doth               = zero
        doth               = -D
        
      endif
      
       h_inter_strain2      = h_inter_strain + doth
       
       h_inter_strainNorm   = zero
       h_inter_strainNorm   = norm2(h_inter_strain2)
       
       rho                  = h_inter_strainNorm/Rmax 
          
       if (rho.gt.1.0d0) then
           doth             =  doth + h_inter_strain2/rho 
     &                      - h_inter_strain2
           h_inter_strain2  =  h_inter_strain2/rho 
           rho              =  1.0d0
                            
           h_inter_strainD  = normalized33(h_inter_strain2)
       endif
      
      return
      end subroutine
      !===============================================================================================  
      ! openGEOFEA
      !===============================================================================================  
      !
      ! SUBROUTINE: PhantomVisco                                                         
      !---------------                                                                    
      !> @author Patrick Staubach, patrick.staubach@yahoo.de                                       
      !                                                                                   
      ! DESCRIPTION:                                                                      
      !> @brief Returns viscous stress and viscous stiffness based on a linear viscosity model
      !                                                                         
      ! REVISION HISTORY                                                        
      !> @date 20.05.2018 - Initial version               
      !---------------                                                          
      !> @param[in] T                : stress in GEOTECHNICAL sign convention and full tensor notation                                             
      !> @param[in] dstran           : strain increment in GEOTECHNICAL sign convention and voigt notation
      !> @param[in] propsVisco       : Props of viscosity as defined in subroutine umat
      !>                               propsVisco(1)   = 0.5d0   minimum of lower 
      !>                               propsVisco(2)   = 2.0d0   maximum of lower 
      !>                               propsVisco(3)   = 2.0d0   maximum of lambda
      !>                               propsVisco(4)   = 0.2d0   minimum of lambda
      !>                               propsVisco(5)   = 2.0d0   maximum of mou
      !>                               propsVisco(6)   = 0.2d0   minimum of mou
      !> @param[in] npropsVisco      : Number of props of viscosity
      !> @param[out] stress_visc     : Viscous stress
      !> @param[out] CMATRIX         : Viscous stiffness
      !===============================================================================================  
      
      subroutine PhantomVisco(T,ntens,dstran,propsVisco,npropsVisco,
     &                        stress_visc,CMATRIX)

      use tools  
      implicit none

      integer,intent(in)    :: ntens,npropsVisco
      real(8),intent(in)    :: T(3,3) 
      real(8),intent(in)    :: propsVisco(npropsVisco),dstran(ntens)
      
      real(8),intent(out)   :: stress_visc(ntens),CMATRIX(ntens,ntens)
      
      real(8)               :: Stressp,lambda,mou,ALAMDA,AMU
      
      real(8), parameter    :: zero=0.0d0
      
      Stressp = tr(T)/3.0d0 !positive
      
      !propsVisco
      ! ------------------------------------
      !propsVisco(1) minimum of lower stress point (about 0.5 kPa, positive!)
      !propsVisco(2) maximum of lower stress point (about 2.0 kPa, positive!)
      !propsVisco(3) maximum of lambda, volumetric part of viscosity
      !propsVisco(4) minimum of lambda, volumetric part of viscosity
      !propsVisco(5) maximum of mou
      !propsVisco(6) minimum of mou

      ALAMDA  = propsVisco(4)
      AMU     = propsVisco(6)
      
      !Check for small stresses that lie between propsVisco(1) and propsVisco(2)
      ! ------------------------------------
      if ((Stressp < propsVisco(2)) .and. (Stressp > propsVisco(1)))then
          
      !Linear interpolation of stiffness
      ! ------------------------------------
          lambda = (propsVisco(3)-propsVisco(4))
     &           / (propsVisco(2)-propsVisco(1))
          
          mou    = (propsVisco(5)-propsVisco(6))
     &           / (propsVisco(2)-propsVisco(1))
                    
          ALAMDA = propsVisco(4) + lambda*(propsVisco(2)-Stressp)   !minimum of lambda + interpolation
          AMU    = propsVisco(6) + mou*(propsVisco(2)-Stressp)      !minimum of mou + interpolation
          
      endif
      
      
      !Check for very small stresses (p<propsVisco(1))
      ! ------------------------------------
      if (Stressp < propsVisco(1)) then  
          
      !Full stiffness
      ! ------------------------------------
      ALAMDA  = propsVisco(3) !maximum of lambda 
      AMU     = propsVisco(5) !maximum of mou 
      
      endif 
      
      !Constitutive stiffness
      ! ------------------------------------
      CMATRIX                   = zero
      
      CMATRIX(1:3,1:3)          = ALAMDA+2.0d0*AMU
      CMATRIX(4:ntens,4:ntens)  = AMU
      CMATRIX(1,2)              = ALAMDA
      CMATRIX(1,3)              = ALAMDA
      CMATRIX(2,3)              = ALAMDA
      CMATRIX(2,1)              = CMATRIX(1,2)
      CMATRIX(3,1)              = CMATRIX(1,3)
      CMATRIX(3,2)              = CMATRIX(2,3)	  
      
      !Compute viscous stress
      ! ------------------------------------
      stress_visc               = zero
      stress_visc               = matmul(CMATRIX,dstran)
          
      
      return
      end subroutine
      
      !===============================================================================================  
      ! openGEOFEA
      !===============================================================================================  
      !
      ! SUBROUTINE: proj_bounding                    
      !---------------                                                                    
      !> @author Patrick Staubach, patrick.staubach@yahoo.de                                       
      !                                                                                   
      ! DESCRIPTION:                                                                      
      !> @brief  Projection of the stress tensor on a bounding surface
      !> @brief First projection: if tr(T)/3 > sig0 with a given sig0 < 0,
      !> @brief a hydrostatic stress is added to T to make tr(T)/3 = sig0.
      !> @brief Second projection: if T is outside the Matsuoka-Nakai surface
      !> @brief with a given phi [rad], T is projected on the surface
      !> @brief in the direction normal to the hydrostatic axis.
      !> @brief The output integer i_proj shows what projection has been made:
      !> @brief 0 - no projection; 
      !> @brief 1 - only the first projection;
      !> @brief 2 - the first and then the second projection;
      !> @brief 3 - only the second projection.                                               
      !> REVISION HISTORY       
      !> 
      !> @date 20.05.2018 - Initial version, taken from niemunis            
      !---------------                                                          
      !> @param[in] T                : stress in MECHANICAL sign convention and full tensor notation                                             
      !> @param[in] phi              : friction angle
      !> @param[in] sig0             : Given minimum stress, defined by user
      !> @param[in] npropsVisco      : Number of props of viscosity
      !> @param[out] i_proj          : integer i_proj shows what projection has been made
      !===============================================================================================  
      
      SUBROUTINE proj_bounding(T,Tcorrected,phi,sig0,i_proj)

      implicit none
      
      real(8), intent(in)    :: phi,sig0
      real(8), intent(in)    :: T(3,3)
      real(8), intent(out)   :: Tcorrected(3,3)
      integer, intent(out)   :: i_proj
      
      
      real(8)                :: tr_sig,rc,w,w1,w2,w3,dummy2
      real(8)                :: s11,s22,s33,s12,s13,s23
      real(8)                :: norm,det,r,dummy,p,q,rr,root
      
      real(8), parameter     :: zero=0.0d0,pi=3.141592653d0
      
      i_proj=0
      

      tr_sig     = T(1,1)+T(2,2)+T(3,3)
      
      Tcorrected = zero
      Tcorrected = T
      
      !First projection, add hydrostatic stress such that p>0
      ! ------------------------------------
      if (tr_sig/3.0d0 .gt.sig0) then
        Tcorrected(1,1) = T(1,1)+sig0-tr_sig/3.0d0
        Tcorrected(2,2) = T(2,2)+sig0-tr_sig/3.0d0
        Tcorrected(3,3) = T(3,3)+sig0-tr_sig/3.0d0
        
        tr_sig          = 3.0d0*sig0
        i_proj          = 1
      endif

      w=dsin(1.4d0*phi)**2
      rc=(1.0d0-w)/(9.0d0-w)
      
      s11  = Tcorrected(1,1)/tr_sig
      s22  = Tcorrected(2,2)/tr_sig
      s33  = Tcorrected(3,3)/tr_sig
      s12  = Tcorrected(1,2)/tr_sig
      s13  = Tcorrected(1,3)/tr_sig
      s23  = Tcorrected(2,3)/tr_sig
      s11  = s11-1.0d0/3.0d0
      s22  = s22-1.0d0/3.0d0
      s33  = s33-1.0d0/3.0d0
      
      norm = dsqrt(s11*s11+s22*s22+s33*s33+
     &        2.0d0*s12*s12+2.0d0*s23*s23+2.0d0*s13*s13)
      
      if (norm/abs(tr_sig).lt.1.0d-3) return
      s11 = s11/norm
      s22 = s22/norm
      s33 = s33/norm
      s12 = s12/norm
      s13 = s13/norm
      s23 = s23/norm
      det = s11*s22*s33+2.0d0*s12*s23*s13-s13*s13*s22-
     &    s23*s23*s11-s12*s12*s33
      r=0.5d0*rc-1.0d0/6.0d0
      dummy=1.0d0/27.0d0-rc/3.0d0
      if (abs(det).lt.1.0d-6*abs(r)) then
          dummy2 = -dummy/r
          root   = 0.0d0
          if (dummy2>0.0d0) root=dsqrt(dummy2)
        go to 10
      END if
      r=r/det
      dummy=dummy/det
      p=-r*r/3.0d0
      q=2.0d0*r*r*r/27.0d0+dummy
      if (q.ge.0.0d0) then 
        rr=dsqrt(abs(p)/3.0d0)
      else 
        rr=-dsqrt(abs(p)/3.0d0)
      END if
      w=dacos(q/(2.0d0*rr*rr*rr))
      w1 = -2.0d0*rr*dcos(w/3.0d0)-r/3.0d0
      w2 = -2.0d0*rr*dcos(w/3.0d0+2.0d0*pi/3.0d0)-r/3.0d0
      w3 = -2.0d0*rr*dcos(w/3.0d0+4.0d0*pi/3.0d0)-r/3.0d0
      w  = zero
      if (w1.gt.0.0d0.and.1.0d0/w1.gt.w) w=1.0d0/w1
      if (w2.gt.0.0d0.and.1.0d0/w2.gt.w) w=1.0d0/w2
      if (w3.gt.0.0d0.and.1.0d0/w3.gt.w) w=1.0d0/w3
      root=1.0d0/w
10    if (root.ge.norm) return
      Tcorrected(1,1) = (s11*root+1.0d0/3.0d0)*tr_sig
      Tcorrected(2,2) = (s22*root+1.0d0/3.0d0)*tr_sig
      Tcorrected(3,3) = (s33*root+1.0d0/3.0d0)*tr_sig
      Tcorrected(1,2) = s12*root*tr_sig
      Tcorrected(1,3) = s13*root*tr_sig
      Tcorrected(2,3) = s23*root*tr_sig
      Tcorrected(2,1) = Tcorrected(1,2)
      Tcorrected(3,1) = Tcorrected(1,3)
      Tcorrected(3,2) = Tcorrected(2,3)
      if (i_proj.eq.0) then
        i_proj=3
      else
        i_proj=2
      END if
	  

      end subroutine
