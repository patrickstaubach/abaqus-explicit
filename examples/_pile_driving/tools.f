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
! Module: tools
!
!> @author Patrick Staubach, patrick.staubach@yahoo.de
!          Bauhaus University Weimar, Ruhr-University Bochum
!
! DESCRIPTION:
! @briefT Tensor operations
! @briefT Originally from A. Niemunis
!
! REVISION HISTORY
!> @date 02.02.2021 - Initial version
!=======================================================================================================
      module tools   !from A. Niemunis  (KIT Karlsruhe)                                  
      implicit none
	  save              
      integer :: ixx,jxx
                                                                        
      real(8),parameter  :: sq2=1.4142135623730950488d0               
      real(8),parameter  :: sq3=1.7320508075688772935d0      
      real(8),parameter  :: sq6=2.4494897427831780982d0      
      real(8),parameter  :: sq23= 0.81649658092772603273d0  
      real(8),parameter  :: pi=3.141592653589793238462643d0  
      real(8),parameter  :: tercja=0.333333333333333333333333d0   
      logical :: ok                                                     
      real(8), parameter,dimension(1:3,1:3) ::
     &                  delta =RESHAPE((/ 1,0,0,0,1,0,0,0,1/),(/3,3/))  
      real(8), parameter,dimension(1:3,1:3,1:3,1:3) ::
     &             Jdelta=RESHAPE( (/ (1, (0, ixx=1,9) ,jxx=1,8),1 /),  
     &                                 (/3,3,3,3/))
      real(8), parameter,dimension(1:3,1:3,1:3,1:3) ::
     &         Idelta= RESHAPE((/ 2,0,0,0,0,0,0,0,0,0,               !  
     &                  1,0,1,0,0,0,0,0,0,0,     1,0,0,0,1,0,0,0,1,0,
     &                  1,0,0,0,0,0,0,0,0,0,     2,0,0,0,0,0,0,0,0,0,
     &                  1,0,1,0,0,0,1,0,0,0,     1,0,0,0,0,0,0,0,1,0,
     &                  1,0,0,0,0,0,0,0,0,0,     2/),(/3,3,3,3/)) /2.0d0
      real(8), parameter,dimension(1:3,1:3,1:3,1:3) ::
     &                 pure_dev_old=RESHAPE((/
     &  2,0,0,0,-1,0,0,0,-1,  0,3,0,0,0,0,0,0, 0,  0,0,3,0, 0,0,0,0,0,
     &  0,0,0,3, 0,0,0,0, 0, -1,0,0,0,2,0,0,0,-1,  0,0,0,0, 0,3,0,0,0,
     &  0,0,0,0, 0,0,3,0, 0,  0,0,0,0,0,0,0,3, 0, -1,0,0,0,-1,0,0,0,2
     &                 /),(/3,3,3,3/))/3.0d0
      real(8), parameter,dimension(1:3,1:3,1:3,1:3) ::
     &                  pure_dev=RESHAPE((/
     & 4,0,0,0,-2,0,0,0,-2,   0,3,0,3,0,0,0,0, 0,   0,0,3,0,0,0,3,0,0,
     & 0,3,0,3, 0,0,0,0, 0,  -2,0,0,0,4,0,0,0,-2,   0,0,0,0,0,3,0,3,0,
     & 0,0,3,0, 0,0,3,0, 0,   0,0,0,0,0,3,0,3,0,   -2,0,0,0,-2,0,0,0,4
     &                 /),(/3,3,3,3/))/6.0d0


      real(8), parameter,dimension(81,81) ::                          
     &            JJdelta=RESHAPE((/ (1,(0, ixx=1,81),jxx=1,80),1 /),
     &                                    (/81,81/))                  
      integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),       
     &                                     j6=(/ 1,2,3,2,3,3/)
      integer, parameter,dimension(1:3,1:3)::                         
     &                   ij= RESHAPE((/1,5,7,4,2,9,6,8,3/),(/3,3/))   
     
      integer, parameter,dimension(1:3,1:3)::                         
     &                   ij6= RESHAPE((/1,4,5,4,2,6,5,6,3/),(/3,3/))  

      integer, parameter,dimension(1:3,1:3,1:3,1:3)::                 
     &                index81= RESHAPE((/ (ixx, ixx=1,81)/),(/3,3,3,3/))
     
      interface operator(.out.)                                   
        module procedure outmal,outmal3
       end interface
       
      interface operator(.xx.)                                        
         module procedure mal,mal2,mal3,mal4,mal5,mal6
      end interface
       
      contains
       
      function map2stran(a,ntens)    !from A. Niemunis                                   
        implicit none                                                
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stran
        integer :: i
        map2stran(1)=a(1,1)
        map2stran(2)=a(2,2)
        map2stran(3)=a(3,3)
        do i=4,ntens
        map2stran(i) = a(i6(i),j6(i)) +   a(j6(i),i6(i))                          
        enddo                                                           
      end function map2stran

      function map2D(a,ntens)   !from A. Niemunis
        implicit none                                                   
        real(8),  dimension(1:3,1:3) :: map2D
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2D=0
        map2D(1,1) = a(1)
        map2D(2,2) = a(2)
        map2D(3,3) = a(3)
        do i=4,ntens
         map2D(i6(i),j6(i))=a(i)/2.0d0                                   
         map2D(j6(i),i6(i))=a(i)/2.0d0                                  
        enddo
      end function map2D

      function map2stress(a,ntens)   !from A. Niemunis                                
        implicit none                                                
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stress
        integer :: i
        do i=1,ntens
          map2stress(i) = a(i6(i),j6(i))                                
        enddo                                                           
      end function map2stress

      function map2T(a,ntens)  !from A. Niemunis                                        
        implicit none                                                  
        real(8),  dimension(1:3,1:3) :: map2T
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2T=0
        map2T(1,1) = a(1)
        map2T(2,2) = a(2)
        map2T(3,3) = a(3)
        do i=4,ntens
        map2T(i6(i),j6(i))=a(i)
        map2T(j6(i),i6(i))=a(i)
        enddo
      end function map2T

      function map2ddsdde(LL,ntens)   !from A. Niemunis                                       
        implicit none                                                   
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: LL
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens,1:ntens) :: map2ddsdde
        integer :: i,j
        do i=1,ntens
        do j=1,ntens
          if (j <= 3) map2ddsdde(i,j) = LL(i6(i),j6(i),i6(j),j6(j))
          if (j >  3) map2ddsdde(i,j) =  0.5d0*
     &     (LL(i6(i),j6(i),i6(j),j6(j))+LL(i6(i),j6(i),j6(j),i6(j)) )
        enddo
        enddo
      end function map2ddsdde

      function outmal(a,b)    !from A. Niemunis                                            
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a,b
        real(8), dimension(1:3,1:3,1:3,1:3) :: outmal
        integer :: i,j,k,l
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
            outmal(i,j,k,l) =  a(i,j)*b(k,l)
        enddo
        enddo
        enddo
        enddo
      end function outmal

      function mal(a,b)    !from A. Niemunis                                             
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a,b  
        real(8) :: mal
         mal          =  a(1,1)*b(1,1)+
     &                    a(1,2)*b(1,2)+
     &                    a(1,3)*b(1,3)+
     &                    a(2,1)*b(2,1)+
     &                    a(2,2)*b(2,2)+
     &                    a(2,3)*b(2,3)+
     &                    a(3,1)*b(3,1)+
     &                    a(3,2)*b(3,2)+
     &                    a(3,3)*b(3,3)
      end function mal

      function mal2(a,b)   !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: a
        real(8), intent(in),  dimension(1:3,1:3) :: b
        real(8), dimension(1:3,1:3):: mal2
        integer :: i,j
        do  i=1,3
        do  j=1,3
        mal2(i,j)    =  a(i,j,1,1)*b(1,1)+
     &                  a(i,j,1,2)*b(1,2)+
     &                  a(i,j,1,3)*b(1,3)+
     &                  a(i,j,2,1)*b(2,1)+
     &                  a(i,j,2,2)*b(2,2)+
     &                  a(i,j,2,3)*b(2,3)+
     &                  a(i,j,3,1)*b(3,1)+
     &                  a(i,j,3,2)*b(3,2)+
     &                  a(i,j,3,3)*b(3,3)
        enddo
        enddo
      end function mal2

      function mal3(a,b)    !from A. Niemunis
        implicit none
        real(8), intent(in),  dimension(1:3,1:3) :: a
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: b
        real(8), dimension(1:3,1:3):: mal3
        integer :: k,l
         do  k=1,3
         do  l=1,3
         mal3(k,l) =   a(1,1)*b(1,1,k,l)+
     &                 a(1,2)*b(1,2,k,l)+
     &                 a(1,3)*b(1,3,k,l)+
     &                 a(2,1)*b(2,1,k,l)+
     &                 a(2,2)*b(2,2,k,l)+
     &                 a(2,3)*b(2,3,k,l)+
     &                 a(3,1)*b(3,1,k,l)+
     &                 a(3,2)*b(3,2,k,l)+
     &                 a(3,3)*b(3,3,k,l)
         enddo
         enddo
      end function mal3

      function mal4(a,b)    !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3):: a,b
        real(8), dimension(1:3,1:3,1:3,1:3):: mal4
        integer :: i,j,k,l
         do  i=1,3
         do  j=1,3
         do  k=1,3
         do  l=1,3
         mal4(i,j,k,l)= a(i,j,1,1)*b(1,1,k,l)+
     &                 a(i,j,1,2)*b(1,2,k,l)+
     &                 a(i,j,1,3)*b(1,3,k,l)+
     &                 a(i,j,2,1)*b(2,1,k,l)+
     &                 a(i,j,2,2)*b(2,2,k,l)+
     &                 a(i,j,2,3)*b(2,3,k,l)+
     &                 a(i,j,3,1)*b(3,1,k,l)+
     &                 a(i,j,3,2)*b(3,2,k,l)+
     &                 a(i,j,3,3)*b(3,3,k,l)
         enddo
         enddo
         enddo
         enddo
      end function mal4

      function mal5(a,b)    !from A. Niemunis
        implicit none
        real(8),intent(in),dimension(1:3,1:3,1:3,1:3,1:3,1:3)::a
        real(8), intent(in), dimension(1:3,1:3):: b
        real(8), dimension(1:3,1:3,1:3,1:3):: mal5
        integer :: i,j,k,l
        do  i=1,3
        do  j=1,3
        do  k=1,3
        do  l=1,3
        mal5(i,j,k,l)=a(i,j,k,l,1,1)*b(1,1)+
     &                    a(i,j,k,l,1,2)*b(1,2)+
     &                    a(i,j,k,l,1,3)*b(1,3)+
     &                    a(i,j,k,l,2,1)*b(2,1)+
     &                    a(i,j,k,l,2,2)*b(2,2)+
     &                    a(i,j,k,l,2,3)*b(2,3)+
     &                    a(i,j,k,l,3,1)*b(3,1)+
     &                    a(i,j,k,l,3,2)*b(3,2)+
     &                    a(i,j,k,l,3,3)*b(3,3)
        enddo
        enddo
        enddo
        enddo
      end function mal5

      function mal6(a,b)    !from A. Niemunis
        implicit none
        real(8),intent(in), dimension(1:3,1:3):: a
        real(8),intent(in),dimension(1:3,1:3,1:3,1:3,1:3,1:3)::b
        real(8), dimension(1:3,1:3,1:3,1:3):: mal6
        integer :: i,j,k,l
        do  i=1,3
        do  j=1,3
        do  k=1,3
        do  l=1,3
        mal6(i,j,k,l) =   a(1,1)*b(1,1,i,j,k,l)+
     &                    a(1,2)*b(1,2,i,j,k,l)+
     &                    a(1,3)*b(1,3,i,j,k,l)+
     &                    a(2,1)*b(2,1,i,j,k,l)+
     &                    a(2,2)*b(2,2,i,j,k,l)+
     &                    a(2,3)*b(2,3,i,j,k,l)+
     &                    a(3,1)*b(3,1,i,j,k,l)+
     &                    a(3,2)*b(3,2,i,j,k,l)+
     &                    a(3,3)*b(3,3,i,j,k,l)
        enddo
        enddo
        enddo
        enddo
      end function mal6

      function outmal3(a,b)    !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: a
        real(8), intent(in), dimension(1:3,1:3)  :: b
        real(8), dimension(1:3,1:3,1:3,1:3,1:3,1:3) :: outmal3
        integer :: i,j,k,l,m,n
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
        do m=1,3
        do n=1,3
            outmal3(i,j,k,l,m,n) =  a(i,j,k,l)*b(m,n)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
      end function outmal3

       function tr(a)     !from A. Niemunis                                                  
          implicit none                                       
          real(8), intent(in), dimension(1:3,1:3)  :: a
          real(8) :: tr
          tr=a(1,1)+a(2,2)+a(3,3)
       end function tr

      function map299(a)      !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: a
        real(8),  dimension(1:9,1:9) :: map299
        integer :: i,j,k,l
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
        map299(ij(i,j),ij(k,l)) = a(i,j,k,l)
        enddo
        enddo
        enddo
        enddo
      end function map299

      function map23333(a)     !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:9,1:9)  :: a
        real(8), dimension(1:3,1:3,1:3,1:3) :: map23333
        integer :: i,j,k,l
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
        map23333(i,j,k,l) = a(ij(i,j),ij(k,l))
        enddo
        enddo
        enddo
        enddo
      end function map23333

      function inv3333(a3333,success)    !from A. Niemunis                                       
        implicit none                                                     
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3)::a3333            
        logical, intent(out)  :: success
        logical, dimension(1:9) :: deleted
        real(8), dimension(1:3,1:3,1:3,1:3) :: inv3333
        real(8), dimension(1:9,1:9) :: a,one
        real(8), dimension(1:9,1:18) :: c
        real(8) :: cc, dd
        integer :: i,j,k
        inv3333 = 0
        deleted = .FALSE.
        success = .FALSE.
        a = map299(a3333)
        do i=5,9,2                                                        
        if (  ALL(  abs(a(i-1,1:9) - a(i,1:9) )<  1.0d-6 ) ) then         
         deleted(i) = .TRUE.
         a(i,1:9)=0
         a(1:9,i)=0
         a(i,i)= max (1.0d0, a(i-1,i-1) )                                 
        endif
        enddo
        one = map299(Jdelta)
        c(1:9,1:9) = a
        c(1:9,10:18)=one
        do i=1,9                                                          
          cc = c(i,i)                                                     
          if (abs(cc).lt.1d-6) return                                     
          c(i,i) = cc-1.0d0
          do k=i+1, 18
             dd=c(i,k)/cc
             do j=1,9
               c(j,k) = c(j,k)-dd*c(j,i)
             enddo
          enddo
        enddo
        a= c(1:9,10:18)
        do i=5,9,2                                                        
        if (deleted(i)) then
         a(i,1:9)  = a(i-1,1:9)/2
         a(i-1,1:9)= a(i-1,1:9)/2
         a(1:9,i)  = a(1:9,i-1)/2
         a(1:9,i-1)= a(1:9,i-1)/2
        endif
        enddo
        inv3333=map23333( a )
        success=.TRUE.
      end function inv3333
      
      function inv33(f,success)    !from A. Niemunis                                        
        implicit none                                                 
        real(8), intent(in), dimension(1:3,1:3) :: f
        logical, intent(out)  :: success
        real(8), dimension(1:3,1:3) :: inv33
        real(8) :: detf
        inv33 = 0
        success = .FALSE.
        detf =   -f(1,3)*f(2,2)*f(3,1)+f(1,2)*f(2,3)*f(3,1)
     &           +f(1,3)*f(2,1)*f(3,2)-f(1,1)*f(2,3)*f(3,2)
     &           -f(1,2)*f(2,1)*f(3,3)+f(1,1)*f(2,2)*f(3,3)
        if(abs(detf) < tiny(detf)) return   
        success = .TRUE.
        inv33(1,1)= (-f(2,3)*f(3,2)+f(2,2)*f(3,3))/detf
        inv33(1,2)= ( f(1,3)*f(3,2)-f(1,2)*f(3,3))/detf
        inv33(1,3)= (-f(1,3)*f(2,2)+f(1,2)*f(2,3))/detf
        inv33(2,1)= ( f(2,3)*f(3,1)-f(2,1)*f(3,3))/detf
        inv33(2,2)= (-f(1,3)*f(3,1)+f(1,1)*f(3,3))/detf
        inv33(2,3)= ( f(1,3)*f(2,1)-f(1,1)*f(2,3))/detf
        inv33(3,1)= (-f(2,2)*f(3,1)+f(2,1)*f(3,2))/detf
        inv33(3,2)= ( f(1,2)*f(3,1)-f(1,1)*f(3,2))/detf
        inv33(3,3)= (-f(1,2)*f(2,1)+f(1,1)*f(2,2))/detf
      end function inv33
      
      function dev(a)    !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a
        real(8), dimension(1:3,1:3) :: dev
        real(8) :: tr3
        tr3 = (a(1,1)+a(2,2)+a(3,3))/3.0d0
        dev = a - delta * tr3
      end function dev

      function hated(a)     !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a
        real(8), dimension(1:3,1:3)  :: hated
        real(8)  :: tr
        tr =a(1,1)+a(2,2) +a(3,3)
         hated = 0.0d0
        if( abs(tr)  >  tiny(tr) )  hated = a/tr
      end function hated

      function normalized33(a)    !from A. Niemunis
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a
        real(8), dimension(1:3,1:3)  :: normalized33
        real(8)  :: sqnorm
        sqnorm =a(1,1)*a(1,1)+a(1,2)*a(1,2)+a(1,3)*a(1,3)+
     &    a(2,1)*a(2,1)+a(2,2)*a(2,2)+a(2,3)*a(2,3)+
     &    a(3,1)*a(3,1)+a(3,2)*a(3,2)+a(3,3)*a(3,3)
        normalized33=0
        if(sqnorm >  tiny(sqnorm)   ) then
            sqnorm=sqrt(sqnorm)
            normalized33=a/sqnorm
        endif
      end function normalized33

      end module
     
