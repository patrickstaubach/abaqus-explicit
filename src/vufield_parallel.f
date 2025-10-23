!=======================================================================================================
!This file is part of VUMAT_HMC_Staubach.
!
!VUMAT_HMC_Staubach is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License !as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!VUMAT_HMC_Staubach is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with VUMAT_HMC_Staubach. If not, see <https://www.gnu.org/licenses/>. 
!=======================================================================================================
!
! SUBROUTINE: VUFIELD
!
!> @author Patrick Staubach, patrick.staubach@yahoo.de
!          Bauhaus University Weimar, Ruhr-University Bochum
!
! DESCRIPTION:
!> @brief Contains the VUFIELD routine used to calculate effective interface friction according to the paper
!> @brief "Hydro-mechanically coupled CEL analyses with effective contact stresses" https://onlinelibrary.wiley.com/doi/10.1002/nag.3725 International Journal for Numerical and Analytical Methods in Geomechanics 
!
! REVISION HISTORY
!> @date 03.03.2024 - Initial version
!> @date 22.03.2024 - Support of up to 8 CPUs
!> @date 23.10.2025 - stylistic cleanup, adds z-range window (min/max)
!=======================================================================================================
      subroutine vufield(field, nblock, nfield, kfield, ncomp,
     1                   kstep, jflags, jnodeuid, time,
     2                   coords, u, v, a)
c
      include 'vaba_param.inc'
c
c ---- indices / parameters
      integer i_ufld_current,   i_ufld_increment, i_ufld_period, i_ufld_total
      integer i_ufld_coordx,    i_ufld_coordy,    i_ufld_coordz
      integer i_ufld_kinc,      i_ufld_kpass
      parameter( i_ufld_current   = 1,
     &           i_ufld_increment = 2,
     &           i_ufld_period    = 3,
     &           i_ufld_total     = 4 )
      parameter( i_ufld_coordx = 1,
     &           i_ufld_coordy = 2,
     &           i_ufld_coordz = 3 )
      parameter( i_ufld_kinc   = 1,
     &           i_ufld_kpass  = 2 )
c
c ---- arguments 
      integer nblock, nfield, kfield, ncomp, kstep
      integer jnodeuid(nblock)
      integer jflags(2)
      dimension time(4)
      dimension coords(3,nblock), u(8,nblock), v(8,nblock), a(8,nblock)
      dimension field(nblock,ncomp,nfield)
c
c ---- locals
      integer           maxrec
      parameter        (maxrec = 200000)
      dimension         coords_field(maxrec,3), field_vusdfld(maxrec)
      real*8            coords_min, field_min, pile_radius, min_coords3, max_coords3
      real*8            dist2, x, y, z, val
      integer           nline, i, k, ios
      integer           kprocessnum
      integer           frequency_update, kinc
      character*256     outdir
      logical           ex
c
c ---- helper vars
      integer           rank, lun
      character*128     fdata, fmark
      integer           fsz
c
!#============================================================================= !# Variables to be edited !#
      pile_radius      = 1.0d0     !! only in the zone 1.4 times the pile radius the field variable is assigned
      frequency_update = 10        !! only in every 10th calculation increment the field variable is updated
      min_coords3      = 11.0d0    !! for z-coords larger than this value the field variable is not assigned, change at end of file if needed!!
      max_coords3      = 0.0d0     !! for z-coords lower than this value the field variable is not assigned, change at end of file if needed!!
c
c ---- rank / outdir
      call vgetrank( kprocessnum )
      call vgetoutdir( outdir, lenoutdir )
      kinc = jflags(i_ufld_kinc)
c
c ---- cadence:
c      every frequency_update-th increment: rank 0 cleans files and returns
      if ( mod(kinc, frequency_update) .eq. 0 ) then
         if (kprocessnum .eq. 0) then
            do rank = 0, 9999
               call make_filenames(outdir, rank, fdata, fmark)
               inquire(file=trim(fdata), exist=ex)
               if (ex) then
                  open(newunit=lun, file=trim(fdata), status='old', iostat=ios)
                  if (ios .eq. 0) close(lun, status='delete')
               endif
               inquire(file=trim(fmark), exist=ex)
               if (ex) then
                  open(newunit=lun, file=trim(fmark), status='old', iostat=ios)
                  if (ios .eq. 0) close(lun, status='delete')
               endif
            end do
         endif
         return
c      only proceed on the increment immediately after cleanup
      else if ( mod(kinc, frequency_update) .ne. 1 ) then
         return
      endif
c
c ---- collect all samples from all rank files 
      nline = 0
      do rank = 0, 9999
         call make_filenames(outdir, rank, fdata, fmark)
c        require marker
         inquire(file=trim(fmark), exist=ex)
         if (.not. ex) cycle
c        require data file with size > 0
         inquire(file=trim(fdata), exist=ex, size=fsz)
         if (.not. ex) cycle
         if (fsz .le. 0) cycle
c        read all lines
         open(newunit=lun, file=trim(fdata), form='formatted',
     &        access='sequential', status='old', action='read', iostat=ios)
         if (ios .ne. 0) cycle
   10    continue
         read(lun,*,iostat=ios) x, y, z, val
         if (ios .ne. 0) goto 11
         nline = nline + 1
         if (nline .le. maxrec) then
            coords_field(nline,1) = x
            coords_field(nline,2) = y
            coords_field(nline,3) = z
            field_vusdfld(nline)  = val
         endif
         goto 10
   11    continue
         close(lun)
         if (nline .ge. maxrec) exit
      end do
c
      if (nline .le. 0) return
c
c ---- nearest-neighbor assignment 
      do 200 k = 1, nblock
         if ( (coords(1,k)**2 + coords(2,k)**2) .lt. (pile_radius*1.2d0)**2
     &        .and. (coords(3,k) .gt. min_coords3)
     &        .and. (coords(3,k) .lt. max_coords3) ) then
c
            coords_min = 1.0d30
            field_min  = 0.0d0
c
            do 210 i = 1, nline
               dist2 = (coords(1,k)-coords_field(i,1))**2
     &               + (coords(2,k)-coords_field(i,2))**2
     &               + (coords(3,k)-coords_field(i,3))**2
               if (dist2 .lt. coords_min) then
                  coords_min = dist2
                  field_min  = field_vusdfld(i)
               endif
  210       continue
c
            field(k,1,1) = field_min
         endif
  200 continue
c
      return
      end

