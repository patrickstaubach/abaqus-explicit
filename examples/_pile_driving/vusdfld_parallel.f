!=======================================================================================================
! this file is part of vumat_hmc_staubach.
!
! vumat_hmc_staubach is free software: you can redistribute it and/or modify it under the terms of the
! gnu general public license as published by the free software foundation, either version 3 of the license,
! or (at your option) any later version.
!
! vumat_hmc_staubach is distributed in the hope that it will be useful, but without any warranty; without
! even the implied warranty of merchantability or fitness for a particular purpose. see the gnu general public
! license for more details.
!
! you should have received a copy of the gnu general public license along with vumat_hmc_staubach.
! if not, see <https://www.gnu.org/licenses/>.
!=======================================================================================================
!
! subroutine: vusdfld
!
!> @author patrick staubach, patrick.staubach@yahoo.de
!          bauhaus university weimar, ruhr-university bochum
!
! description:
!> contains the vusdfld routine used to calculate effective interface friction according to the paper
!> "hydro-mechanically coupled cel analyses with effective contact stresses"
!> international journal for numerical and analytical methods in geomechanics
!> https://onlinelibrary.wiley.com/doi/10.1002/nag.3725
!
! revision history
!> 03.03.2024 - initial version
!> 22.03.2024 - support of up to 8 cpus
!> 23.10.2025 - No CPU limitation, stylistic cleanup
!=======================================================================================================
      subroutine vusdfld(
c read only variables -
     1   nblock, nstatev, nfieldv, nprops, ndir, nshr,
     2   jelem, kintpt, klayer, ksecpt,
     3   steptime, totaltime, dt, cmname,
     4   coordmp, direct, t, charlength, props,
     5   stateold,
c write only variables -
     6   statenew, field )
c
      include 'vaba_param.inc'
c
c ---- arguments (arrays declared with dimension as requested)
      integer nblock, nstatev, nfieldv, nprops, ndir, nshr
      integer kintpt(nblock), klayer(nblock), ksecpt(nblock)
      dimension jelem(nblock), coordmp(nblock,*), direct(nblock,3,3), t(nblock,3,3)
      character*80 cmname
      real*8 steptime, totaltime, dt
      real*8 charlength(nblock), props(nprops)
      dimension stateold(nblock,nstatev), statenew(nblock,nstatev)
      dimension field(nblock,nfieldv)
c
c ---- locals
      integer          kprocessnum, lun, ios, k
      character*256    outdir
      character*128    fdata, fmark
      logical          do_write, ex
      integer          nw
      real*8           normal_vector(2), normal_stress(2)
      real*8           pile_center(2), pile_radius, min_coords3, max_coords3
      real*8           dx, dy, nrm
      
!#============================================================================= !# Variables to be edited !#

      pile_radius   = 1.0d0   !! only in the zone 1.2 times the pile radius the field variable is assigned
      min_coords3   = 0.0d0   !! for z-coords lower than this value the field variable is not assigned, change at end of file if needed!!
      max_coords3   = 11.0d0  !! for z-coords larger than this value the field variable is not assigned, change at end of file if needed!!
      pile_center(1)= 0.0d0 !! define pile center in x
      pile_center(2)= 0.0d0 !! define pile center in y

c
c ---- identify rank / filenames
      call vgetrank(kprocessnum)
      call vgetoutdir(outdir, lenoutdir)
      call make_filenames(outdir, kprocessnum, fdata, fmark)
c
c ---- open per-rank file and truncate (no accumulation)
      open(newunit=lun, file=trim(fdata), form='formatted',
     &     access='sequential', status='replace',
     &     action='write', iostat=ios)
      if (ios .ne. 0) return
c
      nw = 0
c
c ---- compute and write samples (every increment; reader controls cadence)
      do 500 k = 1, nblock
c       initialize field for this point
        field(k,:) = 0.0d0
        do_write = .false.
c       radial filter
        do_write = ( sqrt(coordmp(k,1)**2 + coordmp(k,2)**2)
     &               .lt. (pile_radius*1.2d0) )
c       z window [min_coords3, max_coords3]
        if (do_write) do_write = (coordmp(k,3) .lt. max_coords3)
        if (do_write) do_write = (coordmp(k,3) .gt. min_coords3)
c
        if (do_write) then
c         inward normal (multiply by -1), projected in xy plane
          dx  = coordmp(k,1) - pile_center(1)
          dy  = coordmp(k,2) - pile_center(2)
          nrm = sqrt(dx*dx + dy*dy)
          if (nrm .eq. 0.0d0) goto 500
          normal_vector(1) = -dx / nrm
          normal_vector(2) = -dy / nrm
c
c         inside the pile, flip to outward
          if ( (coordmp(k,1)**2 + coordmp(k,2)**2) .le. pile_radius**2 ) then
            normal_vector(1) = -normal_vector(1)
            normal_vector(2) = -normal_vector(2)
          endif
c
c         stress projection 
          normal_stress(1) = stateold(k,27)*normal_vector(1)
     &                     + stateold(k,30)*normal_vector(2)
          normal_stress(2) = stateold(k,28)*normal_vector(2)
     &                     + stateold(k,30)*normal_vector(1)
c
          field(k,1) = abs( normal_stress(1)*normal_vector(1)
     &                    + normal_stress(2)*normal_vector(2) )
     &               / ( stateold(k,24)
     &                 + abs( normal_stress(1)*normal_vector(1)
     &                      + normal_stress(2)*normal_vector(2) )
     &                 + stateold(k,25) )
c
c         clamp and sanitize
          if (isnan(field(k,1))) field(k,1) = 0.0d0
          if (field(k,1) .gt. 1.0d0) field(k,1) = 1.0d0
          if (field(k,1) .lt. 0.0d0) field(k,1) = 0.0d0
c
c         one line: x y z value
          write(lun,'(4(1x,es24.16))') coordmp(k,1), coordmp(k,2),
     &                                   coordmp(k,3), field(k,1)
          nw = nw + 1
        endif
c
500   continue
c
      close(lun)
c
c ---- create marker only if we actually wrote something (nw > 0)
      if (nw .gt. 0) then
         open(newunit=lun, file=trim(fmark), status='replace',
     &        action='write', iostat=ios)
         if (ios .eq. 0) then
            write(lun,'(a)') 'ok'
            close(lun)
         endif
      else
c        ensure no stale marker exists for empty data
         inquire(file=trim(fmark), exist=ex)
         if (ex) then
            open(newunit=lun, file=trim(fmark), status='old', iostat=ios)
            if (ios .eq. 0) close(lun, status='delete')
         endif
      endif
c
      return
      end
c
c=======================================================================================================
c helper (same signature as in vufield): data and marker names (rbuf = 6 chars)
c=======================================================================================================
      subroutine make_filenames(outdir, rank, fdata, fmark)
      character*(*) outdir
      integer       rank
      character*128 fdata, fmark
      character*6   rbuf
      integer       lod
c
      write(rbuf,'(i6.6)') rank
      fdata = 'vusdfld_out_rank'//rbuf//'.txt'
      fmark = 'vusdfld_out_rank'//rbuf//'.ok'
c
      lod = len_trim(outdir)
      if (lod .gt. 0) then
        if (outdir(lod:lod) .ne. '/' .and. outdir(lod:lod) .ne. '\' ) then
          fdata = outdir(1:lod)//'/'//fdata
          fmark = outdir(1:lod)//'/'//fmark
        else
          fdata = outdir(1:lod)//fdata
          fmark = outdir(1:lod)//fmark
        endif
      endif
c
      return
      end
