      subroutine vusdfld(
c Read only variables -
     1   nblock, nstatev, nfieldv, nprops, ndir, nshr,
     2   jElem, kIntPt, kLayer, kSecPt,
     3   stepTime, totalTime, dt, cmname,
     4   coordMp, direct, T, charLength, props,
     5   stateOld,
c Write only variables -
     6   stateNew, field )
c
      include 'vaba_param.inc'
c
      dimension jElem(nblock), coordMp(nblock,*),
     1 direct(nblock,3,3), T(nblock,3,3),
     2          charLength(nblock), props(nprops),
     3          stateOld(nblock,nstatev),
     4          stateNew(nblock,nstatev),
     5          field(nblock,nfieldv)
      character*80 cmname
c
c     Local arrays from vgetvrm are dimensioned to
c     maximum block size (maxblk)
c
      parameter( nrData=6 )
      character*3 cData(maxblk*nrData)
      dimension rData(maxblk*nrData), jData(maxblk*nrData)
      real(8) normal_vector(2),normal_stress(2), pile_center(2), pile_radius, min_coords3
      integer KPROCESSNUM, myunit

      logical file_exists
      character*256 OUTDIR
      
      CALL VGETOUTDIR( OUTDIR, LENOUTDIR )
      
!# =============================================================================
!# Variables to be edited
!# =============================================================================
      pile_radius      = 0.5d0      !! only in the zone 1.4 times the pile radius the field variable is assigned
      min_coords3      = 85.0d0     !! for z-coords lower than this value the field variable is not assigned, change at end of file if needed!!
      pile_center(1:2) = [0.0d0,0.0d0] !! define pile center
!# =============================================================================
!# Get processes
!# =============================================================================
      !! field(:,2) is filled by vufield
      if(any(field(:,2)<5)) return

      file_exists = .false.

      CALL VGETRANK(KPROCESSNUM) !! This gives the process number 

      if(KPROCESSNUM == 0) then
        myunit =105
        inquire(file = 'vusdfld_out1.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out1.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out1.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 1) then
        myunit =106
        inquire(file = 'vusdfld_out2.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out2.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out2.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 2) then
        myunit =107
        inquire(file = 'vusdfld_out3.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out3.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out3.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 3) then
        myunit =108
        inquire(file = 'vusdfld_out4.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out4.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out4.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 4) then
        myunit =109
        inquire(file = 'vusdfld_out5.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out5.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out5.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 6) then
        myunit =110
        inquire(file = 'vusdfld_out6.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out6.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out6.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 7) then
        myunit =111
        inquire(file = 'vusdfld_out7.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out7.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out7.dat',
     1    status='new', action="write", iostat = ios)
        endif
      elseif(KPROCESSNUM == 8) then
        myunit =112
        inquire(file = 'vusdfld_out8.dat',
     1   exist=file_exists)

        if(file_exists) then
          open(unit = myunit, file = 'vusdfld_out8.dat',
     1    status="old", position="append", action="readwrite", iostat = ios)
        else
          open(unit = myunit, file = 'vusdfld_out8.dat',
     1    status='new', action="write", iostat = ios)
        endif
      else
        write(*,*) 'ERROR SUBROUTINE VUFIELD: Not defined for the given number of processes'
      endif
	  
!# =============================================================================
!# Calculate field variable
!# =============================================================================
      do 100 k = 1, nblock

        if(sqrt(coordMp(k,1)**2+ coordMp(k,2)**2) < pile_radius*1.2 .and.
     1          coordMp(k,3)>min_coords3) then

         !! The normal vector has to point inward (multiply with -1)
         normal_vector(1:2) = -(coordMp(k,1:2)-pile_center(1:2))/norm2(coordMp(k,1:2)-pile_center(1:2))

         !! The pile has a radius of pile_radius, inside of the pile the normal vector has to point outward
         if(sqrt(coordMp(k,1)**2+coordMp(k,2)**2) <= pile_radius)
     1    normal_vector(1:2) = - normal_vector(1:2)

         normal_stress(1) = stateOld(k,27)*normal_vector(1) + stateOld(k,30)*normal_vector(2)
         normal_stress(2) = stateOld(k,28)*normal_vector(2) + stateOld(k,30)*normal_vector(1)

         field(k,1) = abs(dot_product(normal_stress,normal_vector))
     1   /(stateOld(k,24) + abs(dot_product(normal_stress,normal_vector))
     2   +  stateOld(k,25))

         if(field(k,1) > 1.0d0) field(k,1) = 1.0d0
         if(field(k,1) < 0.0d0) field(k,1) = 0.0d0
         if(isnan(field(k,1))) field(k,1) = 0.0d0

          write(myunit,*) coordMp(k,1:3),field(k,1)

        endif

  100 continue

      close(myunit)


c
      return
      end