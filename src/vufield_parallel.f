

      SUBROUTINE VUFIELD(FIELD, NBLOCK, NFIELD, KFIELD, NCOMP,
     1                   KSTEP, JFLAGS, JNODEUID, TIME,
     2                   COORDS, U, V, A)
C
      INCLUDE 'VABA_PARAM.INC'


C     indices for the time array TIME
      PARAMETER( i_ufld_Current   = 1,
     *           i_ufld_Increment = 2,
     *           i_ufld_Period    = 3,
     *           i_ufld_Total     = 4 )

C     indices for the coordinate array COORDS
      PARAMETER( i_ufld_CoordX = 1,
     *           i_ufld_CoordY = 2,
     *           i_ufld_CoordZ = 3 )

C     indices for the displacement array U
      PARAMETER( i_ufld_SpaDisplX = 1,
     *           i_ufld_SpaDisplY = 2,
     *           i_ufld_SpaDisplZ = 3,
     *           i_ufld_RotDisplX = 4,
     *           i_ufld_RotDisplY = 5,
     *           i_ufld_RotDisplZ = 6,
     *           i_ufld_AcoPress  = 7,
     *           i_ufld_Temp      = 8 )

C     indices for the velocity array V
      PARAMETER( i_ufld_SpaVelX   = 1,
     *           i_ufld_SpaVelY   = 2,
     *           i_ufld_SpaVelZ   = 3,
     *           i_ufld_RotVelX   = 4,
     *           i_ufld_RotVelY   = 5,
     *           i_ufld_RotVelZ   = 6,
     *           i_ufld_DAcoPress = 7,
     *           i_ufld_DTemp     = 8 )

C     indices for the acceleration array A
      PARAMETER( i_ufld_SpaAccelX  = 1,
     *           i_ufld_SpaAccelY  = 2,
     *           i_ufld_SpaAccelZ  = 3,
     *           i_ufld_RotAccelX  = 4,
     *           i_ufld_RotAccelY  = 5,
     *           i_ufld_RotAccelZ  = 6,
     *           i_ufld_DDAcoPress = 7,
     *           i_ufld_DDTemp     = 8 )

C     indices for JFLAGS
      PARAMETER( i_ufld_kInc   = 1,
     *           i_ufld_kPass  = 2 )
C
      DIMENSION FIELD(NBLOCK,NCOMP,NFIELD)
      DIMENSION JFLAGS(2), JNODEUID(NBLOCK), TIME(4),
     *          COORDS(3,NBLOCK)
      DIMENSION U(8,NBLOCK), V(8,NBLOCK), A(8,NBLOCK),coords_field(5000,3)
      real(8) coords_min, field_min, field_vusdfld(5000), pile_radius, min_coords3
      integer nline, ios, jj, ii, k
      logical file_exists
      integer frequency_update, bytes_per_line
      character*256 OUTDIR
      
      CALL VGETOUTDIR( OUTDIR, LENOUTDIR )
      
!# =============================================================================
!# Variables to be edited
!# =============================================================================
      pile_radius      = 0.5d0     !! only in the zone 1.4 times the pile radius the field variable is assigned
      frequency_update = 10        !! only in every 10th calculation increment the field variable is updated
      min_coords3      = 85.0d0    !! for z-coords lower than this value the field variable is not assigned, change at end of file if needed!!
      bytes_per_line   = 62        !! bytes for one line of the files to be read, Check if subroutine is reading correct number of lines!!
!# =============================================================================
!# Get processes
!# =============================================================================
      CALL VGETRANK( KPROCESSNUM )

      if(KPROCESSNUM == 0) then
        myunit = 105
      elseif(KPROCESSNUM == 1) then
        myunit = 106
      elseif(KPROCESSNUM == 2) then
        myunit = 107
      elseif(KPROCESSNUM == 3) then
        myunit = 108
      elseif(KPROCESSNUM == 4) then
        myunit = 109
      elseif(KPROCESSNUM == 5) then
        myunit = 110
      elseif(KPROCESSNUM == 6) then
        myunit = 111
      elseif(KPROCESSNUM == 7) then
        myunit = 112
      else
        write(*,*) 'ERROR SUBROUTINE VUFIELD: Not defined for the given number of processes'
      endif

      !! Operations are performed every 10th increment
      if (mod(JFLAGS(i_ufld_kInc),frequency_update) <= 0) then
        field(:,1,2) = 10
        if(KPROCESSNUM == 0) then
          inquire(file = 'vusdfld_out1.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out1.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out2.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out2.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out3.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out3.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out4.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out4.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out5.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out5.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out6.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out6.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out7.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out7.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out8.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out8.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          return
        else
          return
        endif
      endif

      !! Return and delete files for other increments
      if(any(field(:,1,2) <5)) then
        if(KPROCESSNUM == 0) then
          inquire(file = 'vusdfld_out1.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out1.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out2.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out2.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out3.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out3.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out4.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out4.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out5.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out5.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out6.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out6.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          inquire(file = 'vusdfld_out7.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out7.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
            inquire(file = 'vusdfld_out8.dat',
     1    exist=file_exists)
          if(file_exists) then
            open(unit = myunit, file = 'vusdfld_out8.dat',
     1      status='old', position="append", action="readwrite", iostat = ios)
            close(myunit, status='delete')
          endif
          return
        else
          return
        endif
      endif
!# =============================================================================
!# Start reading values to update field variables at nodes
!# =============================================================================
      !! the second field variable indicates that in the next increments no update is performed
      field(:,1,2) = 0

!# =============================================================================
!# For each thread an individual file is written. Currently only up to 4 threads are supported
!# First file from first process
!# =============================================================================
      file_exists = .false.
      inquire(file = 'vusdfld_out1.dat',
     1exist=file_exists,size=isize)
      if(file_exists) then
        open(unit = myunit, file = 'vusdfld_out1.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
      else
        RETURN
      endif
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      nline = 1
      if(isize>0) then
        read1: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read1
          nline = nline + 1
        end do read1
        close(myunit)
      endif

!# =============================================================================
!# Second file from second process
!# =============================================================================
      inquire(file = 'vusdfld_out2.dat',
     1exist=file_exists,size=isize,RECL=iRECL)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0 ) then
        open(unit = myunit, file = 'vusdfld_out2.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read2: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read2
          nline = nline + 1
        end do read2
        close(myunit)
      endif

!# =============================================================================
!# Third file from third process
!# =============================================================================
      inquire(file = 'vusdfld_out3.dat',
     1exist=file_exists,size=isize,RECL=iRECL)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out3.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read3: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios <0) exit read3
          nline = nline + 1
        end do read3
        close(myunit)
      endif

!# =============================================================================
!# Fourth file from fourth process
!# =============================================================================
      inquire(file = 'vusdfld_out4.dat',
     1exist=file_exists,size=isize)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out4.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read4: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read4
          nline = nline + 1
        end do read4
        close(myunit)
      endif
!# =============================================================================
!# Fifth file from fifth process
!# =============================================================================
      inquire(file = 'vusdfld_out5.dat',
     1exist=file_exists,size=isize)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out5.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read5: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read5
          nline = nline + 1
        end do read5
        close(myunit)
      endif
!# =============================================================================
!# Sixth file from sixth process
!# =============================================================================
      inquire(file = 'vusdfld_out6.dat',
     1exist=file_exists,size=isize)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out6.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read6: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read6
          nline = nline + 1
        end do read6
        close(myunit)
      endif
!# =============================================================================
!# Seventh file from seventh process
!# =============================================================================
      inquire(file = 'vusdfld_out7.dat',
     1exist=file_exists,size=isize)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out7.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read7: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read7
          nline = nline + 1
        end do read7
        close(myunit)
      endif
      
!# =============================================================================
!# Eighth file from eighth process
!# =============================================================================
      inquire(file = 'vusdfld_out8.dat',
     1exist=file_exists,size=isize)
      nlines = int(isize/bytes_per_line) !! This solution was required because it appears that ios checking if end of file reached does not work with abaqus
      if(file_exists .and. isize > 0) then
        open(unit = myunit, file = 'vusdfld_out8.dat',
     1  status='unknown', action="read", iostat = ios, form="FORMATTED")
        read8: do ii = 1, nlines
          read(myunit,*,iostat=ios) coords_field(nline,1:3), field_vusdfld(nline)
          if(ios < 0) exit read8
          nline = nline + 1
        end do read8
        close(myunit)
      endif
      
!# =============================================================================
!# Assign field value to nodes using the shortest distance
!# =============================================================================

      do k = 1, nblock

        if(sqrt(COORDS(1,k)**2+ COORDS(2,k)**2) < pile_radius*1.2 .and.
     1         COORDS(3,k)>min_coords3) then

          coords_min = 1d6

          do jj = 1, nline

            if(norm2(COORDS(:,k)-coords_field(jj,:))< coords_min) then
              coords_min = norm2(COORDS(:,k) - coords_field(jj,:))
              field_min = field_vusdfld(jj)
            endif
          enddo
          field(k,1,1)  = field_min
        endif

       enddo

c

      RETURN
      END

