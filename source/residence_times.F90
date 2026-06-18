!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the computation of residence times for reactive
! species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module residence_times
 
  Use atomic_model,     Only: model_type

  Use fileset,          Only: file_type, &
                              FILE_RES_TIMES, &
                              FILE_SET, & 
                              refresh_out

  Use input_types,      Only: in_param,   &
                              in_string
                              
  Use numprec,          Only: wi,& 
                              wp
 
  Use process_data,     Only: set_read_status, &
                              capital_to_lower_case, &
                              check_for_rubbish, &
                              get_word_length

  Use trajectory,       Only: traj_type, &
                              check_time_directive
                              
  Use unit_output,      Only: info, &
                              error_stop                               

  Implicit None
  Private  
  !Type for residence times
  Type, Public :: restimes_type
    Type(in_string)  :: invoke
    Type(in_param)   :: rattling_wait
  End Type

  Public :: read_residence_times, check_residence_times
  Public :: residence_times_reactive_sites
  
Contains

  Subroutine read_residence_times(iunit, restimes_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &residence_times block
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(restimes_type), Intent(InOut) :: restimes_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &residence_times block (within the &reactive_analysis block, SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_residence_times" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_residence_times') Exit
      Call check_for_rubbish(iunit, '&residence_times')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='rattling_wait') Then
         Read (iunit, Fmt=*, iostat=io) restimes_data%rattling_wait%tag, &
                                      & restimes_data%rattling_wait%value,& 
                                      & restimes_data%rattling_wait%units
         Call set_read_status(word, io, restimes_data%rattling_wait%fread,&
                                      & restimes_data%rattling_wait%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_residence_times"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_residence_times
  
  Subroutine check_residence_times(files, restimes_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &residence_times block
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(In   ) :: files(:)
    Type(restimes_type), Intent(InOut) :: restimes_data

    Character(Len=256)  :: error_set

    error_set = '***ERROR in the &residence_times block of file '//Trim(files(FILE_SET)%filename)//' -'

    Call check_time_directive(restimes_data%rattling_wait, 'rattling_wait', error_set, .False.)

    If (.Not. restimes_data%rattling_wait%fread) Then
      restimes_data%rattling_wait%value= 0.0_wp
      restimes_data%rattling_wait%units= 'fs'
    End If
    
  End Subroutine check_residence_times  
 

  Subroutine residence_times_reactive_sites(files, model_data, traj_data, restimes_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence times for the changing species
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(In   ) :: model_data
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(restimes_type), Intent(InOut) :: restimes_data

    Integer(Kind=wi)   :: i, m
    Integer(Kind=wi)   :: iunit
    Real(Kind=wp)      :: time
    Logical            :: set_u0
    Character(Len=256) :: message
    
    Integer(Kind=wi)   :: ref_indx(model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: icount(model_data%reactive_species%N0%value)    
    Real(Kind=wp)      :: values(traj_data%frames, model_data%reactive_species%N0%value)
    Character(Len=8)   :: tag(traj_data%frames, 2, model_data%reactive_species%N0%value)
    Character(Len=8)   :: ref_tag(model_data%reactive_species%N0%value)
    Logical            :: hold(model_data%reactive_species%N0%value)
    
    Real(Kind=wp)      :: rattling
    Real(Kind=wp)      :: tchange, base_time
    
    rattling=restimes_data%rattling_wait%value

    Do m = 1, model_data%reactive_species%N0%value
      base_time=(traj_data%seg_analysis%frame_ini-1)*traj_data%timestep%value
      set_u0=.True.
      icount(m)=0
      i = traj_data%seg_analysis%frame_ini
      hold(m)=.False.
      Do While (i <= traj_data%seg_analysis%frame_last)
        time=(i-1)*traj_data%timestep%value
        If (set_u0) Then
          ref_indx(m)=traj_data%track_chem%config(i,m)%indx
          ref_tag(m)=traj_data%track_chem%config(i,m)%tag
          set_u0=.False.
        Else
          If (traj_data%track_chem%config(i,m)%indx /= ref_indx(m)) Then
            If (.Not. hold(m)) Then
              tchange=time
              hold(m)=.True.
            End If
          Else 
            hold(m)=.False.
          End If  
          If (hold(m)) Then
            If (((time-tchange) > rattling .Or. i==traj_data%seg_analysis%frame_last)) Then
              hold(m)=.False.
              icount(m)=icount(m)+1
              values(icount(m),m)=(tchange-base_time)/1000.0_wp
              tag(icount(m),1,m)=ref_tag(m)
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
              ref_tag(m)=traj_data%track_chem%config(i,m)%tag
              base_time=tchange
            End If
          End If
        End If
        i=i+1
      End Do
    End Do
    
    Open(Newunit=files(FILE_RES_TIMES)%unit_no, File=files(FILE_RES_TIMES)%filename, Status='Replace')
    iunit=files(FILE_RES_TIMES)%unit_no
    
    Do m = 1, model_data%reactive_species%N0%value
      Write (iunit,'(a,i3)') '#  Species', m 
      Write (iunit,'(a)') '#  Residence Time (ps)    Tag for site' 
        Do i =1, icount(m)
          Write(iunit,'(f11.3,15x,a)') values(i,m), Trim(tag(i,1,m))
        End Do
      Write (iunit,'(a)') ' ' 
    End Do 
    
    Write (message,'(1x,a)') 'The Residence Times for each reactive species were&
                             & printed to the "'//Trim(files(FILE_RES_TIMES)%filename)//'" file.'
    Call info(message, 1)
    Close(iunit)
    Call info(' ', 1)
    Call refresh_out(files)

    Call refresh_out(files)
    
  End Subroutine residence_times_reactive_sites

End Module residence_times
