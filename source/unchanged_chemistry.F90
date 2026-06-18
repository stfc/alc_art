!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to track unchanged species along the trajectory
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module unchanged_chemistry

  Use atomic_model,   Only:  model_type, &
                             read_model


  Use constants,      Only : max_components, &
                             max_unchanged_atoms
                      
  Use fileset,        Only : file_type, &
                             FILE_SET, & 
                             FILE_TRAJECTORY, &
                             FILE_UNCHANGED_CHEM, & 
                             refresh_out  

  Use input_types,    Only: in_string
                             
                             
  Use numprec,        Only: wi,& 
                            wp

  Use process_data,   Only: capital_to_lower_case, &
                            check_for_rubbish, &
                            get_word_length, &
                            set_read_status, &
                            check_end
                         
  Use trajectory,     Only: traj_type
                            
  Use unit_output,    Only: info, &
                            error_stop 
                           
  Implicit None
  Private

  !Type to print the position of selected atoms, whose content remain unchanged
  Type, Public :: unchanged_type
    Private
    Type(in_string), Public  :: invoke
    Type(in_string)  :: tag
    Integer(Kind=wi) :: N0
    Integer(Kind=wi) :: indexes(max_unchanged_atoms)
    Type(in_string)  :: list_indexes
  End Type
  
  Public :: read_track_unchanged_chemistry
  Public :: check_initial_unchanged_labels, check_unchanged_chemistry
  Public :: print_unchanged_chemistry
  
Contains

  Subroutine read_track_unchanged_chemistry(iunit, unchanged_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settings to track chemically unchanged species 
    ! along the trajectory
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),     Intent(In   ) :: iunit
    Type(unchanged_type), Intent(InOut) :: unchanged_data 

    Integer(Kind=wi)   :: io, length, j
    Character(Len=256) :: message, messages(2)
    Character(Len=256) :: word
    Character(Len=256) :: set_error
    Logical :: error, fread
    
    set_error = '***ERROR in the &track_unchanged_chemistry block (SETTINGS file).'
    error=.False.
    fread= .True.

    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&track_unchanged_chemistry')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&track_unchanged_chemistry')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, unchanged_data%N0
    If (Trim(word) /= 'number') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "number" is expected to be defined first'
      error=.True.
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number"'
      error=.True.
    Else
      If (unchanged_data%N0<1) Then
        Write (messages(2),'(a)') 'The "number" directive MUST BE >= 1'
        error=.True.
      End If  
      If (unchanged_data%N0>max_components) Then
        Write (messages(2),'(a,i3,a)') 'Directive number: are you sure you want to consider more than ', max_components,&
                                    & '? Please check'
        error=.True.
      End If
      If (unchanged_data%N0>max_unchanged_atoms) Then
        Write (messages(2),'(a,i3,a)') 'Directive "number": the user cannot track more than ', max_unchanged_atoms,&
                                       &' per simulation. In case a larger number is needed, run the code several times'
        error=.True.
      End If
    End If
    ! print erro if any
    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_track_unchanged_chemistry" to close the block.&
                                  & Check if directives "tag" and "list_indexes" are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_track_unchanged_chemistry') Exit
      Call check_for_rubbish(iunit, '&track_unchanged_chemistry')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='list_indexes') Then
        unchanged_data%indexes=-1
        Read (iunit, Fmt=*, iostat=io) unchanged_data%list_indexes%type,&
                                       (unchanged_data%indexes(j), j= 1, unchanged_data%N0)
        Call set_read_status(word, io, unchanged_data%list_indexes%fread, unchanged_data%list_indexes%fail,&
                                     & unchanged_data%list_indexes%type)

      Else If (Trim(word)=='tag') Then
         Read (iunit, Fmt=*, iostat=io) word, unchanged_data%tag%type 
         Call set_read_status(word, io, unchanged_data%tag%fread, unchanged_data%tag%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with&
                                & "&end_track_unchanged_chemistry"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_track_unchanged_chemistry

  Subroutine check_initial_unchanged_labels(files, model_data, traj_data, unchanged_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the atomic tags for each component of the 
    ! list_indexes (&track_unchanged_chemistry block) is the same as the "tag"
    ! directive defined
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),       Intent(InOut) :: files(:)
    Type(model_type),      Intent(InOut) :: model_data
    Type(traj_type),       Intent(InOut) :: traj_data
    Type(unchanged_type),  Intent(InOut) :: unchanged_data

    Character(Len=256)  :: messages(2), word
    Character(Len=256)  :: error_set
    Integer(Kind=wi)    :: j, k
    
    error_set = '***ERROR in the &track_unchanged_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=files(FILE_TRAJECTORY)%filename, Status='old')
    Call read_model(files, model_data, 1, traj_data%ensemble%type)
    Close(files(FILE_TRAJECTORY)%unit_no) 
    
    Do j=1, unchanged_data%N0
      k=unchanged_data%indexes(j) 
      If ((model_data%config%atom(k)%tag)/=(unchanged_data%tag%type)) Then
        Call info(' ', 1)
        Write(word,*) k
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Index "'//Trim(Adjustl(word))//'" (defined in list_indexes)&
                                       & does not correspond to the atomic tag "'//Trim(unchanged_data%tag%type)//'".'  
        Write (messages(2),'((1x,a))') 'According to the &reference_composition block, this index corresponds to atomic&
                                       & tag "'//Trim(model_data%config%atom(k)%tag)//'".&
                                       & Please review the labels and indexes of the atomic model'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 
      
  End Subroutine check_initial_unchanged_labels
  
  Subroutine check_unchanged_chemistry(files, model_data, unchanged_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &track_unchanged_chemistry block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(In   ) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(unchanged_type), Intent(InOut) :: unchanged_data

    Character(Len=256)  :: messages(2), word
    Character(Len=256)  :: error_set
    Integer(Kind=wi)    :: j, k
    Logical             :: flag
    
    error_set = '***ERROR in the &track_unchanged_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (unchanged_data%tag%fread) Then
      If (unchanged_data%tag%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tag" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      ! Check if all tags correspond to the same element (type a)
      j=1
      flag=.True.
      Do While (j <= model_data%reference_composition%atomic_species .And. flag)
        If (model_data%reference_composition%tag(j)==unchanged_data%tag%type) Then
          flag=.False.
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The atomic tag "'//Trim(unchanged_data%tag%type)//&
                                       &'" (defined for the "tag" directive) has not been defined&
                                       & in the &reference_composition block! Please review the settings' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must the "tag" (atomic tag)&
                                    & to track along the trajectory'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    If (unchanged_data%list_indexes%fread) Then
      If (unchanged_data%list_indexes%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "list_indexes" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If

      Do j=1, unchanged_data%N0-1
        Do k=j+1, unchanged_data%N0
          If (unchanged_data%indexes(j)==unchanged_data%indexes(k)) Then
            Write(word,*) unchanged_data%indexes(j)
            Write (messages(1),'(2(1x,a))') Trim(error_set), 'Index "'//Trim(Adjustl(word))//' is repeated in the list!'
            Write (messages(2),'((1x,a))') 'Values in the "list_indexes" must be  different'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End Do
      End Do 
      
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define "list_indexes" for&
                                    & all those atoms that the user wants to print'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine check_unchanged_chemistry

  Subroutine print_unchanged_chemistry(files, traj_data, unchanged_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the positions of those atomic indexes defined in the
    ! &track_unchanged_chemistry block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(traj_type),      Intent(InOut) :: traj_data
    Type(unchanged_type), Intent(InOut) :: unchanged_data
  
    Integer(Kind=wi)   :: iunit, i, k, l 
    Character(Len=256) :: frame
    Character(Len=256) :: num_species, message, messages(3), spec, num1, num2
    Logical            :: flag
  
    flag=.True. 
  
    Open(Newunit=files(FILE_UNCHANGED_CHEM)%unit_no, File=files(FILE_UNCHANGED_CHEM)%filename, Status='Replace')
    iunit=files(FILE_UNCHANGED_CHEM)%unit_no
    
    If (traj_data%seg_analysis%frame_ini==1) Then
      Write(iunit,'(a)') '# Tracking unchanged chemical species over the whole trajectory'    
    Else  
      Write(iunit,'(a,1x,f10.4,1x,a)') '# Tracking the unchanged chemical species ignoring the first',& 
                                   &  traj_data%seg_analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.' 
    End If
    
    If(unchanged_data%N0==1) Then
      Write(iunit,'(a)') '# The label and number for the species is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'
    Else                               
      Write(iunit,'(a)') '# The species labelling, ordering and numbering is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'    
    End If
    
    spec=unchanged_data%tag%type
    Write (num1,*) unchanged_data%indexes(1)
    If(unchanged_data%N0==1) Then
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num1))
    Else
      Write (num2,*) unchanged_data%indexes(unchanged_data%N0)
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num1))//&
                              &'.... XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num2))
    End If 
    
    i=traj_data%seg_analysis%frame_ini
    Do While (i <= traj_data%seg_analysis%frame_last .And. flag)
      l =1
      Do While (l<= unchanged_data%N0 .And. flag)
        k=unchanged_data%indexes(l)
        If (Trim(traj_data%config(i,k)%tag) /= Trim(unchanged_data%tag%type)) Then
          flag=.False.
        End If
        l=l+1
      End Do
      If (flag) Then
        Write(iunit,'(f10.4, 1x, *(f11.3))') (i-traj_data%seg_analysis%frame_ini)*traj_data%timestep%value/1000.0_wp,&
                & (traj_data%config(i,unchanged_data%indexes(l))%r(:), l=1, unchanged_data%N0)
      Else
        Write (messages(1),'(1x,a)') '**********************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: Problems with tracking species defined in the&
                                        & &track_unchanged_chemistry block'
        Write(num_species,*) k
        Write(frame,*)       i 
        Write (messages(2),'(1x,a)') '   Requested index "'//Trim(Adjustl(num_species))//'" has changed&
                                      & its chemistry at frame: '//Trim(Adjustl(frame))
        Write (messages(3),'(1x,a)') '   Please review the settings. The tracking was printed to the "'&
                                    &//Trim(files(FILE_UNCHANGED_CHEM)%filename)//'" file just up to this frame'
        Call info(messages, 3)
        Write (messages(1),'(1x,a)') '**********************************************'
        Call info(messages, 1)
      End If
      i=i+1
    End Do
    
    If (flag) Then
      Write (message,'(1x,a)') 'The tracking of the selected, unchanged chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_UNCHANGED_CHEM)%filename)//'" file'
      Call info(message, 1)
    End If
    
    Close(iunit)
    Call refresh_out(files)

  End Subroutine print_unchanged_chemistry

End Module unchanged_chemistry
