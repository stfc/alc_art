!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to obtain statistics for:
!  * intermolecular angles and distances
!  * intramolecular angles and distances 
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:   -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module nonreact_stat

  Use atomic_model,     Only: nonreactive_type, &
                              check_length_directive, &
                              check_PBC
                               
  Use constants,        Only: Rads_to_degrees, &
                              max_at_species 

  Use fileset,          Only: file_type, &
                              FILE_INTERMOL_DISTANCES_NN1, &
                              FILE_INTERMOL_DISTANCES_NN2, &
                              FILE_INTERMOL_ANGLES, &
                              FILE_INTRAMOL_DISTANCES, &
                              FILE_INTRAMOL_ANGLES, &
                              FILE_SET, &
                              refresh_out

  Use input_types,      Only: in_param, &
                              in_logic, &
                              in_string

  Use numprec,          Only: wi,& 
                              wp

  Use process_data,     Only: set_read_status, &
                              capital_to_lower_case, &
                              check_for_rubbish, &
                              get_word_length                              
                              
  Use trajectory,       Only: traj_type, &
                              within_region  
 

  Use unit_output,      Only: info, &
                              error_stop 
 
  Implicit None
  Private 
  
  ! Type for geometrical paremeter
  Type :: geo_param_type
    Type(in_string)      :: invoke
    Character(Len=256)   :: name
    Type(in_string)      :: tag_species
    Integer(Kind=wi)     :: nspecies
    Character(Len=8)     :: species(max_at_species)
    Integer(Kind=wi)     :: num_spec(max_at_species)
    Type(in_param)       :: lower_bound
    Type(in_param)       :: upper_bound
    Type(in_param)       :: delta 
  End Type
 
  ! Type for computation of the nonreactive statistics
  Type :: geo_spec_type
    Type(in_string)       :: invoke
    Character(Len=256)    :: tag
    Type(geo_param_type)  :: dist
    Type(geo_param_type)  :: angle
    Type(in_logic)        :: only_ref_tags_as_nn
  End Type 
  
  Type, Public :: nonreact_stat_type
    Type(geo_spec_type), Public :: intra_geom
    Type(geo_spec_type), Public :: inter_geom
  End Type nonreact_stat_type
  
  Public :: read_geom_param_nonreactive_species, check_nonreact_stat_settings
  Public :: geometry_statistics_nonreactive_species
  
Contains

  Subroutine read_geom_param(iunit, inblock, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read distance and angle settings
    ! for statistics of the nonreactive species
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),     Intent(In   ) :: iunit
    Character(*),         Intent(In   ) :: inblock
    Type(geo_param_type), Intent(InOut) :: M 
    
    Integer(Kind=wi)   :: io, length, i
    Character(Len=256) :: message, word
    Character(Len=256) :: messages(2)
    Character(Len=256) :: set_error
    
    M%delta%tag='delta'
    
    set_error = '***ERROR in "&'//Trim(M%name)//'" within the "&'//Trim(inblock)//'" block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_'//Trim(M%name)//'" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_'//Trim(M%name)) Exit
      Call check_for_rubbish(iunit, '&'//Trim(M%name))

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word
 
      Else If (Trim(word)=='species') Then
        If (Trim(M%name)=='distance_parameters') Then
          M%nspecies=2
        Else If (Trim(M%name)=='angle_parameters') Then
          M%nspecies=3
        End If
        Read (iunit, Fmt=*, iostat=io) M%tag_species%type, (M%species(i), i=1, M%nspecies) 
        Call set_read_status(word, io, M%tag_species%fread, M%tag_species%fail, M%tag_species%type)

      Else If (Trim(word)=='lower_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%lower_bound%tag, M%lower_bound%value, M%lower_bound%units 
         Call set_read_status(word, io, M%lower_bound%fread, M%lower_bound%fail)

      Else If (Trim(word)=='upper_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%upper_bound%tag, M%upper_bound%value, M%upper_bound%units 
         Call set_read_status(word, io, M%upper_bound%fread, M%upper_bound%fail)

      Else If (Trim(word)=='delta') Then
         Read (iunit, Fmt=*, iostat=io) M%delta%tag, M%delta%value, M%delta%units 
         Call set_read_status(word, io, M%delta%fread, M%delta%fail)

      Else
        Write (messages(1),'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_'//Trim(M%name)//'"? &
                                & Have you defined the directives correctly? See the "use_code.md" file'
        Call info (messages, 2)
        Call error_stop(' ')
      End If
    End Do
  
  End Subroutine read_geom_param

  Subroutine read_geom_param_nonreactive_species(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read blocks with parameters
    ! for the statistical analysis of geometry quantitites 
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),      Intent(In   ) :: iunit
    Type(geo_spec_type),   Intent(InOut) :: T 
    
    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the "&'//Trim(T%tag)//'" block (SETTINGS file).'
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_'//Trim(T%tag)//'" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_'//Trim(T%tag)) Exit
      Call check_for_rubbish(iunit, '&'//Trim(T%tag))

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If ((Trim(word)=='only_ref_tags_as_nn') .And. (Trim(T%tag)=='intermol_statistics')) Then
        Read (iunit, Fmt=*, iostat=io) word, T%only_ref_tags_as_nn%stat
        Call set_read_status(word, io, T%only_ref_tags_as_nn%fread, T%only_ref_tags_as_nn%fail)

      Else If (Trim(word)=='&distance_parameters') Then
        Read (iunit, Fmt=*, iostat=io) T%dist%invoke%type
        Call set_read_status(word, io, T%dist%invoke%fread, T%dist%invoke%fail, T%dist%invoke%type)
        T%dist%name='distance_parameters'
        Call read_geom_param(iunit, T%tag, T%dist)

      Else If (Trim(word)=='&angle_parameters') Then
        Read (iunit, Fmt=*, iostat=io) T%angle%invoke%type
        Call set_read_status(word, io, T%angle%invoke%fread, T%angle%invoke%fail, T%angle%invoke%type)
        T%angle%name='angle_parameters'
        Call read_geom_param(iunit, T%tag, T%angle)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_'//Trim(T%tag)//'"?'
        Call error_stop(message)
      End If

    End Do
  
  End Subroutine read_geom_param_nonreactive_species

  Subroutine check_intramol_stat_settings(error_set, T, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check species defined to 
    ! compute the statistices of intramolecular
    ! geometry
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256),      Intent(In   ) :: error_set
    Type(nonreactive_type),  Intent(In   ) :: T
    Type(geo_param_type),    Intent(InOut) :: M   
  
    Integer(Kind=wi)   ::  k, j, n
    Character(Len=1)   :: num
    Character(Len=256) :: messages(3)
    

    messages(1)=error_set
    Write (messages(2),'(1x,a)') 'Problems in "'//Trim(M%invoke%type)//'" of "&intramol_stat_settings"'

    Do j = 1, M%nspecies
      M%num_spec(j)=0
      Do k = 1, T%num_components
        If (Trim(M%species(j))==Trim(T%element(k)))Then
          M%num_spec(j)=M%num_spec(j)+1
        End If
      End Do      
      If (M%num_spec(j)==0) Then
        Write(num,'(i1)') j 
        Write (messages(3),'(1x,a)') 'Argument '//Trim(num)//' of the "species" directive does not&
                                     & correspond to the elements defined in "&atomic_components"' 
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
    
    
    Do k = 1, T%num_components
      n=0
      Do j = 1, M%nspecies
        If (Trim(M%species(j))==Trim(T%element(k)))Then
          n=n+1
        End If
      End Do
      If (n>T%N0_element(k))Then
        Write(num,'(i1)') n 
        Write (messages(3),'(1x,a)') 'The number of times the element "'//Trim(T%element(k))//'" is listed in the&
                                   & "species" directive ('//Trim(num)//' times) exceeds the value set in "&atomic_components"' 
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
    
  End Subroutine check_intramol_stat_settings
  
  Subroutine check_settings_geom_param(error_set, inblock, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check distance and angle settings
    ! for statistics of the nonreactive species
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256),   Intent(In   ) :: error_set
    Character(Len=256),   Intent(In   ) :: inblock 
    Type(geo_param_type), Intent(InOut) :: M   

    Character(Len=256)  :: messages(3), error
    
    messages(1)=error_set
    
    If (Trim(inblock)=='&intramol_statistics') Then
      If (M%tag_species%fread) Then
         If (M%tag_species%fail) Then
            Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: Problems to read the "'//Trim(M%tag_species%type)//&
                                        &'" directive within "'//Trim(M%invoke%type)//'".' 
            Call info(messages, 2)
            Call error_stop(' ') 
         End If
      Else
        Write (messages(2),'(1x,a)')  'Problems in "'//Trim(inblock)//'": The user must define the species involved&
                                    & inside "'//Trim(M%invoke%type)//'" using the "species" directive, which is missing.' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else If (Trim(inblock)=='&intermol_statistics') Then 
      If (M%tag_species%fread) Then
         Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: the definition of the "'//Trim(M%tag_species%type)//&
                                     &'" directive within "'//Trim(M%invoke%type)//'" is not necessary.' 
         Write (messages(3),'(1x,a)') 'The statistical analysis is carried out using the "reference_tag" defined in&
                                     & "&selected_nonreactive_species". Please remove "'//Trim(M%tag_species%type)//&
                                     &'" from this block.' 

         Call info(messages, 3)
         Call error_stop(' ') 
      End If
      If (M%tag_species%fread) Then
         Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: the definition of the "'//Trim(M%tag_species%type)//&
                                     &'" directive within "'//Trim(M%invoke%type)//'" is not required.' 
         Write (messages(3),'(1x,a)') 'Please remove "'//Trim(M%tag_species%type)//'" from this block.' 

         Call info(messages, 3)
         Call error_stop(' ') 
      End If
    End If 
 
    ! Error message just in case....
    error=Trim(messages(1))//' Check "'//Trim(M%invoke%type)//'" inside "'//Trim(inblock)//'".'
    
    !Check lower_bound, upper_bound and delta
    If (.Not. M%lower_bound%fread) Then
      M%lower_bound%tag='lower_bound'
    End If
    If (.Not. M%upper_bound%fread) Then
      M%upper_bound%tag='upper_bound'
    End If
    If (.Not. M%delta%fread) Then
      M%delta%tag='delta'
    End If
    
    If (Trim(M%invoke%type) == '&distance_parameters') Then
      Call check_length_directive(M%lower_bound, error, .True., 'directive')
      Call check_length_directive(M%upper_bound, error, .True., 'directive')
      Call check_length_directive(M%delta, error, .True., 'directive')
      If (M%lower_bound%value >= M%upper_bound%value) Then
        Write (messages(2),'(1x,a)')  'Problems with "'//Trim(M%invoke%type)//'" in "'//Trim(inblock)//'": The value of&
                                    & "upper_bound" must be larger than "lower_bound" (make sure this is the case if&
                                    & you use different units)' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else If (Trim(M%invoke%type) == '&angle_parameters') Then
      Call check_angle_directive(M%lower_bound, error, .True., 'directive')
      Call check_angle_directive(M%upper_bound, error, .True., 'directive')
      Call check_angle_directive(M%delta, error, .True., 'directive')
      If (M%lower_bound%value >= M%upper_bound%value) Then
        Write (messages(2),'(1x,a)')  'Problems with "'//Trim(M%invoke%type)//'" in "'//Trim(inblock)//'": The value of&
                                    & "upper_bound" must be larger than "lower_bound" (make sure this is the case if&
                                    & you use different units)' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    End If  
     
  End Subroutine check_settings_geom_param
  
  Subroutine check_angle_directive(T, error_set, kill, type_directive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check angle related directives
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill
    Character(Len=*),        Intent(In   )  :: type_directive
  
    Character(Len=256)  :: messages(2)
    
    If (T%fread) Then
      If (T%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "'&
                                      &//Trim(T%tag)//'" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (T%value < epsilon(1.0_wp)) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Input value for "'//Trim(T%tag)//&
                                    &'" MUST be larger than zero'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        Call capital_to_lower_case(T%units)
        If (Trim(T%units) /= 'rads' .And. Trim(T%units) /= 'degrees') Then
           If (Trim(type_directive) /= 'inblock') Then
             Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                      & 'Units for directive "'//Trim(T%tag)//'" must be "Degrees" or "Rads".&
                                      & Have you defined the units? Please review.'
           Else
             Write (messages(1),'(2(1x,a))')  Trim(error_set), '. Units for angles must be "Degrees" or "Rads".&
                                           & Have you defined the units? Please review'
           End If
           Call info(messages, 1)
           Call error_stop(' ')
        End If
        ! Transform to Angstrom
        If (Trim(T%units) == 'rads') Then
           T%value=Rads_to_degrees * T%value
        End If
      End If
    Else 
      If (kill) Then
        If (Trim(type_directive) /= 'inblock') Then
          Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(T%tag)//'" directive'
        Else
          Write (messages(1),'(1x,a)')  Trim(error_set)
        End If  
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_angle_directive  
    

  Subroutine check_nonreact_stat_settings(files, nonreactive_species, nonreact_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings for statistical analysis of nonreactive
    ! species geometries: intra and inter molecular info
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(In   ) :: files(:)
    Type(nonreactive_type),   Intent(In   ) :: nonreactive_species
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data

    Character(Len=256)   :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&selected_nonreactive_species" block.'    
    
    ! Check intramol_stat_settings 
    If (nonreact_stat_data%intra_geom%invoke%fread) Then
      If (nonreact_stat_data%intra_geom%angle%invoke%fread) Then
        If (nonreactive_species%atoms_per_species < 3) Then
        Write (messages(2),'(1x,a)')  'Problems in "'//Trim(nonreact_stat_data%intra_geom%invoke%type)//&
                                     &'": it is not possible to define an internal angle when nonreactive species are diatomic.& 
                                     & Remove "'//Trim(nonreact_stat_data%intra_geom%angle%invoke%type)//'"'
        Call info(messages, 2)
        Call error_stop(' ') 
        End If
      End If
      If (nonreact_stat_data%intra_geom%dist%invoke%fread) Then
        Call check_settings_geom_param(messages(1), nonreact_stat_data%intra_geom%invoke%type, nonreact_stat_data%intra_geom%dist)
        Call check_intramol_stat_settings(messages(1),nonreactive_species, nonreact_stat_data%intra_geom%dist)
      End If
      If (nonreact_stat_data%intra_geom%angle%invoke%fread) Then
        Call check_settings_geom_param(messages(1), nonreact_stat_data%intra_geom%invoke%type, nonreact_stat_data%intra_geom%angle)
        Call check_intramol_stat_settings(messages(1),nonreactive_species, nonreact_stat_data%intra_geom%angle)  
      End If
      If ((.Not. nonreact_stat_data%intra_geom%dist%invoke%fread) .And. &
          (.Not. nonreact_stat_data%intra_geom%angle%invoke%fread)) Then
           Write (messages(2),'(1x,a)')  'Empty "'//Trim(nonreact_stat_data%intra_geom%invoke%type)//&
                                    &'" block! Please define "&distance_parameters" and/or "&angle_parameters",&
                                    & or remove the block.'
           Call info(messages, 2)
           Call error_stop(' ') 
      End If    
    End If

    ! Check intermol_statistics
    If (nonreact_stat_data%inter_geom%invoke%fread) Then
      ! Check if only reference tags will be condired as NNs
      If (nonreact_stat_data%inter_geom%only_ref_tags_as_nn%fread) Then
        If (nonreact_stat_data%inter_geom%only_ref_tags_as_nn%fail) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                    & "only_ref_tags_as_nn" (choose either .True. or .False.)'
          Call info(messages,1)
          Call error_stop(' ')
        End If
      Else
        nonreact_stat_data%inter_geom%only_ref_tags_as_nn%stat=.False.
      End If
    
      If (nonreact_stat_data%inter_geom%dist%invoke%fread) Then
        Call check_settings_geom_param(messages(1), nonreact_stat_data%inter_geom%invoke%type, nonreact_stat_data%inter_geom%dist)
      End If
      If (nonreact_stat_data%inter_geom%angle%invoke%fread) Then
        Call check_settings_geom_param(messages(1), nonreact_stat_data%inter_geom%invoke%type, nonreact_stat_data%inter_geom%angle)
      End If
      If ((.Not. nonreact_stat_data%inter_geom%dist%invoke%fread) .And. &
          (.Not. nonreact_stat_data%inter_geom%angle%invoke%fread)) Then
           Write (messages(2),'(1x,a)')  'Empty "'//Trim(nonreact_stat_data%inter_geom%invoke%type)//&
                                     &'" block! Please define "&distance_parameters" and/or "&angle_parameters",&
                                     & or remove the block.'
           Call info(messages, 2)
           Call error_stop(' ') 
      End If    
    End If
 
  End Subroutine check_nonreact_stat_settings
  
  
  Subroutine geometry_statistics_nonreactive_species(files, traj_data, nonreactive_species, nonreact_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to inter and/or intramolecular distance and angles of the
    ! nonreactive species
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Type(file_type),          Intent(InOut) :: files(:)
    Type(traj_type),          Intent(InOut) :: traj_data
    Type(nonreactive_type),   Intent(In   ) :: nonreactive_species
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data
    
    Character(Len=256)  :: message
    Logical             :: flag_inter_geo_stat

    If (nonreact_stat_data%intra_geom%invoke%fread) Then
      If (nonreact_stat_data%intra_geom%dist%invoke%fread) Then
        Call obtain_intramol_geom_stat(files, traj_data, nonreactive_species%atoms_per_species,&
                                    & nonreact_stat_data%intra_geom%dist)
      End If
      If (nonreact_stat_data%intra_geom%angle%invoke%fread) Then
        Call obtain_intramol_geom_stat(files, traj_data, nonreactive_species%atoms_per_species,&
                                    & nonreact_stat_data%intra_geom%angle)
      End If
    End If

    ! Compute intermolecular properties for nonreactive species
    If (nonreact_stat_data%inter_geom%invoke%fread) Then
      Call find_neighbours_nonreactive_species(traj_data, nonreact_stat_data%inter_geom%only_ref_tags_as_nn%stat,&
                                             & flag_inter_geo_stat)
      If (nonreact_stat_data%inter_geom%dist%invoke%fread .And. flag_inter_geo_stat)  Then
        Call obtain_intermol_geom_stat(files, traj_data, nonreact_stat_data%inter_geom%dist, 1)
        Call obtain_intermol_geom_stat(files, traj_data, nonreact_stat_data%inter_geom%dist, 2)
        Write (message,'(1x,a)') 'The probability distribution of the intermolecular distances were printed to files "'&
                                &//Trim(files(FILE_INTERMOL_DISTANCES_NN1)%filename)//'" and "'&
                                &//Trim(files(FILE_INTERMOL_DISTANCES_NN2)%filename)//'"'
        Call info(message, 1)
        Write (message,'(1x,a)') 'which separetely consider the first and the second nearest nonreactive species, respectively.'
        Call info(message, 1)
        Call info(' ', 1)
      End If
      If (nonreact_stat_data%inter_geom%angle%invoke%fread .And. flag_inter_geo_stat) Then
        Call obtain_intermol_geom_stat(files, traj_data, nonreact_stat_data%inter_geom%angle)
      End If
    End If  
    
    Call refresh_out(files)
  
  End Subroutine geometry_statistics_nonreactive_species
  
  Subroutine obtain_intermol_geom_stat(files, traj_data, M, num_nn)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of geometrical
    ! parameters (distance and angles) between three closest nonreactive
    ! species, with the criteria defined in the &intermol_statistics
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),            Intent(InOut)  :: files(:)
    Type(traj_type),            Intent(InOut)  :: traj_data
    Type(geo_param_type),       Intent(In   )  :: M
    Integer(Kind=wi), Optional, Intent(In   )  :: num_nn

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k, k1, k2, mk
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(2), message
    Logical            :: falloc, flag

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: u(3), u2(3), angle 
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(M%upper_bound%value-M%lower_bound%value)/M%delta%value)
    
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining geometry statistics&
                                & of nonreactive species. Analysis will not be executed.'
      Call info(message, 1)                          
      falloc=.False.
    Else
      falloc=.True.
    End If
 
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0
      
      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
        h=0
        num_var=0
        Do  j= 1, traj_data%Nmax_species
          If (traj_data%region%define%fread) Then
            mk=traj_data%species(i,j)%list(1)
            Call within_region(traj_data, i, mk, flag)
          Else
            flag=.True.
          End If

          If (traj_data%species(i,j)%alive .And. flag) Then
            If (Trim(M%invoke%type) == '&distance_parameters') Then
              k=traj_data%species(i,j)%list(1)
              k2=traj_data%species(i,j)%nn(num_nn)
              u=traj_data%config(i,k2)%r-traj_data%config(i,k)%r          
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              If (norm2(u) >= M%lower_bound%value .And. norm2(u) <= M%upper_bound%value) Then
                mk=Floor((norm2(u)-M%lower_bound%value)/M%delta%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            Else If (Trim(M%invoke%type) == '&angle_parameters') Then
              k=traj_data%species(i,j)%list(1)
              k1=traj_data%species(i,j)%nn(1)
              k2=traj_data%species(i,j)%nn(2)
              u =traj_data%config(i,k1)%r-traj_data%config(i,k)%r
              u2=traj_data%config(i,k2)%r-traj_data%config(i,k)%r
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              Call check_PBC(u2, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              angle=Acos((u(1)*u2(1)+u(2)*u2(2)+u(3)*u2(3))/norm2(u)/norm2(u2))*Rads_to_degrees
              If (angle >= M%lower_bound%value .And. angle <= M%upper_bound%value) Then
                mk=Floor((angle-M%lower_bound%value)/M%delta%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            End If
          End If
        End Do
        
        If(num_var /= 0) Then
          accum=accum+num_var
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do mk= 1, nbins 
            d(mk)= d(mk)+Real(h(mk),Kind=wp)/num_var
          End Do
        End If
      End Do

      ! Print results
      If (accum /= 0) Then
         Do mk=1, nbins 
           d(mk)=d(mk)/net_frames/M%delta%value
         End Do
       
        ! Print File
        If (Trim(M%invoke%type) == '&distance_parameters') Then
          If (num_nn == 1) Then
            Open(Newunit=files(FILE_INTERMOL_DISTANCES_NN1)%unit_no, File=files(FILE_INTERMOL_DISTANCES_NN1)%filename,&
                              &Status='Replace')
            iunit=files(FILE_INTERMOL_DISTANCES_NN1)%unit_no
            Write (iunit,'(a)') '#  Probability distribution of the intermolecular distances&
                               & using only the first nearest nonreactive species and the settings of '//Trim(M%invoke%type)  
            Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          Else If (num_nn == 2) Then
            Open(Newunit=files(FILE_INTERMOL_DISTANCES_NN2)%unit_no, File=files(FILE_INTERMOL_DISTANCES_NN2)%filename,&
                              &Status='Replace')
            iunit=files(FILE_INTERMOL_DISTANCES_NN2)%unit_no
            Write (iunit,'(a)') '#  Probability distribution of the intermolecular distances&
                               & using only the second nearest nonreactive species and the settings of '//Trim(M%invoke%type)
 
            Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          End If
        Else If (Trim(M%invoke%type) == '&angle_parameters') Then
          Open(Newunit=files(FILE_INTERMOL_ANGLES)%unit_no, File=files(FILE_INTERMOL_ANGLES)%filename, Status='Replace')
          iunit=files(FILE_INTERMOL_ANGLES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of intermolecular angles using the first and second nearest&
                            & nonreactive species and the settings of '//Trim(M%invoke%type)
          Write (iunit,'(a)') '#  Value [Degrees]      Probability [1/Degrees]' 
          Write (message,'(1x,a)') 'The probability distribution of the intermolecular angles was printed to the "'&
                                  &//Trim(files(FILE_INTERMOL_ANGLES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        End If
          Do mk=1, nbins
            Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
          End Do
      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)')   '   WARNING: the statistics for the requested intermolecular geometry could not be executed'
        Write (messages(2),'(1x,a)') '   Please verify the settings for the '//Trim(M%invoke%type)//' in &intermol_statistics'
        Call info(messages, 2)
        Write (messages(1),'(1x,a)')   '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Deallocate arrays   
      Deallocate(d,h)
    End If
    
    Call refresh_out(files)
    
  End Subroutine obtain_intermol_geom_stat

  
  Subroutine obtain_intramol_geom_stat(files, traj_data, atoms_per_species, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of internal geometrical
    ! parameters (distance and angles) for nonreactive species 
    ! as defined in the &intramol_statistics
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut)  :: files(:)
    Type(traj_type),      Intent(InOut)  :: traj_data
    Integer(Kind=wi),     Intent(In   )  :: atoms_per_species
    Type(geo_param_type), Intent(In   )  :: M

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k1, k2, k3, mk, l, l1, l2
    Integer(Kind=wi)  :: ni(3), nj(2)
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(2), message
    Logical           :: falloc, flag, flag1, flag2

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: u(3), u2(3), angle 
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(M%upper_bound%value-M%lower_bound%value)/M%delta%value)
    
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining geometry statistics&
                                & of nonreactive species. Analysis will not be executed.'
      Call info(message, 1)                          
      falloc=.False.
    Else
      falloc=.True.
    End If
 
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0
      
      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
        h=0
        num_var=0
        Do  j= 1, traj_data%Nmax_species
          If (traj_data%region%define%fread) Then
            mk=traj_data%species(i,j)%list(1)
            Call within_region(traj_data, i, mk, flag)
          Else
            flag=.True.
          End If

          If (traj_data%species(i,j)%alive .And. flag) Then
            If (Trim(M%invoke%type) == '&distance_parameters') Then 
              Do k1= 1, atoms_per_species
                ni(1)=traj_data%species(i,j)%list(k1)
                Do k2= k1+1, atoms_per_species
                  ni(2)=traj_data%species(i,j)%list(k2)
                  flag1=(traj_data%config(i,ni(1))%element==M%species(1)) .And.&
                        (traj_data%config(i,ni(2))%element==M%species(2))
                  flag2=(traj_data%config(i,ni(1))%element==M%species(2)) .And.&
                        (traj_data%config(i,ni(2))%element==M%species(1))      
                  If (flag1 .Or. flag2) Then
                    u=traj_data%config(i,ni(1))%r-traj_data%config(i,ni(2))%r          
                    Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                    If (norm2(u) >= M%lower_bound%value .And. norm2(u) <= M%upper_bound%value) Then
                      mk=Floor((norm2(u)-M%lower_bound%value)/M%delta%value)+1
                      If (mk <= nbins) Then
                        h(mk)=h(mk)+1
                        num_var=num_var+1
                      End If
                    End If
                  End If
                End Do
              End Do
              
            Else If (Trim(M%invoke%type) == '&angle_parameters') Then
              Do k1= 1, atoms_per_species
                ni(1)=traj_data%species(i,j)%list(k1)
                Do k2= k1+1, atoms_per_species
                  ni(2)=traj_data%species(i,j)%list(k2)
                  Do k3= k2+1, atoms_per_species
                    ni(3)=traj_data%species(i,j)%list(k3)
                    Do l= 1, 3
                      If (traj_data%config(i,ni(l))%element==M%species(2)) Then
                        Do l1= 1, atoms_per_species
                          nj(1)=traj_data%species(i,j)%list(l1)
                          Do l2= l1+1, atoms_per_species
                            nj(2)=traj_data%species(i,j)%list(l2)
                            If (l1 /= l .And. l2 /= l) Then
                              flag1=(traj_data%config(i,nj(1))%element==M%species(1)) .And.&
                                    (traj_data%config(i,nj(2))%element==M%species(3))
                              flag2=(traj_data%config(i,nj(1))%element==M%species(3)) .And.&
                                    (traj_data%config(i,nj(2))%element==M%species(1))
                              If (flag1 .Or. flag2) Then
                                u =traj_data%config(i,nj(1))%r-traj_data%config(i,ni(l))%r
                                u2=traj_data%config(i,nj(2))%r-traj_data%config(i,ni(l))%r
                                Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                                Call check_PBC(u2, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                                angle=Acos((u(1)*u2(1)+u(2)*u2(2)+u(3)*u2(3))/norm2(u)/norm2(u2))*Rads_to_degrees
                                If (angle >= M%lower_bound%value .And. angle <= M%upper_bound%value) Then
                                  mk=Floor((angle-M%lower_bound%value)/M%delta%value)+1
                                  If (mk <= nbins) Then
                                    h(mk)=h(mk)+1
                                    num_var=num_var+1
                                  End If
                                End If
                              End If
                            End If
                          End Do
                        End Do
                      End If
                    End Do  
                  End Do
                End Do
              End Do
            End If
          End If
        End Do

        If(num_var /= 0) Then
          accum=accum+num_var
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do mk= 1, nbins 
            d(mk)= d(mk)+Real(h(mk),Kind=wp)/num_var
          End Do
        End If
      End Do

      ! Print results
      If (accum /= 0) Then
        Do mk=1, nbins 
          d(mk)=d(mk)/net_frames/M%delta%value
        End Do
      
        ! Print File
        If (Trim(M%invoke%type) == '&distance_parameters') Then
          Open(Newunit=files(FILE_INTRAMOL_DISTANCES)%unit_no, File=files(FILE_INTRAMOL_DISTANCES)%filename, Status='Replace')
          iunit=files(FILE_INTRAMOL_DISTANCES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of the intramolecular distances&
                             & using the settings of '//Trim(M%invoke%type)  
          Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          Write (message,'(1x,a)') 'The probability distribution of the intramolecular distances was printed to the "'&
                                  &//Trim(files(FILE_INTRAMOL_DISTANCES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        Else If (Trim(M%invoke%type) == '&angle_parameters') Then
          Open(Newunit=files(FILE_INTRAMOL_ANGLES)%unit_no, File=files(FILE_INTRAMOL_ANGLES)%filename, Status='Replace')
          iunit=files(FILE_INTRAMOL_ANGLES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of the intramolecular angles using the settings of '//Trim(M%invoke%type)
          Write (iunit,'(a)') '#  Value [Degrees]      Probability [1/Degrees]' 
          Write (message,'(1x,a)') 'The probability distribution of the intramolecular angles was printed to the "'&
                                  &//Trim(files(FILE_INTRAMOL_ANGLES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        End If
          Do mk=1, nbins
            Write(iunit,'(2x,f12.4,6x,f13.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
          End Do
      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)')   '   WARNING: the statistics for the requested intramolecular geometry could not be executed'
        Write (messages(2),'(1x,a)') '   Please verify the settings for the '//Trim(M%invoke%type)//' in &intramol_statistics'
        Call info(messages, 2)
        Write (messages(1),'(1x,a)')   '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Deallocate arrays   
      Deallocate(d,h)
    End If

    Call refresh_out(files)
    
  End Subroutine obtain_intramol_geom_stat

  Subroutine find_neighbours_nonreactive_species(traj_data, only_ref_tags_as_nn, flag_exec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the two closest nonreactive species to a
    ! nonreactive species
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),  Intent(InOut) :: traj_data
    Logical,          Intent(In   ) :: only_ref_tags_as_nn
    Logical,          Intent(  Out) :: flag_exec    

    Integer(Kind=wi)  :: accum, net_frames
    Integer(Kind=wi)  :: i, j, k, mk, indx1, indx2, mindx1, mindx2
    
    Character(Len=256) :: messages(3)
    Logical            :: flag, flag1, flag2

    Real(Kind=wp)  :: min_dist1, min_dist2, u(3)
    Logical        :: modified, finclude
    
    net_frames=0
    Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
      accum=0
      Do  j= 1, traj_data%Nmax_species
        min_dist1 = Huge(1.0_wp)
        min_dist2 = Huge(1.0_wp)
        If (traj_data%region%define%fread) Then
          mk=traj_data%species(i,j)%list(1)
          Call within_region(traj_data, i, mk, flag)
        Else
          flag=.True.
        End If

        If (traj_data%species(i,j)%alive .And. flag) Then
          accum=accum+1
          indx1=traj_data%species(i,j)%list(1)
          mindx1=indx1
          mindx2=indx1
          Do  k= 1, traj_data%Nmax_species
            If (only_ref_tags_as_nn) Then
              If (traj_data%species(i,k)%alive) Then
                finclude=.True.
              Else
                finclude=.False.
              End If
            Else
              finclude=.True.
            End If
            If (k /= j .And. finclude) Then  
              indx2= traj_data%species(i,k)%list(1)
              u= traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              flag1= norm2(u) < min_dist1
              flag2= norm2(u) < min_dist2
              If (flag1 .And. flag2) Then
                min_dist2=min_dist1
                mindx2=mindx1
                min_dist1=norm2(u)
                mindx1=indx2
              Else If ((.Not. flag1) .And. flag2) Then
                min_dist2=norm2(u)
                mindx2=indx2
              End If
            End If  
          End Do
          traj_data%species(i,j)%nn(1)=mindx1
          traj_data%species(i,j)%nn(2)=mindx2
        End If
      End Do
      If (accum > 2 ) Then
        net_frames=net_frames+1
      End If
    End Do
    
    If (net_frames==0) Then
      Write (messages(1),'(1x,a)') '*************************************************************************************'
      Call info(messages, 1)
      Write (messages(1),'(1x,a)') '   WARNING: it looks the system has two or less nonreactive species along the trajectory (!?)'
      Write (messages(2),'(1x,a)') '   The intermolecular analysis of geometry parameters will not be executed.' 
      Write (messages(3),'(1x,a)') '   Please review the systems and the settings.'                        
      Call info(messages, 3)
      Write (messages(1),'(1x,a)') '************************************************************************************'
      Call info(messages, 1)
      flag_exec=.False.
    Else
      flag_exec=.True.
    End If
    
  End Subroutine find_neighbours_nonreactive_species 
  
  
End Module nonreact_stat
