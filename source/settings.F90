!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that:
! - reads the SETTINGS file and defines the settings for analysis
! - checks correctness of defined directives
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author      - i.scivetti   2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module settings 

  Use atomic_model,        Only: model_type, &
                                 check_length_directive
                                 
  Use constants,           Only: Bohr_to_A,  &
                                 chemsymbol, & 
                                 NPTE, &
                                 max_at_species,&
                                 max_components,&
                                 max_unchanged_atoms
                           
  Use fileset,             Only: file_type, &
                                 FILE_SET, &  
                                 FILE_OUT, & 
                                 refresh_out
                           
  Use coord_distr,         Only: coord_distr_type, &
                                 read_coord_distrib, &
                                 check_coord_distrib
                                 
  Use nndist_distr,        Only: nndist_distr_type, &
                                 read_selected_nn_distances,&
                                 check_selected_nn_distances                             
                                 
  Use nonreact_stat,       Only: nonreact_stat_type, &
                                 read_geom_param_nonreactive_species, &
                                 check_nonreact_stat_settings
                                 
  Use input_types,         Only: in_integer, &
                                 in_logic,   &
                                 in_scalar,  &
                                 in_param,   & 
                                 in_string                                
                                 
  Use msd,                 Only: msd_type, &
                                 read_msd, &
                                 check_msd
                           
  Use numprec,             Only: wi, &
                                 wp
                           
  Use ocf,                 Only: ocf_type, &
                                 read_ocf_settings, &
                                 check_ocf_nonreactive_species,&
                                 check_ocf_reactive_species
                                  
  Use process_data,        Only: capital_to_lower_case, &
                                 check_for_rubbish, &
                                 get_word_length, &
                                 remove_symbols, &
                                 set_read_status, &
                                 prevent_segmentation, &
                                 check_end
                                                                 
 Use residence_times,      Only: restimes_type, &
                                 read_residence_times, &
                                 check_residence_times
                                 
  Use rdf,                 Only: rdf_type, &
                                 read_rdf, &
                                 check_rdf
                                 
  Use spcf,                Only: spcf_type, &
                                 read_spcf,&
                                 check_spcf 
                                 
  Use tcf,                 Only: tcf_type, &
                                 read_tcf, &
                                 check_tcf
                                 
  Use trajectory,          Only: traj_type, &
                                 check_time_directive

  Use unchanged_chemistry, Only: unchanged_type, &
                                 read_track_unchanged_chemistry, &
                                 check_initial_unchanged_labels, &
                                 check_unchanged_chemistry
                            
  Use unit_output,         Only: error_stop,&
                                 info
                               

  Implicit None
  
  Public :: read_settings, check_settings_for_trajectory_analysis

Contains

  Subroutine read_settings(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data, coord_distr_data,&
                         & nonreact_stat_data, nndist_distr_data, unchanged_data, rdf_data, restimes_data, tcf_data,&
                         & spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings from SETTINGS file.
    ! Lines starting with # are ignored and assumed as comments. 
    ! If a directive is identified during the reading of the file, subroutine "set_read_fail" 
    ! assigns fread=.True. On the contrary, the subroutine assigns fail=.True. (fail=.False.) 
    ! if the format/syntax for the directive is correct (incorrect)
    ! If the directive is repeated the execution is aborted via subroutine duplication 
    ! 
    ! author        - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(InOut) :: files(:)
    Type(model_type),         Intent(InOut) :: model_data 
    Type(traj_type),          Intent(InOut) :: traj_data
    Type(ocf_type),           Intent(InOut) :: ocf_nonreactive
    Type(ocf_type),           Intent(InOut) :: ocf_reactive
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(coord_distr_type),   Intent(InOut) :: coord_distr_data
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data    
    Type(nndist_distr_type),  Intent(InOut) :: nndist_distr_data
    Type(unchanged_type),     Intent(InOut) :: unchanged_data
    Type(rdf_type),           Intent(InOut) :: rdf_data
    Type(restimes_type),      Intent(InOut) :: restimes_data
    Type(tcf_type),           Intent(InOut) :: tcf_data
    Type(spcf_type),          Intent(InOut) :: spcf_data
 
    Logical            :: safe
    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io, iunit
  
    Character(Len=256)  :: message

    Character(Len=32 )  :: set_file
    Character(Len=32 )  :: set_error

    set_file = Trim(files(FILE_SET)%filename)
    set_error = '***ERROR in the '//Trim(set_file)//' file.'

    ! Open the SETTINGS file with settings
    Inquire(File=files(FILE_SET)%filename, Exist=safe)
    
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(set_file), '(settings for analysis) not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_SET)%unit_no, File=Trim(set_file), Status='old')
      iunit=files(FILE_SET)%unit_no 
    End If

     Read (iunit, Fmt=*, iostat=io) word
     ! If nothing is found, complain and abort
     If (is_iostat_end(io)) Then
       Write (message,'(3(1x,a))') Trim(set_error), Trim(set_file), 'file seems to be empty?. Please check'
       Call error_stop(message)
     End If
     ! Check header has "#" as the first character 
     If (word(1:1)/='#') Then
       Write (message,'(4(1x,a))') Trim(set_error), 'Heading comment in file', Trim(set_file), & 
                                  'is required and MUST be preceded with the symbol "#"'
       Call error_stop(message)
     End If

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Exit
      end If
      Call check_for_rubbish(iunit, Trim(set_file)) 
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
        ! Do nothing if line is a comment of we have an empty line
        Read (iunit, Fmt=*, iostat=io) word

      ! Model related variables 
      Else If (word(1:length) == 'reactive_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_chemistry%stat
        Call set_read_status(word, io, model_data%reactive_chemistry%fread, model_data%reactive_chemistry%fail)

      Else If (word(1:length) == 'geometry_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%geometry_format%type
        Call set_read_status(word, io, model_data%geometry_format%fread, model_data%geometry_format%fail, &
                           & model_data%geometry_format%type)

      Else If (word(1:length) == 'cell_units') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%config%cell_units%type
        Call set_read_status(word, io, model_data%config%cell_units%fread, model_data%config%cell_units%fail)
        
      Else If (word(1:length) == 'position_units') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%config%position_units%type
        Call set_read_status(word, io, model_data%config%position_units%fread, model_data%config%position_units%fail)
        
      Else If (word(1:length) == '&reference_composition') Then
        Read (iunit, Fmt=*, iostat=io) model_data%reference_composition%invoke%type
        Call set_read_status(word, io, model_data%reference_composition%invoke%fread, model_data%reference_composition%invoke%fail)
        ! Read information inside the block
        Call read_reference_composition(iunit, model_data)

      Else If (word(1:length) == '&simulation_cell') Then
        Read (iunit, Fmt=*, iostat=io) model_data%config%simulation_cell%type
        Call set_read_status(word, io, model_data%config%simulation_cell%fread, model_data%config%simulation_cell%fail)
        ! Read information inside the block
        Call read_input_cell(iunit, model_data)  

      Else If (word(1:length) == '&reactive_species') Then
        Read (iunit, Fmt=*, iostat=io) model_data%reactive_species%invoke%type
        Call set_read_status(word, io, model_data%reactive_species%invoke%fread, model_data%reactive_species%invoke%fail)
        ! Read information inside the block
        Call read_reactive_species(iunit, model_data)

      Else If (word(1:length) == '&selected_nonreactive_species') Then
        Read (iunit, Fmt=*, iostat=io) model_data%nonreactive_species%invoke%type
        Call set_read_status(word, io, model_data%nonreactive_species%invoke%fread, model_data%nonreactive_species%invoke%fail)
        !Read information inside the block
        Call read_nonreactive_species(iunit, model_data)

      ! Trajectory related variables  
      Else If (word(1:length) == 'print_track_chemistry') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%print_track_chemistry%stat
       Call set_read_status(word, io, traj_data%print_track_chemistry%fread, traj_data%print_track_chemistry%fail)
       
      Else If (word(1:length) == 'print_retagged_trajectory') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%print_retagged_trajectory%stat
       Call set_read_status(word, io, traj_data%print_retagged_trajectory%fread, traj_data%print_retagged_trajectory%fail)       
      
      Else If (word(1:length) == 'ensemble') Then 
        Read (iunit, Fmt=*, iostat=io) word, traj_data%ensemble%type
        Call set_read_status(word, io, traj_data%ensemble%fread, traj_data%ensemble%fail, &
                           & traj_data%ensemble%type)

      Else If (word(1:length) == 'recorded_timestep') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%timestep%tag, traj_data%timestep%value,&
                                      &traj_data%timestep%units
        Call set_read_status(word, io, traj_data%timestep%fread, traj_data%timestep%fail)

      Else If (word(1:length) == '&segment_trajectory') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%seg_analysis%invoke%type
        Call set_read_status(word, io, traj_data%seg_analysis%invoke%fread, traj_data%seg_analysis%invoke%fail)
        Call read_segment_trajectory(iunit, traj_data)

      Else If (word(1:length) == '&region') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%region%define%type
        Call set_read_status(word, io, traj_data%region%define%fread, traj_data%region%define%fail)
        !Read information inside the block
        Call read_region(iunit, traj_data)
        
      Else If (word(1:length) == '&reactive_analysis') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%reactive_analysis%type
        Call set_read_status(word, io, traj_data%reactive_analysis%fread, traj_data%reactive_analysis%fail)
        !Read information inside the block
        Call read_reactive_analysis(iunit, ocf_reactive, restimes_data, tcf_data, spcf_data)
        
      Else If (word(1:length) == '&nonreactive_analysis') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%nonreactive_analysis%type
        Call set_read_status(word, io, traj_data%nonreactive_analysis%fread, traj_data%nonreactive_analysis%fail)
        !Read information inside the block
        Call read_nonreactive_analysis(iunit, model_data, ocf_nonreactive, msd_data, nonreact_stat_data)
        
      Else If (word(1:length) == '&general_analysis') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%general_analysis%type
        Call set_read_status(word, io, traj_data%general_analysis%fread, traj_data%general_analysis%fail)
        !Read information inside the block
        Call read_general_analysis(iunit, coord_distr_data, nndist_distr_data, rdf_data, unchanged_data)
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Directive not recognised. Inform and kill 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Else
        Call check_functionality(word, set_error)
        If (word(1:1)=='&') Then
          Write (message,'(1x,4a)') Trim(set_error), ' Unknown directive found: "', Trim(word),&
                                  &'. Do you use "&" to define a block? If so,&
                                  & make sure the block is valid and has right syntax.'
        Else
          Write (message,'(1x,a)') Trim(set_error)//' Unknown directive found: "'//Trim(word)//'".&
                                  & Have you correctly defined the previous directives? Have you forgotten something maybe?'
        End If 
        Call error_stop(message)
      End If

    End Do
    ! Close file
    Close(files(FILE_SET)%unit_no)

  End Subroutine read_settings
  
   Subroutine check_functionality(word, set_error)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Subroutine to check if a functionality has been defined outside any
     ! of the following blocks:
     ! - &reactive_analysis
     ! - &nonreactive_analysis
     ! - &general_analysis
     !
     ! author    - i.scivetti June 2025
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Character(Len=256), Intent(In   ) :: word
     Character(Len=32 ), Intent(In   ) :: set_error
     
     Character(Len=256)  :: message

     If (Trim(word)=='&ocf_reactive' .Or. &
         Trim(word)=='&tcf'                     .Or. &
         Trim(word)=='&spcf'                    .Or. &
         Trim(word)=='&residence_times') Then
         Write (message,'(1x,a)') Trim(set_error)//' Block "'//Trim(word)//'"&
                                  & must be defined within the "&reactive_analysis" block'
         Call error_stop(message)
     End If

     If (Trim(word)=='&msd'  .Or. &
         Trim(word)=='&ocf_nonreactive'  .Or. &
         Trim(word)=='&intermol_statistics' .Or. &
         Trim(word)=='&intramol_statistics') Then
         Write (message,'(1x,a)') Trim(set_error)//' Block "'//Trim(word)//'"&
                                  & must be defined within the "&nonreactive_analysis" block'
         Call error_stop(message)
     End If

     If (Trim(word)=='compute_amount') Then
         Write (message,'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//'"&
                                  & must be defined within the "&nonreactive_analysis" block'
         Call error_stop(message)
     End If
     
     If (Trim(word)=='&selected_nn_distances'  .Or. &
         Trim(word)=='&rdf'                    .Or. &
         Trim(word)=='&coord_distrib'          .Or. &
         Trim(word)=='&track_unchanged_chemistry') Then
         Write (message,'(1x,a)') Trim(set_error)//' Block "'//Trim(word)//'"&
                                  & must be defined within the "&general_analysis" block'
         Call error_stop(message)
     End If
     
   End Subroutine check_functionality

  Subroutine read_input_cell(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the simulation cell vectors of the input structure
    ! defined in &simulation_cell
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(model_type),   Intent(InOut) :: model_data
   
    Integer(Kind=wi)   :: io, i, j
    Character(Len=64 ) :: error_simulation_cell
    Character(Len=256) :: messages(2), word
    Logical            :: endblock

    error_simulation_cell = '***ERROR in &simulation_cell of SETTINGS file'
    Write (messages(1),'(a)') error_simulation_cell

    i=1
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&simulation_cell')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&simulation_cell') 
          Read (iunit, Fmt=*, iostat=io) (model_data%config%cell(i,j), j=1,3)
          If (io/=0) Then
            Write (messages(2),'(a,i2)') 'Problems with the definition of cell vector', i
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Not all the cell vectors for the&  
                                     & input structure have been defined. Please check.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do
 
    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&simulation_cell')
      Call capital_to_lower_case(word)
      If (word /= '&end_simulation_cell') Then
        If (word(1:1) /= '#') Then
          Write (messages(2),'(a)') 'Block for cell vectors must be closed with&
                                  & sentence &end_simulation_cell.'
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_input_cell

  Subroutine read_reactive_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the definitions to search and identify changes chemisty
    ! Information is read from the &reactive_species block
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
 
    set_error = '***ERROR in &reactive_species of SETTINGS file.'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_reactive_species" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_reactive_species') Exit
      Call check_for_rubbish(iunit, '&reactive_species')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'total_number') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%N0%value
        Call set_read_status(word, io, model_data%reactive_species%N0%fread, model_data%reactive_species%N0%fail)

      Else If (word(1:length) == 'type') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%type%type
        Call set_read_status(word, io, model_data%reactive_species%type%fread,&
                             model_data%reactive_species%type%fail,model_data%reactive_species%type%type)
        
       Else If (word(1:length) == '&search_environment') Then
         Read (iunit, Fmt=*, iostat=io) word
         Call set_read_status(word, io, model_data%reactive_species%search_envr%criteria%fread,&
                            & model_data%reactive_species%search_envr%criteria%fail)
         Call read_environment_settings(iunit, model_data)

       Else If (Trim(word)=='&bonding_criteria') Then
            Read (iunit, Fmt=*, iostat=io) word
            Call set_read_status(word, io, model_data%reactive_species%bonds%criteria%fread,&
                  & model_data%reactive_species%bonds%criteria%fail)
            model_data%reactive_species%bonds%criteria%stat = .True.
         Call read_bonding_criteria(iunit, model_data)
 
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&End_chemistry"?'
        Call error_stop(message)
      End If

    End Do
   
  End Subroutine read_reactive_species

  Subroutine read_bonding_criteria(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settings from the &bonding_criteria block and 
    ! identify the species to be tracked
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
     

    set_error = '***ERROR in the "&bonding_criteria" sub-block (within &reactive_species).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_bonding_criteria" to close the block.&
                                  & Check if directives are set correctly.'         
        Call info(messages, 1)                          
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_bonding_criteria') Exit
      Call check_for_rubbish(iunit, '&bonding_criteria')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'only_element') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%bonds%species%type
        Call set_read_status(word, io, model_data%reactive_species%bonds%species%fread,&
                           & model_data%reactive_species%bonds%species%fail)

      Else If (Trim(word)=='cutoff') Then
          Read (iunit, Fmt=*, iostat=io)  model_data%reactive_species%bonds%cutoff%tag,&
                                        & model_data%reactive_species%bonds%cutoff%value,&
                                        & model_data%reactive_species%bonds%cutoff%units
         Call set_read_status(word, io, model_data%reactive_species%bonds%cutoff%fread,&
                              model_data%reactive_species%bonds%cutoff%fail)

      Else If (Trim(word)=='number_of_bonds') Then
          Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%bonds%N0%value
          Call set_read_status(word, io, model_data%reactive_species%bonds%N0%fread,&
               & model_data%reactive_species%bonds%N0%fail)
        
       Else If (word(1:length) == '&extra_reactive_bonds') Then
          Read (iunit, Fmt=*, iostat=io) model_data%extra_bonds%invoke%type
          Call set_read_status(word, io, model_data%extra_bonds%invoke%fread, model_data%extra_bonds%invoke%fail)
          ! Read information inside the block
          Call read_extra_bonding(iunit, model_data)
        
      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//'" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the sub-block with "&end_bonding_criteria"?'
        Write (messages(3),'(1x,a)') 'Have you included the units for "cutoff"?'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

    End Do
   
  End Subroutine read_bonding_criteria
  
  Subroutine read_environment_settings(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the defined criteria to explore the  environment of 
    ! selected atoms
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, j, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
    

    set_error = '***ERROR in the sub-block &search_environment within &reactive_species (SETTINGS file).'
    Write (messages(2),'(1x,a)') 'Have you properly closed the sub-block with "&end_search_environment"?'
    Write (messages(3),'(1x,a)') 'Have you set the units for directive "cutoff"?'                         

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_search_environment" to close the block.&
                                  & Check if directives are set correctly.' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_search_environment') Exit
      Call check_for_rubbish(iunit, '&search_environment')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='include_tags') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%search_envr%N0_incl
        Call prevent_segmentation(iunit, io, word, model_data%reactive_species%search_envr%N0_incl,&
                                & 'max_components', max_components, set_error)
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%search_envr%N0_incl,&
                                          & (model_data%reactive_species%search_envr%tg_incl(j), j = 1,&
                                          & model_data%reactive_species%search_envr%N0_incl)
        Call set_read_status(word, io, model_data%reactive_species%search_envr%info_include%fread,&
                           & model_data%reactive_species%search_envr%info_include%fail)
        model_data%reactive_species%search_envr%info_include%stat = .True.

      Else If (Trim(word)=='exclude_pairs') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%search_envr%N0_excl 
        Call prevent_segmentation(iunit, io, word, model_data%reactive_species%search_envr%N0_excl,&
                                & 'max_components', max_components, set_error)
        Read (iunit, Fmt=*, iostat=io) word, model_data%reactive_species%search_envr%N0_excl,&
                                          & (model_data%reactive_species%search_envr%tg_excl(j), j = 1,&
                                          & model_data%reactive_species%search_envr%N0_excl)
        Call set_read_status(word, io, model_data%reactive_species%search_envr%info_exclude%fread,&
                           & model_data%reactive_species%search_envr%info_exclude%fail)
        model_data%reactive_species%search_envr%info_exclude%stat = .True.

      Else If (Trim(word)=='cutoff') Then
         Read (iunit, Fmt=*, iostat=io) model_data%reactive_species%search_envr%cutoff%tag, &
                                        model_data%reactive_species%search_envr%cutoff%value,&
                                        model_data%reactive_species%search_envr%cutoff%units 
         Call set_read_status(word, io, model_data%reactive_species%search_envr%cutoff%fread,&
                            & model_data%reactive_species%search_envr%cutoff%fail)

      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//&
                                    &'" is not recognised as a valid settings.'
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
   
  End Subroutine read_environment_settings

  Subroutine read_extra_bonding(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to extra bond settings
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, messages(3)
    Character(Len=256)  :: error_block
    Logical  :: error, endblock, fread 

    error= .False.
    error_block = '***ERROR in &extra_reactive_bonds (inside &bonding_criteria).' 
    Write (messages(1),'(a)') error_block 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&extra_reactive_bonds')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&extra_reactive_bonds')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%extra_bonds%N0

    If (Trim(word) /= 'types_of_bonds') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "types_of_bonds" is expected.'
      error=.True.
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "type_of_bonds"'
      error=.True.
    Else
      If (model_data%extra_bonds%N0<1) Then
        Write (messages(2),'(a)') 'The "type_of_bonds" directive MUST BE >= 1'
        error=.True.
      ElseIf (model_data%extra_bonds%N0>max_components) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider more than ', max_components,&
                                    & ' for "type_of_bonds"? Please check'
        error=.True.
      End If
    End If

    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If

    i=1
    Do While (i <= model_data%extra_bonds%N0)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&extra_reactive_bonds')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&extra_reactive_bonds')
  
        Read (iunit, Fmt=*, iostat=io) model_data%extra_bonds%tg1(i), model_data%extra_bonds%tg2(i),   &
                                       model_data%extra_bonds%bond(i)%value, model_data%extra_bonds%bond(i)%units

        If (io/=0) Then
          If (Trim(model_data%extra_bonds%tg1(i)) == '&end') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for extra bonds. Only ',  i-1,&
                                    &' species set out of ', model_data%extra_bonds%N0, ' (types_of_bonds)'
            Write (messages(3),'(a)') 'Please check. What is the value set for "type_of_bonds"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          Else  
            Write (messages(2),'(a,i3)') 'Problems to read bonding criteria ',  i
            Write (messages(3),'(a)') 'Please check. What is the value set for "type_of_bonds"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          End If
        Else
          Write (messages(1),'(a,i3)') Trim(error_block)//' Problems to read bonding criteria ',  i
          model_data%extra_bonds%bond(i)%fread= .True.
          Call check_length_directive(model_data%extra_bonds%bond(i), messages(1), .True., 'inblock')
        End If

        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&extra_reactive_bonds')
      Call capital_to_lower_case(word)
      If (word /= '&end_extra_reactive_bonds') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%extra_bonds%N0) Then 
            Write (messages(2),'(a)') 'Number of extra bonds specified is larger than&
                                     & the value given by directive "type_of_bonds"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_extra_reactive_bonds. Please check. Is the&
                                     & number of defined bond criteria the same as set for directive "type_of_bonds"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do
    
  End Subroutine read_extra_bonding

  Subroutine read_nonreactive_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the defined criteria to explore the selected species
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut)  :: model_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &selected_nonreactive_species block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_selected_nonreactive_species" to close the block.&
                                  & Check if directives are set correctly.'         
        Call info(messages, 1)                          
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_selected_nonreactive_species') Exit
      Call check_for_rubbish(iunit, '&selected_nonreactive_species')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='name') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%nonreactive_species%name%type
        Call set_read_status(word, io, model_data%nonreactive_species%name%fread, model_data%nonreactive_species%name%fail)

      Else If (Trim(word)=='reference_tag') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%nonreactive_species%reference_tag%type
        Call set_read_status(word, io, model_data%nonreactive_species%reference_tag%fread,&
                           & model_data%nonreactive_species%reference_tag%fail)

      Else If (Trim(word)=='bond_cutoff') Then
        Read (iunit, Fmt=*, iostat=io) model_data%nonreactive_species%bond_cutoff%tag,  &
                                       model_data%nonreactive_species%bond_cutoff%value, &
                                       model_data%nonreactive_species%bond_cutoff%units
                                       
        Call set_read_status(word, io, model_data%nonreactive_species%bond_cutoff%fread,&
                           & model_data%nonreactive_species%bond_cutoff%fail)

      Else If (word(1:length) == '&atomic_components') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, model_data%nonreactive_species%atomic_components%fread,&
                           & model_data%nonreactive_species%atomic_components%fail)
        Call read_components_nonreactive_species(iunit, model_data)

      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//&
                                  &'" is not recognised as a valid settings. See the "use_code.md" file.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_selected_nonreactive_species"?'
        Write (messages(3),'(1x,a)') 'Have you included the units for "bond_cutoff"?'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

    End Do
    
  End Subroutine read_nonreactive_species
  
  Subroutine read_components_nonreactive_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to extra bond settings
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, messages(3)
    Character(Len=256)  :: error_block
    Logical             :: endblock, fread 

    error_block = '***ERROR in the &atomic_components sub-block (wihtin &selected_nonreactive_species)'
    Write (messages(1),'(a)') error_block 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &selected_nonreactive_species)')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&atomic_components')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%nonreactive_species%num_components

    If (Trim(word) /= 'number_components') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "number_components" is expected.'
      Call info(messages, 2) 
      Call error_stop(' ')
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number_components"'
      Call info(messages, 2) 
      Call error_stop(' ')
    Else
      If (model_data%nonreactive_species%num_components<1) Then
        Write (messages(2),'(a)') 'The "number_components" directive MUST BE >= 1'
        Call info(messages, 2) 
        Call error_stop(' ')
      ElseIf (model_data%nonreactive_species%num_components>max_at_species) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider more than ', max_at_species,&
                                    & ' for "number_components"? Please check'
        Write (messages(3),'(a)') 'If you are sure of what you are doing, look for the parameter "max_at_species" in the code,&
                                  & increase its value as needed and recompile.'
        Call info(messages, 3)
        Call error_stop(' ')
      End If
    End If

    i=1
    Do While (i <= model_data%nonreactive_species%num_components)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &selected_nonreactive_species)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&atomic_components')
  
        Read (iunit, Fmt=*, iostat=io) model_data%nonreactive_species%element(i), model_data%nonreactive_species%N0_element(i)

        If (io/=0) Then
          If (Trim(model_data%nonreactive_species%element(i)) == '&end') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for extra bonds. Only ',  i-1,&
                                    &' species set out of ', model_data%nonreactive_species%num_components,&
                                    & ' (number_components)'
            Write (messages(3),'(a)') 'Please check. What is the value set for "number_components"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          Else  
            Write (messages(2),'(a,i3)') 'Problems to read the species component ',  i      
            Write (messages(3),'(a)') 'Please check. What is the value set for "number_components"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          End If 
        End If

        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &selected_nonreactive_species)')
      Call capital_to_lower_case(word)
      If (word /= '&end_atomic_components') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%nonreactive_species%num_components) Then 
            Write (messages(2),'(a)') 'Number of extra bonds specified is larger than&
                                     & the value given by directive "number_components"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_atomic_components. Please check. Is the&
                                     & number of defined bond criteria the same as set for directive "number_components"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do
    
  End Subroutine read_components_nonreactive_species

  Subroutine read_reference_composition(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the amount of atoms and chemical elements for each atomic
    ! species. Information must be defined in the &reference_composition block
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data
 
    Logical  :: endblock, loop, error
    Logical  :: header(3), error_duplication
 
    Integer(Kind=wi)   :: io, i, j, k, ilist, ic
    Character(Len=256) :: messages(10), word
    Character(Len=64 ) :: error_reference_composition
 
    error_reference_composition = '***ERROR in &reference_composition of SETTINGS file'
    Write (messages(1),'(a)') error_reference_composition
 
    header=.False.
    error_duplication=.False.
    error=.False.
    ilist=10
 
    Write (messages(3),'(1x,a)')    'The correct structure for the block must be:'
    Write (messages(4),'(1x,a)')    '&reference_composition'
    Write (messages(5),'(1x,a)')    '  atomic_species    Nsp'
    Write (messages(6),'(1x,a)')    '  tags      tg1    tg2    tg3   .... tgNsp'
    Write (messages(7),'(1x,a)')    '  elements  E_tg1  E_tg2  E_tg3 .... E_tgNsp'
    Write (messages(8),'(1x,a)')    '  amounts   N_tg1  N_tg2  N_tg3 .... N_tgNsp'
    Write (messages(9),'(1x,a)')    '&end_reference_composition'
    Write (messages(10),'(1x,a)')    'See the "use_code.md" file for details'
    
    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%reference_composition%atomic_species
    Call capital_to_lower_case(word) 
    If (Trim(word) /= 'atomic_species') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "atomic_species" is expected.'
      error=.True.
    End If 
    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "atomic_species"'
      error=.True.
    Else
      If (model_data%reference_composition%atomic_species<1) Then
        Write (messages(2),'(a)') 'The "atomic_species" directive MUST BE >= 1'
        error=.True.
      End If  
    End If
   
    If (error) Then
      Call info(messages, ilist) 
      Call error_stop(' ')
    End If

    Call model_data%init_reference_composition()
   
    i=1
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(iunit, '&reference_composition')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&reference_composition') 
              Read (iunit, Fmt=*, iostat=io) word,&
                            & (model_data%reference_composition%tag(j), j = 1, model_data%reference_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='amounts') Then
            If (.Not. header(2)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&reference_composition') 
              Read (iunit, Fmt=*, iostat=io) word,&
                               & (model_data%reference_composition%N0(j), j = 1, model_data%reference_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "amount" directive for each atomic tag'
                Call info(messages, ilist)
                Call error_stop(' ')
              End If  
              header(2)=.True.
            Else
              error_duplication=.True.
            End If  
          Else If (Trim(word)=='elements') Then
            If (.Not. header(3)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&reference_composition') 
              Read (iunit, Fmt=*, iostat=io) word, &
                         & (model_data%reference_composition%element(j), j = 1, model_data%reference_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the chemical element for each atomic tag (directive "element")'
                Call info(messages, ilist)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word), '". Please chek the amount of atomic species&
                                          & defined and the input species.'
            Call info(messages, ilist)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    ' '
          Call info(messages, ilist)
          Call error_stop(' ')
        End If
      End If
 
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within the block.&
                                     & If there is no duplication, then there is an inconsistency&
                                     & in the info provided for value of "atomic_species"'
        Call info(messages, ilist)
        Call error_stop(' ')
      End If
 
    End Do 
 
    endblock=.False.
    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&reference_composition') 
      Call capital_to_lower_case(word)
      If (word /= '&end_reference_composition') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'All info have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_reference_composition.' 
            Call info(messages, ilist)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

    ! Check if species tags contain asterix
    Do i=1, model_data%reference_composition%atomic_species
      ic= Index(Trim(model_data%reference_composition%tag(i)), '*') 
      If (ic > 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%reference_composition%tag(i)), &
                                 '" contains an asterisk. Defined species MUST NOT contain asterisks. Please correct.'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do
    
    ! Check if the number of atoms are correct
    Do i=1, model_data%reference_composition%atomic_species
      If (model_data%reference_composition%N0(i)< 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%reference_composition%tag(i)), '" CANNOT be associated with&
                                 & negative number of atoms within the input structure! Please correct'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do
 
    ! Calculate the number of total atoms set in the block
    model_data%reference_composition%numtot=0
    Do i=1, model_data%reference_composition%atomic_species
      model_data%reference_composition%numtot=model_data%reference_composition%numtot+model_data%reference_composition%N0(i)
    End Do
 
    ! Assing atomic numbers
    Do i=1, model_data%reference_composition%atomic_species
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(model_data%reference_composition%element(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
         Write (messages(2),'(1x,5a)') 'Wrong chemical element "', Trim(model_data%reference_composition%element(i)),&
                                      & '" defined for species tag "', Trim(model_data%reference_composition%tag(i)),&
                                      & '". Please check.' 
         Call info(messages, 2)
         Call error_stop(' ')
      End If
    End Do
 
  End Subroutine read_reference_composition
  
  Subroutine read_region(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the &region block. This block defines the portion of the
    ! system to be analysed.
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

    Integer(Kind=wi)   :: io, length, k
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &region block (SETTINGS file).'
 
    traj_data%region%number=0
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_region" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_region') Exit
      Call check_for_rubbish(iunit, '&region')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='delta_x') Then
        traj_data%region%number(1)=traj_data%region%number(1)+1
        k=traj_data%region%number(1)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(1,k)%type, &
                                    & traj_data%region%domain(1,1,k),     &
                                    & traj_data%region%domain(1,2,k),     &
                                    & traj_data%region%inout(1,k)
        Call set_read_status(word, io, traj_data%region%invoke(1,k)%fread, &
                            & traj_data%region%invoke(1,k)%fail, traj_data%region%invoke(1,k)%type)
         
      Else If (Trim(word)=='delta_y') Then
        traj_data%region%number(2)=traj_data%region%number(2)+1
        k=traj_data%region%number(2)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(2,k)%type, &
                                    & traj_data%region%domain(2,1,k),     &
                                    & traj_data%region%domain(2,2,k),     &
                                    & traj_data%region%inout(2,k)
        Call set_read_status(word, io, traj_data%region%invoke(2,k)%fread, &
                            & traj_data%region%invoke(2,k)%fail, traj_data%region%invoke(2,k)%type)

      Else If (Trim(word)=='delta_z') Then
        traj_data%region%number(3)=traj_data%region%number(3)+1
        k=traj_data%region%number(3)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(3,k)%type, &
                                    & traj_data%region%domain(3,1,k),     &
                                    & traj_data%region%domain(3,2,k),     &
                                    & traj_data%region%inout(3,k)
        Call set_read_status(word, io, traj_data%region%invoke(3,k)%fread, &
                            & traj_data%region%invoke(3,k)%fail, traj_data%region%invoke(3,k)%type)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_region"?'
        Call error_stop(message)
      End If

    End Do
    
    ! Assing to 1 if not read
    Do k = 1, 3
      If (traj_data%region%number(k)==0) Then
        traj_data%region%number(k)=1
      End If
    End Do
    
  End Subroutine read_region

  Subroutine read_segment_trajectory(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the time settings from the &segment_trajectory block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(traj_type),  Intent(InOut) :: traj_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &segment_trajectory block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_segment_trajectory" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_segment_trajectory') Exit
      Call check_for_rubbish(iunit, '&segment_trajectory')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='segment_time') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%seg_analysis%segment_time%tag, &
                                      & traj_data%seg_analysis%segment_time%value,& 
                                      & traj_data%seg_analysis%segment_time%units
        Call set_read_status(word, io, traj_data%seg_analysis%segment_time%fread,&
                                      & traj_data%seg_analysis%segment_time%fail)

      Else If (Trim(word)=='start_time') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%seg_analysis%start_time%tag, &
                                      & traj_data%seg_analysis%start_time%value,& 
                                      & traj_data%seg_analysis%start_time%units
        Call set_read_status(word, io, traj_data%seg_analysis%start_time%fread,&
                                      & traj_data%seg_analysis%start_time%fail)

      Else If (Trim(word)=='restart_every') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%seg_analysis%restart_every%tag, &
                                      & traj_data%seg_analysis%restart_every%value,& 
                                      & traj_data%seg_analysis%restart_every%units
        Call set_read_status(word, io, traj_data%seg_analysis%restart_every%fread,&
                                      & traj_data%seg_analysis%restart_every%fail)

      Else If (Trim(word)=='end_time') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%seg_analysis%end_time%tag, &
                                      & traj_data%seg_analysis%end_time%value,&
                                      & traj_data%seg_analysis%end_time%units
        Call set_read_status(word, io, traj_data%seg_analysis%end_time%fread,&
                                      & traj_data%seg_analysis%end_time%fail)

      Else If (word(1:length) == 'normalise_at_t0') Then
        Read (iunit, Fmt=*, iostat=io) word, traj_data%seg_analysis%normalise_at_t0%stat
       Call set_read_status(word, io, traj_data%seg_analysis%normalise_at_t0%fread, traj_data%seg_analysis%normalise_at_t0%fail)
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_segment_trajectory"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_segment_trajectory

  Subroutine read_general_analysis(iunit, coord_distr_data, nndist_distr_data, rdf_data, unchanged_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read which non-reactive-related quantitities must be computed 
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),        Intent(In   ) :: iunit
    Type(coord_distr_type),  Intent(InOut) :: coord_distr_data
    Type(nndist_distr_type), Intent(InOut) :: nndist_distr_data
    Type(rdf_type),          Intent(InOut) :: rdf_data
    Type(unchanged_type),    Intent(InOut) :: unchanged_data

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message
    Character(Len=256) :: word
    Character(Len=256) :: set_error
    Logical :: error, fread
    
    set_error = '***ERROR in the &general_analysis block (SETTINGS file).'
    error=.False.
    fread= .True.

    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&general_analysis')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&general_analysis')
      End If
    End Do

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_general_analysis" to close the block.&
                                  & Check if directives "tag" and "list_indexes" are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_general_analysis') Exit
      Call check_for_rubbish(iunit, '&general_analysis')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == '&selected_nn_distances') Then
        Read (iunit, Fmt=*, iostat=io) nndist_distr_data%invoke%type
        Call set_read_status(word, io, nndist_distr_data%invoke%fread, nndist_distr_data%invoke%fail)
        !Read information inside the block
        Call read_selected_nn_distances(iunit, nndist_distr_data)

      Else If (word(1:length) == '&rdf') Then
        Read (iunit, Fmt=*, iostat=io) rdf_data%invoke%type
        Call set_read_status(word, io, rdf_data%invoke%fread, rdf_data%invoke%fail)
        !Read information inside the block
        Call read_rdf(iunit, rdf_data)

      Else If (word(1:length) == '&coord_distrib') Then
        Read (iunit, Fmt=*, iostat=io) coord_distr_data%invoke%type
        Call set_read_status(word, io, coord_distr_data%invoke%fread, coord_distr_data%invoke%fail)
        !Read information inside the block
        Call read_coord_distrib(iunit, coord_distr_data)

      Else If (word(1:length) == '&track_unchanged_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) unchanged_data%invoke%type
        Call set_read_status(word, io, unchanged_data%invoke%fread, unchanged_data%invoke%fail)
        !Read information inside the block
        Call read_track_unchanged_chemistry(iunit, unchanged_data)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with&
                                & "&end_general_analysis"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_general_analysis
  
  Subroutine read_nonreactive_analysis(iunit, model_data, ocf_nonreactive, msd_data, nonreact_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read which non-reactive-related quantitities must be computed 
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),         Intent(In   ) :: iunit
    Type(model_type),         Intent(InOut) :: model_data
    Type(ocf_type),           Intent(InOut) :: ocf_nonreactive
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data    
    
    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message
    Character(Len=256) :: word
    Character(Len=256) :: set_error
    Logical :: error, fread
    
    set_error = '***ERROR in the &nonreactive_analysis block (SETTINGS file).'
    error=.False.
    fread= .True.

    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&nonreactive_analysis')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&nonreactive_analysis')
      End If
    End Do

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_nonreactive_analysis" to close the block.&
                                  & Check if directives "tag" and "list_indexes" are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_nonreactive_analysis') Exit
      Call check_for_rubbish(iunit, '&nonreactive_analysis')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == '&ocf_nonreactive') Then
        Read (iunit, Fmt=*, iostat=io) ocf_nonreactive%invoke%type
        Call set_read_status(word, io, ocf_nonreactive%invoke%fread, ocf_nonreactive%invoke%fail)
        !Read information inside the block
        Call read_ocf_settings(iunit, ocf_nonreactive, 'ocf_nonreactive')

      Else If (word(1:length) == '&msd') Then
        Read (iunit, Fmt=*, iostat=io) msd_data%invoke%type
        Call set_read_status(word, io, msd_data%invoke%fread, msd_data%invoke%fail)
        !Read information inside the block
        Call read_msd(iunit, msd_data)

      Else If (word(1:length) == '&intramol_statistics') Then
        Read (iunit, Fmt=*, iostat=io) nonreact_stat_data%intra_geom%invoke%type
        Call set_read_status(word, io, nonreact_stat_data%intra_geom%invoke%fread, &
                            & nonreact_stat_data%intra_geom%invoke%fail, &
                            & nonreact_stat_data%intra_geom%invoke%type)
        nonreact_stat_data%intra_geom%tag='intramol_statistics'
        Call read_geom_param_nonreactive_species(iunit, nonreact_stat_data%intra_geom)

      Else If (word(1:length) == '&intermol_statistics') Then
        Read (iunit, Fmt=*, iostat=io) nonreact_stat_data%inter_geom%invoke%type
        Call set_read_status(word, io, nonreact_stat_data%inter_geom%invoke%fread, &
                            & nonreact_stat_data%inter_geom%invoke%fail, &
                            & nonreact_stat_data%inter_geom%invoke%type)
        nonreact_stat_data%inter_geom%tag='intermol_statistics'
        Call read_geom_param_nonreactive_species(iunit, nonreact_stat_data%inter_geom)

      Else If (word(1:length) == 'compute_amount') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%nonreactive_species%compute_amount%stat
        Call set_read_status(word, io, model_data%nonreactive_species%compute_amount%fread,&
                           & model_data%nonreactive_species%compute_amount%fail)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with&
                                & "&end_nonreactive_analysis"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_nonreactive_analysis

  Subroutine read_reactive_analysis(iunit, ocf_reactive, restimes_data, tcf_data, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read which reactive-related quantitities must be computed 
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(ocf_type),      Intent(InOut) :: ocf_reactive
    Type(restimes_type), Intent(InOut) :: restimes_data
    Type(tcf_type),      Intent(InOut) :: tcf_data
    Type(spcf_type),     Intent(InOut) :: spcf_data


    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message
    Character(Len=256) :: word
    Character(Len=256) :: set_error
    Logical :: error, fread
    
    set_error = '***ERROR in the &reactive_analysis block (SETTINGS file).'
    error=.False.
    fread= .True.

    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&reactive_analysis')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&reactive_analysis')
      End If
    End Do

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_reactive_analysis" to close the block.&
                                  & Check if directives "tag" and "list_indexes" are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_reactive_analysis') Exit
      Call check_for_rubbish(iunit, '&reactive_analysis')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == '&tcf') Then
        Read (iunit, Fmt=*, iostat=io) tcf_data%invoke%type
        Call set_read_status(word, io, tcf_data%invoke%fread, tcf_data%invoke%fail)
        !Read information inside the block
        Call read_tcf(iunit, tcf_data)

      Else If (word(1:length) == '&spcf') Then
        Read (iunit, Fmt=*, iostat=io) spcf_data%invoke%type
        Call set_read_status(word, io, spcf_data%invoke%fread, spcf_data%invoke%fail)
        !Read information inside the block
        Call read_spcf(iunit, spcf_data)

      Else If (word(1:length) == '&residence_times') Then
        Read (iunit, Fmt=*, iostat=io) restimes_data%invoke%type
        Call set_read_status(word, io, restimes_data%invoke%fread, restimes_data%invoke%fail)
        !Read information inside the block
        Call read_residence_times(iunit, restimes_data)

      Else If (word(1:length) == '&ocf_reactive') Then
        Read (iunit, Fmt=*, iostat=io) ocf_reactive%invoke%type
        Call set_read_status(word, io, ocf_reactive%invoke%fread, ocf_reactive%invoke%fail)
        !Read information inside the block
        Call read_ocf_settings(iunit, ocf_reactive, 'ocf_reactive')

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with&
                                & "&end_reactive_analysis"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_reactive_analysis

  Subroutine check_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive,&
                                                  & msd_data, coord_distr_data, nonreact_stat_data, nndist_distr_data,&
                                                  & unchanged_data, rdf_data, restimes_data, tcf_data, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of trajectory-related directives
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(InOut) :: files(:)
    Type(model_type),         Intent(InOut) :: model_data
    Type(traj_type),          Intent(InOut) :: traj_data
    Type(ocf_type),           Intent(InOut) :: ocf_nonreactive
    Type(ocf_type),           Intent(InOut) :: ocf_reactive    
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(coord_distr_type),   Intent(InOut) :: coord_distr_data
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data
    Type(nndist_distr_type),  Intent(InOut) :: nndist_distr_data
    Type(unchanged_type),     Intent(InOut) :: unchanged_data
    Type(rdf_type),           Intent(InOut) :: rdf_data
    Type(restimes_type),      Intent(InOut) :: restimes_data
    Type(tcf_type),           Intent(InOut) :: tcf_data
    Type(spcf_type),          Intent(InOut) :: spcf_data

    
    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%print_retagged_trajectory%fread) Then
      If (traj_data%print_retagged_trajectory%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_retagged_trajectory" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
      If((.Not. model_data%reactive_chemistry%stat) .And. traj_data%print_retagged_trajectory%stat) Then 
        Write (messages(1),'(2(1x,a))') Trim(error_set), ' The user has set "print_retagged_trajectory" to .True. but&
                                      & "reactive_chemistry" is set to .False. Why do you want to retag the trajectory?&
                                      & Please change'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%print_retagged_trajectory%stat=.False.
    End If
    
    If (traj_data%print_track_chemistry%fread) Then
      If (traj_data%print_track_chemistry%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_track_chemistry" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
      If((.Not. model_data%reactive_chemistry%stat) .And. traj_data%print_track_chemistry%stat) Then 
        Write (messages(1),'(2(1x,a))') Trim(error_set), ' The user has set "print_track_chemistry" to .True. but&
                                      & "reactive_chemistry" is set to .False. Please change'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      If (model_data%reactive_chemistry%stat) Then
        traj_data%print_track_chemistry%stat=.True.
      Else
        traj_data%print_track_chemistry%stat=.False.
      End If
    End If    

    ! Check timestep
    Call check_time_directive(traj_data%timestep, 'timestep', error_set, .True.)
    
    ! Check ensemble
    If (traj_data%ensemble%fread) Then
      If (traj_data%ensemble%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "ensemble" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%ensemble%type)/='nve'  .And. &
            Trim(traj_data%ensemble%type)/='nvt'  .And. &
            Trim(traj_data%ensemble%type)/='npt') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "ensemble". Valid options: "NVE", "NVT" and "NPT"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "ensemble" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If

   ! Check infor for data analysis
    Call check_segment_trajectory(files, traj_data)
    
    ! Check info defined in &OCF_NONREACTIVE block 
    If (ocf_nonreactive%invoke%fread) Then
      Call check_ocf_nonreactive_species(files, ocf_nonreactive)
    End If 

    ! Check info defined in &region block 
     Call check_region(files, traj_data)
    
    ! Check info defined in &MSD block 
    If (msd_data%invoke%fread) Then
      Call check_msd(files, msd_data)
    End If 

    ! Check info defined in &RDF block 
    If (rdf_data%invoke%fread) Then
      Call check_rdf(files, model_data, rdf_data)
    End If 

    ! Check settigns for statistical analysis of angles and distances for nonreactive species 
    If (nonreact_stat_data%intra_geom%invoke%fread .Or. nonreact_stat_data%inter_geom%invoke%fread) Then
      Call check_nonreact_stat_settings(files, model_data%nonreactive_species, nonreact_stat_data)
    End If  
  
    ! Check &selected_nn_distances
    If (nndist_distr_data%invoke%fread) Then
      Call check_selected_nn_distances(files, model_data, nndist_distr_data)
    End If
  
    If (unchanged_data%invoke%fread) Then
      ! Check info defined in &track_unchanged_chemistry block 
      Call check_unchanged_chemistry(files, model_data, unchanged_data)
      !Check the labelling against info of the &track_unchanged_chemistry block
      Call check_initial_unchanged_labels(files, model_data, traj_data, unchanged_data)
    End If
  
    ! Check info defined in &coord_distrib block 
    If (coord_distr_data%invoke%fread) Then
      Call check_coord_distrib(files, model_data, coord_distr_data)
    End If 
    
    
    If (unchanged_data%invoke%fread) Then
    End If 

    ! Check info defined in &tcf block 
    If (tcf_data%invoke%fread) Then
      Call check_tcf(files, tcf_data)
    End If 

    ! Check info defined in &spcf block 
    If (spcf_data%invoke%fread) Then
      Call check_spcf(files, spcf_data)
    End If 

    ! Check info defined in &residence_times block 
    If (restimes_data%invoke%fread) Then
      Call check_residence_times(files, restimes_data)
    End If    
    
    ! Check info of the &ocf_reactive block
    If (ocf_reactive%invoke%fread) Then
      Call check_ocf_reactive_species(files, ocf_reactive)
    End If 

    If (Trim(model_data%geometry_format%type)=='xyz') Then
      If (Trim(traj_data%ensemble%type)/='nve' .And. Trim(traj_data%ensemble%type)/='nvt') Then
         Call info(' ', 1)
         Write (messages(1),'(1x,a)') Trim(error_set)//' To date, trajectories in "xyz" format cannot be&
                                    & processed in the "'//Trim(traj_data%ensemble%type)//'" ensemble&
                                    & (this is part of future implementation).'  
         Call info(messages,1)           
         Call error_stop(' ')
      End If
    Else If (Trim(model_data%geometry_format%type)=='vasp') Then      
      If(model_data%config%simulation_cell%fread) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' Trajectories in "vasp" format contain the definition&
                                    & of the simulation cell within the file.'
         Write (messages(2),'(4x,a)') 'Thus, definition of the "&simulation_cell" block is not needed and can cause&
                                    & problems. Please remove/comment "&simulation_cell".'  
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If
   
    If (ocf_nonreactive%invoke%fread) Then
      If (.Not. model_data%nonreactive_species%invoke%fread) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//'The computation of the orientational correlation&
                                     & function (OCF) requires the definition of the "&selected_nonreactive_species" block'
         Call info(messages,1)           
         Call error_stop(' ')
      Else
        If (model_data%nonreactive_species%atoms_per_species == 1) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The computation of the orientational correlation&
                                     & function (OCF) requires that the species defined in the&
                                     & "&selected_nonreactive_species" block is a molecule. Please review the settings'
         Call info(messages,1)           
         Call error_stop(' ')
        End If
      End If
    End If
    
     If (tcf_data%invoke%fread) Then
      If (.Not. model_data%reactive_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &tcf block but&
                                     & the &reactive_species block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of TCFs is&
                                     & only possible for systems with reactive species'
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If

    If (spcf_data%invoke%fread) Then
      If (.Not. model_data%reactive_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &spcf block but&
                                     & the &reactive_species block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of SPCFs is&
                                     & only possible for systems with reactive species'
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If

    If (restimes_data%invoke%fread) Then
      If (.Not. model_data%reactive_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &residence_times block but&
                                     & the &reactive_species block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of residence times is&
                                     & only possible for systems with reactive species'
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If

    If (ocf_reactive%invoke%fread) Then
      If (.Not. model_data%reactive_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &OCF_REACTIVE block but&
                                     & the &reactive_species block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of the orientational chemistry is&
                                     & only possible for systems with reactive species'
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If    
    
  End Subroutine check_settings_for_trajectory_analysis

  Subroutine check_region(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &egion block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer             :: k, m, j
    
    error_set = '***ERROR in the &region block of file '//Trim(files(FILE_SET)%filename)//' -'

    m=0
    If (traj_data%region%define%fread) Then
      Do k=1,3
        Do j = 1, traj_data%region%number(k)
          If (traj_data%region%invoke(k,j)%fread) Then
            If (traj_data%region%invoke(k,j)%fail) Then
              Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the '//&
                                             &Trim(traj_data%region%invoke(k,j)%type)//' directive.'
              Call info(messages, 1)
              Call error_stop(' ')
            Else
              If (traj_data%region%domain(k,1,j) > traj_data%region%domain(k,2,j)) Then
                If (traj_data%region%number(k)==1) Then
                  Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                            &'The lower value of the domain for defined "'&
                                            &//Trim(traj_data%region%invoke(k,j)%type)//&
                                            &'" is larger than the upper value!!! Please change.'
                Else
                  Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                            &'The lower value of the domain for one of the defined "'&
                                            &//Trim(traj_data%region%invoke(k,j)%type)//&
                                            &'" is larger than the upper value!!! Please change.'
                End If
                Call info(messages, 1)
                Call error_stop(' ')
              End If
              If (Abs(traj_data%region%domain(k,1,j) - traj_data%region%domain(k,2,j))<epsilon(1.0_wp)) Then
                If (traj_data%region%number(k)==1) Then
                  Write (messages(1),'(2(1x,a))') Trim(error_set),& 
                                          &'The lower and upper values of the domain for the defined "'&
                                          &//Trim(traj_data%region%invoke(k,j)%type)//&
                                          &'" are exaclty the same! Please change.'
                Else
                  Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                          &'The lower and upper values of the domain for one of the defined "'&
                                          &//Trim(traj_data%region%invoke(k,j)%type)//&
                                          &'" are exaclty the same! Please change.'
                End If
                Call info(messages, 1)
                Call error_stop(' ')
              End If
              Call capital_to_lower_case(traj_data%region%inout(k,j))    
              If (Trim(traj_data%region%inout(k,j)) /= 'inside' .And. Trim(traj_data%region%inout(k,j)) /= 'outside') Then
                 If (traj_data%region%number(k)==1) Then
                   Write (messages(1),'(2(1x,a))')  Trim(error_set),'The last argument of the defined directive "'&
                                                    &//Trim(traj_data%region%invoke(k,j)%type)//&
                                                    &'" must be either "inside" or "outside", referring&
                                                    & to the region defined by the limits. Please change.'
                 Else
                   Write (messages(1),'(2(1x,a))')  Trim(error_set),'The last argument for one of the defined directives "'&
                                                    &//Trim(traj_data%region%invoke(k,j)%type)//&
                                                    &'" must be either "inside" or "outside", referring&
                                                    & to the region defined by the limits. Please change.'
                 End If
                 Call info(messages, 1)
                 Call error_stop(' ')
              Else
                If (Trim(traj_data%region%inout(k,j)) == 'inside') Then
                  traj_data%region%inside(k,j)=.True.
                Else If (Trim(traj_data%region%inout(k,j)) == 'outside') Then
                  traj_data%region%inside(k,j)=.False.
                End If
              End If
            End If
          Else
            m=m+1
            traj_data%region%inside(k,j)=.True.
            traj_data%region%domain(k,1,j)=-Huge(1.0_wp)
            traj_data%region%domain(k,2,j)= Huge(1.0_wp)
          End If
        End Do
      End Do
    End If
    
    If (m==3) Then
       Write (messages(1),'(1x,a)') 'ERROR: the &region block contains no settings!' 
       Call info(messages, 1)
       Call error_stop(' ')
    End If
    
  End Subroutine check_region

  Subroutine check_segment_trajectory(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &segment_trajectory block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set, message

    error_set = '***ERROR in the &segment_trajectory block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%seg_analysis%invoke%fread) Then
      Call check_time_directive(traj_data%seg_analysis%segment_time, 'segment_time',  error_set, .False.)
      Call check_time_directive(traj_data%seg_analysis%end_time, 'end_time',  error_set, .False.)      
      Call check_time_directive(traj_data%seg_analysis%start_time, 'start_time', error_set, .False.)
      Call check_time_directive(traj_data%seg_analysis%restart_every, 'restart_every' ,error_set, .False.)

      If (traj_data%seg_analysis%normalise_at_t0%fread) Then
        If (traj_data%seg_analysis%normalise_at_t0%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                    & "normalise_at_t0" (choose either .True. or .False.)'
          Call info(message,1)
          Call error_stop(' ')
        End If
      Else
        traj_data%seg_analysis%normalise_at_t0%stat=.False.
      End If

    Else
      traj_data%seg_analysis%normalise_at_t0%stat=.False.
    End If
    
  End Subroutine check_segment_trajectory
 
End module settings
