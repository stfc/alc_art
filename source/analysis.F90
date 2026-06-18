!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the analysis of the generated trajectory
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module analysis

  Use atomic_model,        Only: model_type, &
                                 reference_tag_reactivity
                           
  Use coord_distr,         Only: coord_distr_type, &
                                 compute_coordinate_distribution
                                 
  Use fileset,             Only: file_type, &
                                 refresh_out
                                 
  Use nndist_distr,        Only: nndist_distr_type, &
                                 compute_nn_distance_distribution
                                 
  Use nonreact_stat,       Only: nonreact_stat_type, &
                                 geometry_statistics_nonreactive_species
                                
  Use numprec,             Only: wi,& 
                                 wp
                           
  Use msd,                 Only: msd_type, & 
                                 mean_squared_displacement
                           
  Use ocf,                 Only: ocf_type, &                            
                                 compute_ocf_reactive_species, &
                                 compute_ocf_nonreactive_species
                            
  Use rdf,                 Only: rdf_type, &
                                 radial_distribution_function
                                 
  Use residence_times,     Only: restimes_type, &
                                 residence_times_reactive_sites
  
  Use spcf,                Only: spcf_type, &
                                 special_pair_correlation_function
                                 
  Use tcf,                 Only: tcf_type, &
                                 transfer_correlation_function_sites
                                 
  Use trajectory,          Only: traj_type, &
                                 define_trajectory_segments, &
                                 compute_number_nonreactive_species, &
                                 find_active_bonds,&
                                 print_tracking_species
                                 
  Use unchanged_chemistry, Only: unchanged_type, &
                                 print_unchanged_chemistry
                                 
  Use unit_output,         Only: info, &
                                 error_stop 

  Implicit none 
  Private
  
  Public :: trajectory_analysis, print_settings_for_trajectory_analysis
  
Contains

  Subroutine print_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive,&
                                                  & msd_data, coord_distr_data, nonreact_stat_data, nndist_distr_data,&
                                                  & unchanged_data, rdf_data, restimes_data, tcf_data, spcf_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print a summary of the trajectory settings
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(In   ) :: files(:)
    Type(model_type),         Intent(In   ) :: model_data
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
    
    
    Character(Len=256) :: messages(4), word
    Integer(Kind=wi)   :: k, j
    Logical            :: tag_reacts  

    Call info(' ', 1) 
    Call info('Trajectory settings', 1) 
    Call info('===================', 1)     
    
    ! Check settings and define segments along the trajectory
    Call define_trajectory_segments(files, traj_data)
    
    Write (messages(1),'(1x,a)') '- by specification, the trajectory corresponds to the "'&
                                 &//Trim(traj_data%ensemble%type)//'" ensemble and was recorded in "'&
                                 &//Trim(model_data%geometry_format%type)//'" format'
    Call info(messages, 1)
    Write(word,'(f10.2)') traj_data%timestep%value
    Write (messages(1),'(1x,a,f5.2,a)') '- the time step between recorded configurations is '//Trim(Adjustl(word))//' fs'
    Call info(messages, 1)

    If (traj_data%seg_analysis%end_time%fread) Then
      If ((traj_data%seg_analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write(word,'(f8.3)') traj_data%seg_analysis%end_time%value/1000.0_wp
        Write (messages,'(1x,a)') '- the analysis will consider up to '//Trim(Adjustl(word))//' ps of the trajectory' 
        Call info(messages, 1)
      End If
    End If        
    
    If (traj_data%seg_analysis%start_time%fread) Then
      Write(word,'(f8.3)') traj_data%seg_analysis%start_time%value/1000.0_wp
      Write (messages(1),'(1x,a)') '- the initial '//Trim(Adjustl(word))//' ps of the trajectory&
                                   & will be discarded from the analysis' 
      Call info(messages, 1)
    End If
    
    If (tcf_data%invoke%fread       .Or. &
        spcf_data%invoke%fread      .Or. &
        restimes_data%invoke%fread  .Or. &
        ocf_reactive%invoke%fread) Then
        Call info(' ', 1)
        Write (messages(1),'(1x,a)') '=== Settings for analysis of REACTIVE species ==='
        Call info(messages, 1)
    End If    
    
    If (tcf_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&tcf" block will compute the&
                                  & Transfer Correlation Function (TCF) by considering&
                                  & the change of reactive sites& 
                                  & using the method: '//Trim(tcf_data%method%type)
      Call info(messages, 1)                            
    End If

    If (spcf_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&spcf" block will compute the&
                                  & Special Pair Correlation Function (SPCF) for the reactive species&
                                  & using the method: '//Trim(spcf_data%method%type)
      Call info(messages, 1)                            
    End If

    If (restimes_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&residence_times" block will compute the&
                                  & residence times for each changing species (file RES_TIMES):'
      If (restimes_data%rattling_wait%fread) Then
        Write(word,'(f8.3)') restimes_data%rattling_wait%value/1000.0_wp
        Write (messages(2),'(3x,a)') 'Rattling times lower than '//Trim(Adjustl(word))//' ps will&
                                    & be discarded'  
      Else
        Write (messages(2),'(3x,a)') 'Rattling effects will be included'
      End If
      Call info(messages, 2)
    End If

    If (ocf_reactive%invoke%fread) Then
      Write (messages(1),'(1x,a)') 'The definition of the "&OCF_REACTIVE" block will compute the Orientational&
                                  & Correlation Function (OCF) of the reactive species as follows:'
      Write (messages(2),'(1x,a)') '- the unit vector is defined with the method: '//&
                                  & Trim(ocf_reactive%u_definition%type)
      Write (messages(3),'(1x,a,i2)') '- the correlation terms are obtained using the Legendre polynomial of order ',&
                                  & ocf_reactive%legendre_order%value
      Call info(messages, 3)
      If (ocf_reactive%u_definition%type == 'unrattled_special_pair') Then
        Write (messages(1),'(1x,a)') '***WARNING*** The use of "unrattled_special_pair" option is still under revision.&
                                     & We recommend using the "special" pair option instead.'
        Call info(messages, 1)
      End If          
      If (ocf_reactive%legendre_order%value /= 2) Then
        Write (messages(1),'(1x,a)') '***WARNING*** The user is advised to reconsider setting a value of "legendre_order"&
                                    & different than 2 for the computation of OCF_REACTIVE'
        Call info(messages, 1)
      End If
    End If
    
    If (ocf_nonreactive%invoke%fread  .Or. &
        msd_data%invoke%fread  .Or. &
        nonreact_stat_data%intra_geom%invoke%fread .Or. &
        nonreact_stat_data%inter_geom%invoke%fread) Then
        Call info(' ', 1)
        Write (messages(1),'(1x,a)') '=== Settings for analysis of NON-REACTIVE species ==='
        Call info(messages, 1)
    End If    
    
    If (ocf_nonreactive%invoke%fread) Then
      Call info(' ', 1)
      If (model_data%nonreactive_species%atoms_per_species /= 1) Then
        If (model_data%nonreactive_species%atoms_per_species == 2) Then
          If (Trim(ocf_nonreactive%u_definition%type)/='bond_12') Then
            Write (messages(1),'(1x,a)')  '**WARNING: since the nonreactive species is diatomic, the&
                                          & method to compute the rotating unit vector (u_definition) &
                                          & has been reset to "bond_12". Methods "bond_13", "bond_12-13"&
                                          & "bond_123" and "plane" are meaningless for this case!'
            ocf_nonreactive%u_definition%type='bond_12'                             
            Call info(messages, 1)
          End If
        End If
      End If
      Write (messages(1),'(1x,a)') 'The definition of the "&OCF_NONREACTIVE" block will compute the Orientational&
                                  & Correlation Function (OCF) of the nonreactive species "'&
                                  &//Trim(model_data%nonreactive_species%name%type)//&
                                  & '" (defined in the &selected_nonreactive_species block) as follows:'
      Write (messages(2),'(1x,a)') '- the attached rotating unit vector is defined with the method: '//&
                                  & Trim(ocf_nonreactive%u_definition%type)
      Write (messages(3),'(1x,a,i2)') '- the correlation terms are obtained using the Legendre polynomial of order ',&
                                  & ocf_nonreactive%legendre_order%value
      Call info(messages, 3)
      If (model_data%nonreactive_species%atoms_per_species > 2) Then
        If (Trim(ocf_nonreactive%u_definition%type)/='bond_12-13' .And. & 
          Trim(ocf_nonreactive%u_definition%type)/='bond_123') Then 
          Write (messages(1),'(1x,a)')  '*** WARNING: the "'//Trim(ocf_nonreactive%u_definition%type)//'"& 
                                         & option for the rotating unit vector (u_definition)&
                                         & is not recommended for most studies.'
          Write (messages(2),'(1x,a)')  '           Unless the user is fully certain, either the&
                                         & "bond_12-13" or "bond_123" option should be used instead.'
          Call info(messages, 2)
        End If
      End If
      If (ocf_nonreactive%legendre_order%value /= 2) Then
        Write (messages(1),'(1x,a)') '***WARNING*** The user is advised to reconsider setting a value of "legendre_order"&
                                    & different than 2 for the computation of OCF_NONREACTIVE'
        Call info(messages, 1)
      End If
    End If

    If (msd_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&MSD" block will execute a Mean Square&
                                  & Displacement analysis of the species "'&
                                  &//Trim(model_data%nonreactive_species%name%type)//&
                                  & '" (defined in the &selected_nonreactive_species block) as follows:'
      Write (messages(2),'(1x,a)') '- the values will be computed for the coordinates(s): '//&
                                  & Trim(msd_data%select%type)
      Call info(messages, 2)

      If (msd_data%pbc_xyz%fread) Then
        Write (messages(1),'(1x,a)') '- the "pbc_xyz" directive specifies which coordinate uses (or not) periodic&
                                     & boundary conditions' 
        Call info(messages, 1)                             
      End If
    End If

    If (nonreact_stat_data%intra_geom%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&intramol_statistics" block will compute the probability&
                                  & distribution for the intramolecular:'
      Call info(messages, 1)                            
      If (nonreact_stat_data%intra_geom%dist%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- distances, using the settings of "&distance_parameters"'
        Call info(messages, 1)                            
      End If
      If (nonreact_stat_data%intra_geom%angle%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- angles, using the settings of "&angle_parameters"'
        Call info(messages, 1)                            
      End If
      Write (messages(1),'(1x,a)') 'corresponding to the species "'//Trim(model_data%nonreactive_species%name%type)//&
                                  & '" (defined in the &selected_nonreactive_species block).'
      Call info(messages, 1)                          
    End If

    
    If (nonreact_stat_data%inter_geom%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&intermol_statistics" block will compute the probability&
                                  & distribution for the intermolecular:'
      Call info(messages, 1)                            
      If (nonreact_stat_data%intra_geom%dist%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- distances, using the settings of "&distance_parameters"'
        Call info(messages, 1)                            
      End If
      If (nonreact_stat_data%intra_geom%angle%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- angles, using the settings of "&angle_parameters"'
        Call info(messages, 1)                            
      End If
      tag_reacts=.False.
      Write (messages(1),'(1x,a)') 'by considering the two closest "'//Trim(model_data%nonreactive_species%name%type)//& 
                                  &'" species to each "'//Trim(model_data%nonreactive_species%name%type)//'" species&
                                  & (see the &selected_nonreactive_species block).'
      If (model_data%reactive_species%invoke%fread) Then
        Call reference_tag_reactivity(model_data, tag_reacts)
      End If
    
      If ((.Not. nonreact_stat_data%inter_geom%only_ref_tags_as_nn%stat) .And. tag_reacts) Then
        Write (messages(1),'(1x,a)') 'by considering each "'//Trim(model_data%nonreactive_species%name%type)//'" species&
                                    & and the two closest species containing "'&
                                    &//Trim(model_data%nonreactive_species%reference_tag%type)//'" or "'&
                                    &//Trim(model_data%nonreactive_species%reference_tag%type)//'*" as reference tags'
      End If
      Call info(messages, 1)
    End If    
    
    If (rdf_data%invoke%fread          .Or. &
        coord_distr_data%invoke%fread  .Or. &
        nndist_distr_data%invoke%fread .Or. &
        unchanged_data%invoke%fread) Then
        Call info(' ', 1)
        Write (messages(1),'(1x,a)') '=== Settings for GENERAL analysis ==='
        Call info(messages, 1)
    End If    
    
    If (rdf_data%invoke%fread) Then
      Call info(' ', 1)
      Write(word,'(f10.2)') rdf_data%dr%value
      Write (messages(1),'(1x,a)') 'The definition of the "&RDF" block will compute the Radial&
                                  & Distribution Function (RDF) and the Coordination Numbers (CN) using:'
      Write (messages(2),'(1x,a)') '- the tags defined in "tags_species_a"  and "tags_species_b"'
      Write (messages(3),'(1x,a)') '- a discretization of '//Trim(Adjustl(word))//' Angstrom'
      Call info(messages, 3)
    End If

    If (coord_distr_data%invoke%fread) Then
      Call info(' ', 1)
      Write(word,'(f10.2)') coord_distr_data%delta%value
      Write (messages(1),'(1x,a)') 'The definition of the "&coord_distrib" block will compute the distribution&
                                  & of the '//Trim(coord_distr_data%coordinate%type)//'-values for all&
                                  & the "'//Trim(coord_distr_data%species)//'" species in the whole system'
      Write (messages(2),'(1x,a)') 'with a selected discretization of '//Trim(Adjustl(word))//' Angstrom for the coordinate.'
      Call info(messages, 2)
    End If
    
    If (nndist_distr_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&selected_nn_distances" block will compute probability distribution&
                               & of the shortest distance of the selected pair of species (this is not RDF)'
      Call info(messages, 1)                             
    End If
    
    If (unchanged_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&track_unchanged_chemistry" block will print the positions&
                                  & of the selected atomic indexes with unchanged chemistry along the trajectory.'
      Call info(messages, 1)
    End If
    
    If (ocf_nonreactive%invoke%fread      .Or. &
        ocf_reactive%invoke%fread .Or. &
        msd_data%invoke%fread      .Or. &
        tcf_data%invoke%fread      .Or. &
        spcf_data%invoke%fread)    Then
      If (traj_data%seg_analysis%segment_time%fread .Or. &
          (traj_data%seg_analysis%restart_every%fread  .And. (traj_data%seg_analysis%N_seg /=1))) Then
        Call info(' ', 1)  
        Write (messages(1),'(1x,a)') 'Instructions for analysis in time segments (&segment_trajectory)&
                                    & only applies to:' 
        Call info(messages, 1)
        If (tcf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* TCF'
          Call info(messages, 1)
        End If  
        If (spcf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* SPCF'
          Call info(messages, 1)
        End If  
        If (ocf_reactive%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* OCF_REACTIVE'
          Call info(messages, 1)
        End If
        If (ocf_nonreactive%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* OCF_NONREACTIVE'
          Call info(messages, 1)
        End If  
        If (msd_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* MSD'
          Call info(messages, 1)
        End If
        If (traj_data%seg_analysis%segment_time%fread) Then
          Write(word,'(f8.3)') traj_data%seg_analysis%segment_time%value/1000.0_wp
          Write (messages(1),'(2x,a)') '- using segments of '//Trim(Adjustl(word))//' ps' 
          Call info(messages, 1)
        End If
        If (traj_data%seg_analysis%restart_every%fread  .And. (traj_data%seg_analysis%N_seg /=1)) Then
          Write(word,'(f8.3)') traj_data%seg_analysis%restart_every%value/1000.0_wp
          Write (messages(1),'(2x,a)') '- the starting points of segments for analysis are separated&
                                         & by '//Trim(Adjustl(word))//' ps'
          Call info(messages, 1)
        End If
      End If
    End If

    If (traj_data%region%define%fread) Then
      Call info(' ', 1)
      If (ocf_nonreactive%invoke%fread .Or. &
          msd_data%invoke%fread .Or. &
          rdf_data%invoke%fread .Or. &
          ocf_reactive%invoke%fread .Or. &
          nndist_distr_data%invoke%fread .Or. &
          nonreact_stat_data%intra_geom%invoke%fread .Or. &
          nonreact_stat_data%inter_geom%invoke%fread) Then
        Write (messages(1),'(1x,a)') 'From the definition of the "&region" block, the computation of'
        Call info(messages, 1)
        If (tcf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* TCF'
          Call info(messages, 1)
        End If  
        If (spcf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* SPCF'
          Call info(messages, 1)
        End If  
        If (ocf_reactive%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* OCF_REACTIVE'
          Call info(messages, 1)
        End If
        If (ocf_nonreactive%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* OCF_NONREACTIVE'
          Call info(messages, 1)
        End If  
        If (msd_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* MSD'
          Call info(messages, 1)
        End If
        If (nonreact_stat_data%intra_geom%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* Intramolecular parameters (nonreactive species)'
          Call info(messages, 1)
        End If
        If (nonreact_stat_data%inter_geom%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* Intermolecular parameters (nonreactive species)'
          Call info(messages, 1)
        End If
        If (rdf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* RDF'
          Call info(messages, 1)
        End If  
        If (nndist_distr_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '* shortest distance distribution for the selected pair'
          Call info(messages, 1)
        End If
        Write (messages(1),'(1x,a)') 'will be only carried out:'
        Call info(messages, 1)
        Do k = 1, 3
          Do j = 1, traj_data%region%number(k)
            If (traj_data%region%invoke(k,j)%fread) Then
              Write (messages(1),'(3x,a,2f9.2)') '* '//Trim(traj_data%region%inout(k,j))//' the "'&
                                        //Trim(traj_data%region%invoke(k,j)%type)//'" region with&
                                        & lower and upper value: ', traj_data%region%domain(k,1,j), &
                                        & traj_data%region%domain(k,2,j)
              Call info(messages, 1)
            End If
          End Do 
        End Do
        If (rdf_data%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the RDF analysis, the definition of the &region block&
                                      & only applies to the species listed in "tags_species_a" (&rdf block)'
          Call info(messages, 1)
        End If  
        If (nndist_distr_data%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the analysis of the shortest distance distribution,&
                                      & the definition of the &region block only applies to the species& 
                                      & listed in "reference_species" (&selected_nn_distances block)'
          Call info(messages, 1)
        End If  
      End If 
    End If
    
    ! Refresh 
    Call refresh_out(files)
    
  End Subroutine print_settings_for_trajectory_analysis
  
  Subroutine trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data,&
                               & coord_distr_data, nonreact_stat_data, nndist_distr_data, unchanged_data,&
                               & rdf_data, restimes_data, tcf_data, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the trajectory depending on the options of the
    ! SETTINGS file
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(InOut) :: files(:)
    Type(model_type),         Intent(In   ) :: model_data
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
  
    Character(Len=256)  :: message

    traj_data%active_bonds_computed=.False.
    
    If(model_data%reactive_chemistry%stat) Then
      If(traj_data%print_track_chemistry%stat) Then 
        Call print_tracking_species(files, traj_data, model_data)
      End If
    ! Compute quantities for reactive species
      If (traj_data%reactive_analysis%fread) Then
        Call compute_reactive_quantities(files, model_data, traj_data, ocf_reactive, restimes_data, tcf_data, spcf_data)
      End If
    End If

    ! Compute quantities for non-reactive species
    If (traj_data%nonreactive_analysis%fread) Then
      Call compute_nonreactive_quantities(files, model_data, traj_data, ocf_nonreactive, msd_data, nonreact_stat_data)
    End If
    
    ! Compute quantities for reactive and/or non-reactive species, depending on the requested calculation
    If (traj_data%general_analysis%fread) Then
      Call compute_general_quantities(files, model_data, traj_data, coord_distr_data, nndist_distr_data,&
                                    & unchanged_data, rdf_data)
    End If
    
    If (Trim(traj_data%ensemble%type)/='nve') Then
      If (tcf_data%invoke%fread       .Or.&
          spcf_data%invoke%fread      .Or.&
          ocf_nonreactive%invoke%fread       .Or.&
          ocf_reactive%invoke%fread  .Or.&
          msd_data%invoke%fread) Then   
        Call info(' ', 1)
        Call info(' ****************************************************************************************', 1)
        Write (message,'(1x,a)') 'IMPORTANT: The user should bear in mind that the computed properties&
                                 & might be influenced'
        Call info(message, 1)                         
        If (Trim(traj_data%ensemble%type)=='nvt') Then
          Write (message,'(12x,a)') 'by the "thermostat" used to generate the trajectory'  
        Else If (Trim(traj_data%ensemble%type)=='npt') Then
          Write (message,'(12x,a)') 'by the "thermostat" and "barostat" used to generate the trajectory'  
        End If
        Call info(message, 1)
        Call info(' ****************************************************************************************', 1)
      End If
    End If
    
  End Subroutine trajectory_analysis
 
  Subroutine compute_general_quantities(files, model_data, traj_data, coord_distr_data, nndist_distr_data, unchanged_data, rdf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute general quantities 
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),         Intent(InOut) :: files(:)
    Type(model_type),        Intent(In   ) :: model_data
    Type(traj_type),         Intent(InOut) :: traj_data
    Type(coord_distr_type),  Intent(InOut) :: coord_distr_data
    Type(nndist_distr_type), Intent(InOut) :: nndist_distr_data
    Type(unchanged_type),    Intent(InOut) :: unchanged_data
    Type(rdf_type),          Intent(InOut) :: rdf_data
  
    Character(Len=256)  :: message

    If (rdf_data%invoke%fread            .Or. &
        coord_distr_data%invoke%fread  .Or. &
        nndist_distr_data%invoke%fread        .Or. &
        unchanged_data%invoke%fread) Then
        Call info(' ', 1)
        Write (message,'(1x,a)') '=== Generated information from general analysis ==='
        Call info(message, 1)
        Call info(' ', 1)
    End If    

    ! Compute coordinate distribution
    If (coord_distr_data%invoke%fread) Then
      Call compute_coordinate_distribution(files, model_data, traj_data, coord_distr_data)
    End If

    ! Compute RDF 
    If (rdf_data%invoke%fread) Then
      Call radial_distribution_function(files, model_data, traj_data, rdf_data)
    End If

    ! Print coordinates for selected unchanged species along the MD trajectory
    If (unchanged_data%invoke%fread) Then
      Call print_unchanged_chemistry(files, traj_data,unchanged_data)
    End If

    ! Compute the distribution between the shortest distance of a selected pair
    If (nndist_distr_data%invoke%fread) Then
      Call compute_nn_distance_distribution(files, traj_data, model_data, nndist_distr_data)
    End If
    
  End Subroutine compute_general_quantities
    
  Subroutine compute_nonreactive_quantities(files, model_data, traj_data, ocf_nonreactive, msd_data, nonreact_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute quantities related to non-reactive species
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),          Intent(InOut) :: files(:)
    Type(model_type),         Intent(In   ) :: model_data
    Type(traj_type),          Intent(InOut) :: traj_data
    Type(ocf_type),           Intent(InOut) :: ocf_nonreactive
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(nonreact_stat_type), Intent(InOut) :: nonreact_stat_data
  
    Character(Len=256)  :: message

    If (ocf_nonreactive%invoke%fread  .Or. &
        msd_data%invoke%fread  .Or. &
        nonreact_stat_data%intra_geom%invoke%fread .Or. &
        nonreact_stat_data%inter_geom%invoke%fread) Then
        Call info(' ', 1)
        Write (message,'(1x,a)') '=== Generated information for NON-REACTIVE species ==='
        Call info(message, 1)
        Call info(' ', 1)
    End If    
    
    ! Compute OCF for nonreactive species
    If (ocf_nonreactive%invoke%fread) Then
      Call compute_ocf_nonreactive_species(files, traj_data, ocf_nonreactive)
    End If

    ! Compute MSD for nonreactive species
    If (msd_data%invoke%fread) Then
      Call mean_squared_displacement(files, model_data, traj_data, msd_data)
    End If

    ! Compute the total average of nonreactive species within the system
    If (model_data%nonreactive_species%invoke%fread .And. model_data%nonreactive_species%compute_amount%stat) Then
      Call compute_number_nonreactive_species(traj_data, model_data)
    End If

    ! Compute intramolecular and/or intermolecular statistics of nonreactive species
    If (nonreact_stat_data%intra_geom%invoke%fread .Or. nonreact_stat_data%inter_geom%invoke%fread) Then
      Call geometry_statistics_nonreactive_species(files, traj_data, model_data%nonreactive_species, nonreact_stat_data)
    End If
    
  End Subroutine compute_nonreactive_quantities
    
  Subroutine compute_reactive_quantities(files, model_data, traj_data, ocf_reactive, restimes_data, tcf_data, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute quantities related to reactive species
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(In   ) :: model_data
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(ocf_type),      Intent(InOut) :: ocf_reactive
    Type(restimes_type), Intent(InOut) :: restimes_data
    Type(tcf_type),      Intent(InOut) :: tcf_data
    Type(spcf_type),     Intent(InOut) :: spcf_data
  
    Character(Len=256)  :: message

    If (tcf_data%invoke%fread  .Or. &
        spcf_data%invoke%fread .Or. &
        restimes_data%invoke%fread .Or. &
        ocf_reactive%invoke%fread) Then
        Call info(' ', 1)
        Write (message,'(1x,a)') '=== Generated information for REACTIVE species ==='
        Call info(message, 1)
        Call info(' ', 1)
    End If    
    
    ! Compute tcf
    If (tcf_data%invoke%fread) Then
      Call transfer_correlation_function_sites(files, model_data, traj_data, tcf_data)
    End If
    
    ! Compute residence_times
    If (restimes_data%invoke%fread) Then
      Call residence_times_reactive_sites(files, model_data, traj_data, restimes_data)
    End If
    
    ! Compute SPCF
    If (spcf_data%invoke%fread) Then
      If (.Not. traj_data%active_bonds_computed) Then
        Call find_active_bonds(traj_data, model_data)
        traj_data%active_bonds_computed=.True.
      End If
      Call special_pair_correlation_function(files, model_data, traj_data, spcf_data)
    End If
    
    ! Compute OCF for the reactive species
    If (ocf_reactive%invoke%fread) Then
       If (.Not. traj_data%active_bonds_computed) Then
         Call find_active_bonds(traj_data, model_data)
         traj_data%active_bonds_computed=.True.
       End If
       Call compute_ocf_reactive_species(files, model_data, traj_data, ocf_reactive)
    End If

  End Subroutine compute_reactive_quantities
   
End Module analysis
