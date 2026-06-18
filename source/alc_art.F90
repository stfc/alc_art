!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Welcome to ALC_ART: a ALC software to analyse MD trajectories
! both for reactive and non-reactive systems. ART stands for 
! "Analysis of Reactive Trajectories"
! while MDS refers to "MD Simulatios".
! The main purpose of this code is to offer the posibility to 
! compute simultaneously:
!
! * Radial Distribution Functions
! * Transfer Correlation Functions (only for reactive systems)
! * Mean Square Dislpacements
! * Orientational Correlation Functions
! * Special Pair Orientational and Transfer Correlation Functions
!
! Please refer to file "use_code.md" for a detailed explanation of the capabilities
! Example cases can be found in the examples folder
!
! This code is available under the BSD 3-Clause License.
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)  
!               
! Author:            Ivan Scivetti (i.scivetti)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_art

  Use analysis,            Only: trajectory_analysis, &
                                 print_settings_for_trajectory_analysis 
                           
  Use atomic_model,        Only: model_type, &
                                 atomistic_model, &
                                 read_model,&
                                 print_model_settings,&
                                 check_model_settings
                           
  Use fileset,             Only: file_type, &
                                 NUM_FILES, &
                                 print_header_out, &
                                 set_system_files, &
                                 wrapping_up
                                 
  Use coord_distr,         Only: coord_distr_type                          
                                 
  Use msd,                 Only: msd_type
                           
  Use nndist_distr,        Only: nndist_distr_type
                           
  Use nonreact_stat,       Only: nonreact_stat_type
                           
  Use numprec,             Only: wi,& 
                                 wp
                                 
  Use ocf,                 Only: ocf_type
                           
  Use rdf,                 Only: rdf_type
                          
  Use residence_times,     Only: restimes_type
                               
  Use settings,            Only: read_settings, &
                                 check_settings_for_trajectory_analysis
                          
  Use spcf,                Only: spcf_type
                               
  Use tcf,                 Only: tcf_type                           
                           
  Use trajectory,          Only: traj_type, &
                                 extract_trajectory,&
                                 trajectory_setup
                            
  Use unchanged_chemistry, Only: unchanged_type 
                          
  Use unit_output,         Only: info
                       

Implicit None

! Definition of variables
  Type(file_type)           :: files(NUM_FILES)
  Type(model_type)          :: model_data
  Type(traj_type)           :: traj_data
  Type(ocf_type)            :: ocf_nonreactive
  Type(ocf_type)            :: ocf_reactive
  Type(msd_type)            :: msd_data
  Type(coord_distr_type)    :: coord_distr_data
  Type(nonreact_stat_type)  :: nonreact_stat_data
  Type(nndist_distr_type)   :: nndist_distr_data
  Type(rdf_type)            :: rdf_data
  Type(restimes_type)       :: restimes_data
  Type(tcf_type)            :: tcf_data
  Type(spcf_type)           :: spcf_data
  Type(unchanged_type)      :: unchanged_data
  
  !Time related variables
  Integer(kind=wi)   :: start,finish,rate

  ! Array to print information
  Character(Len=256) :: message

  ! Start of the code 
  !!!!!!!!!!!!!!!!!!!
  ! Record initial time
  Call system_clock(count_rate=rate)
  Call system_clock(start)
  ! Initialise settings for input/output files
  Call set_system_files(files)
  ! Print header of OUTPUT
  Call print_header_out(files) 
  ! Read settings from SET
  Call read_settings(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data,&
                   & coord_distr_data, nonreact_stat_data, nndist_distr_data, unchanged_data,&
                   & rdf_data, restimes_data, tcf_data, spcf_data)
  ! Check the specification of directives from the SETTINGS file
  Call check_model_settings(files, model_data)
  Call check_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data,&
                                            & coord_distr_data, nonreact_stat_data, nndist_distr_data, unchanged_data,&
                                            & rdf_data, restimes_data, tcf_data, spcf_data)    
  ! Print model related settings
  Call print_model_settings(files, model_data)
  ! Prepare the trajectory
  Call trajectory_setup(files, model_data, traj_data)
  ! Print trajectory relatred information according to the definition in SETTINGS
  Call print_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data,&
                                            & coord_distr_data, nonreact_stat_data, nndist_distr_data, unchanged_data, &
                                            & rdf_data, restimes_data, tcf_data, spcf_data)
  ! Read and define trajectory
  Call extract_trajectory(files, model_data, traj_data)
  
  ! Analyse trajectory
  Call trajectory_analysis(files, model_data, traj_data, ocf_nonreactive, ocf_reactive, msd_data,&
                         & coord_distr_data, nonreact_stat_data, nndist_distr_data, unchanged_data, & 
                         & rdf_data, restimes_data, tcf_data, spcf_data)

  ! Record final time
  Call system_clock(finish)

  ! Print execution time
  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUT_EQCM file
  Call wrapping_up(files)

End Program alc_art
