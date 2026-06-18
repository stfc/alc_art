!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to analyse the trajectory. If the directive reactive_chemistry
! is set to .True., the algorithm searches and tracks changes of 
! chemical species based on the information of the &reactive_species
! block. 
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Sept 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module trajectory 
  Use atomic_model, Only : model_type, &
                           about_cell, &
                           min_intra, &
                           atomistic_model, & 
                           check_definition_bonds, &
                           check_PBC,&
                           check_cell_consistency,&
                           check_length_directive, &
                           check_orthorhombic_cell,&
                           compute_distance_pbc,&
                           read_model,&
                           obtain_maximum_number_species,&
                           identify_monitored_indexes

  Use constants,    Only : max_components, &
                           max_at_species, &
                           initial_tolerance

  Use fileset,      Only : file_type, &
                           FILE_SET, & 
                           FILE_TRAJECTORY,&
                           FILE_TRACK_CHEMISTRY,&
                           FILE_TAGGED_TRAJ, &
                           refresh_out

  Use input_types,  Only : in_integer, &
                           in_logic,   &
                           in_scalar,  &
                           in_param,   & 
                           in_string

  Use numprec,      Only : wi,&
                           wp 
  
  Use process_data, Only : capital_to_lower_case, &
                           detect_rubbish,        &
                           remove_symbols,        &
                           remove_front_tabs
                           
  Use unit_output,  Only : error_stop, &
                           info

  Implicit None
  Private

  Type :: atom_type
     Real(Kind=wp)    :: r(3)
     Integer(Kind=wi) :: indx
     Character(Len=8) :: tag
     Character(Len=2) :: element
    Integer(Kind=wi)  :: nn_indx(3)
  End Type 

  Type :: box_type
     Real(Kind=wp)    :: cell(3,3)
     Real(Kind=wp)    :: invcell(3,3)
     Real(Kind=wp)    :: volume
     Real(Kind=wp)    :: cell_length(3)
  End Type 

  Type :: seg_analysis_type
    Type(in_string)   :: invoke
    Type(in_logic)    :: normalise_at_t0
    Type(in_param)    :: segment_time
    Type(in_param)    :: end_time
    Type(in_param)    :: start_time
    Type(in_param)    :: restart_every
    Integer(Kind=wi)  :: N_seg
    Integer(Kind=wi)  :: Np_segment
    Integer(Kind=wi)  :: frame_ini
    Integer(Kind=wi)  :: frame_last
    Logical           :: normalised
    Integer(Kind=wi), Allocatable :: seg_indx(:,:)
    Real(Kind=wp),    Allocatable :: variable(:,:)
    Integer(Kind=wi), Allocatable :: max_points(:) 
  End Type
 
  ! Type to describe species
  Type :: species_type
    Integer(Kind=wi) :: list(max_at_species)
    Integer(Kind=wi) :: nn(2)
    Logical          :: alive
    Real(Kind=wp)    :: u(3,2)
    Real(Kind=wp)    :: u0(3,2)
  End Type
  
  !Type to describe the region where to constrain the analysis
  Type :: region_type
    Type(in_string)  :: define
    Logical          :: belong(3,max_components)
    Type(in_string)  :: invoke(3,max_components)
    Character(Len=8) :: inout(3,max_components)
    Logical          :: inside(3,max_components)
    Real(Kind=wp)    :: domain(3,2,max_components)
    Integer(Kind=wi) :: number(3)
  End Type
  
  !Type for tracking species
  Type :: track_type
    Type(atom_type),    Allocatable :: config(:,:)
  End Type

  ! Type for eqcm data and analysis
  Type, Public :: traj_type
    Private
    Type(species_type),      Public,  Allocatable :: species(:,:)
    Type(atom_type),         Public,  Allocatable :: config(:,:)
    Type(box_type),          Public,  Allocatable :: box(:)
    Type(region_type),       Public   :: region
    Type(seg_analysis_type), Public   :: seg_analysis
    Type(track_type),        Public   :: track_chem              
    Type(in_param),          Public   :: timestep
    Type(in_logic),          Public   :: print_retagged_trajectory
    Type(in_logic),          Public   :: print_track_chemistry
    Integer(Kind=wi),        Public   :: frames
    Integer(Kind=wi),        Public   :: Nmax_species
    Integer(Kind=wi),        Public   :: N_species
    Logical                           :: reload_trajectory
    Logical,                 Public   :: active_bonds_computed
    Type(in_string),         Public   :: ensemble
    Type(in_string),         Public   :: reactive_analysis
    Type(in_string),         Public   :: nonreactive_analysis
    Type(in_string),         Public   :: general_analysis
  Contains
    Private
      Procedure         :: alloc_trajectory  => allocate_trajectory_arrays
      Procedure         :: alloc_analysis    => allocate_analysis_arrays
      Final             :: cleanup
  End Type traj_type

  Public :: extract_trajectory, within_region, check_time_directive  
  Public :: trajectory_setup, compute_number_nonreactive_species
  Public :: define_trajectory_segments, find_active_bonds, average_segments 
  Public :: print_tracking_species

Contains

  Subroutine allocate_trajectory_arrays(T, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate trajectory arrays
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(traj_type),   Intent(InOut)  :: T
    Type(model_type),   Intent(In   )  :: model_data
    
    Integer(Kind=wi)    :: fail(2)
    Character(Len=256)  :: message
    Logical             :: error_alloc

    error_alloc=.False.
    fail=0
    
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for the trajectory&
                                & (subroutine allocate_trajectory_arrays). It is likely that the&
                                & trajectory and/or the system is exceedingly large.'

    Allocate(T%config(T%frames, model_data%config%num_atoms), Stat=fail(1))
    Allocate(T%box(T%frames),                                 Stat=fail(2))
    If (Any(fail > 0)) Then
       error_alloc=.True.
    End If

    If(model_data%reactive_chemistry%stat) Then 
      Allocate(T%track_chem%config(T%frames, max_components), Stat=fail(1))
      If (Any(fail > 0)) Then
        error_alloc=.True.
      End If
    End If

    If (model_data%nonreactive_species%invoke%fread) Then
      Allocate(T%species(T%frames, model_data%config%Nmax_species), Stat=fail(1))
      If (Any(fail > 0)) Then
        error_alloc=.True.
      End If
      T%Nmax_species=model_data%config%Nmax_species
    End If

    If (error_alloc) Then
      Call error_stop(message)
    End If
      
  End Subroutine allocate_trajectory_arrays

  Subroutine allocate_analysis_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for date analysis
    !
    ! author    - i.scivetti Mrch 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(traj_type),   Intent(InOut)  :: T
    
    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    fail=0
    Allocate(T%seg_analysis%seg_indx(2,T%seg_analysis%N_seg), Stat=fail(1))
    Allocate(T%seg_analysis%max_points(T%seg_analysis%N_seg), Stat=fail(2))
    Allocate(T%seg_analysis%variable(T%seg_analysis%Np_segment, T%seg_analysis%N_seg), Stat=fail(3))

    If (Any(fail > 0)) Then
       Write (message,'(1x,1a)') '***ERROR: Allocation problems for trajectory analysis &
                               & (subroutine allocate_analysis_arrays). Please review&
                               & settings of the &segment_trajectory block'
       Call info(message, 1)
       Call error_stop(' ')
    End If

  End Subroutine allocate_analysis_arrays
  
  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type) :: T

    If (Allocated(T%config)) Then
      Deallocate(T%config)
    End If 

    If (Allocated(T%track_chem%config)) Then
      Deallocate(T%track_chem%config)
    End If 
    
    If (Allocated(T%species)) Then
      Deallocate(T%species)
    End If     

    If (Allocated(T%seg_analysis%seg_indx)) Then
      Deallocate(T%seg_analysis%seg_indx)
    End If 

    If (Allocated(T%seg_analysis%variable)) Then
      Deallocate(T%seg_analysis%variable)
    End If 

    If (Allocated(T%seg_analysis%max_points)) Then
      Deallocate(T%seg_analysis%max_points)
    End If 
   
  End Subroutine cleanup
  
  Subroutine trajectory_setup(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to setup variables and check variables against the trajectory
    !
    ! refact    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    
    Logical            :: safe, fortho
    Character(Len=256) :: message
    Character(Len=256) :: input_file, set_error
    Integer(Kind=wi)   :: j
    
    input_file=(files(FILE_TRAJECTORY)%filename)
    set_error = '***ERROR -'

    Inquire(File=input_file, Exist=safe)
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(input_file), 'not found'
      Call error_stop(message)
    Else
      traj_data%reload_trajectory=.False. 
    End If

    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=input_file,Status='old')
    Call read_model(files, model_data, 1, traj_data%ensemble%type)
    Close(files(FILE_TRAJECTORY)%unit_no) 
    ! Scale simulation cell 
    model_data%config%cell=model_data%config%cell_scaling * model_data%config%cell
    If (Trim(model_data%config%position_units%type) == 'bohr') Then  
      Do j=1,3 
        model_data%config%atom(:)%r(j)=model_data%config%position_scaling* model_data%config%atom(:)%r(j) 
      End Do
    End If
    ! Compute cell related quantities for first checks 
    Call about_cell(model_data%config%cell,model_data%config%invcell,&
                  & model_data%config%cell_length, model_data%config%volume)
                    
    ! Check consistency between the system and the input cell
    Call check_cell_consistency(model_data)
    Call check_orthorhombic_cell(model_data%config%cell, fortho) 
      If (.Not.fortho) Then
        Call info(' ', 1)
        Write (message,'(1x,1a)') '***WARNING: the atomic model is not orthorhombic. This code has only been tested&
                                  & for orthorhombic cells. We do not guarantee a correct functioning.'
        Call info(message, 1)
     End If        
    
    If (model_data%nonreactive_species%invoke%fread) Then
      ! Compute the maximum amount of nonreactive species to be monitored and allocate them      
      Call obtain_maximum_number_species(model_data)
      Call model_data%init_species() 
      Call identify_monitored_indexes(model_data)
      model_data%config%species(:)%alive=.False.
    End If

    ! Check bonds againts the size of the simulation cell
    If(model_data%reactive_chemistry%stat) Then 
      Call check_definition_bonds(model_data, 1)
    End If
    
    ! Check if the defined region for analysis is within the simulation cell
    If (traj_data%region%define%fread) Then
      Call check_region_domain(model_data, traj_data, 1)
    End If 

    ! Search for the number of frames
    Call obtain_number_frames(files, model_data, traj_data)
    ! Allocate trajectory arrays
    Call traj_data%alloc_trajectory(model_data)
    
  End Subroutine trajectory_setup  
    
  Subroutine extract_trajectory(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to extract the trajectory. Atomic positions and elements must
    ! be provided in the TRAJECTORY file. If there are changes in the chemistry,
    ! the initial atomic tags asociated to each atoms are redefined according to
    ! the settings of &reactive_species (and &extra_reactive_bonds, if defined) 
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    
    Logical            :: loop_traj
    Character(Len=256) :: message, messages(4), nframes, md_length, net_md_length, nsegments
    Character(Len=256) :: input_file
    Integer(Kind=wi)   :: i, j
    
    If ((.Not. model_data%reactive_species%invoke%fread) .And. &
        (.Not. model_data%nonreactive_species%invoke%fread)) Then
        Call info(' ', 1)
        Write (messages(1),'(1x,a)') '****************************************************************************************'
        Write (messages(2),'(1x,a)') '*** WARNING: Neither "&reactive_species" nor "&selected_nonreactive_species" blocks are '
        Write (messages(3),'(1x,a)') '             defined. Are you sure about what you want to do? Please review the settings'
        Write (messages(4),'(1x,a)') '****************************************************************************************'
        Call info(messages, 4)
        If (.Not. traj_data%general_analysis%fread) Then
            Call info(' ', 1)
            Write (messages(1),'(1x,a)') '***************************************************************'
            Write (messages(2),'(1x,a)') '*** ERROR: No analysis requested. There is no point to continue'
            Write (messages(3),'(1x,a)') '***************************************************************'
            Call info(messages, 3)
            Call error_stop(' ')
        End If        
    Else
      If (.Not. traj_data%general_analysis%fread) Then
          Call info(' ', 1)
          Write (messages(1),'(1x,a)') '*********************************************'
          Write (messages(2),'(1x,a)') '*** WARNING: No general analysis is requested'
          Write (messages(3),'(1x,a)') '*********************************************'
          Call info(messages, 3)
      End If
    End If
    
    input_file=(files(FILE_TRAJECTORY)%filename)

    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=(input_file),Status='old')
    Write(md_length,'(f12.3)') (traj_data%frames-1)*traj_data%timestep%value/1000.0_wp
    Write(nframes,'(i8)') traj_data%frames
    Call info(' ', 1)
    Call info('Start of the analysis', 1)
    Call info('=====================', 1)
    Write (message,'(1x,a)') 'The code has identified a total of '//Trim(Adjustl(nframes))//' frames. From&
                                 & the setting of the "recorded_timestep" directive, the recorded MD trajectory is '&
                                 &//Trim(Adjustl(md_length))//' ps long.'
    Call info(message, 1)

    If (traj_data%seg_analysis%N_seg > 1) Then
      Write(nsegments,'(i8)') traj_data%seg_analysis%N_seg
      Write (message,'(1x,a)') 'From the information of the &segment_trajectory block, a total of '&
                               &//Trim(Adjustl(nsegments))//' time segments will be considered.'
      Call info(message, 1)
    End If
    
    If (traj_data%seg_analysis%end_time%fread) Then
      If ((traj_data%seg_analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write(net_md_length,'(f8.3)') traj_data%seg_analysis%end_time%value/1000.0_wp
        Write (message,'(1x,a)') 'Nevertheless, from the set value of the "'//Trim(traj_data%seg_analysis%end_time%tag)//&
                                    &'" directive (&segment_trajectory block), the analysis will consider up to '&
                                    &//Trim(Adjustl(net_md_length))//' ps of the MD trajectory.' 
        Call info(message, 1)
      End If
    End If    
    
    Call info(' ', 1)
    Call info(' Reading trajectory from the "'//Trim(input_file)//'" file...', 1)
    i=1
    loop_traj=.True.
    Do While (i <= traj_data%frames)
      Call read_model(files, model_data, i, traj_data%ensemble%type)
      If (Trim(model_data%config%position_units%type) == 'bohr') Then
        Do j=1,3 
          model_data%config%atom(:)%r(j)=model_data%config%position_scaling* model_data%config%atom(:)%r(j) 
        End Do
      End If
      If (Trim(traj_data%ensemble%type) == 'npt') Then
        model_data%config%cell=model_data%config%cell_scaling * model_data%config%cell
        Call about_cell(model_data%config%cell,model_data%config%invcell,&
                        model_data%config%cell_length, model_data%config%volume)
        Call check_definition_bonds(model_data, i)
      End If
      ! Identify the components of the model
      Call atomistic_model(model_data, i)
      ! Copy to trajectory arrays for later analysis
      Call copy_to_trajectory(traj_data, model_data, i)
      i=i+1
    End Do

    Close(files(FILE_TRAJECTORY)%unit_no) 
    Call info(' The trajectory has been defined successfully!', 1)
    Call info(' ', 1)
    Call refresh_out(files)

    If (traj_data%frames==1) Then
      Call info(' **********************************************************', 1)
      Call info(' ** WARNING: ONLY ONE FRAME WAS DETECTED IN THE TRAJECTORY!', 1)
      Call info(' **********************************************************', 1)
      Call info(' ', 1)
    End If

    If(model_data%reactive_chemistry%stat) Then
     Call residence_percentage(traj_data, model_data)
     Call refresh_out(files) 
    End If
    
    If(traj_data%print_retagged_trajectory%stat) Then 
      input_file=(files(FILE_TAGGED_TRAJ)%filename)
      Call print_retagged_trajectory(files, model_data, traj_data)
      Write (message,'(1x,a)') 'A copy of the trajectory with modified tags for the atomic species was printed&
                              & to the "'//Trim(input_file)//'" file'
      Call info(message, 1)
      Call refresh_out(files)
    End If 

    If((.Not. traj_data%print_track_chemistry%stat) .And. model_data%reactive_chemistry%stat) Then 
      input_file=(files(FILE_TRACK_CHEMISTRY)%filename)
      Write (message,'(1x,a)') 'The user has instructed not to print the "'//Trim(input_file)//'" file'
      Call info(message, 1)
      Call refresh_out(files)
    End If 
    
  End Subroutine extract_trajectory
  
  Subroutine within_region(traj_data, i, m, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the atom under consideration is within the defined region 
    ! for analysis
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer,           Intent(In   ) :: m
    Logical,           Intent(  Out) :: flag        

    Integer(Kind=wi)  :: k, j

    Logical :: fpass(3)
   
    Do k = 1, 3
      Do j = 1, traj_data%region%number(k)
        If (traj_data%region%inside(k,j)) Then
          If (traj_data%region%domain(k,1,j) <= traj_data%config(i,m)%r(k) .And. &
              traj_data%region%domain(k,2,j) >= traj_data%config(i,m)%r(k)) Then
            traj_data%region%belong(k,j) = .True.
          Else
            traj_data%region%belong(k,j) = .False.
          End If       
        Else
          If (traj_data%region%domain(k,1,j) >  traj_data%config(i,m)%r(k) .Or. &
              traj_data%region%domain(k,2,j) <  traj_data%config(i,m)%r(k)) Then
            traj_data%region%belong(k,j) = .True.
          Else
            traj_data%region%belong(k,j) = .False.
          End If       
        End If
        If (j==1) Then
          fpass(k)=traj_data%region%belong(k,j)
        Else
          fpass(k)=fpass(k) .Or. traj_data%region%belong(k,j)
        End If
      End Do
    End Do 

    flag=fpass(1) .And. fpass(2) .And. fpass(3)
    
  End Subroutine within_region
    

  Subroutine residence_percentage(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence percentage of reactive 
    ! species along the trajectory
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, j, k, l, m
    Integer(Kind=wi)   :: counts(model_data%reactive_species%search_envr%N0_incl)
    Real(Kind=wp)      :: amount
    Character(Len=8)   :: word
    Character(Len=256) :: messages(5)
    
    l=0
    counts=0
    Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
      l=l+1
      k=0
      Do j = 1, model_data%reactive_species%N0%value
        Do m = 1, model_data%reactive_species%search_envr%N0_incl
          word=Trim(model_data%reactive_species%search_envr%tg_incl(m))//'*'
          If (traj_data%track_chem%config(i,j)%tag==word) Then
            counts(m)=counts(m)+1
            k=k+1 
          End If
        End Do
      End Do
    End Do

    Write (messages(1),'(1x,a)') 'Population probabilities of "'//Trim(model_data%reactive_species%type%type)//&
                                &'" species along MD trajectory'
    Write (messages(2),'(1x,a)') '------------------------------------'
    Write (messages(3),'(1x,a)') 'Fraction (%)    Reference Atomic Tag'
    Write (messages(4),'(1x,a)') '------------------------------------'
    Call info(messages, 4)                            
    Do m = 1, model_data%reactive_species%search_envr%N0_incl
      word=Trim(model_data%reactive_species%search_envr%tg_incl(m))//'*' 
      amount= 100.0_wp * Real(counts(m),Kind=wp)/(l*model_data%reactive_species%N0%value)
      Write (messages(1),'(6x,f7.3,4x,a)') amount, Trim(word)
      Call info(messages, 1)                                                                    
    End Do                    
    Write (messages(1),'(1x,a)') '------------------------------------'
    Call info(messages, 1)
    If (model_data%reactive_species%search_envr%N0_incl==1) Then
      Write (messages(1),'(1x,a)') 'NOTE: 100% population is consistent with having defined a single species in "include_tags"'
      Call info(messages, 1)
    End If
    Call info(' ', 1)
    
  End Subroutine residence_percentage

  Subroutine find_active_bonds(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to identify the active bond for the changing sites along the
    ! trajectory
    !
    ! author    - i.scivetti April 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, m
    Integer(Kind=wi)   :: s1, s2, s3
    
    Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
      Do  m= 1, model_data%reactive_species%N0%value
        Call compute_closest_pairs(traj_data, model_data, i, m, s1, s2, s3)
        traj_data%track_chem%config(i,m)%nn_indx(1)=s1
        traj_data%track_chem%config(i,m)%nn_indx(2)=s2
        traj_data%track_chem%config(i,m)%nn_indx(3)=s3
      End Do
    End Do
  
  End Subroutine find_active_bonds
  
  Subroutine compute_closest_pairs(traj_data, model_data, frame, nchem, s1, s2, s3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the closest possible acceptor
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(In   ) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: frame
    Integer(Kind=wi),  Intent(In   ) :: nchem
    Integer(Kind=wi),  Intent(Out  ) :: s1
    Integer(Kind=wi),  Intent(Out  ) :: s2
    Integer(Kind=wi),  Intent(Out  ) :: s3
  
    Integer(Kind=wi)   :: i, j, k, l, mindx(3)
    Real(Kind=wp)      :: dist, min_dist(3)
    Logical            :: match_j, fexcl, flag(3)
    Character(Len=8)   :: tgexcl 

    i=traj_data%track_chem%config(frame,nchem)%indx
    min_dist(1)=Huge(1.0_wp) 
    min_dist(2)=Huge(1.0_wp)
    min_dist(3)=Huge(1.0_wp)
    
    If (model_data%reactive_species%search_envr%info_exclude%fread) Then
      tgexcl=traj_data%config(frame,i)%tag  
      Call remove_symbols(tgexcl, '*')
      fexcl=.False.
      l=1
      Do While (l <= model_data%reactive_species%search_envr%N0_excl .And. (.Not. fexcl))
        If (tgexcl==model_data%reactive_species%search_envr%tg_excl(l)) Then
           fexcl=.True.
        End If
        l=l+1
      End Do  
    Else 
      fexcl=.False. 
    End If

    mindx(1)=i
    mindx(2)=i
    mindx(3)=i
    
    j=1
    Do While (j <= model_data%config%num_atoms)
      If (i/=j) Then
        match_j=.False.
        k=1
        Do While (k <= model_data%reactive_species%search_envr%N0_incl .And. (.Not. match_j))
          If (traj_data%config(frame,j)%tag==model_data%reactive_species%search_envr%tg_incl(k)) Then
            match_j=.True.
            If (fexcl) Then
              l=1
              Do While (l <= model_data%reactive_species%search_envr%N0_excl .And. match_j)
                If (traj_data%config(frame,j)%tag==model_data%reactive_species%search_envr%tg_excl(l)) Then
                   match_j=.False.
                End If
                l=l+1
              End Do
            End If
          End If
          k=k+1
        End Do

        If(match_j) Then
          Call compute_distance_PBC(traj_data%config(frame,i)%r, traj_data%config(frame,j)%r,&
                                  & traj_data%box(frame)%cell, traj_data%box(frame)%invcell, dist)
          flag(1)= dist < min_dist(1)
          flag(2)= dist < min_dist(2)
          flag(3)= dist < min_dist(3)
          If (flag(1) .And. flag(2) .And. flag(3)) Then
            min_dist(3)=min_dist(2)
            min_dist(2)=min_dist(1)
            min_dist(1)=dist
            mindx(3)=mindx(2)
            mindx(2)=mindx(1)
            mindx(1)=j
          Else If ((.Not. flag(1)) .And. flag(2) .And. flag(3)) Then
            min_dist(3)=min_dist(2)
            min_dist(2)=dist
            mindx(3)=mindx(2)
            mindx(2)=j
          Else If ((.Not. flag(1)) .And. (.Not. flag(2)) .And. flag(3)) Then
            min_dist(3)=dist
            mindx(3)=j
          End If
        End If
      End If
      j=j+1
    End Do
    
    s1=mindx(1)
    s2=mindx(2)
    s3=mindx(3)
    
    If (s1==s2) Then
      call error_stop('ERROR')
    End If
    
    If (s1==s3) Then
      call error_stop('ERROR')
    End If
    
    If (s2==s3) Then
      call error_stop('ERROR')
    End If
  
  End Subroutine compute_closest_pairs
 
  Subroutine average_segments(files, traj_data, file_number, what, msd_coord)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to average physical quantities computed for each time segment
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: file_number
    Character(Len=*),  Intent(In   ) :: what
    Character(Len=*),  Optional, Intent(In   ) :: msd_coord
  
    Integer(Kind=wi) :: i, k, iunit, Nnet 
    Real(Kind=wp)    :: sum_i, average, average_net, std, norma
    Logical          :: flag, fcompute
    Character(Len=256) :: message
    Character(Len=256) :: quantity

    If (Trim(what) == 'OCF_REACTIVE') Then
      quantity='OCF'
    Else
      quantity=what
    End If
    
    fcompute=.True.
    traj_data%seg_analysis%normalised=.False.
    
    ! Check the value at t=0
    If (Trim(what) /= 'MSD') Then
      sum_i=0.0_wp
      Nnet=0
      Do k= 1, traj_data%seg_analysis%N_seg
        Nnet=Nnet+1
        sum_i=sum_i+traj_data%seg_analysis%variable(1,k)
      End Do
      
      If (Nnet > 0) Then
        norma=sum_i/Nnet
        If (Abs(norma-1.0_wp)>initial_tolerance)then
          If (Abs(norma)< initial_tolerance) Then
            Write (message,'(1x,a)') '*** WARNING: Problems with the computation of the average '//Trim(what)
            Call info(message, 1)
            Write (message,'(1x,a)') '             This is likely due to poor statistics'
            Call info(message, 1)
            If (traj_data%region%define%fread) Then
              Write (message,'(1x,a)') '             Please check the settings: the &region block might be too small.'
            Else
              Write (message,'(1x,a)') '             Please check the settings'
            End If
            Call info(message, 1)
            Call info('***', 1)
          Else
            If (traj_data%seg_analysis%normalise_at_t0%stat) Then
              Write (message,'(1x,a)') '*** INFO: The average '//Trim(what)//' has been normalised at t=0.'
              Call info(message, 1)
              traj_data%seg_analysis%variable=traj_data%seg_analysis%variable/norma
              traj_data%seg_analysis%normalised=.True.
            Else
              Write (message,'(1x,a)') '*** WARNING: The average '//Trim(what)//' is NOT normalised at t=0.'
              Call info(message, 1)
              If (traj_data%region%define%fread) Then
                Write (message,'(1x,a)') '             Please check the settings: the region defined in the &region&
                                         & block for analysis might be too small.'
              End If
              
              If (traj_data%seg_analysis%invoke%fread) Then
                Write (message,'(1x,a)') '    To normalise, set the "normalise_at_t0" directive to .True. in the&
                                        & &segment_trajectory block.'
              Else 
                Write (message,'(1x,a)') '    To normalise, use the &segment_trajectory block and set the "normalise_at_t0"&
                                        & directive to .True.'
              End If
              Call info(message, 1)                       
              Call info(' ***', 1) 
            End If
          End If
        End If
      Else
        Write (message,'(1x,a)') '**** PROBLEMS: The average '//Trim(what)//' could not be computed.'
        Call info(message, 1)
        Write (message,'(1x,a)') '             This is likely due to poor statistics'
        Call info(message, 1)
         If (traj_data%region%define%fread) Then
           Write (message,'(1x,a)') '             Please check the settings: the region defined in the &region&
                                   & block for analysis might be too small.'
         End If
        Call info(message, 1)
        Call info(' ***', 1)
      End If
    End If
    
    If (fcompute) Then             
      ! Print header
      Open(Newunit=files(file_number)%unit_no, File=files(file_number)%filename, Status='Replace')
      iunit=files(file_number)%unit_no
      If (Trim(what) == 'MSD') Then
        Write (iunit,'(a)') '#  Average MSD and STD (in Angstrom^2) for the coordinate(s) "'//&
                        &Trim(msd_coord)//'" of the selected nonreactive species'
      Else If (Trim(what) == 'OCF') Then
        Write (iunit,'(a)') '#  Average OCF and STD (dimensionless)&
                        & for the selected nonreactive species'
      Else If (Trim(what) == 'OCF_REACTIVE') Then
        Write (iunit,'(a)') '#  Average OCF and STD (dimensionless)&
                        & for the "reactive" species'
      Else If (Trim(what) == 'TCF') Then
        Write (iunit,'(a)') '#  Average TCF and STD (dimensionless)&
                        & for the "reactive" species'
      Else If (Trim(what) == 'SPCF') Then
        Write (iunit,'(a)') '#  Average SPCF and STD (dimensionless)&
                        & for the special pairs associated to the "reactive" species'                        
      End If
      
      Write (iunit,'(a)') '#  Time (ps)      '//Trim(quantity)//'          STD' 
      
      i=1
      flag=.True.    
      Do While ((i<=traj_data%seg_analysis%Np_segment) .And. flag)
        sum_i=0.0_wp
        Nnet=0
        Do k= 1, traj_data%seg_analysis%N_seg
          If(i<=traj_data%seg_analysis%max_points(k)) Then
            Nnet=Nnet+1
            sum_i=sum_i+traj_data%seg_analysis%variable(i,k)
          End If
        End Do
      
        ! Compute average
        If (Nnet > 0) Then
          average=sum_i/Nnet
          If (Nnet>1) Then
            sum_i=0.0_wp
            Do k= 1, traj_data%seg_analysis%N_seg
              If(i<=traj_data%seg_analysis%max_points(k)) Then
                sum_i=sum_i+(traj_data%seg_analysis%variable(i,k)-average)**2
              End If
            End Do
            std=sqrt(sum_i/(Nnet-1))
          Else
            std=0.0_wp
          End If
        Else
          flag=.False.
        End If
        
        If (flag) Then
          If (Trim(what)=='MSD') Then
            average_net=average
          Else
            If(average > 1.0_wp) Then
              average_net=1.0_wp
            Else
              average_net=average
            End If
          End If  
          Write(iunit,'(3(f10.3, 3x))') (i-1)*traj_data%timestep%value/1000.0_wp, average_net, std
        End If
        i=i+1
        
      End Do
      Write (message,'(1x,a)') 'The average '//Trim(what)//' was printed to the "'//&
                               &Trim(files(file_number)%filename)//'" file.'
      Call info(message, 1)
      Close(iunit)
    End If 
    
    
  End Subroutine average_segments
  
  
  Subroutine obtain_number_frames(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to obtain the number of frames recorded in the TRAJECTORY file 
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data

    Logical            :: loop_traj
    Character(Len=256) :: check
    Integer(Kind=wi)   :: i, stat

    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=files(FILE_TRAJECTORY)%filename,Status='old')

    i=1
    loop_traj=.True.
    Do While (loop_traj)
      Call read_model(files, model_data, i, traj_data%ensemble%type)
      Read(files(FILE_TRAJECTORY)%unit_no, Fmt= *, iostat=stat) check
      If (is_iostat_end(stat)) Then
        loop_traj=.False.
      Else 
        backspace files(FILE_TRAJECTORY)%unit_no
      End If
      i=i+1
    End Do
    traj_data%frames=i-1
    
    Close(files(FILE_TRAJECTORY)%unit_no) 
    
  End Subroutine obtain_number_frames

  Subroutine define_trajectory_segments(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &segment_trajectory block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set, warn_set, message
    Integer(Kind=wi)    :: i, k, j, l, kref, kini, frames
    Logical             :: flag, one_seg
    Real(Kind=wp)       :: time, tend, teff, tini, tref
    
    error_set = '***ERROR in the &segment_trajectory block of file '//Trim(files(FILE_SET)%filename)//' -'
    warn_set = '***WARNING from the &segment_trajectory block of file '//Trim(files(FILE_SET)%filename)//' -'

   
    If (traj_data%seg_analysis%end_time%fread) Then
      If (traj_data%seg_analysis%end_time%value > (traj_data%frames-1)*traj_data%timestep%value) Then
        Write (message,'(2(1x,a))') Trim(warn_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%end_time%tag)//&
                                 &'" is larger than the total time for the trajectory. The analysis will be performed&
                                 & up to largest recorded time.'
        Call info(message, 1)          
        tend=(traj_data%frames-1)*traj_data%timestep%value
        frames=traj_data%frames
      Else

        If (traj_data%timestep%value>=traj_data%seg_analysis%end_time%value) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%end_time%tag)//&
                               &'" must be larger that the timestep for the trajectory. Please check the directives.'
          Call info(message, 1) 
          Call error_stop(' ')
        End If      
        
        If (traj_data%seg_analysis%end_time%value > (traj_data%frames-2)*traj_data%timestep%value) Then
          Write (message,'(2(1x,a))') Trim(warn_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%end_time%tag)//&
                                   &'" is in between the last two recorded times. The analysis will be performed&
                                   & up to largest recorded time.'
          Call info(message, 1)          
          tend=(traj_data%frames-1)*traj_data%timestep%value
          frames=traj_data%frames
        Else
          tend=traj_data%seg_analysis%end_time%value
          i=1
          flag=.True.
          Do While (i <= traj_data%frames .And. flag)
            time=(i-1)*traj_data%timestep%value
            If (time >= tend) Then
              frames=i
              flag=.False.
            End If
            i=i+1
          End do           
        End If
      End If
    Else
      tend=(traj_data%frames-1)*traj_data%timestep%value
      frames=traj_data%frames
    End If
    
    ! Set the net number of frames
    traj_data%seg_analysis%frame_last=frames
    
    If (.Not. traj_data%seg_analysis%start_time%fread) Then
       traj_data%seg_analysis%start_time%value=-traj_data%timestep%value
       traj_data%seg_analysis%frame_ini = 1
       tini=0.0_wp
    Else
      If (tend <= traj_data%seg_analysis%start_time%value) Then
        Call info(' ', 1)
        If (.Not. traj_data%seg_analysis%end_time%fread) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%start_time%tag)//&
                                 &'" is larger than (or equal) the total time for the trajectory. Please check&
                                 & the settings and the value for the "recorded_timestep" directive.'
        Else
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%start_time%tag)//&
                                 &'" is larger than (or equal) the value set for "'//Trim(traj_data%seg_analysis%end_time%tag)//&
                                 &'". Please check settings'        
        End If
        Call info(message, 1)
        Call error_stop(' ')
      Else  
        i=1
        flag=.True.
        Do While (i <= traj_data%frames .And. flag)
          time=(i-1)*traj_data%timestep%value
          If (time >= traj_data%seg_analysis%start_time%value) Then
            traj_data%seg_analysis%frame_ini = i
            tini=(i-1)*traj_data%timestep%value
            flag=.False.
          End If
          i=i+1
        End do 
      End If   
    End If
  
  
    teff=tend-tini 

    ! Compare timestep with other time settings of &segment_trajectory
    If (traj_data%seg_analysis%segment_time%fread) Then
      If (traj_data%timestep%value>=traj_data%seg_analysis%segment_time%value) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%segment_time%tag)//&
                               &'" must be larger that the timestep for the trajectory.&
                               & Please check the "recorded_timestep" directive.'
        Call info(message, 1) 
        Call error_stop(' ')
      End If
      If (teff<traj_data%seg_analysis%segment_time%value) Then
           Write (message,'(2(1x,a))') Trim(warn_set), 'The input value for the "'&
                                   &//Trim(traj_data%seg_analysis%segment_time%tag)//&
                                   &'" directive was too large and has been redefined to comply with the rest&
                                   & of the settings and the length of the trajectory.'
          Call info(message, 1)
          traj_data%seg_analysis%segment_time%value=teff
      End If
    Else
      traj_data%seg_analysis%segment_time%value=teff
    End If
 

    If (traj_data%seg_analysis%restart_every%fread) Then
      If (traj_data%timestep%value > traj_data%seg_analysis%restart_every%value) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%seg_analysis%restart_every%tag)//&
                                 &'" must be larger that the timestep for the trajectory. Please check values (and units).'
        Call info(message, 1) 
        Call error_stop(' ')
      End If
      If (teff<traj_data%seg_analysis%restart_every%value) Then
           Write (message,'(2(1x,a))') Trim(warn_set), 'The input value for the "'&
                                   &//Trim(traj_data%seg_analysis%restart_every%tag)//&
                                   &'" directive was too large and has been redefined to comply with the rest of the&
                                   & directive and the length of the trajectory.'
          Call info(message, 1)
          traj_data%seg_analysis%segment_time%value=teff+traj_data%timestep%value
      End If      
      
    Else
      traj_data%seg_analysis%restart_every%value=teff+traj_data%timestep%value      
    End If
    
    ! Calculate the number of segments
    i=0; j=0; l=0
    one_seg=.True.
    k=traj_data%seg_analysis%frame_ini
    tref=tini; kini=k; flag=.True.
    Do While (k <= frames)
      time=(k-1)*traj_data%timestep%value
      If (time>=(tref+traj_data%seg_analysis%segment_time%value)) Then
        i=i+1  
        l=k-kini+1 
        j=0
        If (time>=(tref+traj_data%seg_analysis%restart_every%value)) Then
          If (flag) Then
            kref=k
          End If
        Else
           If (traj_data%seg_analysis%restart_every%fread) Then
             kref=Nint((tref+traj_data%seg_analysis%restart_every%value)/traj_data%timestep%value)+1
           Else
             kref=k
           End If
        End If
        tref=(kref-1)*traj_data%timestep%value
        k=kref
        kini=k
        one_seg=.False.
        flag=.True.
      Else
        If (time>=(tref+traj_data%seg_analysis%restart_every%value) .And. flag) Then
          kref=k
          flag=.False.
        End If
        j=j+1
      End If
      k=k+1
    End Do

    If(one_seg) Then
      traj_data%seg_analysis%N_seg=1
      traj_data%seg_analysis%Np_segment=j
    Else
      traj_data%seg_analysis%N_seg=i
      traj_data%seg_analysis%Np_segment=l
    End If
    
    ! Allocate arrays
    Call traj_data%alloc_analysis()
    
    ! Calculate the number of segments
    If (traj_data%seg_analysis%N_seg /= 1) Then
      i=0; j=0; l=0
      k=traj_data%seg_analysis%frame_ini
      tref=tini; kini=k; flag=.True.
      Do While (k <= frames)
        time=(k-1)*traj_data%timestep%value
        If (time>=(tref+traj_data%seg_analysis%segment_time%value)) Then
          i=i+1
          l=k-kini+1
          j=0
          traj_data%seg_analysis%seg_indx(1,i)=kini
          traj_data%seg_analysis%seg_indx(2,i)=k
          If (time>=(tref+traj_data%seg_analysis%restart_every%value)) Then
            If (flag) Then
              kref=k
            End If
          Else
            If (traj_data%seg_analysis%restart_every%fread) Then
              kref=Nint((tref+traj_data%seg_analysis%restart_every%value)/traj_data%timestep%value)+1
            Else
              kref=k
            End If
          End If
          tref=(kref-1)*traj_data%timestep%value
          k=kref
          kini=k
          flag=.True.
        Else
          If (time>=(tref+traj_data%seg_analysis%restart_every%value) .And. flag) Then
            kref=k
            flag=.False.
          End If
          j=j+1
        End If
        k=k+1
      End Do
    Else
      traj_data%seg_analysis%seg_indx(1,1)=traj_data%seg_analysis%frame_ini
      traj_data%seg_analysis%seg_indx(2,1)=traj_data%seg_analysis%Np_segment+traj_data%seg_analysis%frame_ini-1
    End If
    
  End Subroutine define_trajectory_segments   

  Subroutine compute_number_nonreactive_species(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the average number (and STD) of the nonreactive species 
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)  :: i, j
    Integer(Kind=wi)  :: num_at_a, net_frames
    Integer(Kind=wi)  :: accum_a
    
    Character(Len=256) :: messages(3)
    Logical            :: flag

    Real(Kind=wp) :: average, std, sum_i
    
    ! counting
    Real(Kind=wp), Allocatable  :: nat(:)
   
    ! In case &region is defined
    flag=.True.
    net_frames=0
    accum_a=0
    
    Allocate(nat(traj_data%seg_analysis%frame_last-traj_data%seg_analysis%frame_ini+1))
    
    ! Compute the histogram for atoms of type a and b
    Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
      ! Define the number and list of indexes for type of species "a"
      num_at_a=0
      net_frames=net_frames+1 
      Do j = 1, model_data%config%num_atoms
        If (model_data%nonreactive_species%reference_tag%type==traj_data%config(i,j)%tag) Then
          If (traj_data%region%define%fread) Then
             Call within_region(traj_data, i, j, flag)
          End If
          If (flag) Then
            num_at_a=num_at_a+1
          End If
        End If
      End Do
      ! Accummulators
      nat(net_frames)=num_at_a 
      accum_a=accum_a+num_at_a
    End Do
      
    average= Real(accum_a,Kind=wp)/net_frames
    
    ! Compute average
    If (net_frames > 1) Then
      sum_i=0.0_wp
      j=0 
      Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
        j=j+1 
        sum_i=sum_i+(Real(nat(j),Kind=wp)-average)**2
      End Do
      std=sqrt(sum_i/(net_frames-1))
    Else
      std=0.0_wp
    End If

    If (traj_data%region%define%fread) Then
       Write (messages(1),'(1x,a)') 'Amount of nonreactive species "'//Trim(model_data%nonreactive_species%name%type)//&
                                   &'" within the selected region as specified in the &region block'
    Else
       Write (messages(1),'(1x,a)') 'Amount of nonreactive species "'//Trim(model_data%nonreactive_species%name%type)//&
                                   &'" within the simulation cell'
    End If
    Write (messages(2),'(1x,f8.2,5x,a,f8.2)')  average, '+/-', STD
    Call info(messages, 2)
    Call info(' ', 1)
    
  End Subroutine compute_number_nonreactive_species

  Subroutine check_region_domain(model_data, traj_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks the definition of the &region block against the size of the 
    ! simulation cell
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),  Intent(InOut) :: traj_data
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: frame
  
    Character(Len=256) :: messages(2)  
    Integer(Kind=wi)   :: k, j
    Real(Kind=wp)      :: min_cell, max_cell, vector(3)
    Logical            :: flag1, flag2

    Write (messages(1),'(1x,a,i6)') '***ERROR: inconsistency between the size of the simulation cell and&
                                & the sepecifications of the &region block for frame: ', frame     

    Do k = 1, 3
      Do j = 1, traj_data%region%number(k)
        If (traj_data%region%invoke(k,j)%fread) Then
          vector(:)=model_data%config%cell(:,k)
          min_cell=Minval(vector)
          max_cell=Maxval(vector)
          If (traj_data%region%inside(k,j)) Then
            flag1 = (traj_data%region%domain(k,1,j) <= min_cell) .And.&
                    (traj_data%region%domain(k,2,j) <= min_cell)
            flag2 = (traj_data%region%domain(k,1,j) >= max_cell) .And.&
                    (traj_data%region%domain(k,2,j) >= max_cell)
            If (flag1 .Or. flag2) Then
               Write (messages(2),'(1x,a)') 'There are NO atoms inside the domain range defined for "'//&
                                          &Trim(traj_data%region%invoke(k,j)%type)//'". Please change' 
               Call info(messages,2)
               Call error_stop(' ')
            End If
          Else  
            flag1 = (traj_data%region%domain(k,1,j) <= min_cell) .And.&
                    (traj_data%region%domain(k,2,j) >= max_cell)
            If (flag1) Then
               Write (messages(2),'(1x,a)') 'There are NO atoms outside the domain range defined for "'//&
                                          &Trim(traj_data%region%invoke(k,j)%type)//'". Please change' 
               Call info(messages,2)
               Call error_stop(' ')
            End If
          End If
        End If
      End Do
    End Do
  
  End Subroutine check_region_domain
 
  Subroutine check_time_directive(T, tag, error_set, kill)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check time related directivesd
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: tag 
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill

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
        If (Trim(T%units) /= 'fs' .And. &
           Trim(T%units) /= 'ps') Then
           Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                    & 'Units for directive "'//Trim(T%tag)//&
                                    &'" must be "fs" or "ps". Have you included the units?'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        ! Transform to fs
        If (Trim(T%units) == 'ps') Then
           T%value=1000_wp* T%value
        End If
      End If
    Else 
      If (kill)then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(tag)//'" directive'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_time_directive  
  
  Subroutine print_retagged_trajectory(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the trajectory
    !
    ! author    - i.scivetti July 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
  
    Integer(Kind=wi)   :: iunit, i, l, k 
  
    ! Print tracked species
      If(model_data%reactive_chemistry%stat) Then
        Open(Newunit=files(FILE_TAGGED_TRAJ)%unit_no, File=files(FILE_TAGGED_TRAJ)%filename, Status='Replace')
        iunit=files(FILE_TAGGED_TRAJ)%unit_no
        Write(iunit,'(2i8,1x,a)') model_data%config%num_atoms, traj_data%frames, ' # number of total atoms and trajectory frames'
        Do l = 1, traj_data%frames
          Write(iunit,'(a,1x,i8)') 'Frame=', l 
          Do i = 1, model_data%config%num_atoms 
            Write(iunit,'(a, 4x, 3(f11.3), 4x, a)') Trim(traj_data%config(l,i)%element),&
                                                 & (traj_data%config(l,i)%r(k), k=1, 3),&
                                                 &  Trim(traj_data%config(l,i)%tag) 
          End Do               
        End Do
        Close(iunit)
      End If
  
  End Subroutine print_retagged_trajectory

  Subroutine print_tracking_species(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the positions of those species that change their 
    ! chemsitry along the trajectory
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
  
    Integer(Kind=wi)   :: iunit, i, l 
    Character(Len=256) :: num_species
    Character(Len=256) :: message
  
    ! Print tracked species
    Open(Newunit=files(FILE_TRACK_CHEMISTRY)%unit_no, File=files(FILE_TRACK_CHEMISTRY)%filename, Status='Replace')
    iunit=files(FILE_TRACK_CHEMISTRY)%unit_no
    If (traj_data%seg_analysis%frame_ini==1) Then
      Write(iunit,'(a)') '# Tracking the change of chemical species over the whole trajectory'    
    Else
      Write(iunit,'(a,1x,f10.4,1x,a)') '# Tracking the change of chemical species ignoring the first',& 
                                   &  traj_data%seg_analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.'
    End If

    If (model_data%reactive_species%N0%value==1) Then
      Write (iunit,'(a,9x,a)') '# Time (ps)', 'XYZ_Species_1' 
    Else
      Write(num_species,*) model_data%reactive_species%N0%value
      Write (iunit,'(a,9x,2a)') '# Time (ps)', 'XYZ_Species_1 .... XYZ_Species_', Trim(Adjustl(num_species)) 
    End If
    
    Do i = traj_data%seg_analysis%frame_ini, traj_data%seg_analysis%frame_last
       Write(iunit,'(f10.4, 4x, *(f11.3))') (i-traj_data%seg_analysis%frame_ini)*traj_data%timestep%value/1000_wp,&
                                       & (traj_data%track_chem%config(i,l)%r(:), l=1, model_data%reactive_species%N0%value)
    End Do
    Write (message,'(1x,a)') 'The tracking of the reactive species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_TRACK_CHEMISTRY)%filename)//'" file'
    Call info(message, 1)
    Call refresh_out(files)
    Close(iunit)
    Call info(' ', 1)
  
  End Subroutine print_tracking_species

  Subroutine copy_to_trajectory(traj_data, model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to copy model arrays to each 
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data
    Integer(Kind=wi),   Intent(In   ) :: frame

    Integer(Kind=wi) :: l
    
    traj_data%config(frame,:)%tag=model_data%config%atom(:)%tag
    traj_data%config(frame,:)%element=model_data%config%atom(:)%element
    traj_data%box(frame)%cell=model_data%config%cell
    traj_data%box(frame)%invcell=model_data%config%invcell
    traj_data%box(frame)%volume=model_data%config%volume
    traj_data%box(frame)%cell_length=model_data%config%cell_length
    Do l = 1,3
      traj_data%config(frame,:)%r(l)=model_data%config%atom(:)%r(l)
    End Do
    
    ! Copy tracked species only if reactive_chemistry is set to True
    If(model_data%reactive_chemistry%stat) Then 
      Do l = 1, model_data%reactive_species%N0%value
        traj_data%track_chem%config(frame,l)%r=model_data%track_chem(l)%r
        traj_data%track_chem%config(frame,l)%indx=model_data%track_chem(l)%indx
        traj_data%track_chem%config(frame,l)%tag=model_data%track_chem(l)%tag
      End Do
    End If

    ! Copy to species arrays
    If (model_data%nonreactive_species%invoke%fread) Then
      Do l = 1, model_data%config%Nmax_species
         traj_data%species(frame,l)%alive=model_data%config%species(l)%alive
         traj_data%species(frame,l)%list=model_data%config%species(l)%list
      End Do
    End If
      
  End Subroutine copy_to_trajectory
 
End Module trajectory
