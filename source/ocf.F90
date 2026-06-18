!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the Orientational Correlation Function (OCF)
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module ocf

  Use atomic_model,     Only: model_type, &
                              check_PBC

  Use fileset,          Only: file_type, &
                              FILE_OCF_NONREACTIVE_ALL, &
                              FILE_OCF_NONREACTIVE_AVG, &
                              FILE_OCF_REACTIVE_ALL, &
                              FILE_OCF_REACTIVE_AVG, &
                              FILE_SET, &
                              refresh_out
                              
  Use input_types,      Only: in_integer, &
                              in_logic,   &
                              in_string

  Use numprec,          Only: wi,& 
                              wp

  Use process_data,     Only: capital_to_lower_case, &
                              check_for_rubbish, &
                              get_word_length, &
                              remove_symbols, &
                              set_read_status                              
                           
  Use trajectory,       Only: traj_type, &
                              average_segments, &
                              within_region   
                              
  Use unit_output,      Only: info, &
                              error_stop 
      
  Implicit None
  Private
  
  Type, Public :: ocf_type
    Private
    Type(in_string),  Public  :: invoke
    Type(in_integer), Public  :: legendre_order
    Type(in_string),  Public  :: u_definition
    Type(in_logic),   Public  :: print_all_segments
  End Type ocf_type

  Public :: read_ocf_settings
  Public :: check_ocf_nonreactive_species, check_ocf_reactive_species
  Public :: compute_ocf_nonreactive_species, compute_ocf_reactive_species
  
Contains

  Subroutine read_ocf_settings(iunit, ocf, label)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &ocf_reactive 
    ! and ocf_nonreactive blocks
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(ocf_type),    Intent(InOut) :: ocf
    Character(Len=*),  Intent(In   ) :: label 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &'//Trim(label)//' block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_'//Trim(label)//'" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_'//Trim(label)) Exit
      Call check_for_rubbish(iunit, Trim(label))

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='u_definition') Then
        Read (iunit, Fmt=*, iostat=io) word, ocf%u_definition%type
        Call set_read_status(word, io, ocf%u_definition%fread,&
                                     & ocf%u_definition%fail,&
                                     & ocf%u_definition%type)

      Else If (Trim(word)=='legendre_order') Then
         Read (iunit, Fmt=*, iostat=io) word, ocf%legendre_order%value
         Call set_read_status(word, io, ocf%legendre_order%fread,&
                            & ocf%legendre_order%fail)
                            
      Else If (word(1:length) == 'print_all_segments') Then
       Read (iunit, Fmt=*, iostat=io) word, ocf%print_all_segments%stat
       Call set_read_status(word, io, ocf%print_all_segments%fread,&
                         & ocf%print_all_segments%fail)
                            
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings. See the "use_code.md" file.&
                                & Have you properly closed the block with "&end_'//Trim(label)//'"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_ocf_settings

  Subroutine check_ocf_reactive_species(files, ocf_reactive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &ocf_reactive block
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(In   ) :: files(:)
    Type(ocf_type),   Intent(InOut) :: ocf_reactive

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &ocf_reactive block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (ocf_reactive%legendre_order%fread) Then
      If (ocf_reactive%legendre_order%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "legendre_order" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (ocf_reactive%legendre_order%value < 1 .Or. ocf_reactive%legendre_order%value>4) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                &'Input value for "legendre_order" must be a value between 1 and 4 (polynomial order).&
                                & We recommend setting this value to 2.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "legendre_order" directive.&
                                      & We recommend setting this value to 2.'
       Call info(messages, 1)
       Call error_stop(' ')
    End If    
    
    If (ocf_reactive%u_definition%fread) Then
      If (ocf_reactive%u_definition%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "variable" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(ocf_reactive%u_definition%type)/='special_pair'     .And. &
            Trim(ocf_reactive%u_definition%type)/='unrattled_special_pair')  Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    & 'Wrong input for "variable". Valid options:&
                                    & "special_pair" or "unrattled_special_pair"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "variable" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "special_pair" or "unrattled_special_pair"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    If (ocf_reactive%print_all_segments%fread) Then
      If (ocf_reactive%print_all_segments%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_segments" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      ocf_reactive%print_all_segments%stat=.False.
    End If
    
  End Subroutine check_ocf_reactive_species

  Subroutine compute_ocf_reactive_species(files, model_data, traj_data, ocf_reactive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational chemistry along the trajectory
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(In   ) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(InOut) :: ocf_reactive

    Integer(Kind=wi)   :: i, j, k, l, m
    Integer(Kind=wi)   :: iunit, ini_indx, ncouples
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Real(Kind=wp)      :: u(3,model_data%reactive_species%N0%value), u0(3,model_data%reactive_species%N0%value)

    Logical            :: set_s0, modified
    Logical            :: first_change(model_data%reactive_species%N0%value)
    Integer(Kind=wp)   :: first_index(model_data%reactive_species%N0%value)

    Character(Len=256) :: message
    
    Integer(Kind=wi)   :: indexes(2,model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: ref_indx(model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: sites(4,traj_data%seg_analysis%Np_segment,model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: s1, s2
 
    If (ocf_reactive%print_all_segments%stat) Then
      ! Print header
      Open(Newunit=files(FILE_OCF_REACTIVE_ALL)%unit_no, File=files(FILE_OCF_REACTIVE_ALL)%filename,  Status='Replace')
      iunit=files(FILE_OCF_REACTIVE_ALL)%unit_no
      Write (iunit,'(a)') '#  Orientational Correlation Function for "reactive" species (OCF_REACTIVE)' 
      Write (iunit,'(a)') '#  Results for all the time segments' 
      Write (iunit,'(a)') '#  Time (ps)         OCF' 
    End If

    !Set max_points to beyond the segment
    traj_data%seg_analysis%max_points=traj_data%seg_analysis%Np_segment+1

    Do k= 1, traj_data%seg_analysis%N_seg
      set_s0=.True.
      first_change=.True.
      l=0
      ini_indx=traj_data%seg_analysis%seg_indx(1,k)
      Do i = traj_data%seg_analysis%seg_indx(1,k), traj_data%seg_analysis%seg_indx(2,k)
        l=l+1
        If (Trim(ocf_reactive%u_definition%type)=='special_pair') Then
        ! Obtain the index of the closest acceptor
          Do m = 1, model_data%reactive_species%N0%value
            sites(1,l,m)=traj_data%track_chem%config(i,m)%indx
            sites(2,l,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
          End Do
        Else If (Trim(ocf_reactive%u_definition%type)=='unrattled_special_pair') Then
          If (set_s0) Then
            Do m = 1, model_data%reactive_species%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%indx
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_s0=.False.
          Else
            Do m = 1, model_data%reactive_species%N0%value
              If (traj_data%track_chem%config(i,m)%indx/=indexes(2,m)) Then
                indexes(1,m)=indexes(2,m)
                indexes(2,m)=traj_data%track_chem%config(i,m)%indx
                If (first_change(m)) Then
                  Do j=1, l
                    sites(1,j,m)=ref_indx(m)
                    sites(2,j,m)=indexes(2,m)
                  End Do
                  first_change(m)=.False.
                  first_index(m)=l
                Else
                  sites(1,l,m)=indexes(1,m)
                  sites(2,l,m)=indexes(2,m)
                End If
              Else
                If (.Not. first_change(m)) Then
                  sites(1,l,m)=indexes(1,m)
                  sites(2,l,m)=indexes(2,m)
                End If
              End If
            End Do
          End If
        End If  
      End Do
      
      ! Set initial vector for the transfer couple at the start of the time segment
      Do m=1, model_data%reactive_species%N0%value
        s1=sites(1,1,m)
        s2=sites(2,1,m)
        
        u0(:,m)=traj_data%config(ini_indx,s2)%r-traj_data%config(ini_indx,s1)%r
        Call check_PBC(u0(:,m), traj_data%box(ini_indx)%cell, traj_data%box(ini_indx)%invcell, 0.5_wp, modified)
        u0(:,m)=u0(:,m)/norm2(u0(:,m))
      End Do

      l=0
      ! Compute the orientational correlation function from the chaning chemistry
      base_time=(traj_data%seg_analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%seg_analysis%seg_indx(1,k), traj_data%seg_analysis%seg_indx(2,k)
        l=l+1
        time=(i-1)*traj_data%timestep%value
        suma_i=0.0_wp
        ncouples=0
        Do m=1, model_data%reactive_species%N0%value
          s1=sites(1,l,m)
          s2=sites(2,l,m)

          u(:,m)=traj_data%config(i,s2)%r-traj_data%config(i,s1)%r
          Call check_PBC(u(:,m), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
          u(:,m)=u(:,m)/norm2(u(:,m))

          If (Trim(ocf_reactive%u_definition%type)=='special_pair') Then
            Call orientational_correlation_term_transfer_couple(traj_data, i, s1, u(:,m), u0(:,m), suma_i, ncouples)
          Else If (Trim(ocf_reactive%u_definition%type)=='unrattled_special_pair') Then
            If (l<first_index(m)) Then
              Call orientational_correlation_term_transfer_couple(traj_data, i, s1, u(:,m), u0(:,m), suma_i, ncouples)
            Else
              Call orientational_correlation_term_transfer_couple(traj_data, i, s2, u(:,m), u0(:,m), suma_i, ncouples) 
            End If
          End If
        End Do 
        
        If (i==traj_data%seg_analysis%seg_indx(2,k)) Then
          If (ncouples /= 0) Then
            suma_i=suma_i/ncouples
            traj_data%seg_analysis%variable(l,k)=suma_i
            If (ocf_reactive%print_all_segments%stat) Then
              Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If
          If (ocf_reactive%print_all_segments%stat) Then
            If ((traj_data%seg_analysis%N_seg /=1) .And. (k /= traj_data%seg_analysis%N_seg)) Then
              If (k /= traj_data%seg_analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)         OCF' 
              End If
            End If
           End If 
        Else
          If (ncouples /= 0) Then
            suma_i=suma_i/ncouples
          Else
            suma_i=0.0_wp
          End If
          traj_data%seg_analysis%variable(l,k)=suma_i
          If (ocf_reactive%print_all_segments%stat) Then
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If  
        End If
        
      End Do
      
    End Do  

    If (ocf_reactive%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The OCF_REACTIVE analysis for the multiple time segments was printed to the "'&
                                 &//Trim(files(FILE_OCF_REACTIVE_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The OCF_REACTIVE analysis was printed to the "'&
                                 &//Trim(files(FILE_OCF_REACTIVE_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time segment.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    ! Compute average
    Call average_segments(files, traj_data, FILE_OCF_REACTIVE_AVG, 'OCF_REACTIVE')
    If (traj_data%seg_analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time segment was used to compute the average OCF_REACTIVE! The computed STD&
                              & is zero. Use/Check the &segment_trajectory block to improve the statistics.'
      Call info(message, 1)                        
    End If
    
    If (.Not. ocf_reactive%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the OCF_REACTIVE analysis for all time segments,&
                                & the "print_all_segments" directive (within the &ocf_reactive) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%seg_analysis%N_seg ==1 .And. (.Not. traj_data%seg_analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_OCF_REACTIVE_ALL)%filename)//'" and "'&
                               &//Trim(files(FILE_OCF_REACTIVE_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If

    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine compute_ocf_reactive_species 

  Subroutine orientational_correlation_term_transfer_couple(traj_data, i, s2, u, u0, suma_i, ncouples)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation for the relevant 
    ! transfer couple at the MD frame i (cij)
    !
    ! author    - i.scivetti Sept 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: s2
    Real(Kind=wp),     Intent(In   ) :: u(3)
    Real(Kind=wp),     Intent(In   ) :: u0(3)
    Real(Kind=wp),     Intent(InOut) :: suma_i
    Integer(Kind=wi),  Intent(InOut) :: ncouples
    
    Real(Kind=wp)     :: x, cij
    Logical           :: flag
    
    If (traj_data%region%define%fread) Then
      Call within_region(traj_data, i, s2, flag)
    Else
      flag=.True.
    End If

    If (flag) Then
      x=Dot_product(u,u0)
      ncouples=ncouples+1
      cij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
      suma_i=suma_i+cij
    End If

  End Subroutine orientational_correlation_term_transfer_couple  
  
  Subroutine check_ocf_nonreactive_species(files, ocf_nonreactive)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    ! Subroutine to check the settings of the &OCF_NONREACTIVE block  
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(ocf_type),     Intent(InOut) :: ocf_nonreactive

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &OCF_NONREACTIVE block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (ocf_nonreactive%legendre_order%fread) Then
      If (ocf_nonreactive%legendre_order%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "legendre_order" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (ocf_nonreactive%legendre_order%value < 1 .Or. ocf_nonreactive%legendre_order%value>4) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                &'Input value for "legendre_order" must be a value between 1 and 4 (polynomial order).&
                                & We recommend setting this value to 2.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "legendre_order" directive.&
                                      & We recommend setting this value to 2.'
       Call info(messages, 1)
       Call error_stop(' ')
    End If

    If (ocf_nonreactive%u_definition%fread) Then
      If (ocf_nonreactive%u_definition%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "u_definition" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(ocf_nonreactive%u_definition%type)/='bond_12'  .And. &
            Trim(ocf_nonreactive%u_definition%type)/='bond_13'  .And. &
            Trim(ocf_nonreactive%u_definition%type)/='bond_123' .And. &
            Trim(ocf_nonreactive%u_definition%type)/='bond_12-13'  .And. &
            Trim(ocf_nonreactive%u_definition%type)/='plane') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "u_definition". Valid options: "bond_12", "bond_13",&
                                    & "bond_12-13", "bond_123" or "plane"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "u_definition" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If
    
    If (ocf_nonreactive%print_all_segments%fread) Then
      If (ocf_nonreactive%print_all_segments%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_segments" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      ocf_nonreactive%print_all_segments%stat=.False.
    End If
    
  End Subroutine check_ocf_nonreactive_species
 
  Subroutine compute_ocf_nonreactive_species(files, traj_data, ocf_nonreactive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational correlation function (OCF)
    ! of nonreactive species.
    ! Different possible flavours are available depending on the settings
    ! of the &OCF_NONREACTIVE block
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(InOut) :: ocf_nonreactive

    Integer(Kind=wi)  :: i, j, k, l
    Integer(Kind=wi)  :: Nini_species, iunit
    Real(Kind=wp)     :: suma_i 
    Real(Kind=wp)     :: time
    Real(Kind=wp)     :: base_time
    Logical           :: set_u0
    Character(Len=256) :: message
    
    Logical           :: terminated(traj_data%Nmax_species)

    If (ocf_nonreactive%print_all_segments%stat) Then
      ! Print header
      Open(Newunit=files(FILE_OCF_NONREACTIVE_ALL)%unit_no, File=files(FILE_OCF_NONREACTIVE_ALL)%filename, Status='Replace')
      iunit=files(FILE_OCF_NONREACTIVE_ALL)%unit_no
      Write (iunit,'(a)') '#  Orientational Correlation Function (OCF) for the selected non-reactive species' 
      Write (iunit,'(a)') '#  Results for all the time segments' 
      Write (iunit,'(a)') '#  Time (ps)         OCF'
    End If
    
    !Set max_points to beyond the segment
    traj_data%seg_analysis%max_points=traj_data%seg_analysis%Np_segment+1

    Do k= 1, traj_data%seg_analysis%N_seg
      set_u0=.True.
      l=0
      ! Initialise terminated tag
      Do j = 1, traj_data%Nmax_species
        terminated(j)=.False.
      End Do
      base_time=(traj_data%seg_analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%seg_analysis%seg_indx(1,k), traj_data%seg_analysis%seg_indx(2,k)
        l=l+1
        time=(i-1)*traj_data%timestep%value
        If (.Not. set_u0) Then
          Do j=1,3
            traj_data%species(i,:)%u0(j,1)=traj_data%species(i-1,:)%u0(j,1)
            If (Trim(ocf_nonreactive%u_definition%type) == 'bond_12-13') Then
              traj_data%species(i,:)%u0(j,2)=traj_data%species(i-1,:)%u0(j,2)
            End If
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_nonreactive_species(traj_data, ocf_nonreactive, i, j)
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              If (Trim(ocf_nonreactive%u_definition%type) == 'bond_12-13') Then
                traj_data%species(i,j)%u0(:,2)=traj_data%species(i,j)%u(:,2)
              End If
              Nini_species=Nini_species+1
            Else
              terminated(j)=.True.
            End If
          End Do
          set_u0=.False.
          If (Nini_species==0) Then
            Write (message,'(1x,a,2x,i6,a)') '***PROBLEMS: the code could not identify a single nonreactive species for frame ', i,&
                                            & '. Please review the settings for the &selected_nonreactive_species block'
            Call info(message, 1)
            Call error_stop(' ')
          End If
        End If
      
        suma_i=0.0_wp
        traj_data%N_species=0
        Do j = 1, traj_data%Nmax_species
          If(.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_nonreactive_species(traj_data, ocf_nonreactive, i, j)
              Call orientational_correlation_term_nonreactive_species(traj_data, ocf_nonreactive, i, j, suma_i)  
            Else
              !terminated(j)=.True.
            End If
          End If  
        End Do
    
        If (i==traj_data%seg_analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
            suma_i=suma_i/traj_data%N_species
            traj_data%seg_analysis%variable(l,k)=suma_i
            If (ocf_nonreactive%print_all_segments%stat) Then
              Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If  
          terminated=.False.
          set_u0=.True.
          If (ocf_nonreactive%print_all_segments%stat) Then
            If ((traj_data%seg_analysis%N_seg /=1) .And. (k /= traj_data%seg_analysis%N_seg)) Then
              If (k /= traj_data%seg_analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)         OCF' 
              End If
            End If  
          End If                      
        Else
          If ((traj_data%N_species) /= 0) Then
            suma_i=suma_i/traj_data%N_species
          Else  
            suma_i=0.0_wp
          End If
          traj_data%seg_analysis%variable(l,k)=suma_i
          If (ocf_nonreactive%print_all_segments%stat) Then
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If
        End If
      End Do
    End Do
    
    If (ocf_nonreactive%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The OCF analysis for the multiple time segments was printed to the "'&
                                 &//Trim(files(FILE_OCF_NONREACTIVE_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The OCF analysis was printed to the "'//Trim(files(FILE_OCF_NONREACTIVE_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time segment.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    
    Call average_segments(files, traj_data, FILE_OCF_NONREACTIVE_AVG, 'OCF')
    If (traj_data%seg_analysis%N_seg ==1 ) Then
      Write (message,'(1x,a)') 'WARNING: A single time segment was used to compute the average OCF! The computed STD&
                                & is zero. Use/Check the &segment_trajectory block to improve the statistics.'
      Call info(message, 1)
    End If
    
    If (.Not. ocf_nonreactive%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the OCF analysis for all time segments,&
                                & the "print_all_segments" directive (within the &ocf_nonreactive block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%seg_analysis%N_seg ==1 .And. (.Not. traj_data%seg_analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_OCF_NONREACTIVE_ALL)%filename)//'" and "'&
                               &//Trim(files(FILE_OCF_NONREACTIVE_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If
    
    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine compute_ocf_nonreactive_species

  Subroutine rotation_vector_nonreactive_species(traj_data, ocf_nonreactive, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the rotation vector of the nonreactive species
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),   Intent(InOut)  :: ocf_nonreactive
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j

    Integer(Kind=wi) :: indx1, indx2, indx3, k
    Logical          :: modified
    Real(Kind=wp), Dimension(3)  :: u12, u13


    indx1=traj_data%species(i,j)%list(1)
    indx2=traj_data%species(i,j)%list(2)
    indx3=traj_data%species(i,j)%list(3)
    
    If (Trim(ocf_nonreactive%u_definition%type) == 'bond_12') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_nonreactive%u_definition%type) == 'bond_13') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_nonreactive%u_definition%type) == 'bond_12-13') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Do k=1,3
        traj_data%species(i,j)%u(k,1)=u12(k)
        traj_data%species(i,j)%u(k,2)=u13(k)
      End Do
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
      traj_data%species(i,j)%u(:,2)=traj_data%species(i,j)%u(:,2)/norm2(traj_data%species(i,j)%u(:,2))
    Else If (Trim(ocf_nonreactive%u_definition%type) == 'bond_123') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Do k=1, 3
        traj_data%species(i,j)%u(k,1)=u12(k)+u13(k)
      End Do
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_nonreactive%u_definition%type) == 'plane') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call cross_product(u12, u13, traj_data%species(i,j)%u(:,1))
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    End If    
    
  End Subroutine rotation_vector_nonreactive_species  


  Subroutine orientational_correlation_term_nonreactive_species(traj_data, ocf_nonreactive, i, j, suma_i)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation for the relevant 
    ! species j at the MD frame i (cij)
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(InOut) :: ocf_nonreactive
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j
    Real(Kind=wp),     Intent(InOut) :: suma_i
    
    Real(Kind=wp)     :: x, cij, c2ij, c0ij
    Logical           :: flag
    Integer(Kind=wi)  :: m
    
    If (traj_data%region%define%fread) Then
      m=traj_data%species(i,j)%list(1)
      Call within_region(traj_data, i, m, flag)
    Else
      flag=.True.
    End If
    
    If (flag) Then
      x=Dot_product(traj_data%species(i,j)%u(:,1),traj_data%species(i,j)%u0(:,1))
      traj_data%N_species=traj_data%N_species+1
      Select Case (ocf_nonreactive%legendre_order%value)  
        Case (1)
          cij=x
        Case (2)
          cij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
        Case (3)
          cij=(5.0_wp*(x)**3-3.0_wp*x)/2.0_wp
        Case (4)
          cij=(35.0_wp*(x)**4-30.0_wp*x**2+3.0_wp)/8.0_wp
      End Select  
      
      If (Trim(ocf_nonreactive%u_definition%type) == 'bond_12-13') Then
        x=Dot_product(traj_data%species(i,j)%u(:,2),traj_data%species(i,j)%u0(:,2))
        c0ij=cij
        Select Case (ocf_nonreactive%legendre_order%value)  
          Case (1)
            c2ij=x
          Case (2)
            c2ij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
          Case (3)
            c2ij=(5.0_wp*(x)**3-3.0_wp*x)/2.0_wp
          Case (4)
            c2ij=(35.0_wp*(x)**4-30.0_wp*x**2+3.0_wp)/8.0_wp
        End Select  
        cij=(c0ij+c2ij)/2.0_wp
      End If
      
      suma_i=suma_i+cij

    End If

  End Subroutine orientational_correlation_term_nonreactive_species  

  Subroutine cross_product(a, b, cross)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the cross_product 
    !
    ! author    - i.scivetti May 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3)
    Real(Kind=wp), Intent(In   ) :: b(3)
    Real(Kind=wp), Intent(  Out) :: cross(3) 

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  End Subroutine cross_product
  
End Module  
