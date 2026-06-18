!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the Special Pair Correlation Function (TCF) for the
! resctive species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module spcf

  Use atomic_model,     Only: model_type
                        
                              
  Use fileset,          Only: file_type, &
                              FILE_SPCF_ALL, &
                              FILE_SPCF_AVG, &
                              FILE_SET, &
                              refresh_out

  Use input_types,      Only: in_logic,   &
                              in_string

  Use numprec,          Only: wi,& 
                              wp

  Use process_data,     Only: set_read_status, &
                              capital_to_lower_case, &
                              check_for_rubbish, &
                              get_word_length
                              
  Use trajectory,       Only: traj_type, &
                              average_segments, &
                              within_region

  Use unit_output,      Only: info, &
                              error_stop                                                             
 
  Implicit None
  
  Private
  !Type for spcf
  Type, Public :: spcf_type
    Type(in_string)  :: invoke
    Type(in_string)  :: method
    Type(in_logic)   :: print_all_segments
  End Type
  
  Public :: read_spcf, check_spcf
  Public :: special_pair_correlation_function

Contains

  Subroutine read_spcf(iunit, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &spcf block
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(spcf_type), Intent(InOut)  ::  spcf_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &spcf block (within the &reactive_analysis block, SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_spcf" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_spcf') Exit
      Call check_for_rubbish(iunit, '&spcf')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='method') Then
        Read (iunit, Fmt=*, iostat=io) word, spcf_data%method%type
        Call set_read_status(word, io, spcf_data%method%fread,&
                                     & spcf_data%method%fail,&
                                     & spcf_data%method%type)

      Else If (word(1:length) == 'print_all_segments') Then
       Read (iunit, Fmt=*, iostat=io) word, spcf_data%print_all_segments%stat
       Call set_read_status(word, io, spcf_data%print_all_segments%fread, spcf_data%print_all_segments%fail)
                            
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_spcf"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_spcf
  
  Subroutine check_spcf(files, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &spcf block
    !
    ! author    - i.scivetti June 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(spcf_type),    Intent(InOut) :: spcf_data

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &spcf block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (spcf_data%method%fread) Then
      If (spcf_data%method%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "method" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(spcf_data%method%type)/='hicf'     .And. &
            Trim(spcf_data%method%type)/='hdcf')  Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    & 'Wrong input for "method". Valid options:&
                                    & "HICF" and "HDCF"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "method" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "HICF" and "HDCF"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    If (spcf_data%print_all_segments%fread) Then
      If (spcf_data%print_all_segments%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_segments" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      spcf_data%print_all_segments%stat=.False.
    End If
    
  End Subroutine check_spcf    
  
  Subroutine special_pair_correlation_function(files, model_data, traj_data, spcf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the correlation function of the special pair (SPCF),
    ! which is defined by the atomic pair of the chemical site and the closest NN, 
    ! either donor or acceptor
    !
    ! author    - i.scivetti April 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(In   ) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(spcf_type),   Intent(InOut) :: spcf_data

    Integer(Kind=wi)   :: i, j, k, l, m, n, ntop
    Integer(Kind=wi)   :: iunit, ini_indx
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Logical            :: set_t0
    Character(Len=256) :: message
    
    Logical            :: terminated(traj_data%seg_analysis%Np_segment, model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: indexes(2,model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: ref_indx(2,model_data%reactive_species%N0%value)
    
    Integer(Kind=wi)   :: time_indx(model_data%reactive_species%N0%value)
    Integer(Kind=wi)   :: Nnet
    
    Logical            :: follow(model_data%reactive_species%N0%value)
    Logical            :: first(model_data%reactive_species%N0%value)
    Logical            :: flag
    Logical            :: hold(model_data%reactive_species%N0%value)

    Real(Kind=wp)      :: tchange(model_data%reactive_species%N0%value)
    Character(Len=256) :: method
    
    method=spcf_data%method%type
    
    If (spcf_data%print_all_segments%stat) Then
     ! Print header
     Open(Newunit=files(FILE_SPCF_ALL)%unit_no, File=files(FILE_SPCF_ALL)%filename, Status='Replace')
     iunit=files(FILE_SPCF_ALL)%unit_no
     Write (iunit,'(a)') '#  Special Pair Correlation Function (SPCF) associated to reactive species' 
     Write (iunit,'(a)') '#  Results for all the time segments' 
     Write (iunit,'(a)') '#  Time (ps)         SPCF'
    End If

    !Set max_points to beyond the segment
    traj_data%seg_analysis%max_points=traj_data%seg_analysis%Np_segment+1

    Do k= 1, traj_data%seg_analysis%N_seg
      set_t0=.True.
      follow=.True.
      first=.True.
      tchange=0.0_wp
      hold=.False.
      ! Initialise terminated tag
      Do m = 1, model_data%reactive_species%N0%value
        terminated(:,m)=.False.
      End Do
      
      ini_indx=traj_data%seg_analysis%seg_indx(1,k)
      Do i = traj_data%seg_analysis%seg_indx(1,k), traj_data%seg_analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        If (Trim(method)=='hicf' .Or. Trim(method)=='hdcf') Then
          If (set_t0) Then
            Do m = 1, model_data%reactive_species%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
              ref_indx(1,m)=indexes(1,m)
              ref_indx(2,m)=indexes(2,m)
            End Do  
            set_t0=.False.
          Else
    
            Do m = 1, model_data%reactive_species%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%nn_indx(1)/=indexes(2,m)) Then
                   indexes(1,m)=traj_data%track_chem%config(i,m)%indx
                   indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
                End If
                
                If ((indexes(1,m) /= ref_indx(1,m) .Or. indexes(2,m) /= ref_indx(2,m)) .And. &
                    (indexes(1,m) /= ref_indx(2,m) .Or. indexes(2,m) /= ref_indx(1,m))) Then
                  If (.Not. hold(m)) Then
                   tchange(m)=time
                   hold(m)=.True.
                   time_indx(m)=i
                  End If
                Else
                  hold(m)=.False.
                End If
    
                If (hold(m)) Then
                  If (Trim(method)=='hicf') Then
                    ntop=i
                  Else If (Trim(method)=='hdcf') Then
                    ntop=traj_data%seg_analysis%seg_indx(2,k)
                    follow(m)=.False.
                  End If
                  
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  hold(m)=.False.
                End If
              End If
            End Do
          End If
        End If

        If (Trim(method)=='hicf*' .Or. Trim(method)=='hdcf*') Then
          If (set_t0) Then
            Do m = 1, model_data%reactive_species%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
              ref_indx(1,m)=indexes(1,m)
            End Do  
            set_t0=.False.
          Else
    
            Do m = 1, model_data%reactive_species%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%nn_indx(1)/=indexes(2,m)) Then
                   indexes(1,m)=traj_data%track_chem%config(i,m)%indx
                   indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
                End If
                
                If (indexes(1,m) /= ref_indx(1,m) .And. indexes(2,m) /= ref_indx(1,m)) Then
                  If (.Not. hold(m)) Then
                   tchange(m)=time
                   hold(m)=.True.
                   time_indx(m)=i
                  End If
                Else
                  hold(m)=.False.
                End If
    
                If (hold(m)) Then
                  If (Trim(method)=='hicf*') Then
                    ntop=i
                  Else If (Trim(method)=='hdcf*') Then
                    ntop=traj_data%seg_analysis%seg_indx(2,k)
                    follow(m)=.False.
                  End If
                  
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  hold(m)=.False.
                End If
              End If
            End Do
          End If
        End If
        
      End Do
      
      l=0
      base_time=(traj_data%seg_analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%seg_analysis%seg_indx(1,k), traj_data%seg_analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        l=l+1
        Nnet=0
        Do j = 1, model_data%reactive_species%N0%value
          If (traj_data%region%define%fread) Then
            Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
          Else
            flag=.True.
          End If
          If (flag) Then
            Nnet=Nnet+1
          End If
        End Do
    
        suma_i=0.0_wp
        If (.Not. All(terminated(l,:))) Then
          suma_i=0.0_wp
          Do j = 1, model_data%reactive_species%N0%value
            If(.Not. terminated(l,j)) Then
              If (traj_data%region%define%fread) Then
                Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
              Else
                flag=.True.
              End If
              If (flag) Then
                suma_i=suma_i+1.0_wp    
              End If
            End If  
          End Do
          If (Nnet > 0) Then
            suma_i=suma_i/Nnet
          End If
        End If
       
        traj_data%seg_analysis%variable(l,k)=suma_i
        If (spcf_data%print_all_segments%stat) Then
          Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          If (i==traj_data%seg_analysis%seg_indx(2,k) .And. (traj_data%seg_analysis%N_seg /=1)) Then
             If (k /= traj_data%seg_analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)        SPCF' 
             End If
          End If   
        End If
      End Do
    End Do
    
    If (spcf_data%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The SPCF analysis for the multiple time segments was printed to the "'&
                                 &//Trim(files(FILE_SPCF_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The SPCF analysis was printed to the "'//Trim(files(FILE_SPCF_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time segment.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    
    ! Compute average
    Call average_segments(files, traj_data, FILE_SPCF_AVG, 'SPCF')
    If (traj_data%seg_analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time segment was used to compute the average SPCF! The computed STD&
                              & is zero. Use/Check the &segment_trajectory block to improve the statistics.'
      Call info(message, 1)
    End If
    
    If (.Not. spcf_data%print_all_segments%stat) Then
      If (traj_data%seg_analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the SPCF analysis for all time segments,&
                                & the "print_all_segments" directive (within the &spcf block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%seg_analysis%N_seg ==1 .And. (.Not. traj_data%seg_analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_SPCF_ALL)%filename)//'" and "'//Trim(files(FILE_SPCF_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If

    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine special_pair_correlation_function

End Module spcf
