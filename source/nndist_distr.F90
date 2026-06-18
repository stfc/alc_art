!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to obtain statistics for distribution of the NN 
! distances between selected species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:   -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module nndist_distr

  Use atomic_model,     Only: model_type, &
                              check_length_directive,& 
                              check_PBC

  Use constants,        Only: max_at_species, &
                              max_components
  
  Use fileset,          Only: file_type, &
                              FILE_SELECTED_NN_DISTANCES, &
                              FILE_SET, &
                              refresh_out
  
  Use input_types,      Only: in_param, &
                              in_string

  Use numprec,          Only: wi,& 
                              wp
 
  Use process_data,     Only: set_read_status, &
                              capital_to_lower_case, &
                              check_for_rubbish, &
                              prevent_segmentation, &
                              remove_symbols, &
                              get_word_length

  Use trajectory,       Only: traj_type, &
                              within_region

  Use unit_output,      Only: info, &
                              error_stop 
 
  Implicit None
  Private 
  
  ! Type for shortest distance
  Type, Public :: nndist_distr_type
    Private
    Type(in_string),  Public  :: invoke
    Type(in_string)    :: tag_reference_species
    Type(in_string)    :: tag_nn_species
    Character(Len=8)   :: reference_species
    Integer(Kind=wi)   :: num_nn_species
    Character(Len=8)   :: nn_species(max_at_species)
    Type(in_param)     :: lower_bound
    Type(in_param)     :: upper_bound
    Type(in_param)     :: dr 
  End Type nndist_distr_type
  
  Public :: read_selected_nn_distances, check_selected_nn_distances
  Public :: compute_nn_distance_distribution

Contains

  Subroutine read_selected_nn_distances(iunit, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read parameters from the
    ! &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),        Intent(In   ) :: iunit
    Type(nndist_distr_type), Intent(InOut) :: M 
    
    Integer(Kind=wi)   :: io, length, i
    Character(Len=256) :: message, word
    Character(Len=256) :: messages(2)
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in "&selected_nn_distances" block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_selected_nn_distances" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_selected_nn_distances') Exit
      Call check_for_rubbish(iunit, '&end_selected_nn_distances')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word
 
      Else If (Trim(word)=='reference_species') Then
        Read (iunit, Fmt=*, iostat=io) M%tag_reference_species%type, M%reference_species
        Call set_read_status(word, io, M%tag_reference_species%fread, M%tag_reference_species%fail, M%tag_reference_species%type)

      Else If (Trim(word)=='nn_species') Then
        Read (iunit, Fmt=*, iostat=io) M%tag_nn_species%type, M%num_nn_species
        Call prevent_segmentation(iunit, io, M%tag_nn_species%type, M%num_nn_species,&
                                & 'max_components', max_components, set_error)
        M%nn_species=' '
        Read (iunit, Fmt=*, iostat=io) M%tag_nn_species%type, M%num_nn_species, (M%nn_species(i), i=1, M%num_nn_species) 
        Call set_read_status(word, io, M%tag_nn_species%fread, M%tag_nn_species%fail, M%tag_nn_species%type)
        
      Else If (Trim(word)=='lower_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%lower_bound%tag, M%lower_bound%value, M%lower_bound%units 
         Call set_read_status(word, io, M%lower_bound%fread, M%lower_bound%fail)

      Else If (Trim(word)=='upper_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%upper_bound%tag, M%upper_bound%value, M%upper_bound%units 
         Call set_read_status(word, io, M%upper_bound%fread, M%upper_bound%fail)

      Else If (Trim(word)=='dr') Then
         Read (iunit, Fmt=*, iostat=io) M%dr%tag, M%dr%value, M%dr%units 
         Call set_read_status(word, io, M%dr%fread, M%dr%fail)

      Else
        Write (messages(1),'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_selected_nn_distances"? &
                                & Have you defined the directives correctly? See the "use_code.md" file'
        Call info (messages, 2)
        Call error_stop(' ')
      End If
    End Do
  
  End Subroutine read_selected_nn_distances

  Subroutine check_selected_nn_distances(files, model_data, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the definition of the
    ! parameters defined in the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),         Intent(In   ) :: files(:)
    Type(model_type),        Intent(In   ) :: model_data
    Type(nndist_distr_type), Intent(InOut) :: M 

    Character(Len=256)  :: messages(3), error_set
    Character(Len=8)    :: tagj, tagk
    Logical             :: flag
    Integer(Kind=wi)    :: j, k
    
    ! Error message just in case....
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&selected_nn_distances" block.'

    If(M%tag_nn_species%fread) Then
      If (M%tag_nn_species%fail) Then
        Write (messages(2),'(1x,a)')  'Problems to define the "nn_species" directive'  
        call info(messages, 2)
        call error_stop(' ')
      End If
    Else
      Write (messages(2),'(1x,a)')  'The user must define the "nn_species" directive'  
      call info(messages, 2)
      call error_stop(' ')
    End If

    If(M%tag_reference_species%fread) Then
      If (M%tag_reference_species%fail) Then
        Write (messages(2),'(1x,a)')  'Problems to define the "reference_species" directive'  
        call info(messages, 2)
        call error_stop(' ')
      End If
    Else
      Write (messages(2),'(1x,a)')  'The user must define the "reference_species" directive'  
      call info(messages, 2)
      call error_stop(' ')
    End If

    !Check if the reference_species is defined in the &reference_composition block  
    tagk=Trim(M%reference_species)
    Call remove_symbols(tagk,'*')
    flag=.True.
    j=1
    Do While (j <= model_data%reference_composition%atomic_species .And. flag)
      If (Trim(model_data%reference_composition%tag(j))==Trim(tagk)) Then
        flag=.False.
      End If  
      j=j+1
    End Do
    If (flag) Then
      Write (messages(2),'(1x,a)')   'The tag "'//Trim(M%reference_species)//&
                                     &'" defined in the "reference_species" directive is not a valid option.&
                                     & Please check the definition of the &reference_composition block' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If 
    
    !Check if tags in nn_species are defined in the &reference_composition block  
    Do k=1, M%num_nn_species
      tagk=Trim(M%nn_species(k))
      Call remove_symbols(tagk,'*')
      flag=.True.
      j=1
      Do While (j <= model_data%reference_composition%atomic_species .And. flag)
        If (Trim(model_data%reference_composition%tag(j))==Trim(tagk)) Then
          flag=.False.
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(2),'(1x,a)')   'The tag "'//Trim(M%nn_species(k))//&
                                       &'" defined in the "nn_species" directive is not a valid option.&
                                       & Please check the definition of the &reference_composition block' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If 
    End Do

    !Check if tags defined in nn_species are repeated
    Do j=1, M%num_nn_species-1
      tagj=Trim(M%nn_species(j))
      Do k=j+1, M%num_nn_species 
        tagk=Trim(M%nn_species(k))
        If (Trim(tagj)==Trim(tagk)) Then
          Write (messages(2),'(1x,a)')   'The tag "'//Trim(tagj)//&
                                         &'" is repeated in the specification of the "nn_species" directive.&
                                         & Please remove this duplication.' 
          Call info(messages, 2)
          Call error_stop(' ') 
        End If 
      End Do
    End Do

    !Check lower_bound, upper_bound and delta
    If (.Not. M%lower_bound%fread) Then
      M%lower_bound%tag='lower_bound'
    End If
    If (.Not. M%upper_bound%fread) Then
      M%upper_bound%tag='upper_bound'
    End If
    If (.Not. M%dr%fread) Then
      M%dr%tag='delta'
    End If
    
    Call check_length_directive(M%lower_bound, messages(1), .True., 'directive')
    Call check_length_directive(M%upper_bound, messages(1), .True., 'directive')
    Call check_length_directive(M%dr,          messages(1), .True., 'directive')
    If (M%lower_bound%value >= M%upper_bound%value) Then
      Write (messages(2),'(1x,a)')  'The value of "upper_bound" must be larger than "lower_bound"&
                                  & (make sure this is the case if you use different units)' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
     
  End Subroutine check_selected_nn_distances

  Subroutine compute_nn_distance_distribution(files, traj_data, model_data, nndist_distr_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of the shortest distance between
    ! the "reference_species" and those defined by the "nn_species" directive
    ! in the distance domain defined by the lower_bound and upper_bound
    ! Analysis is performed using the definitions of the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),         Intent(InOut) :: files(:)
    Type(traj_type),         Intent(InOut) :: traj_data
    Type(model_type),        Intent(In   ) :: model_data
    Type(nndist_distr_type), Intent(InOut) :: nndist_distr_data 

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k, mk
    Integer(Kind=wi)  :: num_at(2)
    Integer(Kind=wi)  :: list_indx(model_data%config%num_atoms,2)
    
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(5), message
    Logical            :: falloc, flag, flag1, found

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: rj(3), rk(3), rjk(3)
    Real(Kind=wp)  :: rmin
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(nndist_distr_data%upper_bound%value-nndist_distr_data%lower_bound%value)/nndist_distr_data%dr%value)
   
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining the statistics&
                                & of the shortest distances. Analysis will not be executed.'
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
        num_at=0
        list_indx=0
        Do j = 1, model_data%config%num_atoms
          If (nndist_distr_data%reference_species==traj_data%config(i,j)%tag) Then
              num_at(1)=num_at(1)+1
              list_indx(num_at(1),1)=j
          End If 
          k=1
          flag=.False.
          Do While (k <= nndist_distr_data%num_nn_species .And. (.Not. flag))
            If (nndist_distr_data%nn_species(k)==traj_data%config(i,j)%tag) Then
              num_at(2)=num_at(2)+1
              list_indx(num_at(2),2)=j
              flag=.True.  
            End If
            k=k+1
          End Do
        End Do
        
        If (num_at(1) /= 0 .And. num_at(2) /=0) Then
          Do j = 1, num_at(1)
            rj=traj_data%config(i,list_indx(j,1))%r
            rmin=Huge(1.0_wp)
            found=.False.
            Do k= 1, num_at(2)
              If (list_indx(k,2) /= list_indx(j,1)) Then
                rk=traj_data%config(i,list_indx(k,2))%r
                If (traj_data%region%define%fread) Then
                  Call within_region(traj_data, i, list_indx(j,1), flag1)
                  If (flag1) Then
                    flag=.True.
                  Else
                    flag=.False. 
                  End If
                Else
                   flag=.True.
                End If   
                If (flag) Then
                  rjk=rj-rk
                  Call check_PBC(rjk, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                  If (norm2(rjk) < rmin) Then
                    found=.True.
                    rmin=norm2(rjk)
                  End If
                End If
              End If
            End Do
            If (found) Then
              If (rmin >= nndist_distr_data%lower_bound%value .And. rmin <= nndist_distr_data%upper_bound%value) Then
                mk=Floor((rmin-nndist_distr_data%lower_bound%value)/nndist_distr_data%dr%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            End If
          End Do
        End If
        
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
           d(mk)=d(mk)/net_frames/nndist_distr_data%dr%value
         End Do
       
        ! Print File
        Open(Newunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no, File=files(FILE_SELECTED_NN_DISTANCES)%filename,&
                          &Status='Replace')
        iunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no
        Write (iunit,'(a)') '#  Probability distribution of the shortest distances between'
        Write (iunit,'(a)') '#  the reference species "'//Trim(nndist_distr_data%reference_species)//&
                             &'" and the species defined in the "nn_species" directive'    
        Write (iunit,'(a)') '#  Distance [Angstrom]     Probability [1/Angstrom]' 
        Do mk=1, nbins
          Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5_wp)*nndist_distr_data%dr%value+&
                                            & nndist_distr_data%lower_bound%value, d(mk)
        End Do
        Write (message,'(1x,a)') 'The probability distribution of shortest distances between the "reference_species"&
                                  & and the species defined in the "nn_species" directive was printed to the "'&
                                  &//Trim(files(FILE_SELECTED_NN_DISTANCES)%filename)//'" file.'
        Call info(message, 1)

      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Write (messages(2),'(1x,a)')   '   WARNING: the statistics for the shortest distances between "reference_species" and'
        Write (messages(3),'(1x,a)')   '   "nn_species" could not be executed. File "'&
                                         &//Trim(files(FILE_SELECTED_NN_DISTANCES)%filename)//'" was not generated.'
        Write (messages(4),'(1x,a)')   '   Please check the settings of the "&selected_nn_distances" block.       '
        Write (messages(5),'(1x,a)')   '************************************************************************************'
        Call info(messages, 5)
       End If

      Deallocate(d,h)
    End If
    
    Call refresh_out(files)
    
  End Subroutine compute_nn_distance_distribution
 
End Module nndist_distr
