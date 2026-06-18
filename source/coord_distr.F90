!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to obtain statistics for the coordinate (x, y, or z) 
! distribution of a selected species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:   -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module coord_distr

  Use atomic_model,     Only: model_type, &
                              check_length_directive

  Use fileset,          Only: file_type, &
                              FILE_COORD_DISTRIB, &
                              FILE_SET, &
                              refresh_out

  Use input_types,      Only: in_param, &
                              in_string

  Use numprec,          Only: wi,& 
                              wp

  Use process_data,     Only: capital_to_lower_case, &
                              check_for_rubbish, &
                              get_word_length, &
                              remove_symbols, &
                              set_read_status
                              
  Use trajectory,       Only: traj_type

  Use unit_output,      Only: info, &
                              error_stop 

  Implicit None
  Private
  
  !Type to describe the coordinate distribution
  Type, Public :: coord_distr_type
    Private
    Type(in_string),  Public :: invoke
    Type(in_string),  Public :: species_dir
    Character(Len=8), Public :: species
    Type(in_string),  Public :: coordinate
    Type(in_param),   Public :: delta
    Integer(Kind=wi) :: indx
  End Type
  
  Public :: read_coord_distrib, check_coord_distrib, compute_coordinate_distribution
  
Contains

  Subroutine read_coord_distrib(iunit, coord_distr_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for computing the coordinate distribution
    ! of selective species. Information must be provided in the 
    ! &coord_distrib block 
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),     Intent(In   ) :: iunit
    Type(coord_distr_type), Intent(InOut) :: coord_distr_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &coord_distrib block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_coord_distrib" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_coord_distrib') Exit
      Call check_for_rubbish(iunit, '&coord_distrib')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='species') Then
        Read (iunit, Fmt=*, iostat=io) coord_distr_data%species_dir%type, coord_distr_data%species
        Call set_read_status(word, io, coord_distr_data%species_dir%fread,&
                           & coord_distr_data%species_dir%fail,coord_distr_data%species_dir%type)

      Else If (Trim(word)=='delta') Then
         Read (iunit, Fmt=*, iostat=io) coord_distr_data%delta%tag, &
                                      & coord_distr_data%delta%value,&
                                      & coord_distr_data%delta%units 
         Call set_read_status(word, io, coord_distr_data%delta%fread, coord_distr_data%delta%fail)

      Else If (Trim(word)=='coordinate') Then
        Read (iunit, Fmt=*, iostat=io) word, coord_distr_data%coordinate%type
        Call set_read_status(word, io, coord_distr_data%coordinate%fread,& 
                           & coord_distr_data%coordinate%fail,&
                           & coord_distr_data%coordinate%type)
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_coord_distrib"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_coord_distrib

  Subroutine check_coord_distrib(files, model_data, coord_distr_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),         Intent(In   ) :: files(:)
    Type(model_type),        Intent(In   ) :: model_data
    Type(coord_distr_type),  Intent(InOut) :: coord_distr_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j
    Logical             :: flag

    Character(Len=8)  :: tg
    Character(Len=8)  :: coord(3)
    
    ! Define coordinates to check directive "coordinate"
    coord(1)='x'
    coord(2)='y'
    coord(3)='z'
    
    error_set = '***ERROR in the &coord_distrib block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (.Not. coord_distr_data%delta%fread) Then
      coord_distr_data%delta%tag='delta'
    End If
    Call check_length_directive(coord_distr_data%delta, error_set, .True., 'directive')

    ! Check definition of "species" directive
    If (coord_distr_data%species_dir%fread) Then
      If (coord_distr_data%species_dir%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "species" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "species" directive to&
                                    & compute the coordinate distribution.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
   ! Check if the definition of "species" is valid
    tg=Trim(coord_distr_data%species)
    Call remove_symbols(tg,'*')
    flag=.True.
    j=1
    Do While (j <= model_data%reference_composition%atomic_species .And. flag)
      If (Trim(model_data%reference_composition%tag(j))==Trim(tg)) Then
        flag=.False.
      End If  
      j=j+1
    End Do
    If (flag) Then
      Write (messages(1),'(2(1x,a))') Trim(error_set), '"'//Trim(coord_distr_data%species)//'"&
                                     & defined for the "species" directive is not a valid species.&
                                     & Please review the definition of the &reference_composition block' 
      Call info(messages, 1)
      Call error_stop(' ') 
    End If 

    ! Check definition of "species" directive
    If (coord_distr_data%species_dir%fread) Then
      If (coord_distr_data%species_dir%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "species" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "species" directive to&
                                    & compute the coordinate distribution.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    ! Check definition of "coordinate" directive
    If (coord_distr_data%coordinate%fread) Then
      If (coord_distr_data%coordinate%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "coordinate" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        flag=.True. 
        Do j=1, 3
          If (Trim(coord_distr_data%coordinate%type)==Trim(coord(j))) Then
            flag=.False.
            coord_distr_data%indx=j
          End If
        End Do
        If (flag) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Definition for the "coordinate" directive&
                                    & is not valid. Valid options: "x", "y" or "z".&
                                    & Check correctness of the directives within the block.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
      
      
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "coordinate" directive to&
                                    & compute the coordinate distribution. Valid options: "x", "y" or "z".&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine check_coord_distrib
  
  Subroutine compute_coordinate_distribution(files, model_data, traj_data, coord_distr_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distribution of the coordinates (x, y or z) of
    ! the species selected in the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2025
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(traj_type),      Intent(InOut) :: traj_data
    Type(coord_distr_type), Intent(InOut) :: coord_distr_data

    Integer(Kind=wi)  :: i, j, m, iunit, indx_ini
    Integer(Kind=wi)  :: num_at, nbins, net_frames
    Integer(Kind=wi)  :: accum
    
    Real(Kind=wp)     :: clim(2), coord_value
    Real(Kind=wp)     :: vector(3) 

    Integer(Kind=wi)  :: list_indx(model_data%config%num_atoms)
    
    Character(Len=256) :: ctap
    
    Character(Len=256) :: messages(3), message
    Logical            :: falloc

    Character(Len=256) :: type_error
    Integer(Kind=wi)   :: fail(2)  
    
    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
    
    ! Search for the value of cmax and cmin 
    indx_ini=traj_data%seg_analysis%frame_ini
    vector=0.0_wp
    Do i = 1, 3
      vector(:)=vector(:)+traj_data%box(indx_ini)%cell(i,:)
    End Do
    If (vector(coord_distr_data%indx) > 0.0_wp) Then
      clim(2)=vector(coord_distr_data%indx)
      clim(1)=0.0_wp
      ctap='top'
    Else
      clim(1)=vector(coord_distr_data%indx) 
      clim(2)=0.0_wp 
      ctap='bottom'
    End If
    
    Do i = traj_data%seg_analysis%frame_ini+1, traj_data%seg_analysis%frame_last
      vector=0.0_wp
      Do j = 1, 3
        vector(:)=vector(:)+traj_data%box(i)%cell(j,:)
      End Do
      If (Trim(ctap)=='bottom') Then
        If (vector(coord_distr_data%indx) < clim(1)) Then
          clim(1)=vector(coord_distr_data%indx)
        End If
      Else If (Trim(ctap)=='top') Then
        If (vector(coord_distr_data%indx) > clim(2)) Then
          clim(2)=vector(coord_distr_data%indx)
        End If
      End If
    End Do
     
    ! Define number of bins
    nbins=Floor(Abs(clim(1)-clim(2))/coord_distr_data%delta%value)

    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for coordinate distribution.&
                                & Analysis will not be executed.'
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
        ! Define the number and list of indexes
        num_at=0
        list_indx=0
        Do j = 1, model_data%config%num_atoms
          If (coord_distr_data%species==traj_data%config(i,j)%tag) Then
            num_at=num_at+1
            list_indx(num_at)=j
          End If
        End Do
      
        ! Calculate the histogram for this particular frame of the trajectory
        If (num_at/=0) Then
          h=0
          Do j=1, num_at
            coord_value=Abs(traj_data%config(i,list_indx(j))%r(coord_distr_data%indx))
            If (Trim(ctap)=='top') Then
              m=Floor(coord_value/coord_distr_data%delta%value)+1
            Else
              m=nbins+1-(Floor(coord_value/coord_distr_data%delta%value)+1)
            End If
            If (m <= nbins) Then
              h(m)=h(m)+1
            End If
          End Do 
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do m=1, nbins 
            d(m)= d(m)+Real(h(m),Kind=wp)/num_at
          End Do
        End If
        accum=accum+num_at
      End Do

      Do m=1, nbins 
        d(m)=d(m)/net_frames/coord_distr_data%delta%value      
      End Do
      
      ! Print results
      If (accum /= 0) Then
        ! Print File
        Open(Newunit=files(FILE_COORD_DISTRIB)%unit_no, File=files(FILE_COORD_DISTRIB)%filename, Status='Replace')
        iunit=files(FILE_COORD_DISTRIB)%unit_no
        Write (iunit,'(a)') '#  Distribution of the '//Trim(coord_distr_data%coordinate%type)//'-coordinate&
                           & for the "'//Trim(coord_distr_data%species)//'" species'
        Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
        
        Do m=1, nbins
          Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(m,Kind=wp)-0.5)*coord_distr_data%delta%value, d(m)
        End Do
        Write (message,'(1x,a)') 'The distribution of the '//Trim(coord_distr_data%coordinate%type)//&
                                &'-coordinate for the "'//Trim(coord_distr_data%species)//'" species was&
                                & printed to the "'//Trim(files(FILE_COORD_DISTRIB)%filename)//'" file.'
        Call info(message, 1)
      Else
        type_error=Trim(coord_distr_data%species)
        Write (messages(1),'(1x,a)') '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: coordinate distribution analysis could not be executed'
          Write (messages(2),'(1x,a)') '   Requested species '//Trim(type_error)//' as specified in the &coord_distrib&
                                  & block could not be identified along the trajectory.'
          Write (messages(3),'(1x,a)') '   Please verify the settings for the &coord_distrib block.' 
        Call info(messages, 3)
        Write (messages(1),'(1x,a)') '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Close file
      Close(iunit)
      ! Deallocate arrays   
      Deallocate(d,h)
   End If
    
   Call refresh_out(files) 
    
  End Subroutine compute_coordinate_distribution
  
End Module coord_distr
