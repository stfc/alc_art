# Instructions for use
ALC_ART requires two files:  

* TRAJECTORY: must contain the atomic positions recorded along the MD trajectory. To date, formats "xyz" and "vasp" (XDATCAR file) are the only two implemented options so far.
* SETTINGS: contains the instructions for data analysis.

Both SETTINGS and TRAJECTORY files must be present in the folder where ALC_ART is executed from, otherwise the program will print an error message and abort. Comments to the SETTINGS file can be added using the symbol "#". It is recommended to add a descriptive header in the first lines of the file for revision purposes. Different files will be printed depending on the selected options. The code is flexible enough to recognize directives independently of capitalization. For example, directive ***Ensemble***, ***eNsemBle***, ***ENSEMBLE***, etc, will all be interpreted as ***ensemble***. Together with the individual directives, it is also required to specify blocks and sub-blocks, which are declared with the character "***&***", followed by the ***name*** of the block (i.e. ***&name***), and must be closed with ***&end_name***. ALC_ART generates the OUTPUT file, which details the input settings (based on the content of the SETTINGS file) and the generated information from the data analysis.

Upon execution, the code first checks the correctness of the syntax and format for the defined settings. If a problem is found, an error message is printed to OUTPUT file and the execution is aborted. The code then reads the TRAJECTORY file, and the collected information is compared against the directives of the SETTINGS file. In case there is an inconsistency but the calculation can still proceed, ALC_ART will print a warning message (hopefully!). If the specified settings are incorrect, the program will print an error message instructing the user what to fix. 

The structure of the SETTINGS file can be divided in model and trajectory related directives/blocks. In the following sections, we will use the example case of nano-confined water in Nafion to explain the different functionalities. This model system constitutes an example of a highly reactive system, where chemical species change along the trajectory. The SETTINGS file can be found in the ***examples/Nafion-hydrated/*** folder. The TRAJECTORY file is also attached, so the user can run the code and change the settings. This trajectory is rather short and provided for the purpose of explaining the code and its functionalities.  

The files for the analysis of bulk liquid water (modelled by 64 water molecules in a cubic box) can be found in the ***examples/bulk-water/*** folder. In contrast to the Nafion, the model for bulk liquid water constitutes an example of a **non-reactive** system. We shall not discuss this case in the following notes.  The definition of the ***&reactive_species*** block  (see the section "Identification of chemistry changes" below) is NOT NEEDED for non-reactive trajectories.   

Accepted units are fs (femtoseconds) and ps (picoseconds) for time-related directives, and Angstrom and Bohr for distance-related directives. 

## Model settings
To describe the implemented functionalities for analysis, we will use the example case of nano-confined water in Nafion. The SETTINGS file can be found in the ***examples/Nafion-hydrated/*** folder.

### Format for the atomic positions
The format of the the TRAJECTORY file will be determined by the computational code used to run the MD simulation. To date, ALC_ART only accepts "xyz" and "vasp" formats. The format for the trajectory is specified as follows:  

***geometry_format***&nbsp;&nbsp;xyz

### Reference composition
To specify the atomic details of the model the code requires a reference composition, for which is necessary to define the ***&reference_composition*** block. This block allows tagging the atoms of the model according to the structure of the system under consideration. Tags are needed to classify elements with different chemical environments. For example, a model can have H atoms as part of water molecules and H atoms as part of a backbone membrane: the chemical element is the same (H), but the chemical environment around each hydrogen is completely different. The structure of this block for the present example is defined as follows:  

***&reference_composition***  
&nbsp;&nbsp;&nbsp;atomic_species&nbsp;&nbsp;  7  
&nbsp;&nbsp;&nbsp;tags     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cch&nbsp; Hch&nbsp;  Hs&nbsp;  Och&nbsp;  Sch&nbsp; Hw&nbsp;&nbsp; Ow  
&nbsp;&nbsp;&nbsp;amounts  &nbsp; 56&nbsp;&nbsp;&nbsp;       60&nbsp;&nbsp;&nbsp;  4&nbsp;&nbsp;&nbsp;&nbsp;  12 &nbsp;&nbsp;           4 &nbsp;&nbsp;       112 &nbsp; 56  
&nbsp;&nbsp;      elements &nbsp;  C&nbsp;&nbsp;&nbsp;&nbsp;  H &nbsp;&nbsp;&nbsp;       H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  O&nbsp;&nbsp;&nbsp;&nbsp; S &nbsp;&nbsp;&nbsp;   H &nbsp;&nbsp;&nbsp;&nbsp; O  
***&end_reference_composition***

The first directive must be ***atomic_species***, in this case equal to 7. This setting indicates that 7 different atomic tags will be used to label the atoms. The next task is to define the ***tags***, ***amounts*** and ***elements*** directives. The order for the definition of such directives is irrelevant. These settings describe what type of atoms and how many of them are part of the model. Atomic tags are selected based on a reference configuration (see the [DL_FIELD](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00323) code, for example). In reactive systems, the tags of certain atoms are expected to change with respect to the tagging that corresponds to the reference configuration. We shall return to this point later when we describe the ***&reactive_species*** block. Consistent with the definition of the ***&reference_composition*** block, each frame (atomic configuration) of the trajectory must have the following structure/order:  

* 56&nbsp;&nbsp; C atoms with tag Cch, followed by  
* 60&nbsp;&nbsp; H atoms with tag Hch, followed by  
* 4&nbsp;&nbsp;&nbsp;&nbsp;  H atoms with tag Hs,  followed by  
* 12&nbsp;&nbsp; O atoms with tag Och, followed by  
* 4&nbsp;&nbsp;&nbsp;&nbsp; S atoms with tag Sch, followed by 
* 112 H atoms with tag Hw, followed by
* 56&nbsp;&nbsp; O atoms with  tag Ow.

The suffixes "ch", "s" and "w" here indicate "chain", "sulphonic" and "water", respectively, but this labelling is arbitrary and subject to user's criteria. Ideally, the initial configuration of the trajectory should arrange atoms in sinlge groups to facilitate the definition of the block. Unfortunately this is not always possible. Still, the block is sufficiently flexible to allow multiple definitions of the same tags. For example, let's assume that the atomic configuration of 30 Hch atoms (out of 60) have been grouped at the end. The structure for the block would be:
 
***&reference_composition***  
&nbsp;&nbsp;&nbsp;atomic_species&nbsp;&nbsp;  8  
&nbsp;&nbsp;&nbsp;tags     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cch&nbsp; Hch&nbsp;  Hs&nbsp;  Och&nbsp;  Sch&nbsp; Hw&nbsp;&nbsp; Ow &nbsp;&nbsp; Hch  
&nbsp;&nbsp;&nbsp;amounts  &nbsp; 56&nbsp;&nbsp;&nbsp; 30&nbsp;&nbsp;&nbsp;  4&nbsp;&nbsp;&nbsp;&nbsp;  12 &nbsp;&nbsp;           4 &nbsp;&nbsp;       112 &nbsp; 56&nbsp;&nbsp; 30  
&nbsp;&nbsp;      elements &nbsp;  C&nbsp;&nbsp;&nbsp;&nbsp;  H &nbsp;&nbsp;&nbsp;       H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  O&nbsp;&nbsp;&nbsp;&nbsp; S &nbsp;&nbsp;&nbsp;   H &nbsp;&nbsp;&nbsp;&nbsp; O&nbsp;&nbsp;&nbsp; H  
***&end_reference_composition***

It is important to clarify that although the last group of 30 Hch atoms are equivalent to the previous 30 Hch, it is needed to increase "atomic_species" by 1, from 7 to 8. If there is an inconsistency between the settings in the block and the structure of the TRAJECTORY file, the code will inform the user and abort the execution.

### Simulation cell
Only for those trajectories recorded in "xyz" format (see the ***geometry_format*** directive above), the code will ask to define the simulation cell, which remains unchanged along the trajectory. In fact, trajectories in "xyz" format are only compatible with NVE and NVT simulations (see the ***ensemble*** directive below). Information for the simulation cell for the current example is set as follows:  

***cell_units*** &nbsp;&nbsp;  Angstrom  

***&simulation_cell***  
&nbsp;&nbsp;             10.064000&nbsp;&nbsp;  0.0000000&nbsp;&nbsp;    0.0000000  
&nbsp;&nbsp;&nbsp;&nbsp;  0.000000             13.0735195&nbsp;&nbsp;    0.0000000  
&nbsp;&nbsp;&nbsp;&nbsp;  0.000000&nbsp;&nbsp;  0.0000000               20.0000000  
***&end_simulation_cell***

Both ***cell_units*** and ***&simulation_cell*** are NOT NEEDED if the TRAJECTORY was recorded in "vasp" format. If there is an inconsistency between this block and the atomic positions of the trajectory, the code will inform the user and abort the execution.

### Identification of chemistry changes
Reactive systems are characterized by the breaking and formation of chemical bonds, which leads to the changes of the constituents species. To compute reactive systems, it is first necessary to set the following directive 

***reactive_chemistry***&nbsp;&nbsp;&nbsp;.True.

which instructs ALC_ART to search for changes of chemical species alonf the trajectory. Changes are monitored through the identification of donor and acceptor sites. The user needs to define the ***&reactive_species*** block, which in this example is set as follows: 

***&reactive_species***  
&nbsp;&nbsp;&nbsp;type &nbsp;&nbsp;&nbsp;&nbsp;donor  
&nbsp;&nbsp;&nbsp;total_number &nbsp;&nbsp;4  

&nbsp;&nbsp;&nbsp;&bonding_criteria  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number_of_bonds&nbsp;&nbsp;   3  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;only_element   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   H  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3  &nbsp;&nbsp;Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &extra_reactive_bonds  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;types_of_bonds&nbsp;&nbsp; 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Och&nbsp;&nbsp;Hw&nbsp;&nbsp;1.25&nbsp;&nbsp;  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Och&nbsp;&nbsp;Hs&nbsp;&nbsp;&nbsp;1.25&nbsp;&nbsp;&nbsp;Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&end_extra_reactive_bonds  
&nbsp;&nbsp;&nbsp;&end_bonding_criteria  

&nbsp;&nbsp;&nbsp;&search_environment  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;include_tags&nbsp;&nbsp;&nbsp; 2&nbsp;&nbsp;Ow&nbsp; Och  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;exclude_pairs&nbsp;            1&nbsp; Och  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cutoff       &nbsp;&nbsp; 3.20 &nbsp;&nbsp; Angstrom  
&nbsp;&nbsp;&nbsp;&end_search_environment  

***&end_reactive_species***

The ***type*** directive defines the type of reactive species, in this case a donor type, as the H<sub>3</sub>O<sup>+</sup> species donate protons to their surroundings (*note: this directive can also be defined as acceptor when dealing with systems like anion-exchange membranes, whose dynamics is driven by the diffusion of HO<sup>-</sup>*). The ***total_number*** directive sets to identify a total of 4 donor sites at each frame of the trajectory. In this example, the total number of 4 follows from the modelling of 4 SO<sub>3</sub>H groups (with tag Hs in the ***&reference_composition*** block). Consequently, we set to identify H<sub>3</sub>O<sup>+</sup> species using the ***&bonding_criteria*** sub-block, which specifies that each reference sites (O atoms in this case, see ***&search_environment***) forms 3 bonds only with H atoms, independently of theirs tag. A bond between the reference site and the H atoms is subject to a cutoff distance criterion of 1.3 Angstrom, beyond which the bond is considered to be broken.  
This information is enough to track all H<sub>3</sub>O<sup>+</sup> species. For Nafion, however, H<sub>2</sub>O and H<sub>3</sub>O<sup>+</sup> species are not part of bulk liquid water but constitute an interface with the backbone structure of the membrane. Consequently, protons could also form bonds with the oxygen atoms of the SO<sub>3</sub><sup>-</sup> groups. Such oxygens are labelled as Och in the ***&reference_composition*** block. To account for the formation of these bonds, the user must define the ***&extra_reactive_bonds*** sub-block and specify the number of type of bonds (2 in this case). Together with the atomic tags to define the bond, the user must also define the cutoff value and the units as shown. This ends the specifications for the ***&bonding_criteria*** block. To indicate that the possible tags (Ow and Och) become the reference site of a newly formed chemical species, the implemented algorithm retags the sites by adding an asterisk symbol `*` as apex, so the tags become Ow`*` and Och`*`. The same is done with the corresponding H atoms.   

So far, we have implicitly discussed about the donors (Ow and Och) based on chemical knowledge for the system under consideration. This is not enough for executing the code. In fact, we need to specify the criteria not only to identify the donor sites but also to track the transferring of atomic species to neighbouring acceptors. This information must be added using the ***&search_environment*** block. The ***include_tags*** directive sets the number of atomic tags (2, Ow and Och) that are possible candidates to become the reference site of the new reactive species. The implemented algorithm first identifies the reference sites for the donors. For the subsequent frames of the MD trajectory, ALC_ART determines the bonding pattern for all those Ow and Och neighbouring sites, only within the environment region around each donor site. The environment is defined using a distance cutoff criterion, the ***cutoff*** directive, of 3.2 Angstrom (note: this cutoff is different from the bonding cutoff. A reasonable good value for this directive can be obtained from the RDF analysis). Evaluating the bonding criteria for donors and possible acceptors allows detecting (proton) transfer events and tracking changes of the chemical species. The settings for the ***include_tags*** directive also indicate that when the site is Ow`*` or Och`*`, both Ow and Och can be considered as acceptors within the distance cutoff. For this system, we know that Och* -> Och transitions are unphysical. To accelerate the search, we can (optionally) opt to discard the calculation of Och-Och pairs in the search. This is specified with the ***exclude_pairs***.  
It is important to emphasise that the setting of the ***&reactive_species*** block ONLY requires atomic labels, defined in ***&reference_composition***, and the definition of few cutoff values. This structure for the settings is rather flexible and general. Most importantly, it avoids the complexity of setting system-specific bond lengths and angle parameters, as required by other analysis tools.  
The ***&reactive_species*** block prints by default the TRACK_CHEMISTRY file with the xyz position for the donor sites along the trajectory. As explained, such donors can be interpreted as the location of the chemical species, in this case H<sub>3</sub>O<sup>+</sup> or SO<sub>3</sub>H. The population probabilities for the involved species are printed (in percentages) in the OUTPUT file. 
### Definition of nonreactive species to be monitored
The model related settings of the previous sections are already sufficient to evaluate chemistry changes along the trajectory.  However, it is also important to account for species that are not the chemical reactive, but part of the environment. As shown later, such identification allows  computing Mean Square Displacement (MSD) and Orientational Correlation Function (OCF) analyses for these non-reactive part of the system. We shall define and refer to such species as "nonreactive species". In the present example, nonreactive species correspond to the nano-confined water molecules. The following block defines all the directives related to nonreactive species definition, whici is trivial in this case (water), but illustrates the flexibility of the code if more complex nonreactive species need to be defined:  

***&selected_nonreactive_species***  
&nbsp;&nbsp;&nbsp;name H2O  
&nbsp;&nbsp;&nbsp;reference_tag Ow  
&nbsp;&nbsp;&nbsp; ***&atomic_components***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number_components 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;H &nbsp;&nbsp; 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;O &nbsp;&nbsp; 1  
&nbsp;&nbsp;&nbsp;***&end_atomic_components***  
&nbsp;&nbsp;&nbsp;bond_cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  1.2 &nbsp;&nbsp;Angstrom  
***&end_selected_nonreactive_species***

In this block we name the species as" H2xO" and use the Ow as the reference atomic tag. The atomic composition is set using the ***&atomic_components*** block. The cutoff sets that the maximum bonding distance criterion between the O and the H atoms for the set to be considered as a "nonreactive species" to be monitored alomng the trajectorye. 

## Settings for trajectory analysis
Having defined the model related settings, the next step is to specify trajectory related directives and blocks for MD analysis. The first two compulsory directives must be the ensemble and timestep. It is important to clarify that this timestep is not the time step used for the numerical integration of the equation of motion, but the time step used to record the trajectory. For this example we have:  

***ensemble***&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;       NVE  
***recorded_timestep***&nbsp;&nbsp;&nbsp;  fs  

These two settings are important to evaluate and process the data. It is user's responsibility to ensure the correctness of these directives. Options for ***ensemble*** are NVE, NVT and NPT. Analysis for NPT ensembles is NOT compatible for "xyz" formats because the simulation cell is NOT fixed along the trajectory. In contrast, NPT is a valid option for trajectories in "vasp" format.  
For reactive systems and only for debugging purposes, the user can request printing a trajectory that contains the labels used to identify the reactive species:
  
***print_retagged_trajectory***&nbsp;&nbsp;&nbsp; .True. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (only if ***reactive_chemistry***&nbsp;&nbsp;&nbsp;.True.)

This setting prints the relabelled trajectory to the TAGGED_TRAJECTORY file. 
  
***print_track_chemistry***&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .True.

By default, the TRACK_CHEMISTRY file is generated for reactive systems, so the above line is redundant. Clearly, for those cases where several runs are needed (keeping the same values of "start_time" and "end_time" in the ***&segment_trajectory*** below), it is not convenient to print the same TRACK_CHEMISTRY file repetitively. After the first run, the user can prevent the printing of the TRACK_CHEMISTRY file by setting the ***print_track_chemistry*** directive to .False.  

### Definition of time segments for data analysis
In order to make the most of the collected statistics, the computation of TCF, SPCF, OCF, OCF_REACTIVE and MSD (see below) can be optimised by dividing of the whole trajectory in time segments. In ALC_ART, this is done through the definition of the ***&segment_trajectory*** block, which in this case reads:

***&segment_trajectory***  
&nbsp;&nbsp;&nbsp;segment_time&nbsp;&nbsp;   10.0&nbsp;&nbsp;  ps  
&nbsp;&nbsp;&nbsp;start_time&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  2.00&nbsp;&nbsp;ps  
&nbsp;&nbsp;&nbsp;restart_every &nbsp;&nbsp;   0.50&nbsp; ps  
&nbsp;&nbsp;&nbsp;end_time &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  25.0&nbsp; ps  
***&end_segment_trajectory***  

The ***segment_time*** directive specifies the duration of the time segment used to compute the properties. In this case, the time segment is set to 10 ps, and this duratio is used for the computation of all time-related quantities. It might well happen that the computation of quantities require different time segments, for which the user should compute each quantity separately in different runs. The user can delay the start of the analysis by using the ***start_time*** directive, which ignores the first part of the trajectory (2 ps here).  
Generally, the size for systems compatible with standard DFT simulations are often not large enough statistically-wise, while the length of computed trajectories are limited to tens of picoseconds. Consequently, computed properties are subject to significantly large errors. To reduce this effect, we follow the strategy of S. Kim et al. [\[1\]](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.202000498) and use multiple time starting points, which are separated by the settings of the ***restart_every*** directive (0.50 ps). In this example the analysis for the first segment starts at 2 ps and it lasts for 10 ps (up to the first 12 ps of the trajectory). At 2.50 ps, a new analysis starts, which also lasts for 10 ps and extends to the 12.50 ps of the trajectory.  Likewise, a third analysis starts at 3.0 ps and finishes at 13 ps. The process is repeated until the starting time of the last segment, which is the total time recorded for the trajectory minus 10 ps. One can also truncate the length of the analysis with the ***end_time*** directive (here 25 ps, of the total 30 ps, will be considered for the analysis). At the end of the cycle, we will have multiple analyses, each of 10 ps length. ALC_ART will compute the average quantity (AVG) using the multiple analyses, together with the standard deviation (STD) along the selected time segment. Results are printed to files TCF_AVG, SPCF_AVG, OCF_NONREACTIVE_AVG, OCF_REACTIVE_AVG and MSD_AVG, depending on which block are defined (see below). Directives ***start_time*** and ***end_time*** are optional. If ***restart_every*** is ommited, time segments follow one after another which, in this example, is equivalent to setting the ***restart_every*** directive to 10 ps.  

**IMPORTANT**: by default, the values computed for all the multiple time segments are not printed. In case the user wants to print the information for all time segments, the ***print_all_segments*** directive must be set to .True. inside the relevant block. Results will be printed to files whose names contain the relevant quantity with the "_ALL" appex. WE ADVISE ITS USE ONLY FOR DEBUGGING PURPOSES

## <span style="color: black">Analysis for the reactive part of the system
Having defined the composition of the reactive part of the system with the ***reactive_species*** block, the code allows computing transfer and orientational related quantities involving such reactive species. To this purpose, one needs to define the ***&reactive_analysis*** block which contains sub-blocks with the settings of the quantities to be computed. For the present case we have:

***&reactive_analysis***  
&nbsp;&nbsp;&nbsp;# Transfer Correlation Function (TCF)  
&nbsp;&nbsp;&nbsp;***&tcf***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;method &nbsp;&nbsp; HiCF  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print_all_segments&nbsp;&nbsp;  .True.  
&nbsp;&nbsp;&nbsp;***&end_tcf***  

&nbsp;&nbsp;&nbsp;# Special Pair Correlation Function (TCF)  
&nbsp;&nbsp;&nbsp;***&spcf***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;method &nbsp;&nbsp; HiCF  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print_all_segments&nbsp;&nbsp;  .True.  
&nbsp;&nbsp;&nbsp;***&end_spcf***  

&nbsp;&nbsp;&nbsp;# Residence times  
&nbsp;&nbsp;&nbsp;***&residence_times***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rattling_wait&nbsp;&nbsp;  0.2&nbsp;&nbsp; ps  
&nbsp;&nbsp;&nbsp;***&end_residence_times***  

&nbsp;&nbsp;&nbsp;# Orientational correlation function for the reactive species  
&nbsp;&nbsp;&nbsp;***&ocf_reactive***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;legendre_order&nbsp;&nbsp;  2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;u_definition &nbsp;&nbsp;&nbsp;&nbsp; special_pair  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print_all_segments&nbsp;&nbsp;  .True.  
&nbsp;&nbsp;&nbsp;***&end_ocf_reactive***  

***&end_reactive_analysis***  

The functionalities that correspond to each block are detailed separately in the following sections

### Transfer Correlation Function (TCF)  
TCFs can be computed by defining the ***&tcf*** sub-block. The ***method*** directive sets the approximation to compute the TCF. Both History-Independent and History-Dependent Correlation Functions (HiCF and HDCF, respectively) have been implemented, following the work of M. Tuckerman et al. [\[2\]](https://aip.scitation.org/doi/10.1063/1.3474625) and T. Zelovich et al. [\[3\]](https://pubs.acs.org/doi/10.1021/acs.jpcc.8b10298). This block generates the TCF_AVG file with the computed average (AVG) and standard deviation (STD) as a function of time, in this case 10 ps long, according to the specification of the ***&segment_trajectory*** block. By setting the ***print_all_segments*** directive to .True. inside this block, the computed correlation for all the 10 ps segments will be printed to the TCF_ALL.  

### Special Pair Correlation Function (SPCF)  
SPCFs can be computed by defining the ***&spcf*** sub-block. The ***method*** directive sets the approximation to compute the SPCF using the special pair. Both History-Independent and History-Dependent Correlation Functions (HiCF and HDCF, respectively) have been implemented. This block generates the SPC_AVG file with the computed average (AVG) and standard deviation (STD) as a function of time, in this case 10 ps long, according to the specification of the ***&segment_trajectory*** block. By setting the ***print_all_segments*** directive to .True. inside this block, the computed correlation for all the 10 ps segments will be printed to the SPCF_ALL.  

### Residence times  
The definition of the ***&residence_times*** sub-block computes the residence times for each species before being transferred. It is important the user understands that residence times and TCF (and SPCF) are not the same concept. To account for the problem of the "rattling", the user can set the ***rattling_wait*** directive, which handles the analysis as follows: once the transfer of the atomic species (in this case proton) from donor to acceptor has occurred, the algorithm waits/holds for some time (0.2 ps in this case) to accept the transfer. If the proton returns back to its donor, the transfer is discarded. Results are printed to the RES_TIMES file for each of the species (either Ow`*` or Och`*`). If the ***rattling_wait*** directive is not declared, rattling effects are included in the analysis.

### Orientational chemistry 
For reactive systems it is also possible to compute the orientional anisotropy related to the reactive species. Such information can be related to measured anisotropies of proton complexes. We refer to  Refs. [\[4\]](https://doi.org/10.1063/1.5108907) and [\[5\]](https://doi.org/10.1021/jacs.1c08552) for a detailed explanation of these quantities. The orientational correlation function for the changing chemical species is computed using the ***&ocf_reactive*** sub-block.  

The ***u_definition*** directive specifies the definition used to define orientation of each reactive species, which is defined the position of the reference tag and its neighbours. In this case, we instruct the code to use the "special pair", defined as the closest acceptor(donor) to the donor(acceptor) taking the reference atomic species (either Ow\* or Och\*) and the closest oxygen neighbour (either Ow or Och). The other implemented option for ***u_definition*** is "unrattled_special_pair", defined as the pair between the reference atom of the reactive species and the neighbouring oxygen, while neglecting any rattling. At this stage of development, however, we only recommend the use of "special_pair" option.      
The correlation is computed using the order of the Legendre polynomial, defined by the ***legendre_order*** directive, equal to 2 in this case.  

This sub-block generates the OCF_REACTIVE_AVG file with the computed average (AVG) and standard deviation (STD) as a function of time, in this case 10 ps long, according to the specification of the ***&segment_trajectory*** block. The user could set the ***print_all_segments*** directive to .True. inside this block, and the computed correlation for all the 10 ps segments will be printed to the OCF_REACTIVE_ALL file.  

## <span style="color: black">Analysis for the selected, nonreactive part of the system
The ***&selected_nonreactive_species*** block above defines the nonreactive species that will be monitored along the trajectory. To compute quantities involving these nonreactive species the user must define the ***&nonreactive_analysis*** block. In this example we have,  

***&nonreactive_analysis***  
&nbsp;&nbsp;&nbsp;compute_amount&nbsp;&nbsp;  .True.  

&nbsp;&nbsp;&nbsp;***&intramol_statistics***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&distance_parameters  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;species&nbsp;&nbsp;&nbsp;&nbsp;  H O  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;delta&nbsp;&nbsp;&nbsp;&nbsp;    0.005  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lower_bound&nbsp;&nbsp;&nbsp;&nbsp;    0.8 Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;upper_bound&nbsp;&nbsp;&nbsp;&nbsp;    1.4 Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&end_distance_parameters  
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&angle_parameters  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;species&nbsp;&nbsp;&nbsp;&nbsp;    H O H  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lower_bound&nbsp;&nbsp;&nbsp;&nbsp;     80   degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;upper_bound&nbsp;&nbsp;&nbsp;&nbsp;     130  degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;delta&nbsp;&nbsp;&nbsp;&nbsp;       0.5    degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&end_angle_parameters  
&nbsp;&nbsp;&nbsp;***&end_intramol_statistics***  

&nbsp;&nbsp;&nbsp;***&intermol_statistics***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;only_ref_tags_as_nn     .True.  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  &distance_parameters  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    delta    0.01  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    lower_bound    2.3  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    upper_bound    4.0 Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&end_distance_parameters  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&angle_parameters  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    lower_bound    30   degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    upper_bound    180  degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    delta          2    degrees  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&end_angle_parameters  
&nbsp;&nbsp;&nbsp;***&end_intermol_statistics***  

&nbsp;&nbsp;&nbsp;***&msd***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;select &nbsp;&nbsp;xy  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pbc_xyz&nbsp;&nbsp;   T&nbsp;&nbsp;  T&nbsp;&nbsp;  T  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print_all_segments&nbsp;&nbsp; .True.  
&nbsp;&nbsp;&nbsp;***&end_msd***  

&nbsp;&nbsp;&nbsp;***&ocf_nonreactive***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;legendre_order&nbsp;&nbsp;  2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;u_definition&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    bond_12-13  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print_all_segments&nbsp;&nbsp; .True.  
***&end_ocf_nonreactive***  

***&end_nonreactive_analysis***  

The functionalities defined in the block above are detailed separately in the following sections.

### Average amount of selected, nonreactive species 
When setting the optional ***compute_amount***  directive to .True. (.False. by default), the code will compute the average number (and the standard deviation) of nonreactive species along the trajectory, and the result is printed to the OUTPUT file. This functionality is particularly useful for the analysis of reactive systems, where the number of non-reactive species is expected to change along the trajectory.  

### Intra and Intermolecular statistics for distance and angles 
ALC_ART also offers the possibility to compute probability densities for ***intra***molecular and ***inter***molecular distances and angles involving the nonreactive species, for which the user must define the ***&intramol_statistics*** and the ***&intermol_statistics*** sub-blocks, respectively. Likewise, information of distances and angles are derived from the definition of the &distance_parameters and &angle_parameters, respectively, as shown above.  

***&intramol_statistics***  
For intramolecular information, the user must define the &intramol_statistics*** block. 
For distances, the user must define the &distance_parameters with the pair of atomic elements involved, which must agree with the definitions of the  ***&atomic_components*** block. To set the range of distances to be considered the user must also define the lower and upper bounds as well as the delta directive, which corresponds to the distance discretization. Allowed units are Bohr and Angstrom. The definition of this block generates the INTRAMOL_DISTANCES file.  
For angles, the &distance_parameters sub-block require three atomic elements, which must also be in agreement with the definition of  ***&atomic_components***. The order of the elements is important. In this case, only the H-O-H angle is computed, which is the internal angle for water. To set the range of angles to be considered the user must also define ***lower_bound*** and ***upper_bounds*** as well as the delta directive, which corresponds to the discretization of the selected range. The definition of this block generates the INTRAMOL_ANGLES file.  

***&intermol_statistics***  
The definition of this block identifies the first and the second nearest nonreactive species around each nonreactive species. In contrast to the ***&intramol_statistics*** block, the user must not define the involved elements, as the algorithm will use the reference atomic species as defined in the ***reference_tag*** directive. Instead, the ***only_ref_tags_as_nn*** directive can be defined: if set to .True., the statistics will only include those species whose reference tag is exactly the same as the tag defined by the ***reference_tag*** directive; if set to .False., all nonreactive species will be considered for the analysis, even if they change their chemistry. The information provided is used to computed the probability density for of the first and second nearest nonreactive species if the &distance_parameters sub-block is defined, which will generate the INTERMOL_DISTANCES_NN1 and INTERMOL_DISTANCES_NN2 files. The description for the required directives of the &distance_parameters sub-block is the same as the  &distance_parameters sub-block within ***&intramol_statistics***. 
If the &angle_parameters sub-block is defined, ALC_ART computes the angle that each nonreactive species form with the first and second nearest nonreactive species. This information is used to obtain the probability density of this angle, considering all the nonreactive species of the system along trajectory. The description for the required directives of the &angle_parameters sub-block is the same as the  &distance_parameters sub-block within ***&intramol_statistics***. The result is printed to the INTERMOL_ANGLES_NN file.  

### Mean Square Displacement (MSD)
The computation of MSDs is a common feature in most available software for MD data analysis. In ALC_ART, the MSD is only computed for the nonreactive species, for which the definition of the ***&selected_nonreactive_species*** is compulsory. Inside the ***&msd*** sub-block, the ***select*** directive instructs the code to consider selected components of the atomic positions to compute the MSD. For those users that are not familiar with the formula to compute the MSD, we refer to the Supporting Information of Ref. [\[6\]](https://pubs.acs.org/doi/10.1021/acs.jpclett.1c04071?ref=PDF). In this example, the option "xy" computes the MSD only in the plane parallel to the Nafion membrane. Available options for ***select*** are: x, y, z, xy, xz, yz and xyz.  
As an optional directive, ***pbc_xyz*** allows including (or not) the effect of periodic boundary conditions (PBCs) for each Cartesian coordinate. PBCs are set by default.  
Activation of the ***&msd*** block generates the MSD_AVG file with the computed average (AVG) and standard deviation (STD) as a function of time, in this case 10 ps long, according to the specification of the ***&segment_trajectory*** block. The user could optionally set the ***print_all_segments*** directive to .True., which will print the results for all the 10 ps segments to the MSD_ALL file.

### Orientational Correlation Function (OCF)
To compute the OCF of the nonreactive species, the ***&ocf_nonreactive*** sub-block must be defined. The nonreactive species must be a molecule (water in this case), otherwise the code will abort the execution. At time zero, the vector unit **u**(0) is computed using a chosen geometry criterion for each the nonreactive species (excluding the chemical species). Implemented options for geometry criterion are discussed below. Unit vectors define the orientations of the molecules. Likewise, at trajectory time **t** the unit vector **u**(t) is also computed for each nonreactive species. The OCT at time **t** is computed as `<`P<sub>l</sub> [**u**(t)**.u**(0)]`>`, where l is the order of the Legendre polynomial, the dot is the inner product and the average < > is done over the number of nonreactive species [\[7\]](https://pubs.acs.org/doi/10.1021/acs.jpcb.5b02936). A remarkable feature of the implementation is the possibility to compute OCF for reactive systems. The implemented algorithm identifies all the nonreactive species at time zero. If one of the species become a reactive site (from having accepted an atomic species), such site (formerly Ow and now Ow* )  is discarded. Likewise, if the same site donates the species and becomes a nonreactive species again, it will be discarded for the whole time segment, as the correlation between **u**(t) and **u**(0) is already lost when the nonreactive species became an reactive site. 
The order of the polynomial "l" is set with the ***legendre_order*** directive. Although the OCF analysis is often performed by setting l=2, the user is free to select orders between 1 and 4.
To define the vector unit **u**, different geometrical criteria are implemented. The **u_definition** directive can be set to:

* bond_12: the unit vector between atoms 1 and 2 of the nonreactive species. This is the only option for diatomic molecules.
* bond_13: the unit vector between atoms 1 and 3 of the nonreactive species. 
* plane: the unit vector from the cross product between vector_12 and vector_13.
* bond_123: the unit vector from the **sum** of vector_12 and vector_13. 
* bond_12-13: the unit vector_12 and vector_13 are evaluated together within the average. 

It is users responsibility to test all these possible settings for the interpretation of the computed OCF. IMPORTANT: For nonreactive species with more than 2 atoms we do not recommend the use of options ***bond_12***, ***bond_13*** and ***plane***. 

Activation of the ***&ocf_nonreactive*** block generates the OCF_NONREACTIVE_AVG file with the computed average (AVG) and standard deviation (STD) as a function of time, in this case 10 ps long, according to the specification of the ***&segment_trajectory*** block. By setting the ***print_all_segments*** directive to .True. inside this block, the computed correlation for all the 10 ps segments will be printed to the OCF_NONREACTIVE_ALL file.

## Functionalities for general analysis
Regardless the system is reactive or not, the code also is capable of computing additional, general quantities from the MD trajectory. To this purpose, the block ***&general_analysis*** must be defined. In this example, we have:  

***&general_analysis***  
&nbsp;&nbsp;&nbsp;***&rdf***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dr &nbsp;&nbsp;0.02&nbsp;&nbsp;  Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;tags_species_a&nbsp;&nbsp;  1&nbsp;&nbsp;  Ow  
&nbsp;&nbsp;&nbsp;tags_species_b&nbsp;&nbsp;  2&nbsp;&nbsp;  Ow\* &nbsp;&nbsp; Och\*  
&nbsp;&nbsp;&nbsp;***&end_rdf***  

&nbsp;&nbsp;&nbsp;***&track_unchanged_chemistry***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; number &nbsp;&nbsp;&nbsp;4  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; tag  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sch  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; list_indexes &nbsp;&nbsp;&nbsp; 133&nbsp;&nbsp;&nbsp;  134&nbsp;&nbsp;&nbsp; 135 &nbsp;&nbsp;&nbsp; 136  
&nbsp;&nbsp;&nbsp;***&end_track_unchanged_chemistry***  

&nbsp;&nbsp;&nbsp;***&coord_distrib***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; species &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Ow  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; coordinate  &nbsp;&nbsp;&nbsp;z  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; delta  &nbsp;&nbsp;0.1&nbsp;&nbsp;  Angstrom  
&nbsp;&nbsp;&nbsp;***&end_coord_distrib***  

&nbsp;&nbsp;&nbsp;***&selected_nn_distances***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;reference_species&nbsp;&nbsp;&nbsp;&nbsp;  Ow  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;nn_species&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  4  Ow Ow\*  Och Och\*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lower_bound&nbsp;&nbsp;&nbsp;&nbsp;    2.3 Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;upper_bound&nbsp;&nbsp;&nbsp;&nbsp;    3.0 Angstrom  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dr&nbsp;&nbsp;&nbsp;&nbsp;    0.01  Angstrom  
&nbsp;&nbsp;&nbsp;***&end_selected_nn_distances***  

***&end_general_analysis***  

### Radial Distribution Function (RDF)
This capability is a common feature in the majority of available software for MD data analysis. Inside the ***&rdf*** sub-block, the user must define the ***tags_species_a*** and ***tags_species_b*** directives, with the amount for species for each type. The chemical elements defined in type **a** (and **b**) can be different. In addition, the user needs to specify the discretization for the distance, ***dr***, set to 0.02 Angstrom in this case.

In this block, ALC_ART is instructed to compute the RDF between the water oxygens (Ow) and all the oxygen sites of the chemically formed species (Ow\* and Och\*). Note the values of 1 and 2, preceding the definition of the tags for type a and b, define the amount of atomic tags to be considered. Results are printed to the RDF file.

### Tracking non-reactive species (unchanged chemistry) 
In addition to the tracking of the changing chemical species, the user can also choose to track atomic sites that do not change their chemistry. This can be important to compare how the location of reactive and non-reactive sites are distributed along the trajectory. The user must define the ***&track_unchanged_chemistry*** sub-block. The first directive inside this sub-block must be ***number***, which indicates how many sites the user wants to track. A maximum of 10 sites are allowed. The directive ***tag*** specifies the atomic species to be tracked. Finally, the user must define the 4 atomic indexes with the ***indexes*** directive. In case that the declared indexes do not correspond to the defined ***tag***, the code will abort the execution. The positions of the tracked indexes are printed to the UNCHANGED_CHEMISTRY file.  
In the hypothetical scenario the users aims to track Ow atoms for this example, the code would track the Ow species as long as they do not change their chemistry. If the chemistry changes, results would be printed up to the frame where the  change has occurred.

### Spatial probability distribution of selected species 
For anisotropic models such as membranes or layered materials it might be convenient to compute the probability distribution of selected species along a particular coordinate (x, y or z). This information can be important to identify the role of confinement for example. To compute this quantity, the user must define the ***&coord_distrib*** sub-block. In this example, ALC_ART will compute the probability distribution of all Ow species (directive ***species***) along the z coordinate (directive ***coordinate***), which is the coordinate perpendicular to the Nafion backbone membrane. The discretization for the z coordinate is defined with the ***delta*** directive, set to 0.1 Angstrom in this case. Results are printed to the COORD_DISTRIBUTION file. 

### Computing the probability distribution for the shortest distance between defined species 
To withdraw further information about the interactions between species, ALC_ART offers the possibility to compute the probability distribution for the distance between reference species and nearest neighbour species selected by the user. For further fundamental detail about this quantity, we refer to Ref. [\[8\]](https://pubs.acs.org/doi/10.1021/acsaenm.4c00836).
To compute this quantity, the user must define the ***&selected_nn_distances*** sub-block. Inside this sub-block the user must define the ***reference_species*** using the tag for the species under consideration. In this example, the reference tag is set to Ow. The ***nn_species*** directive refer for all the species for which the algorithm will search as nearest neighbours to each ***reference_species***. One must first specify the number of species and the list of atomic tags to be considered. We remind that the asterix is used to identify a relevant atomic tag that has become part of a new chemical species from the definition of the ***&reactive_species*** block. Here the tag Ow* corresponds to the hydronium oxygen. The specification for the rest for the remaining directives define the distance range and the discretization to generate the probability distribution. Results are printed to the SELECTED_NN_DISTANCES file.

## Constraining the analysis to a selected region within the system  
By default the trajectory analysis is carried out for the whole simulation cell. ALC_ART also offers the possibility to compute quantities for a particular region. Only rectangular regions within the simulation cell are allowed. To defined a region within the modelled system, the user must set the ***&region*** block:  

***&region***  
&nbsp;&nbsp;&nbsp;Delta_x &nbsp;&nbsp;&nbsp;&nbsp;        -1.6 &nbsp;&nbsp; 11.6 &nbsp;&nbsp; inside  
&nbsp;&nbsp;&nbsp;Delta_y &nbsp;&nbsp;&nbsp;&nbsp;        -1.0 &nbsp;&nbsp; 14.0&nbsp;&nbsp;  inside  
&nbsp;&nbsp;&nbsp;Delta_z &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   8.0 &nbsp;&nbsp; 12.0&nbsp;&nbsp;  outside  
***&end_region***  

The choice of the settings for ***Delta_x***, ***Delta_y*** and ***Delta_x*** is subject to the dimensions of the simulation cell. The first and second values are the minimum and maximum values for the domain along the corresponding coordinate in units of Angstrom. The third argument indicates if the region of interest for analysis are within (inside) or outside the given range. By comparison with the definition of the ***simulation_cell*** block above, we realise that ***Delta_x*** and ***Delta_y*** are redundant for this case, and they can be omitted. In fact, when the **Delta** setting is omitted, ALC_ART will use the whole spatial domain of the simulation cell for that coordinate. In this example, the analysis is focused on the interface regions (up and down) near the SO<sub>3</sub><sup>-</sup> groups. The central region between 8 and 12 Angstrom is not considered. If any of the defined settings is not consistent with the simulation cell (for example, replacing inside by outside in the specification of the ***Delta_x*** directive), the code will abort the execution. The user can double check the definition of the region in the generated OUTPUT file. Finally, multiple definitions for ***Delta_x***, ***Delta_y*** and ***Delta_z*** are allowed.  

IMPORTANT: the ***&region*** block affects the computation of:  

* **TCF**: Transfer Correlation Function  
* **SPCF**: Special Pair Correlation Function  
* **OCF**: Orientational Correlation Function  
* **OCF_REACTIVE**: OCF for the changing chemical species  
* **RDF**: Radial Distribution Function  
* **MSD**: Mean Square Displacement  
* **selected_nn_distances**: probability distribution for the shortest distance between selected species  
* **Intramolecular parameters**:  probability distribution for **intra**molecular angles and distances 
* **Intermolecular parameters**:  probability distribution for **inter**molecular angles and distances  

It is user's responsability to set a correct and meaningful region where to constrain the analysis. Of course, this depends on the type of analysis and the system under consideration.
