## About the code
**ALC_ART** (**A**da **L**ovelace **C**entre **A**nalysis of **R**eactive **T**rajectories) offers an open-source tool to extract orientational anisotropies and transfer correlations from Molecular Dynamics (MD) simulations of reactive systems. By "reactive" we refer to systems in gas and condensed phase where the constituent atomic species change their chemistry composition along the trajectory, forming and breaking bonds, as it occurs in anion [\[1\]](https://pubs.acs.org/doi/10.1021/acs.jpcc.8b10298) and proton exchange membranes [\[2\]](https://pubs.acs.org/doi/10.1021/acs.jpclett.1c04071?ref=PDF), for example. The implemented capabilities allow computing: 

* the location of the changing chemical species along the trajectory,
* residence times, 
* Transfer Correlation Functions (TCFs),
* Special Pair Correlation Functions (SPCFs),
* Orientational Correlation Functions (OCFs),
* Radial Distribution Functions (RDFs),  
* Mean Square Displacements (MSDs)  

***all at once***, thus offering a novel platform to analyse reactive systems where the changing chemical species play a decisive role. During the development of the code, our efforts were focused on designing a simple and flexible input structure to facilitate the definition and tracking of changing chemical species for different systems and chemical environments. Although the applicability of **ALC_ART** is general, its development and optimization have been focused to analyse systems within the size range of Density Functional Theory (DFT) simulations. It is also important to remark that **ALC_ART** also allows the computation of OCFs, RDFs and MSDs of non-reactive systems.  

**ALC_ART** is a serial code written in Fortran, and the result of a collaboration with the [**CLF-ULTRA**](https://www.clf.stfc.ac.uk/Pages/Ultra_Facility.aspx) at the [**Science and Technology Facilities Council (STFC)**](https://www.ukri.org/councils/stfc/) funded by the [**ALC**](https://www.adalovelacecentre.ac.uk/) of the [**Scientific Computing Department (SCD)**](https://www.sc.stfc.ac.uk/). Its structure for development (and maintenance) follows the Continuous Integration (CI) practice and it is integrated within the GitLab DevOps of the STFC. In the root folder, the user will find several Markdown files, which are intended to provide help with the compilation and execution as well as guidance with the multiple available functionalities.  

## Disclaimer
The ALC does not fully guarantee the code is free of errors and assumes no legal responsibility for any incorrect outcome or loss of data.

## Contributors
**Original author:** Ivan Scivetti (SCD, STFC)  

## Structure of files and folders
ALC_ART contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains the tests files (in .tar format) needed for CI purposes. The user should execute the available scripts of the [***tools***](./tools) folder to run the test automatically and verify the code has been installed properly (see the [build_code.md](./build_code.md) file for instructions).
* [***examples***](./examples): example cases to help the user to become familiarised with the code. The SETTINGS files are described in detail.  
* [***scripts***](./scripts): contains scripts for data processing.
* [***source***](./source): contains the source code. Files have the *.F90* extension
* [***tools***](./tools): shell files for building, compiling and testing the code automatically.
* [.gitignore](./.gitignore): instructs Git which file to ignore.
* [CMakeList.txt](./CMakeList.txt): sets the framework for code building and testing with CMake.
* [LICENSE](./LICENSE): BSD 3-Clause License for ALC_ART. 
* README.md: this file.
* [build_code.md](./build_code.md): steps to build, compile and run tests using the CMake platform.
* [use_code.md](./use_code.md): provides instructions for use together with a detailed description of the implemented capabilities. 

## Dependencies
The user must have access to the following software (locally):

* GNU-Fortran (11.2.0) or Intel-Fortran (ifx 2023.1.0)
* Cmake (3.16)
* Make (4.2.1)
* Git (2.34.1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work.

## Getting started

### Obtaining the code
The user can clone the code locally by executing the following command with the SSH protocol
```sh
$ git clone git@github.com:stfc/alc_art.git
```
Instead, if the user wants to use the HTTPS protocol it must execute
```sh
$ git clone https://github.com/stfc/alc_art.git
```
Both ways generate the ***alc_art*** folder as the root directory. Alternatively, the code can be downloaded from any of the available assets.

### Building and testing the code with CMake
Details can be found in file [build_code.md](./build_code.md)

### Making use of the software
Once the code has been installed and tested, the user should create a folder where to run the code from. In such folder, the MD trajectory must be copied to the TRAJECTORY file. The user also needs to provide the SETTINGS file with instructions for the type of analysis to execute. Instructions of the implemented capabilities can be found in the [use_code.md](./use_code.md) file. The SETTINGS files in folder [***examples***](./examples) offer explanatory templates, which are intended to help new users in the setting of input directives for execution. In each of the directories, the user will also find the corresponding TRAJECTORY files.

## Contributing 
Contributions from STFC staff are welcome as long as they comply with the coding protocols, as described in the [coding_protocol.md](./coding_protocol.md) file. Instructions for CI practices are provided in the [CI_instructions.md](./CI_instructions.md) file. New implementation must compulsory include test(s) for CI purposes. 

Contributors should first create a fork using the [**Gitlab STFC**](https://gitlab.stfc.ac.uk/) web user-interface from the main [**ALC_ART**](https://gitlab.stfc.ac.uk/alc_clf-ultra/alc_art) repository. For instructions of how to create a fork, please refer to the following [**link**](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html#creating-a-fork). Access to the generated fork will be available from the *Project* tab of [**Gitlab STFC**](https://gitlab.stfc.ac.uk/), which will be located at the address <span style="color:blue">https://gitlab.stfc.ac.uk/user_id/alc_art</span>. In this address, *user_id* will be the user identification in [**Gitlab STFC**](https://gitlab.stfc.ac.uk/). Once the fork is created, the user can clone the *main* branch in the account *"username"* of the local machine *"wherever"* by executing the following command
```sh
username@wherever:/home/username/codes$ git clone -b main git@gitlab.stfc.ac.uk:user_id/alc_art.git alc_art
```
where ***alc_art*** is set as the root directory. This folder will have been created following the execution of the command above. Any other name can be chosen.  

Contributors should also request to become part of the ALC_ART project by contacting the individual in charge (owner) of the repository.  

## Acknowledgements
* ALC for funding.  
* Gilberto Teobaldi (SCD,STFC) and Paul Donaldson (CLF-ULTRA,STFC) for scientific discussions and support.  
* The Innovation Department of the STFC for assistance with the licensing process.  
* Lesley Mansfield for project management support.
