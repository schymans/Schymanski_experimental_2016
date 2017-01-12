# Schymanski_experimental_2016

Supporting information for article: 
[**Technical note: An experimental setup to measure latent and sensible heat fluxes from (artificial) plant leaves.**](http://www.hydrol-earth-syst-sci-discuss.net/hess-2016-643)
Hydrol. Earth Syst. Sci. Discuss., doi:10.5194/hess-2016-643, 2017

Authors: Stan Schymanski (stan.schymanski@env.ethz.ch), Daniel Breitenstein, Dani Or


## Instructions on how to use these files
- Download newest version of files as [zip file](https://github.com/schymans/Schymanski_experimental_2016/archive/master.zip) and extract to a local directory on your computer.
- Install [sage](http://www.sagemath.org) on your local computer and start the ipython notebook to work with the files, 

or:

- Create an account at [sagemath cloud](https://cloud.sagemath.com)
- Create a new project
- Start newly created project
- Click on "Create or upload files..."
- Drag and drop the downloaded zip file onto the field "Drop files to upload"
- Click on the zip file and extract it when prompted.
- Click on "Files" at the top left of the browser window to see the files. Click on any of them to execute and edit.

### To reproduce computational results and plots:
Open `Worksheet_update.ipynb` in SageMath and click on `Cell` -> `Run All`. Here you can also see in which order the different worksheets need to be run, as their results are re-used by one another.

## Description of folders
### Main folder
The main folder contains the various worksheets (.ipynb files), the README.md and licence information. 
### data
The data folder contains data files accessed by the worksheets
### figures
Figures produced by the worksheets are placed in the figures worksheet. They are equivalent to those in the paper.
### temp
The temp folder contains temporary files created by the worksheets for exchange of data and variables between worksheets, such as .sage and .sobj files.



## Description of files
### .sage files
These files consist of input code extracted from the .ipynb worksheets with the same name. They are used internally to re-use variables and equations defined in various worksheets.
### .sobj files
These files contain snapshots of worksheet data for quick import into other worksheets.
### .ipynb files
These files are the actual worksheets containing code and documentation. Here the actual calculations and generation of figures is performed. 

## Description of worksheets
### Worksheet_update.ipynb
Converts all other worksheets to .sage files, so that they can be more easily imported by other worksheets. It also contains code to execute all other worksheets in the correct order to reproduce the results and figures. Just open this worksheet in SageMath and click on `Cell` -> `Run all`.
### Worksheet_setup.ipynb
Contains code to setup the other worksheets. Only relevant if you are interested in the inner workings of the code.
### leaf_enbalance_eqs.ipynb
Contains variable definitions and equations related to the numerical leaf energy balance model.
### leaf_enbalance_eqs2s.ipynb
Based on leaf_enbalance_eqs.ipynb, but contains additional variable definitions and equations related to the two-sided numerical leaf energy balance model (separate leaf temperatures on each side).
### E_PM_eqs.ipynb
Contains variable definitions and equations related to the analytical solutions of the leaf energy balance. It builds on definitions and equations provided in leaf_enbalance_eqs.ipynb
### leaf_chamber_eqs.ipynb
Equations to compute wind tunnel energy balance etc. It builds on definitions and equations provided in leaf_enbalance_eqs.ipynb
### perforated_foils_eqs.ipynb
Equations to compute stomatal conductance of perforated foils
### perforated_foils_data.ipynb
Calculation of conductance values of perforated foils based on confocal laser scanning microscopy images.
### leaf_chamber_data.ipynb
Evaluation and plotting of experimental data from leaf wind tunnel. 
