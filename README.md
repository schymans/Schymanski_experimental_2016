# Schymanski_leaf-scale_2016
Supporting information for article: http://www.hydrol-earth-syst-sci-discuss.net/hess-2016-363/

You can either access the files directly on the cloud.sagemath.com server [here]( https://cloud.sagemath.com/projects/e66470cf-1fa4-48bc-8a49-2513f177cf8a/files/Schymanski_leaf-scale_2016/), or follow the below instructions.

## Instructions on how to use these files
- Download all files as [zip file](https://github.com/schymans/Schymanski_leaf-scale_2016/archive/master.zip) and extract to a local directory on your computer.
- Install [sage](http://www.sagemath.org) on your local computer and start the ipython notebook to work with the files, 

or:

- Create an account at [sagemath cloud](https://cloud.sagemath.com)
- Create a new project
- Start newly created project
- Click on "Create or upload files..."
- Drag and drop the downloaded zip file onto the field "Drop files to upload"
- Click on the zip file and extract it when prompted.
- Click on "Files" at the top left of the browser window to see the files. Click on any of them to execute and edit.

or:

- Access the published worksheets directly on the cloud.sagemath.com server [here]( https://cloud.sagemath.com/projects/e66470cf-1fa4-48bc-8a49-2513f177cf8a/files/Schymanski_leaf-scale_2016/)

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
### Worksheet_setup.ipynb
Contains code to setup the other worksheets. Only relevant if you are interested in the inner workings of the code.
### leaf_enbalance_eqs.ipynb
Contains variable definitions and equations related to the numerical leaf energy balance model.
### leaf_enbalance_eqs2s.ipynb
Based on leaf_enbalance_eqs.ipynb, but contains additional variable definitions and equations related to the two-sided numerical leaf energy balance model (separate leaf temperatures on each side).
### E_PM_eqs.ipynb
Contains variable definitions and equations related to the analytical solutions of the leaf energy balance. It builds on definitions and equations provided in leaf_enbalance_eqs.ipynb
