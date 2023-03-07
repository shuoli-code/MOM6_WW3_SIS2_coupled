# WW3-MOM6-SIS2 coupled model

The coupled model WW3-MOM6-SIS2 which is adapted from 
Geophysical Fluid Dynamics Laboratory (GFDL) global ocean and sea ice model OM4.0(https://github.com/breichl/FMS_Wave_Coupling/tree/Dec2020)
The wave model WaveWatch 3 is coupled to MOM6-SIS2 considering non-breaking 
wave-induced mixing at ocean surface (Babanin2006) and wave generated momentum 
flux to ice (Boutin2020). 

### Clone this repository, update submodules.

To clone this repository from github use the following:

> Using https:
> git clone https://github.com/shuoli-code/MOM6_WW3_SIS2_coupled.git

Download and update the submodules:

> cd MOM6_WW3_SIS2_coupled
> git submodule update --init --recursive

Now you should be ready to compile.


### Compiling on Gaea

1. First step is to set up the wave model.  The wave model source code must be pre-processed using WAVEWATCH provided programs to convert to standard FORTRAN 90 files.  A sample script to complete this step is provided in "tools/SetUpWW3.csh", which works on Gaea and GFDL workstations to use a particular switch file and extract FORTRAN 90 code from the WW3 ftn code.  This script sets up code to compile WW3 preprocessor routines for building the grid binary (ww3_grid), the forcing binaries (ww3_prnc, needed for standalone WW3), and the stand alone WW3 driver (ww3_multi).  It also sets up code to compile WW3 postprocessor routines for converting the output binary into NetCDF (ww3_ounf).  Note that the wave model needs to know a valid compiler to unpack its "ftn" files into "f90" files, but you shouldn't need to use the same fortran compiler that you will use to create executables.

> NCI/GADI Instructions:
>
> cd MOM6_WW3_SIS2_coupled
> ./Set_Up_WW3.csh

2. The second step is to compile.  Again, a script to do this is provided.  This script will compile 1) FMS library, 2) ww3_grid, 3) ww3_prnc, 4) ww3_multi, 5) WW3 library (for linking within coupled model), 6) the coupled model , and 7) ww3_ounf.

> Gaea Instructions:
>
> cd MOM6_WW3_SIS2_coupled
> ./Wave_Compile.csh

If working on NCI gadi HPC system, these steps should successfully compile libraries and executables needed to set-up and run the WW3 coupled system with FMS.

### Running - OM4_025.JRA example

1. First follow the instructions to download the MOM6-examples input data (https://github.com/NOAA-GFDL/MOM6-examples/wiki/Getting-started#downloading-input-data).  Link this directory into the main directory for this repository as ".datasets", exactly as you would in MOM6-examples to use those test cases.  On gaea we simply execute "ln -sf /lustre/f2/pdata/gfdl/gfdl_O/datasets .datasets".  You would replace the source file location with the location you have put the datasets file you download.

2. Change to the OM4_025.JRA test case directory

> cd examples/OM4_025.JRA

2. We next have to set-up the grid.

a. (Optional) If you are on a system that has access to the FMS/FRE tools, you can build your own grids and exchange grids (these tools are available on NOAA/GFDL github: FRE-NCtools).  See the GRID directory for some ideas of how these grids can be created to work with the wave model.  You would need to run the BuildExchangeGrid.csh and WaveMosaics.py (cshell and python) scripts both to create all necesssary files (there are some pseudo-hacks needed to set-up the grid_spec.nc and the wave related exchange grids, if you are having trouble getting the wind into the wave model this grid_spec step is critical).  But, this is all optional, you can simply run with the example here to have a 0p25 degree resolution case on the MOM6 native 1/4 degree grid.  

b. You will need to update the WW3 grid files in WW3/PreProc, see WW3/PreProc/GenGrid.py for an example (this script should create the basic files for you without modification). This generates the WW3 Depth grid for a curvilinear tripolar grid on the MOM6 native 1/4 degree grid. Note that on gadi, you need netCDF4 in python3.

> cd WW3/PreProc  
> python3 GenGrid.py

c. Next you need to create the mod_def.ww3 file for WW3's internal grid.  Navigate to WW3/PreProc and execute the ww3_grid (note on GFDL/Gaea you need to load the NetCDF libraries used to compile, e.g., on GFDL system: "module load netcdf/4.2" and on Gaea "module load cray-netcdf" before this will work):

> cd WW3/PreProc  
> ../../../../build/intel/wave_ice_ocean/ww3_grid/ww3_grid

2. Don't forget to create a RESTART directory within the OM4_025.JRA example directory:

> mkdir RESTART

4.  Now we should be ready to run.  To run on Gaea we can either submit a job through the batch system, or run the job interactively.  In this example we are going to run it interactively.

> Gaea Instructions:
>
> cd examples/OM4_025.JRA  
> salloc --clusters=c3 --qos=normal --nodes=10 --x11  
> srun -n320 ../../build/intel/wave_ice_ocean/repro/MOM6

On the GFDL workstation I am using custom installed MPI libraries.  This will have to be reproduced for your own set-up.

From the examples/Idealized_Hurricane/0p5 directory:
>  /net2/bgr/Software/build/bin/mpirun -n 4 ../../../build/intel/wave_ice_ocean/repro/MOM6

On NCI/gadi, you can run with a simple script which is provided in examples/OM4_025.JRA/job.sh

5.  If this works, then congratulations, you have successful set-up, compiled, and run this example.  There is much more we can do with this including customizing set-ups, and processing and manipulating output.  Adding more instructions to this README.md is always welcome!

6.  Basic instructions for creating WW3 output files:  MOM6 creates output in an easy to use netCDF format, but WW3 output is stored as binary files which need a further processing step to convert to netCDF.  You should have all the tools to do this already in place.  All you need to do is navigate to WW3/PostProc, and find the build/intel/wave_ice_ocean/ww3_ounf/ww3_ounf executable.  If you run it in this place, the proper files should already be linked to create netCDF from the mod_def and out_grd files.  Note that the ww3_ounf.inp can be edited (along with ww3_multi.inp) to specify diagnostics to save and frequency.  More information can be found in the WW3 directory. On gadi, a simple script is provided in examples/OM4_025.JRA/WW3/PostProc/output.sh
