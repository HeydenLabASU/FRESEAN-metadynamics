# FREquency-SElective ANharmonic (FRESEAN) Mode Analysis of Thermally Excited Vibrations in Proteins (with coarse-graining) 
This codebase allows the user to run MD simulations of proteins and generate vibrational motions utilizing FRESEAN, a method used to extract low-frequency anharmonic vibrational modes of proteins. 
Please read through each of the following sections to understand how to use this repository. 
<details>
  
<summary> FRESEAN Installation </summary>

## Tested Dependencies
FFTW 3.3.10
GNU make 3.8.1
gcc 14.2
python 3.12
gromacs 2022.5
plumed 2.8.2

## FRESEAN Installation
Please follow the following instruction to install our suite of tools.
```
git clone https://github.com/HeydenLabASU/FRESEANCOARSE.git
cd FRESEANCOARSE
make
make install
make clean
source ~/.bashrc
```
If you have already set up GROMACS 2022.5 with Plumed 2.8.2, please proceed to [FRESEAN Toolbox](#FRESEAN-Toolbox) to get an overview of the provided tools.
</details>

<details>
  
<summary> Setting up the MD Engine </summary>


# Setting up the MD Engine

Gromacs 2022.5 is used as the MD engine and Plumed 2.8.2 is used as a plugin to run metadynamics.

## Step 1: Compiling PLUMED 2.8.2
Download PLUMED 2.8.2 from here: https://www.plumed.org/download
06.20.2023: PlUMED updated version on downloads section of website to 2.8.3
```
interactive
tar xfz plumed-2.8.2.tgz
cd plumed-2.8.2
./configure --prefix=$HOME/plumed-2.8.2
make -j 4
make install
```

Make sure that these paths are included in your `.bashrc` file otherwise `plumed` won't be found.
## Step 1B: BASHRC FILE FOR SOURCING
```
export PATH=$PATH:$HOME/plumed-2.8.2/bin
export PLUMED_VIMPATH=$PLUMED_VIMPATH:$HOME/plumed-2.8.2/lib/plumed/vim
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HOME/plumed-2.8.2/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/plumed-2.8.2/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/plumed-2.8.2/lib/pkgconfig
export PLUMED_KERNEL="$HOME/plumed-2.8.2/lib/libplumedKernel.so"
```

Once plumed is installed, gromacs can be installed with the plumed patch. If patching returns `--runtime not found`, make sure that there are two dashes in front of `runtime`.
## Step 2: Patching PLUMED 2.8.2 in Gromacs 2022.5
Download GROMACS 2022.5 from here: https://manual.gromacs.org/documentation/2022.5/download.html
```
cd ..
tar xfz gromacs-2022.5.tar.gz
mv gromacs-2022.5 gromacs-2022.5-plumed-2.8.2
cd gromacs-2022.5-plumed-2.8.2
plumed patch -p --runtime
cd ..
```

## Step 3: Compiling Gromacs 2022.5 with PLUMED 2.8.2 (example for the SOL cluster at Arizona State University)
```
module load gcc-11.2.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
cd gromacs-2022.5-plumed-2.8.2
mkdir build
cd build
cmake .. -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2022.5-plumed-2.8.2 -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_plumed -DGMX_LIBS_SUFFIX=_plumed -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make -j 4
make check
make install
```
You can now use GROMACS with plumed with the command `gmx_plumed`. I also have just normal gromacs installed, which is why I have the weird name change. If you just have gromacs with plumed, the cmake command should be modified to `cmake .. -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2022.5-plumed-2.8.2 -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON`.

## Configuring your BASHRC
Add the following line to your `~/.bashrc` file. Don't forget to run `source ~/.bashrc`!
```
source '$HOME/gromacs-2022.5-plumed-2.8.2/bin/GMXRC.bash'
```
</details>

<details>
<summary> FRESEANCOARSE V2 </summary>
# FRESEANCOARSE V2
There are example scripts provided at `scripts/metad_workflow`. This workflow starts with a pdb file and runs 250 ns well-tempered metadynamics with 0 THz FRESEAN modes. There is a `run.sh` script in each folder that runs the respective step. Each step is briefly described below.

## 00-prep/run.sh 
Prepare your simulation by adding a box around the protein, adding solvent, and generating ions. Keep in mind that the pdb filename, force field, water model and box size will have to be set manually. The default is force field used in our scripts is AMBER99sb-ILDN and tip3p.

## 01-em+equi/run.sh 

Energy minimization (`em.mdp`) and 100ps NPT equilibration (`equi.mdp`). 

## 02-MD/run.sh

20 ns NPT sampling simulation (`sample-NPT.mdp`) with 20 fs output frequency.

## 03-CG/run.sh

Coarse-grain simulation using `fresean coarse`.

## 04-FRESEAN/run.sh

Generate velocity cross-correlation matrix (`fresean covar`), diagonalize the matrix (`fresean eigen`) and extract 0 THz Modes 7 and 8 (`fresean extract`) to `.xyz` format.

## 05-ModeProj/run.sh

Displacement projection of 20 ns trajectory onto FRESEAN modes.

## 06-metadyn/run.sh

Run 250 ns NPT WT-MetaD simulation with plumed input file `plumed-mode-metadyn.dat`. Hills file will be `plumed-mode-metadyn.hills`. Calculate 2D free energy surface in FRESEAN space (`plumed-mode-metadyn.fes`).

## 07-reweight/run.sh

Reweight 2D free energy surface in FRESEAN space to new collective variable space. This will require a plumed input file (`07-reweight/plumed-reweight-CV.dat`) where you can define the space you are reweighting into.

as collective variables
</details>

<details>

<summary> FRESEAN Toolbox Workflow </summary>

# FRESEAN Toolbox Workflow
This is a flowchart of the process used to generate a 2D Free Energy Surface using FRESEAN Modes as CVs. The scripts used to do so are provided in the `scripts/metad_workflow` folder. <br>

<p align="center">
  <img src="imgs/workflow.jpg" width="500" />
</p>

Here is a summarized version of the workflow described above. First, a high-frequency short simulation (10 ns with 20 fs output frequency is what was tested) is run and cross-correlation matrices are generated from the simulation. Then, these matrices are diagonalized and the vibrational modes are extracted into PLUMED format and used as collective variables in matedynamics simulations.

</details>

<details>
  
<summary> Generating FRESEAN Modes </summary>

# FRESEAN Toolbox Programs
Information about the available tools are also accessible by running `fresean`.
## What frequencies can I analyze?
`fresean freqs` will provide output on the current frequency resolution as well as a file containing all frequenceis analyzed given a certain correlation function window. The following command will calculate the frequency resolution and analyzed frequency if we use 500 correlation points spaced at 0.004 ps (for a total correlation time of 2 ps) and output the result to `freqs.txt`.
```
fresean freqs -n 500 -t 0.004 -o freqs.txt
```
## What is an mtop file?
`fresean mtop` is __required__ to convert the topology provided by gromacs into a `.mtop` topology that is recognized by `fresean`. The following command will prompt the user for required information needed to generate the `.mtop` file.
```
fresean mtop -p complex.top
```
## How do I spatially coarsen my system?
`fresean coarse` is used to convert the all-atom trajectory provided by gromacs into a coarsened trajectory containing sidechain and backbone "beads". 

> **Warning**
> Currently, spatial coarse graining is only suported for proteins only containing canonical amino acids and no cofactors. Functionality to be added for future releases.

> **Warning**
> Coarse grain trajectory is output in `.gro` format and must be converted to `.trr` manually. This can be done with `gmx trjconv` (see example script `scripts/metad_workflow/03-CG/run.sh` for more info on how this can be done)
>

```
fresean coarse -f coarse.inp
```

## How can I generate my FRESEAN modes?
`fresean covar` is used to generate the frequency dependent cross-correlation matrices using information from a required input file. Once this code is ran, there will be a `.mmat` binary file containing these matrices. Immediately, `fresean eigen` should be used to read in this binary matrix and generate a new `.mmat` file containing the FRESEAN modes.

> **Note**
> The length of the correlation function must remain consistent between `fresean covar` and `fresean eigen`.

> **Note**
> It is not recommended to perform this analysis with an all-atom high-frequency trajectory. Coarse-graining the trajectory with `fresean coarse` reduces the computation time __significantly__ without loss of important information and is essentially required for performing this analysis on proteins. If you have long trajectories that have been previously coarse grained into 1ns chunks, please see next section for instructions.

```
fresean covar -f gen-modes.inp
fresean eigen -m covar_fresean.mmat -n 500
```
</details>
<details>
<summary> FRESEAN WT-Metadynamics Overview </summary>

## Now I have my modes, now what?
To use the FRESEAN modes as collective variables in an enhanced sampling simulation (such as metadynamics), the vibrational modes must be converted to a format that `Plumed` can understand. <br>
`fresean extract` can be used to convert the binary vibrational modes generated from `fresean eigen` into human-readable xyz format. 
```
fresean extract -f extract.inp
```
Once the modes are converted to `.xyz` format, they must be uncoarsened and converted again into the `Plumed` format.

> **Warning**
> Plumed requires __one__ input file that contains the pdb reference structure and all collective variables. Please refer to the PLUMED documentation (https://www.plumed.org/doc-v2.8/user-doc/html/master-_i_s_d_d-2.html) or utilize `notebooks/PLUMED_Formatting.ipynb` to put this input file together. The python notebook contains further documentation to compile the metadynamics input file.

## Selecting the Correct Metadynamics Parameters

 Selecting the Gaussian Width can be done using short vanilla MD fluctuations along the collective variables. Utilize `notebooks/proj.ipynb`, which contains futher docs on how to select this parameter. Generally, a gaussian width of 0.001 seems to be safe based on collected data. <br>
 
 Selecting the Gaussian height and deposition frequency cannot be reliably chosen a priori. However, a good rule of thumb is to begin by depositing smaller gaussians (=< 0.2 kJ/mol) every 1000 time points. Larger gaussian height will sample the collective variable space much faster, but constantly depositing these large gaussians (> 1 kJ/mol) can cause you to leap over unintended large barriers. <br>

 The bias factor should remain 10. The biasing factor effects the speed of the gaussian height rescaling. However, the gaussian height rescaling can also be controlled by the initial gaussian heigh value. If you are having difficulties overcoming energy barriers, jsut increase the gaussian height by ~0.2 kJ/mol. <br>

## Running WT-Metadynamics
You can run metadynamics using `scripts/template_metad/06-metadyn`. When running metadynamics, the plumed parameter file is provided at `plumed-files/metad/plumed.dat` and is required in addition to the previously prepared pdb file.

</details>


<details>
<summary> Analyzing WT-METAD </summary>

## Analyzing WT-METAD
  The output trajectory of metadynamics includes a hills file of all deposited gaussians (HILLS_PCAVARS1) which can either be directly converted to a free energy surface or reweighted into a different collective variable space. Direct conversion can be done using `plumed-files/metad/sum_hills.sh`. Reweighting can be done using `plumed-files/metad/reweight.sh` with `plumed-files/metad/plumed_reweight_metad.dat` as an input file.

> **Warning**
> These analysis files have to either be manually moved to the directory containing the metadynamics reuslt **OR** the relative path of the input file and all files mentioned in the plumed .dat file have to be changed. I recommend just moving these files. You wont need them anywhere else.
  
</details>
<details>
<summary> FRESEAN Toolbox Input Files </summary>

# FRESEAN Toolbox Input Files
For most of the toolbox, input files are utilized to set parameters. These input files can be found in the `inp-files` folder. Lines starting with a hash are ignored but serve as a header for the variable on the line below. In the following example, the first line is ignored and `complex.mtop` is read in as the parameter value for the `fnTop` variable. 
```
#fnTop
complex.mtop
```
Since this file format is not standard, I will explain each parameter below.
> **Note**
> Parameters with a green ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) box are system-specific and should be modified depending on your protein and filenames. Parameters with a red ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) box are essentially static and the recommended values are provided. 

## 01-coarse.inp
This input file is utilized by `fresean coarse`.<br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnTop`: `.mtop` file name generated from `fresean mtop` <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnCrd`: All-atom trajectory (`.trr` recommended) <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnVel`: If `fnCrd` is `.trr` format, this does not need to be set. For `.xyz` format, `fnCrd` should be the position file and `fnVel` should be the velocity data. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `fnJob`: `.job` file used to define atom groups. Does not need to be changed and can be found in `inp-files` folder. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `nRead`: Number of frames to read from trajectory. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `nSample`: How often we should read from trajectory. Should be set to 1. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `boxSize`: Should be set to the box size used to generate original trajectory. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOutRef`: Name of output all-atom reference file. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOutTraj`: Name of output coarse trajectory in `.gro` format. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOutTopol`: Name of output `.mtop` coarsened topology in `.mtop` format. <br>

## 02-gen-modes.inp and 02-cgen-modes.inp
This file is utilized by `fresean covar`. There are two files provided in `inp-files`. `02-gen-modes.inp` is used for all-atom analysis and `02-cgen-modes` is used if `fresean coarse` was run first (which is recommended!!). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnTop`: `.mtop` file name generated from `fresean mtop` <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnCrd`: All-atom trajectory (`.trr` recommended) <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`fnVel`: If `fnCrd` is `.trr` format, this does not need to be set. For `.xyz` format, `fnCrd` should be the position file and `fnVel` should be the velocity data. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `nRead`: Number of frames to read from trajectory. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`analysisInterval`: How often we should read from trajectory. Used for time coarse graining. Please see [Time Coarse Graining](#Time-Coarse-Graining) for more information on how to vary this parameter. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`fnRef`: Reference file used for translational and rotational fitting <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `alignGrp`: Group of atoms to align to. Should be 0 if provided job file is being used. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `analyzeGrp`: Group of atoms to analyze. Should be 0 if provided job file is being used. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `wrap`: Boundary Conditions. Only 0 supported. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`nCorr`: Number of correlation points. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `winSigma`: Length of Gaussian smoothing function. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `binaryMatrix`: Format of matrix output. 0 for ASCII, 1 for binary `.mmat`. Recommend 1 due to storage cost of ASCII format. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `doGenModes`: Perform Jacobi diagonalization. Should be 0 if using [FRESEAN Toolbox Workflow](#FRESEAN-Toolbox-Workflow).
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `convergence`: Convergence criteria of Jacobi diagonalization. 
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `maxIter`: Maximum number of swaps for Jacobi diagonalization.
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOut`: Name of output `.mmat` file.


## 03-extract.inp
This input file is utilized by `fresean extract`.
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnEigVec`: `.mmat` file containing eigenvectors (output of `fresean eigen`) <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `extractMode`: If set to 0, `freqSel` will extract based on nearest frequency in wavenumbers (cm^-1). If set to 1, `freqSel` will extract based on matrix index. Mode 1 is recommended. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `freqSel`: If `extractMode` set to 0, this is the frequency (in wavenumbers) that you want to extract. If `extractMode` set to 1, this is the matrix index that you want to extract. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `timestep`: Timestep of high frequency simulation. Used to calculate length of correlation function. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `modeStart`: First mode to extract. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `modeEnd`: Last mode to extract. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOut`: String to append to output file. <br>
</details>

<details>
  
<summary> Generating FRESEAN Modes from Longer Trajectories (Dont use yet) </summary>

# Generating FRESEAN Modes from Longer Trajectories

## Averaging Over Cross-Correlation Matrices
All of the required example scripts can be located in `scripts/FRESEAN_average_scripts`. __These scripts may not work for your particular analyte and are only meant to serve as an example of my workflow.__ You are encouraged to write your own scripts depending on your particular needs. <br><br> To generate FRESEAN Modes on longer trajectories (longer than 1 ns), we must first generate individual cross-correlation matrices in 1 ns chunks and then average over all correlation matrices. <br><br>
For example, instead of running `fresean covar` on a 10 ns trajectory, the trajectory should be split using `gmx trjconv` and then `fresean covar` should be run __10 times__ on each trajectory fragment. <br><br>
`fresean avg` can be used to average over all of the generated cross-correlation matrices as follows.
```
fresean avg -f avg.inp
```
After running `fresean avg`, only one cross-correlation matrix will remain and can be diagonalized like normal using `fresean eigen`. <br>

## 05-average.inp
This input file is utilized by `fresean avg`.<br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`nFiles`: Number of cross-correlation matrices that you wish to average over <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`fnList`: A line-delimited list containing the names of individual cross-correlation matrices <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png)`fnOut`: Output file name <br>
