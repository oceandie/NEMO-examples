# SEAMOUNT demonstration case

This notes describe how to install, compile and run the SEAMOUNT test case on the Met Office HPC (XCE/XCF). 

## CHECKING OUT THE CODE 

In your local linux VDI:

1. `svn co https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/UKMO/NEMO_4.0-TRUNK_r14960_HPG` 

2. `cd NEMO_4.0-TRUNK_r14960_HPG/` 

3. `svn update makenemo -r 14946` 

4. `svn update tools -r 14946` 

5. `svn update mk -r 14946` 

6. `git clone https://github.com/oceandie/NEMO-examples.git tests` 

7. `rsync -ave ssh  /YOUR/PATH/NEMO_4.0-TRUNK_r14960_HPG <username>@<HPC-platform>:/desired/path/`

The code for the NEMO HPG schemes to test and the setup of the SEAMOUNT test case can be found at:

`/YOUR/PATH/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/MY_SRC`


## COMPILING THE SEAMOUNT TESTCASE ON THE HPC (XCE/XCF)

1. cd to your SEAMOUNT path (e.g `cd /YOUR/PATH/HPG/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT`)

2. Switch to the development branch: `git checkout dev_seamount`

3. cd the main NEMO path (e.g `cd /YOUR/PATH/HPG/NEMO_4.0-TRUNK_r14960_HPG`) 

4. `source ukmo_utils/use_intel_hpc.sh` 

5. `./makenemo -a SEAMOUNT -m XC40_METO_IFORT`



## RUNNING SEAMOUNT EXPERIMENTS 

1. `cd /YOUR-PATH/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/run_exp`

2. Adapt run_exp_seamount.sh to your needs.

3. Open run_job.sh and adapt the number of cores you want to use 

4. `./run_exp_seamount.sh`

