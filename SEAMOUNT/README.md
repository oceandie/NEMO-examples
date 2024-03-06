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

`cd /YOUR/PATH/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/MY_SRC`

By default, when cloning (aka checking-out) a repository, we clone the entire git repository, including all the branches.
Also, once cloned a repository, by default we are in the `master` (or `main`, depending on the version of git) branch, meaning that we "can" see the code of the `master` (`main`) branch. If you want to see/work on the `dev_seamount` development branch, you need to:

`git checkout dev_seamount`

After this, the code of the `dev_seamount` branch will be visible. 

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


## COMMITTING THE MODIFICATIONS TO GIT

1) Since XCE/XCF are not connected with the web, we need to copy the modified files in the copy of the repository checked out in the VDI. On your local linux machine, run:
```
scp  <username>@<HPC-platform>:/YOUR/PATH/HPG/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/your_files /YOUR/PATH/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/
```

2) Once you copied the files you want to commit to your VDI:
   ```
   cd /YOUR/PATH/NEMO_4.0-TRUNK_r14960_HPG/tests/SEAMOUNT/
   git checkout dev_seamount
   ``` 

