#!/bin/bash

# INPUT section --------------------------------------------------------------------
basedir=${PWD}
nemodir="/projects/jmmp/dbruciaf/NEMO/CHECKOUTS_4.0/NEMO_4.0-TRUNK_r14960_HPG"
testdir=${nemodir}"/tests/SEAMOUNT"
exp_ref=${testdir}"/EXPREF"
nemoexe=${testdir}"/BLD/bin/nemo.exe"

# SETTING the environment for python
isloaded=`module -t list 2> >(grep scitools/default-current)`
[ -z "$isloaded" ] && module load scitools/default-current

for seamount in steep moderate ; do

    nam_tmp=${testdir}"/EXPREF/nam_cfg/namelist_cfg."${seamount}".template"

    for hpg in sco prj djc ; do
        for vcoord in sig sh94 vqs ; do

            exp_dir=${testdir}"/"${hpg}"_"${vcoord}"_"${seamount}
            nam_cfg=${exp_dir}"/namelist_cfg"
            
            if [ ! -d "${exp_dir}" ]; then 
               mkdir ${exp_dir}
               cp ${exp_ref}/*{xml,ref} ${exp_dir}
               ln -s ${nemoexe} ${exp_dir}"/nemo"
               ln -s ${basedir}"/run_job.sh" ${exp_dir}"/run_job.sh"  
               #cp ${nam_tmp} ${nam_cfg}

               # UPDATING NAMELIST -----------
               # 1. HPG scheme
               f90nml -g namdyn_hpg -v ln_hpg_${hpg}=.true. -p ${nam_tmp} ${nam_cfg}".tmp"

               # 2. Vertical coordinates
               case ${vcoord} in
                 sig) 
                   f90nml -g namusr_def -v ln_sco=.true. -p ${nam_cfg}".tmp" ${nam_cfg}
                   ;;
                 sh94)
                   f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -p ${nam_cfg}".tmp" ${nam_cfg}
                   ;;
                 vqs)
                   f90nml -g namusr_def -v ln_sco=.true. -v ln_vqs=.true. -p ${nam_cfg}".tmp" ${nam_cfg}
                   ;;
               esac
               rm ${nam_cfg}".tmp"

               # UPDATING file_def_nemo-oce.xml
               xml_exp=${exp_dir}"/file_def_nemo-oce.xml"
               sed -i 's%<field field_ref="rhd_hpg" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="u_force_west" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="u_force_upper" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="v_force_south" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="v_force_upper" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="gdepw_hpg" />%%g' ${xml_exp}
               sed -i 's%<field field_ref="pressure" />%%g' ${xml_exp}

               # LAUNCHING SIMULATION ------------
               cd ${exp_dir}
               echo "  ${hpg} ${vcoord}"
               qsub run_job.sh
               cd ${basedir}
            else
               "... experiment already exists"
            fi
        done
    done
done
