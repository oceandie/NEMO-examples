#!/bin/bash

# INPUT section --------------------------------------------------------------------
debug=false
basedir=${PWD}
nemodir="/projects/jmmp/dbruciaf/NEMO/CHECKOUTS_4.0/NEMO_4.0-TRUNK_r14960_HPG"
testdir=${nemodir}"/tests/SEAMOUNT"
exp_ref=${testdir}"/EXPREF"
nemoexe=${testdir}"/BLD/bin/nemo.exe"
nam_tmp=${testdir}"/EXPREF/nam_cfg/namelist_cfg.steep.template"

# SETTING the environment for python
isloaded=`module -t list 2> >(grep scitools/default-current)`
[ -z "$isloaded" ] && module load scitools/default-current

for vco in s94 vqs zco; do
#vco=vqs

    for hpg in sco prj djc ffl ffq fflr; do
    #hpg=fflr        

        #for ini in pnt ave; do
        ini=pnt

            #for cor in fp4 fp5; do
            cor=fp4
    
                exp_dir=${testdir}"/"${vco}"_"${hpg}"_"${ini}"_"${cor}
                nam_cfg=${exp_dir}"/namelist_cfg"
            
                if [ ! -d "${exp_dir}" ]; then 
                   mkdir ${exp_dir}
                   cp ${exp_ref}/*{xml,ref} ${exp_dir}
                   ln -s ${nemoexe} ${exp_dir}"/nemo"
                   ln -s ${basedir}"/run_job.sh" ${exp_dir}"/run_job.sh"  
                   if [[ "$debug" == true ]]; then
                      ln -s ${exp_dir}"/file_def_nemo-oce_1ts.xml" ${exp_dir}"/file_def_nemo-oce.xml"
                   else
                      ln -s ${exp_dir}"/file_def_nemo-oce_1h.xml" ${exp_dir}"/file_def_nemo-oce.xml"
                   fi
                   # UPDATING NAMELIST -----------
                   # 1. Vertical coordinates
                   case ${vco} in
                     s94)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -p ${nam_tmp} ${nam_cfg}".tmp0"
                       ;;
                     vqs)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -v ln_vqs=.true. -v rn_rmax=0.2 -p ${nam_tmp} ${nam_cfg}".tmp0"
                       ;;
                     zco)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -v ln_vqs=.true. -v rn_rmax=0.0 -p ${nam_tmp} ${nam_cfg}".tmp0"
                   esac

                   # 2. HPG scheme
                   case ${hpg} in
                     sco|prj|djc)
                       f90nml -g namdyn_hpg -v ln_hpg_${hpg}=.true. -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"
                       ;;
                     ffl)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"
                       ;;
                     ffq)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -v ln_hpg_ffr_vrt_quad=.true. -v ln_hpg_ffr_hor_cub=.true. -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"
                       ;;
                     fflr)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -v ln_hpg_ref=.true. -v ln_hpg_ref_str=.true. -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"
                   esac

                   # 3. Initial condition formulation
                   case ${ini} in
                     pnt)
                       f90nml -g namusr_def -v ln_init_pt_val=.true. -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     ave)
                       f90nml -g namusr_def -v ln_init_pt_val=.false. -v ln_init_curved=.false. -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                   esac

                   # 4. Coriolis parameter
                   case ${cor} in
                     fp4)
                       f90nml -g namusr_def -v rn_fplane=1.0e-4 -p ${nam_cfg}".tmp2" ${nam_cfg}
                       ;;
                     fp5)
                       f90nml -g namusr_def -v rn_fplane=1.0e-5 -p ${nam_cfg}".tmp2" ${nam_cfg}
                       ;;
                   esac

                   rm ${nam_cfg}".tmp"?

                   case ${hpg} in
                     sco|prj|djc)
                       # UPDATING file_def_nemo-oce.xml
                       xml_exp=${exp_dir}"/file_def_nemo-oce.xml"
                       sed -i 's%<field field_ref="rhd_hpg" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_west" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_south" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="gdepw_hpg" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="pressure" />%%g' ${xml_exp}
                       ;;
                   esac

                   # LAUNCHING SIMULATION ------------
                   cd ${exp_dir}
                   echo "  ${vco}-${hpg}-${ini}-${cor}"
                   qsub run_job.sh
                   cd ${basedir}
                else
                   "... experiment already exists"
                fi
            #done
        #done
    done
done
