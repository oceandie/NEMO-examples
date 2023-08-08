#!/bin/bash

#-------------------------------------------------------------------------------------
# INPUT section 
exp_cfg="4"  # 0: Shchepetkin and McWilliams (2002)
             # 1: Ezer, Arango and Shchepetkin (2002)
             # 2: Amy Young
             # 3: Domain: Ezer et al 2002 
             #    Initialisation: Amy Young
             # 4: Domain: Amy Young
             #    Initialisation: Ezer et al. 2022 
outfreq="3h" # "1ts" "1h" "3h" "1d"
expdays=180  # number of days we want to run our simulations
basedir=${PWD}
nemodir="/projects/jmmp/dbruciaf/NEMO/CHECKOUTS_4.0/NEMO_4.0-TRUNK_r14960_HPG"
testdir=${nemodir}"/tests/SEAMOUNT"
exp_ref=${testdir}"/EXPREF"
nemoexe=${testdir}"/BLD/bin/nemo.exe"
nam_tmp=${testdir}"/EXPREF/nam_cfg/namelist_cfg.steep.template"

#-------------------------------------------------------------------------------------
# SETTING the environment for python
isloaded=`module -t list 2> >(grep scitools/default-current)`
[ -z "$isloaded" ] && module load scitools/default-current

#-------------------------------------------------------------------------------------
# SETTING the general experimental setup in the namelist

case ${exp_cfg} in
  0) # Shchepetkin and McWilliams (2002)
    # 1. INITIAL CONDITION
    f90nml -g namusr_def \
           -v nn_ini_cond=${exp_cfg} \
           -p ${nam_tmp} "namcfg.tmpA"
    # 2. EOS
    f90nml -g nameos \
           -v ln_seos=.true. \
           -p "namcfg.tmpA" "namcfg.cfg"    
    ;;
  1) # Ezer, Arango and Shchepetkin (2002)
    # 1. INITIAL CONDITION
    f90nml -g namusr_def \
           -v nn_ini_cond=${exp_cfg} \
           -p ${nam_tmp} "namcfg.tmpA"
    # 2. EOS
    f90nml -g nameos \
           -v ln_eos80=.true. \
           -p "namcfg.tmpA" "namcfg.cfg"
    ;; 
  2) # Amy Young
    # 0.  Forcing model output with freq of 3h since xios on-line averaging requires that 
    #     the period of averaging is a multiple of the timestep
    if [ "${outfreq}" != "3h" ] && [ "${outfreq}" != "1d" ]; then 
       outfreq="3h"
    fi
    # 1. INITIAL CONDITION
    f90nml -g namusr_def \
           -v nn_ini_cond=${exp_cfg} \
           -p ${nam_tmp} "namcfg.tmpA"
    # 2.  Model domain and mesh
    f90nml -g namusr_def \
           -v rn_xdim=380.0 \
           -v rn_ydim=288.0 \
           -v rn_dx=4000.0 \
           -v rn_dz=450. \
           -p "namcfg.tmpA" "namcfg.tmpB"
    # 3.  Initial condition
    f90nml -g namusr_def \
           -v rn_initrho=0.1 \
           -v rn_s=2.0 \
           -p "namcfg.tmpB" "namcfg.tmpC"
    # 4.  Baroclinic timestep
    f90nml -g namdom \
           -v rn_Dt=432. \
           -p "namcfg.tmpC" "namcfg.tmpD"
    # 5.  EOS
    f90nml -g nameos \
           -v ln_eeos=.true. \
           -p "namcfg.tmpD" "namcfg.tmpE"
    # 6.  Momentum advection scheme
    f90nml -g namdyn_adv \
           -v ln_dynadv_vec=.false. \
           -v ln_dynadv_cen2=.true. \
           -p "namcfg.tmpE" "namcfg.tmpF"
    # 7.  Vorticity / Coriolis scheme
    f90nml -g namdyn_vor \
           -v ln_dynvor_een=.false. \
           -v ln_dynvor_ene=.true. \
           -p "namcfg.tmpF" "namcfg.tmpG"
    # 8.  Barotropic timestep
    f90nml -g namdyn_spg \
           -v nn_e=36 \
           -p "namcfg.tmpG" "namcfg.tmpH"
    # 9.  Lateral viscosity
    f90nml -g namdyn_ldf \
           -v rn_Uv=2.0 \
           -p "namcfg.tmpH" "namcfg.tmpI"
    # 10.  Vertical diff/visc
    f90nml -g namzdf \
           -v rn_avm0=0.0 \
           -v rn_avt0=0.0 \
           -p "namcfg.tmpI" "namcfg.tmpJ"
    # 11. MPI halos
    f90nml -g nammpp \
           -v nn_hls=2 \
           -p "namcfg.tmpJ" "namcfg.cfg"
    ;;
  3) # Domain: Ezer et al 2002, Initialisation: Amy Young
    # 1.  Initial condition
    f90nml -g namusr_def \
           -v nn_ini_cond=2 \
           -p ${nam_tmp} "namcfg.tmpA"
    f90nml -g namusr_def \
           -v rn_initrho=0.1 \
           -v rn_s=2.0 \
           -p "namcfg.tmpA" "namcfg.tmpB"
    # 2.  EOS
    f90nml -g nameos \
           -v ln_eeos=.true. \
           -p "namcfg.tmpB" "namcfg.cfg"   
    ;;
  4) # Domain: Amy Young, Initialisation: Ezer et al 2002
    # 1. INITIAL CONDITION
    f90nml -g namusr_def \
           -v nn_ini_cond=1 \
           -p ${nam_tmp} "namcfg.tmpA"
    # 2.  Model domain and mesh
    f90nml -g namusr_def \
           -v rn_xdim=380.0 \
           -v rn_ydim=288.0 \
           -v rn_dx=4000.0 \
           -v rn_dz=450. \
           -p "namcfg.tmpA" "namcfg.tmpB"
    # 3. EOS
    f90nml -g nameos \
           -v ln_eos80=.true. \
           -p "namcfg.tmpB" "namcfg.cfg"
    ;;
esac

rm namcfg.tmp?
nam_tmp="${basedir}/namcfg.cfg"

#-------------------------------------------------------------------------------------
#for vco in sig s94 vqs zco; do
vco=sig

    for hpg in sco prj djc djcr ffl ffq_cub ffq_ccs fflr ffq_cubr ffq_ccsr; do

        for ini in pnt ave; do

            #for cor in fp4 fp5; do
            cor=fp4
    
                expname="SEAMOUNT_cfg_${exp_cfg}_${vco}_${hpg}_${ini}_${cor}"
                exp_dir=${testdir}"/cfg_${exp_cfg}_${vco}_${hpg}_${ini}_${cor}_${outfreq}"
                nam_cfg=${exp_dir}"/namelist_cfg"
            
                if [ ! -d "${exp_dir}" ]; then 
                   mkdir ${exp_dir}
                   cp ${exp_ref}/*{xml,ref} ${exp_dir}
                   ln -s ${nemoexe} ${exp_dir}"/nemo"
                   ln -s ${basedir}"/run_job.sh" ${exp_dir}"/run_job.sh"  
                   ln -s ${exp_dir}"/file_def_nemo-oce_"${outfreq}".xml" ${exp_dir}"/file_def_nemo-oce.xml"

                   # UPDATING NAMELIST -----------

                   # 0. Experience  name
                   f90nml -g namrun \
                          -v cn_exp=${expname} \
                          -p ${nam_tmp} ${nam_cfg}".tmp0"

                   # 1. Run length
                   if [ "${outfreq}" == "1ts" ]; then
                      nn_itend=1 # 1 timesteps
                   else
                      secxday=86400
                      if [ "${exp_cfg}" == "2" ]; then
                         rn_rDT=432
                      else
                         rn_rDT=360
                      fi
                      nn_itend=$((secxday/rn_rDT*expdays))
                   fi
                   f90nml -g namrun \
                          -v nn_itend=${nn_itend} \
                          -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"  

                   # 2. Vertical coordinates
                   case ${vco} in
                     sig)
                       f90nml -g namusr_def \
                              -v ln_sco=.true. \
                              -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     s94)
                       f90nml -g namusr_def \
                              -v ln_sco=.true. \
                              -v ln_s_sh94=.true. \
                              -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     vqs)
                       f90nml -g namusr_def \
                              -v ln_sco=.true. \
                              -v ln_s_sh94=.true. \
                              -v ln_vqs=.true. \
                              -v rn_rmax=0.2 \
                              -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     zco)
                       f90nml -g namusr_def \
                              -v ln_sco=.true. \
                              -v ln_s_sh94=.true. \
                              -v ln_vqs=.true. \
                              -v rn_rmax=0.0 \
                              -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                   esac

                   # 3. HPG scheme
                   case ${hpg} in
                     sco|prj|djc)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_${hpg}=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     djcr)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_djr=.true. \
                              -v ln_hpg_ref=.true. \
                              -v ln_hpg_ref_str=.true. \
                              -v ln_hpg_ref_ccs=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffl)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffq_cub)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -v ln_hpg_ffr_vrt_quad=.true. \
                              -v ln_hpg_ffr_hor_cub=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffq_ccs)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -v ln_hpg_ffr_vrt_quad=.true. \
                              -v ln_hpg_ffr_hor_ccs=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     fflr)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -v ln_hpg_ref=.true. \
                              -v ln_hpg_ref_str=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffq_cubr)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -v ln_hpg_ffr_vrt_quad=.true. \
                              -v ln_hpg_ffr_hor_cub=.true. \
                              -v ln_hpg_ref=.true. \
                              -v ln_hpg_ref_str=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffq_ccsr)
                       f90nml -g namdyn_hpg \
                              -v ln_hpg_ffr=.true. \
                              -v ln_hpg_ffr_vrt_quad=.true. \
                              -v ln_hpg_ffr_hor_ccs=.true. \
                              -v ln_hpg_ref=.true. \
                              -v ln_hpg_ref_str=.true. \
                              -v ln_hpg_ref_ccs=.true. \
                              -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3" 
                   esac
 
                   # Computational halo
                   f90nml -g nammpp \
                          -v nn_hls=2 \
                          -p ${nam_cfg}".tmp3" ${nam_cfg}".tmp4"

                   # 4. Initial condition formulation
                   case ${ini} in
                     pnt)
                       f90nml -g namusr_def \
                              -v ln_init_pt_val=.true. \
                              -p ${nam_cfg}".tmp4" ${nam_cfg}".tmp5"
                       ;;
                     ave)
                       f90nml -g namusr_def \
                              -v ln_init_pt_val=.false. \
                              -v ln_init_curved=.false. \
                              -p ${nam_cfg}".tmp4" ${nam_cfg}".tmp5"
                       ;;
                   esac

                   # 5. Coriolis parameter
                   case ${cor} in
                     fp4)
                       f90nml -g namusr_def \
                              -v rn_fplane=1.0e-4 \
                              -p ${nam_cfg}".tmp5" ${nam_cfg}
                       ;;
                     fp5)
                       f90nml -g namusr_def \
                              -v rn_fplane=1.0e-5 \
                              -p ${nam_cfg}".tmp5" ${nam_cfg}
                       ;;
                   esac

                   rm ${nam_cfg}".tmp"?

                   xml_exp=${exp_dir}"/file_def_nemo-oce.xml"
                   case ${hpg} in
                     sco|prj|djc)
                       # UPDATING file_def_nemo-oce.xml
                       sed -i 's%<field field_ref="rhd_hpg" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_west" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_south" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="e3w" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="pressure" />%%g' ${xml_exp}
                       ;;
                     djcr)
                       # UPDATING file_def_nemo-oce.xml
                       sed -i 's%<field field_ref="u_force_west" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_south" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_upper" />%%g' ${xml_exp}
                       ;;
                   esac
                   if [ "${outfreq}" != "1ts" ] && [ "${hpg}" != fflr ] && [ "${hpg}" != djcr ]; then
                      sed -i 's%<field field_ref="rhd_ref" />%%g' ${xml_exp}
                      sed -i 's%<field field_ref="jk_bot_ref" />%%g' ${xml_exp}
                      sed -i 's%<field field_ref="jk_ref_for_tgt" />%%g' ${xml_exp}
                      sed -i 's%<field field_ref="zfld_ref_on_tgt" />%%g' ${xml_exp}
                      sed -i 's%<field field_ref="p_fld_tgt" />%%g' ${xml_exp}
                   fi

                   # LAUNCHING SIMULATION ------------
                   cd ${exp_dir}
                   echo "  ${expname}  ${outfreq}"
                   qsub run_job.sh
                   cd ${basedir}
                else
                   echo "   ... experiment already exists"
                fi
            #done
        done
    done
#done

rm ${nam_tmp}
