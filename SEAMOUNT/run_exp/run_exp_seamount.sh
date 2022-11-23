#!/bin/bash

#-------------------------------------------------------------------------------------
# INPUT section 
exp_cfg="2"  # 0: Shchepetkin and McWilliams (2002)
             # 1: Ezer, Arango and Shchepetkin (2002)
             # 2: Amy Young
outfreq="1h" # "1ts"
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

f90nml -g namusr_def -v nn_ini_cond=${exp_cfg} -p ${nam_tmp} "namcfg.tmpA"

case ${exp_cfg} in
  0)
    # 1. EOS
    f90nml -g nameos -v ln_seos=.true. -p "namcfg.tmpA" "namcfg.cfg"
    ;;
  1)
    # 1. EOS
    f90nml -g nameos -v ln_eos80=.true. -p "namcfg.tmpA" "namcfg.cfg"
    ;; 
  2) 
    # 1.  Model domain and mesh
    f90nml -g namusr_def -v rn_xdim=380.0 -v rn_ydim=288.0 -v rn_dx=4000.0 -v rn_dz=450. -p "namcfg.tmpA" "namcfg.tmpB"
    # 2.  Initial condition
    f90nml -g namusr_def -v rn_initrho=0.1 -v rn_s=2.0 -p "namcfg.tmpB" "namcfg.tmpC"
    # 3.  Baroclinic timestep
    f90nml -g namdom -v rn_Dt=432. -p "namcfg.tmpC" "namcfg.tmpD"
    # 4.  EOS
    f90nml -g nameos -v ln_eeos=.true. -p "namcfg.tmpD" "namcfg.tmpE"
    # 5.  Momentum advection scheme
    f90nml -g namdyn_adv -v ln_dynadv_vec=.false. -v ln_dynadv_cen2=.true. -p "namcfg.tmpE" "namcfg.tmpF"
    # 6.  Vorticity / Coriolis scheme
    f90nml -g namdyn_vor -v ln_dynvor_een=.false. -v ln_dynvor_ene=.true. -p "namcfg.tmpF" "namcfg.tmpG"
    # 7.  Barotropic timestep
    f90nml -g namdyn_spg -v nn_e=36 -p "namcfg.tmpG" "namcfg.tmpH"
    # 8.  Lateral viscosity
    f90nml -g namdyn_ldf -v rn_Uv=2.0 -p "namcfg.tmpH" "namcfg.tmpI"
    # 9.  Vertical diff/visc
    f90nml -g namzdf -v rn_avm0=0.0 -v rn_avt0=0.0 -p "namcfg.tmpI" "namcfg.tmpJ"
    # 10. MPI halos
    f90nml -g nammpp -v nn_hls=2 -p "namcfg.tmpJ" "namcfg.cfg"
esac

rm namcfg.tmp?
nam_tmp="${basedir}/namcfg.cfg"

#-------------------------------------------------------------------------------------
#for vco in sig s94 vqs zco; do
#for vco in sig s94 vqs; do
vco=sig

    for hpg in sco prj djc ffl ffq fflr; do
    #for hpg in djc djr ffl fflr; do
    #hpg=ffq    

        #for ini in pnt ave; do
        ini=pnt

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
                   f90nml -g namrun -v cn_exp=${expname} -p ${nam_tmp} ${nam_cfg}".tmp0"

                   # 1. Run length
                   case ${outfreq} in
                     1ts)
                       nn_itend=1 # 1 timesteps
                       ;;
                     1h) # 180 days
                       if [ "${exp_cfg}" == 2 ]; then
                          nn_itend=36000 # rn_Dt=432 => 200 time-steps per day
                       else
                          nn_itend=43200 # rn_Dt=360 => 240 time-steps per day
                       fi
                       ;;
                   esac
                   f90nml -g namrun -v nn_itend=${nn_itend} -p ${nam_cfg}".tmp0" ${nam_cfg}".tmp1"  

                   # 2. Vertical coordinates
                   case ${vco} in
                     sig)
                       f90nml -g namusr_def -v ln_sco=.true. -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     s94)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     vqs)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -v ln_vqs=.true. -v rn_rmax=0.2 -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                       ;;
                     zco)
                       f90nml -g namusr_def -v ln_sco=.true. -v ln_s_sh94=.true. -v ln_vqs=.true. -v rn_rmax=0.0 -p ${nam_cfg}".tmp1" ${nam_cfg}".tmp2"
                   esac

                   # 3. HPG scheme
                   case ${hpg} in
                     sco|prj|djc)
                       f90nml -g namdyn_hpg -v ln_hpg_${hpg}=.true. -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     djr)
                       f90nml -g namdyn_hpg -v ln_hpg_djr=.true. -v ln_hpg_ref=.true. -v ln_hpg_ref_str=.true. -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffl)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     ffq)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -v ln_hpg_ffr_vrt_quad=.true. -v ln_hpg_ffr_hor_cub=.true. -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                       ;;
                     fflr)
                       f90nml -g namdyn_hpg -v ln_hpg_ffr=.true. -v ln_hpg_ref=.true. -v ln_hpg_ref_str=.true. -p ${nam_cfg}".tmp2" ${nam_cfg}".tmp3"
                   esac
 
                   # Computational halo
                   #if [ "${hpg}" == djr ] || [ "${hpg}" == ffq ]; then
                      f90nml -g nammpp -v nn_hls=2 -p ${nam_cfg}".tmp3" ${nam_cfg}".tmp4"
                   #else
                   #   cp ${nam_cfg}".tmp3" ${nam_cfg}".tmp4"
                   #fi

                   # 4. Initial condition formulation
                   case ${ini} in
                     pnt)
                       f90nml -g namusr_def -v ln_init_pt_val=.true. -p ${nam_cfg}".tmp4" ${nam_cfg}".tmp5"
                       ;;
                     ave)
                       f90nml -g namusr_def -v ln_init_pt_val=.false. -v ln_init_curved=.false. -p ${nam_cfg}".tmp4" ${nam_cfg}".tmp5"
                       ;;
                   esac

                   # 5. Coriolis parameter
                   case ${cor} in
                     fp4)
                       f90nml -g namusr_def -v rn_fplane=1.0e-4 -p ${nam_cfg}".tmp5" ${nam_cfg}
                       ;;
                     fp5)
                       f90nml -g namusr_def -v rn_fplane=1.0e-5 -p ${nam_cfg}".tmp5" ${nam_cfg}
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
                     djr)
                       # UPDATING file_def_nemo-oce.xml
                       sed -i 's%<field field_ref="u_force_west" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="u_force_upper" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_south" />%%g' ${xml_exp}
                       sed -i 's%<field field_ref="v_force_upper" />%%g' ${xml_exp}
                       ;;
                   esac
                   if [ "${outfreq}" != "1ts" ] && [ "${hpg}" != fflr ] && [ "${hpg}" != djr ]; then
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
        #done
    done
#done

rm ${nam_tmp}
