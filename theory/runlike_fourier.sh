# rm shear_power_tomo_like
# gcc -Wall -I/usr/local/include -L/usr/local/lib -o shear_power_tomo_like shear_power_tomo_like.c -lfftw3 -lgsl -lgslcblas -lm 
# ./shear_power_tomo_like

# rm cl_power_tomo_like
# gcc -Wall -I/usr/local/include -L/usr/local/lib -o cl_power_tomo_like cl_power_tomo_like.c -lfftw3 -lgsl -lgslcblas -lm 
# ./cl_power_tomo_like

# gcc -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -shared -o like.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 

rm like
gcc -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -o like like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 
./like


#zodiac
# setenv LD_LIBRARY_PATH /home/teifler/lib
# gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o like like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 

#gprof on zodiac
# gcc -Wno-missing-braces -Wno-missing-field-initializers -pg -I/home/teifler/include -L/home/teifler/lib -o like like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 
# ./like
#gprof test_gprof gmon.out > analysis.txt



# gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -shared -o like.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D1500 ./runDES_2pt_clusterN_clusterWL_no_sys_1500.py > & /aurora_nobackup/sunglass/teifler/myoutput1.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D2000 ./runDES_2pt_clusterN_clusterWL_no_sys_2000.py > & /aurora_nobackup/sunglass/teifler/myoutput2.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D3000 ./runDES_2pt_clusterN_clusterWL_no_sys_3000.py > & /aurora_nobackup/sunglass/teifler/myoutput3.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D4000 ./runDES_2pt_clusterN_clusterWL_no_sys_4000.py > & /aurora_nobackup/sunglass/teifler/myoutput4.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D4500 ./runDES_2pt_clusterN_clusterWL_no_sys_4500.py > & /aurora_nobackup/sunglass/teifler/myoutput5.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N D5000 ./runDES_2pt_clusterN_clusterWL_no_sys_5000.py > & /aurora_nobackup/sunglass/teifler/myoutput6.log



# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=12:00:00 -q mediumq -o /home/teifler/output/ -e /home/teifler/output/ -N Lpp ./runLSST_pos_pos_bias_phot_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput1.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LclN ./runDES > & /aurora_nobackup/sunglass/teifler/myoutput2.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N Lss ./runLSST_shear_shear_phot_shear_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput3.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2pt ./runLSST_2pt_phot_shear_bias_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput4.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N LclN+WL ./runLSST_clusterN_clusterWL_shear_Mobs_phot_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput5.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WL ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput6.log


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2pt20 ./runLSST_2pt_phot_shear_bias_Rmin20_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput7.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2pt50 ./runLSST_2pt_phot_shear_bias_Rmin50_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput8.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptHO ./runLSST_2pt_phot_shear_bias_HOD_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput9.log 


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WL20 ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_Rmin20_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput10.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WL50 ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_Rmin50_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput11.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLHO ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_HOD_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput12.log


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LssIAim ./runLSST_shear_shear_phot_shear_IA_impact_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput13.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=96:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N LssIAma ./runLSST_shear_shear_phot_shear_IA_marg_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput14.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptim ./runLSST_2pt_phot_shear_bias_IA_impact_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput15.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=96:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptma ./runLSST_2pt_phot_shear_bias_IA_marg_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput16.log


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LssNB ./runLSST_shear_shear_phot_shear_sys_NOBLEND.py > & /aurora_nobackup/sunglass/teifler/myoutput17.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptNB ./runLSST_2pt_phot_shear_bias_sys_NOBLEND.py > & /aurora_nobackup/sunglass/teifler/myoutput18.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLNB ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_NOBLEND.py > & /aurora_nobackup/sunglass/teifler/myoutput19.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LssIAim ./runLSST_shear_shear_IA_impact_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput20.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N Lssc2 ./runLSST_shear_shear_phot_shear_sys_cosmomax2.py > & /aurora_nobackup/sunglass/teifler/myoutput121.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptc2 ./runLSST_2pt_phot_shear_bias_sys_cosmomax2.py > & /aurora_nobackup/sunglass/teifler/myoutput22.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLc2 ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_cosmomax2.py > & /aurora_nobackup/sunglass/teifler/myoutput23.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LppHL ./runLSST_pos_pos_bias_phot_sys_HIGH_NLENS.py > & /aurora_nobackup/sunglass/teifler/myoutput24.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptHL ./runLSST_2pt_phot_shear_bias_sys_HIGH_NLENS.py > & /aurora_nobackup/sunglass/teifler/myoutput25.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLHL ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_HIGH_NLENS.py > & /aurora_nobackup/sunglass/teifler/myoutput26.log


# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptcm ./runLSST_2pt_phot_shear_bias_sys_cosmomin.py > & /aurora_nobackup/sunglass/teifler/myoutput25.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLcm ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_cosmomin.py > & /aurora_nobackup/sunglass/teifler/myoutput26.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptim ./runLSST_2pt_IA_impact_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput27.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N Lssim ./runLSST_shear_shear_phot_shear_IA_impact_sys_offset.py > & /aurora_nobackup/sunglass/teifler/myoutput28.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=12:00:00 -q mediumq -o /home/teifler/output/ -e /home/teifler/output/ -N Lpp ./runLSST_pos_pos_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput30.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LclN ./runLSST_clusterN_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput31.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N Lss ./runLSST_shear_shear_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput32.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2pt ./runLSST_2pt_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput33.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N LclN+WL ./runLSST_clusterN_clusterWL_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput34.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=72:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WL ./runLSST_2pt_clusterN_clusterWL_no_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput35.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=96:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N LssIAma ./runLSST_shear_shear_phot_shear_IA_marg_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput29.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=96:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptIAma ./runLSST_2pt_phot_shear_bias_IA_marg_sys.py > & /aurora_nobackup/sunglass/teifler/myoutput29.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptcovOS ./runLSST_2pt_no_sys_covvaryOS.py > & /aurora_nobackup/sunglass/teifler/myoutput36.log 

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptcovW ./runLSST_2pt_no_sys_covvaryW.py > & /aurora_nobackup/sunglass/teifler/myoutput37.log 











# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2pt ./runLSST_2pt_nosys.py > & /aurora_nobackup/sunglass/teifler/myoutput7.log

# qsub -S /bin/bash -V -l select=1:ncpus=16:mem=1GB -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptcos ./runLSST_2pt_nosys_cosmomax.py > & /aurora_nobackup/sunglass/teifler/myoutput8.log





# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLCOS ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_cosmomax.py > & /nobackup0/sunglass/teifler/myoutput.log 

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLrmin ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_Rmin20_sys.py > & /nobackup0/sunglass/teifler/myoutput.log 

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptrmin ./runLSST_2pt_phot_shear_bias_Rmin20_sys.py > & /nobackup0/sunglass/teifler/myoutput.log 

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptrmin ./runLSST_2pt_phot_shear_bias_Rmin50_sys.py > & /nobackup0/sunglass/teifler/myoutput3.log 


# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLrmin ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_Rmin50_sys.py > & /nobackup0/sunglass/teifler/myoutput2.log 


# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLim ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_IA_impact_sys.py

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WLTm ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_IA_marg_sys.py



# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N L2ptclN+WL2 ./runLSST_2pt_clusterN_clusterWL_phot_shear_bias_Mobs_sys_cosmomax.py


# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=24:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N LssIAma ./runLSST_shear_shear_phot_shear_IA_marg_sys.py

# #!/bin/bash -xv

# ##to create shared library

# # gcc -Wall -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -shared -o combi_power_tomo_like.so -fPIC combi_power_tomo_like.c -lfftw3 -lgsl -lgslcblas -lm 

# # rm combi_power_tomo_like
# # gcc -Wall -I/usr/local/include -L/usr/local/lib -o combi_power_tomo_like combi_power_tomo_like.c -lfftw3 -lgsl -lgslcblas -lm
# # ./combi_power_tomo_like


# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB -l place=free -l walltime=24:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N DES_Y1 ./runDES_Y1_cosmolike.py

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=24:00:00 -q longq -o /home/teifler/output/ -e /home/teifler/output/ -N Euclid ./runEuclid_cosmolike.py

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=12:00:00 -q mediumq -o /home/teifler/output/ -e /home/teifler/output/ -N SmULDB ./runSmall_ULDB_cosmolike.py

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N MedULDB ./runMedium_ULDB_cosmolike.py

# qsub -S /bin/bash -V -l select=1:ncpus=12:mem=1GB:rack1=false -l place=free -l walltime=48:00:00 -q verylongq -o /home/teifler/output/ -e /home/teifler/output/ -N LaULRDB ./runLarge_ULDB_cosmolike.py



# # gcc -Wall -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o combi_power_tomo_like combi_power_tomo_like.c -lfftw3 -lgsl -lgslcblas -lm 








