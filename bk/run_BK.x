qsub -cwd -V -N BK -l h_data=4G,h_rt=10:00:00,highp -pe shared 4 $PWD/run.bash
