qsub -cwd -V -N auto -l h_data=4G,h_rt=8:00:00,highp -pe shared 5 $PWD/run_auto.bash
