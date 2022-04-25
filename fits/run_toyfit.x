qsub -cwd -V -N toyfit -l h_data=4G,h_rt=12:00:00,highp -pe shared 5 $PWD/run_toyfit.bash
