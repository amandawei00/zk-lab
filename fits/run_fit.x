qsub -cwd -V -N fit -l h_data=4G,h_rt=8:00:00,highp -pe shared 50 $PWD/run.bash
