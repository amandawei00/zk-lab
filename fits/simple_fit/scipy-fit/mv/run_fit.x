qsub -cwd -V -N fit -l h_data=4G,h_rt=24:00:00,highp -pe shared 5 $PWD/run.bash
