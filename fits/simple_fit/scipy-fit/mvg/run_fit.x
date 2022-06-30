qsub -cwd -V -N fit -l h_data=8G,h_rt=100:00:00,highp -pe shared 5 $PWD/run.bash
