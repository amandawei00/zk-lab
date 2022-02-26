qsub -cwd -V -N pp-pA -l h_data=4G,h_rt=4:00:00,highp $PWD/run.bash
