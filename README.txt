Hello, welcome! Here you'll find the sum of all the work I've done with Professor Zhongob Kang in his nuclear physics lab. My research primarily involves producing numerical solutions to evolution equations, structure functions, and differential cross sections and analyzing the results to study properties of the color glass condensate.
Below, you'll find a brief description of each of the folders and what's inside. Check the README in each of the folders for more information!
1. bk
	- Solution to Balitsky-Kovchegov solution (main file: bk_solver.py)
	- results: contains results of various parameters and run-time parameters in 'directory.csv'.
	- Reference papers (papers)
2. fits
	- Program to run global fitting over DIS, pp and pA data to fit four parameters in the BK initial conditions and running-coupling equations
3. dis
	- Solutions to structure function and reduced cross section equations
	- Results for various initial conditions of the BK solution
	- Code to plot solutions
	- Reference papers (papers)
4. pp-pA
	- Solutions to differential cross section describing single inclusive hadron production in the forward region (solver.py)
	- Results for various initial conditions of the BK solution
	- Code to plot solutions
	- Reference papers (papers)
5. dihadron
	- Solutions to transverse momentum distributions (still under construction, check back later!)
	- Reference papers and notes (papers)
6. data
	- Data gathered from hepdata.net that I use to compare all my solutions to
7. python_tools
	- Scripts I wrote to do stuff I needed.

Everything here I wrote myself, except the contents of the 'LHAPDF-6.4.0', 'bin', 'include', 'lib', 'lib64' and 'share' folders, which were installed from lhapdf.hepforge.org. However, I could not have done all this without the generous help, support and advice of Professor Kang, Jared Reiten, John Terry, Farid Salazar and Daniel Callos. 


If you have tips or suggestions for improvement, please email me at amandawei00@gmail.com!
