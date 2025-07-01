This is an example of rsRMC program in the simulation of alpha-MAPbI3 structure with a 4*4*4 supercell. 

The rsRMC program is coded with C++.  
recip-mpi-444.c is the main body of rsRMC, it is programed in parallel. 
recip-single-444.c is responsible for the the post-procesing of the structure obtained by recip-mpi-444.c, it is programed in sequence. 
recip-mpi-444.c has a dependency on Eigen library to calculate rotation matrix. To compile the recip-mpi-444.c, you will need to install Eigen in the first place. 
Then run the command: mpic++ -I ~/home/eigen-xx recip-mpi-444.c -o recip-mpi-444.out. 
Also run the command: c++ recip-single-444.c -o recip-single-444.out. 

"supercell" contains the test structure information (elements and coordinates), "exp-pdf" is the measured PDF, "order" contains a serial of numbers in random order, it guides the rsRMC program to displace which atom, and will be updated after each RMC cycle. "str-pdf.py" is a pre-processing program that collects PDF and structure of previous cycles for the rsRMC program running in the next cycle.

After the compilation of recip-mpi-444.c and recip-single-444.c, one need only type: ./pdf.script to run the simulation. The number of RMC cycles and the cores used for the simulation can be modified in the text of "pdf.script".
