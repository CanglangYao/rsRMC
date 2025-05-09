The rsRMC program is coded with C++. 
recip-mpi-444.c is the main body of rsRMC, it is programed in parallel.
recip-single-444.c is responsible for the the post-procesing of the structure obtained by recip-mpi-444.c, it is programed in sequence.
recip-mpi-444.c has a dependency on Eigen library to calculate rotation matrix. To compile the recip-mpi-444.c, you will need to install Eigen in the first place..
Then run the command: mpic++ -I ~/home/eigen-xx recip-mpi-444.c -o recip-mpi-444.out
Run the command: c++ recip-single-444.c -o recip-single-444.out
