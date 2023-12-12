The C file is provided with (MPI.c) which get compiled using the following command:

--    mpicc -g MPI.c

Once the code will be compiled it, it would be saved as and executable called a.out.
This executable can be used for running the program, following is the command to do so:
 mpirun -np <number> ./a.out [N] [random Seed]
here, 
number: number of processors to run the program on.
N: dimension of the input matrix, if N = 100, matrix = 100x100
randomseed = random seed for initial values of the matrix.

To measure the correctness of algorithm, N>=5 can be given, this will output the input 
vector and the processed vector along with the X vector values. The scaling of the algorithm can 
be checked by given different number and N, this is also done in the pdf report.
