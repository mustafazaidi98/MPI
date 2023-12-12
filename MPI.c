#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

// Max matrix size
#define MAX_N 2500
//matrix Dimension
int N;
float* AB;
// Create random seed using time
int time_seed() {
    struct timeval t;
    struct timezone tzdummy;

    gettimeofday(&t, &tzdummy);
    return (int)(t.tv_usec);
}
/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
        int seed = 0;  /* Random seed */
        // char uid[32]; /*User name */

        /* Read command-line arguments */
        srand(time_seed());  /* Randomize */

        if (argc == 3) {
                seed = atoi(argv[2]);
                srand(seed);
                printf("Random seed = %i\n", seed);
        }
        if (argc >= 2) {
                N = atoi(argv[1]);
                if (N < 1 || N > MAX_N) {
                        printf("N = %i is out of range.\n", N);
                        exit(0);
                }
        }
        else {
                printf("Usage: %s <matrix_dimension> [random seed]\n",
                                argv[0]);
                exit(0);
        }

        /* Print parameters */

}
void printXvector(float vector[], int len) {
    printf("Result: [");
    for (int i = 0; i < len; i++){
        printf("%.6f ", vector[i]);
        }
    printf("]\n");
}
void print2DMatrix(float matrix[], int rows) {
    for (int row = 0; row < rows; row++){
        for (int column = 0; column < N + 1; column++){
            printf("%.6f ", matrix[row * (N + 1) + column]);
            if(column >= N)
                 printf("\n\t");
        }
    }
        printf("\n");
}
int main(int argc, char** argv) {
    int processID, NoOfPrc;
    double startTime;
    // Initialize mpi
    MPI_Init(&argc, &argv);
    // Get total size of process
    MPI_Comm_size(MPI_COMM_WORLD, &NoOfPrc);
    // Get Current process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    // Read command line args
    parameters(argc, argv);
    //Last column for B
    int columns = N + 1;
    int localRowSize;
    // Assigning Rows to Process
    if((N % NoOfPrc) > processID)
        localRowSize = N / NoOfPrc + 1;
    else
        localRowSize = N / NoOfPrc;
    float* localMatrix = (float*) malloc((localRowSize + 1) * columns  * sizeof(float));
        // first processor responsible for initializing matrix
    if (processID == 0 ) {
        AB = malloc(N * (columns) * sizeof(float));
        for (int n = 0; n <  N * (columns); n++){
                AB[n] = (float)rand() / 32768.0;;
            }
        if(N<=5){
            printf("Input Vector:\n\t");
            print2DMatrix(AB, N);
        }
    }
    // first processor starts the timer
    if (processID == 0)
        startTime = MPI_Wtime();
    // Distribute the data across the proccessors through first processor
    for (int row = 0; row < N; row++) {
        int loc_row = row / NoOfPrc;
        int rowProcessor = row % NoOfPrc;
        MPI_Status status;
        // if row belongs to first processor, copy it to memory
        if (processID == 0 && rowProcessor == 0){
            memcpy(localMatrix + loc_row * columns, AB + row * columns, columns * sizeof(float));
        }
                // distribute data among processors
        else if (processID == 0){
            MPI_Send(AB + row * columns,columns, MPI_FLOAT, rowProcessor, row, MPI_COMM_WORLD);
        }
        // receive data from first processor
        else if (processID == rowProcessor){
            MPI_Recv(localMatrix + loc_row * columns,columns, MPI_FLOAT,0, row, MPI_COMM_WORLD, &status);
        }
        // continue if row doesnt belong to the prcsr
    }
    // rows relevant to processor running are shared
    // now the elimination row will be shared through the processor working on it
    float eliminationRow[columns];
    for (int currentRow = 0; currentRow < N; currentRow++) {
        int rowIndex = currentRow / NoOfPrc;
        int processIndex = currentRow % NoOfPrc;
        if (processID == processIndex) {
            // Loading to buffer
            memcpy(eliminationRow, localMatrix + rowIndex * columns, columns * sizeof(float));
            // Broadcasting
            MPI_Bcast(localMatrix + rowIndex * columns, columns, MPI_FLOAT, processIndex, MPI_COMM_WORLD);
            for (int i = rowIndex + 1; i < localRowSize; i++) {
                float multiplier = localMatrix[i * columns + currentRow]/localMatrix[rowIndex * columns + currentRow];
                 for (int j = currentRow; j < columns; j++)
                    localMatrix[i * columns + j] -= multiplier * eliminationRow[j];
            }
        } else {
            MPI_Bcast(eliminationRow, columns, MPI_FLOAT, processIndex, MPI_COMM_WORLD);
            for (int i = rowIndex; i < localRowSize; i++)
                if (processIndex < processID || i > rowIndex)
                {
                    float multiplier = localMatrix[i * columns + currentRow]/eliminationRow[currentRow];
                    for (int j = currentRow ; j < columns; j++)
                        localMatrix[i * columns + j] -= multiplier * eliminationRow[j];
                }
        }
    }
    // barrier to ensure all processors are complete with gauss elimination
    MPI_Barrier(MPI_COMM_WORLD);
    // Now we recompile the rows in a matrix
    for (int currentRow = 0; currentRow < N; currentRow++) {
        int rowIndex = currentRow / NoOfPrc;
        int rowProcessor = currentRow % NoOfPrc;
        MPI_Status status;
        // if row belongs to first prcsr copy it to memory
        if (processID == 0 && rowProcessor == 0)
            memcpy(AB + currentRow * columns, localMatrix + rowIndex * columns, columns * sizeof(float));
        else if (processID == 0)
            MPI_Recv(AB + currentRow * columns,columns, MPI_FLOAT,rowProcessor, currentRow, MPI_COMM_WORLD, &status);
        else if (processID == rowProcessor)
            MPI_Send(localMatrix + rowIndex * columns,columns, MPI_FLOAT,0, currentRow, MPI_COMM_WORLD);
        }
    if (processID == 0) {
        // Backsubstitution
        float X[N];
        // finding vector X
        for (int row = N -1; row >= 0; row--) {
            X[row] = AB[row * columns + columns - 1];
            for (int col = N - 1; col > row; col--)
                X[row] = X[row] - AB[row * columns + col] * X[col];
            X[row] = X[row] / AB[row * columns + row];
        }
        double currentTime = MPI_Wtime();
        double time = currentTime - startTime;
        printf("elapsed time: %f seconds", time);

        if(N<=5){
            printf("\nProcessed Matrix:\n\t");
            print2DMatrix(AB, N);
            printXvector(X, N);
        }
    }
    MPI_Finalize();
}

