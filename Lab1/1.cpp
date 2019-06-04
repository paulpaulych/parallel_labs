#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

const int N = 10;
const double t = 10e-6;
const double eps = 10e-9;

bool answerIsCorrect(double * answer){
    for(int i = 0; i < N; ++i) {
        if(fabs(fabs(answer[i]) - 1) >= eps) {
            return false;    
        }
    }
    return true;
}

int initPerThread(int * perThread, int size, int rank){
    int startLine = 0;
    for(int i = 0, tmp = size - (N % size); i < size; ++i) {
        perThread[i] = i < tmp ? (N / size) : (N / size + 1);
        if(i < rank) {
            startLine += perThread[i];
        }
    }
    return startLine;
}

int initData(int *perThread, double *matrix, int startLine, double *b, double *x, int rank) {
    for(int i = 0; i < perThread[rank]; ++i) {
        for(int j = 0; j < N; ++j) {
            matrix[i * N + j] = (startLine + i) == j ? 2 : 1;
        }
        b[i] = N + 1;
        x[i] = 0;
    }
    for(int i = 0; i < N; ++i) {
        b[i] = N + 1;
        x[i] = 0;
    }
}

int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(size > N && rank == 0){
        std::cout << "Incorrect number of processes" << std::endl;
        return 0;
    }

    int *perThread = new int[size];
    int startLine = initPerThread(perThread, size, rank);    

    double *matrix = new double[perThread[rank] * N];
    double *b = new double[N];
    double *x = new double[N];
    
    initData(perThread, matrix, startLine, b, x, rank);
    
    double normB = 0;
    double start = 0;
    
    if(rank == 0) {
        start = MPI_Wtime();
        for(int i = 0; i < N; ++i) {
            normB += b[i] * b[i];
        }
        normB = sqrt(normB);
    }
    
    double *partX = new double[perThread[rank]];
    
    int iterations = 0;
    int isReady = 0;
    while(!isReady) {
        iterations++;
        if(rank == 0){
            printf("it %d\n", iterations);
        }
        double partAnswer = 0;
        for(int i = 0; i < perThread[rank]; ++i) {
            double sum = 0;
            for(int j = 0; j < N; ++j) {
                sum += matrix[i * N + j] * x[j];
            }
            sum -= b[i];
            partX[i] = x[i + startLine] - sum * t;
            partAnswer += sum * sum;
        }
        if(rank == 0) {
            for(int i = startLine, c = startLine + perThread[rank]; i < c; ++i) {
                x[i] = partX[i - startLine];
            }
            double sum = partAnswer;
            int currentLine = perThread[rank];
            for(int i  = 1; i < size; ++i) {
                MPI_Recv(x + currentLine, perThread[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                currentLine += perThread[i];

                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sum += tmp;
            }
            sum = sqrt(sum);
            isReady = (sum / normB) <= eps;
            printf("%.10f\n", sum / normB);
        
        }else{    
            MPI_Send(partX, perThread[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&partAnswer, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }
        MPI_Bcast(&isReady, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if(rank == 0) {
        double finish = MPI_Wtime();
        std::cout << "Answer correct: " << answerIsCorrect(x)<< std::endl;
        std::cout << "Time taken: " << finish - start << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
    }

    delete[] partX;
    delete[] x;
    delete[] b;
    delete[] matrix;
    delete[] perThread;

    MPI_Finalize();
    return 0;
}