#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

const int N = 1000;
const double t = 2e-6;
const double eps = 10e-9;

bool answerIsCorrect(double * answer){
    for(int i = 0; i < N; ++i) {
        if(fabs(fabs(answer[i]) - 1) >= eps) {
            return false;    
        }
    }
    return true;
}

int initPerThread(int * perThread, int * offsets, int size, int rank){
    int startLine = 0;
    offsets[0] = 0;
    for(int i = 0, tmp = size - (N % size); i < size; ++i) {
        perThread[i] = i < tmp ? (N / size) : (N / size + 1);
        if(i < rank) {
            startLine += perThread[i];
        }
    }
    for(int i = 1; i < size; ++i){
        offsets[i] = offsets[i-1]+perThread[i-1];
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
    int *offsets = new int[size];
    int startLine = initPerThread(perThread, offsets, size, rank);

    if(rank == 0) {
        printf("perThread\n");
        for(int i = 0; i < size; ++i) {
            printf("%d ", perThread[i]);
        }
        printf("\n");
    }

    double *matrix = new double[perThread[rank] * N];
    double *b = new double[perThread[rank]];
    double *x = new double[perThread[rank]];
    
    initData(perThread, matrix, startLine, b, x, rank);

    double start = 0;
    double normB = 0;
    if(rank == 0) {
        start = MPI_Wtime();
        for(int i = 0; i < perThread[rank]; ++i) {
            normB += b[i] * b[i];
        }
        for(int i = 1; i < size; ++i) {
            double tmp;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //получаем столбец
            normB += tmp;
        }
        normB = sqrt(normB);
    } else {
        for(int i = 0; i < perThread[rank]; ++i) {
            normB += b[i] * b[i];
        }
        MPI_Send(&normB, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD); //отправляем b        
    }

    double *partSum = new double[perThread[rank]];
    double *tmpXtoSend = new double[N / size + 1];
    double *tmpXtoRecv = new double[N / size + 1];
    int iterations = 0;
    int isReady = 0;
    while(!isReady) {
        iterations++;

        for(int i = 0; i < perThread[rank]; ++i) {
            partSum[i] = 0;
            tmpXtoSend[i] = x[i];
        }

        for(int i = 0; i < size; ++i) {
            int offset = offsets[(rank + i + size) % size];
            for(int j = 0; j < perThread[rank]; ++j) {
                for(int k = 0; k < perThread[(rank + i + size) % size]; ++k) {
                    partSum[j] += matrix[j * N + offset + k] * tmpXtoSend[k];
                }
            }
            MPI_Sendrecv(tmpXtoSend, N / size + 1, MPI_DOUBLE, (rank - 1 + size) % size, 0,
            tmpXtoRecv, N / size + 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::swap(tmpXtoSend, tmpXtoRecv);
        }

        double partAnswer = 0;

        for(int i = 0; i < perThread[rank]; ++i) {
            partSum[i] -= b[i];
            x[i] = x[i] - partSum[i] * t;
            partAnswer += partSum[i] * partSum[i];
        }

        if(rank == 0) {
            double sum = partAnswer;
            for(int i  = 1; i < size; ++i) {
                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sum += tmp;
            }
            sum = sqrt(sum);
            isReady = sum / normB < eps;
        } else {
            MPI_Send(&partAnswer, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        }
        MPI_Bcast(&isReady, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if(rank == 0) {
        double *fullX = new double[N];
        for(int i = 0; i < perThread[rank]; ++i) {
            fullX[i] = x[i];
        }

        for(int i = 1, currentLine = perThread[rank]; i < size; ++i) {
            MPI_Recv(fullX+currentLine, perThread[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            currentLine += perThread[i];
        }
    
        double finish = MPI_Wtime();
        std::cout << "Answer correct: " << answerIsCorrect(fullX) << std::endl;
        std::cout << "Time taken: " << finish - start << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
        
        delete[] fullX;
    } else {
        MPI_Send(x, perThread[rank], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    delete[] x;
    delete[] b;
    delete[] matrix;
    delete[] perThread;
    
    MPI_Finalize();
    return 0;
}