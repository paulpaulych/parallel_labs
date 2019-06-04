#include <stdio.h>
#include <mpi.h>
#include <algorithm>

const int N1 = 1200;
const int N2 = 1200;
const int N3 = 2560;

double * initA(){
    double * matr = new double[N1 * N2];
    for(int i = 0; i < N1; i++){
        for(int j = 0; j < N2; j++){
            matr[i * N2 + j] = (i == j) ? i : 0;
        }        
    }
    return matr;
}
 
double * initB(){
    double * matr = new double[N2 * N3];
    for(int i = 0; i < N2; ++i){
        for(int j = 0; j < N3; ++j){
            matr[i * N3 + j] = (i == j) ? 2 : 0;
        }
    }
    return matr;
}

bool checkAnswer(const double * answer){
    for(int i = 0; i < std::min(N1, N2); ++i){
        for(int j = 0; j < std::min(N2, N3); ++j){
            if(i == 1 && j == 1){
                printf("%f\n", answer[i * N3 + j]);
            }
            if( (answer[i*N3 + j] - ((i==j)?(i*2):0)) > 0.0001 ){
                printf("i = %d, j = %d\n", i, j);
                return false;
            }
        }
    }
    return true;
}

double *multiply(double *a, double *b, int n2, int n1, int n3) {
    double *res = new double [n1 * n3];
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n3; ++j){
            for(int k = 0; k < n2; ++k) {
                res[i * n3 + j] += a[i * n2 + k] * b[j * n2 + k];
            }
        }
    }
    return res;
}
 
void transpose(double *a, int h, int w) {
    double * tmp = new double[h * w];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            tmp[j*h + i] = a[i*w + j];
        }
    }
    std::copy(tmp, tmp + h * w, a);
    delete[] tmp;
}
 
int main(int argc, char **argv) {
    int rank, size, rows, cols;
    double* matrC;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    int dims[2] = { 0, 0 };
    MPI_Dims_create(size, 2, dims);
    rows = dims[0];
    cols = dims[1];
   
    if(rank == 0){
        printf("rows = %d, cols = %d\n", rows, cols);
    }
    
    if((N1 % rows != 0 || N3 % cols != 0) && rank == 0) {
        printf("Wrong parameters!\n");
        return 1;
    }
    
    MPI_Comm rowComm, colComm;
    int row_key = rank / cols;
    int col_key = rank % cols;

    MPI_Comm_split(MPI_COMM_WORLD, row_key, rank, &rowComm);
    MPI_Comm_split(MPI_COMM_WORLD, col_key, rank, &colComm);
    
    int rowsPerThreadA = N1 / rows;
    int colsPerThreadB = N3 / cols;
    
    double *partA = new double[rowsPerThreadA * N2];
    double *partB = new double[colsPerThreadB * N2];
    double *matrA, *matrB;

    if(rank == 0) {
        matrA = initA();
        matrB = initB();
        matrC = new double[N1 * N3];
        transpose(matrB, N2, N3);
    }

    double start = MPI_Wtime();
    if(col_key == 0){
        MPI_Scatter(matrA, rowsPerThreadA * N2, MPI_DOUBLE, partA, rowsPerThreadA * N2, MPI_DOUBLE, 0, colComm);
    }
    if(row_key == 0) {
        MPI_Scatter(matrB, colsPerThreadB * N2, MPI_DOUBLE, partB, colsPerThreadB * N2, MPI_DOUBLE, 0, rowComm);
    }

    MPI_Bcast(partA, rowsPerThreadA * N2, MPI_DOUBLE, 0, rowComm);
    MPI_Bcast(partB, colsPerThreadB * N2, MPI_DOUBLE, 0, colComm);

    double *partC = multiply(partA, partB, N2, rowsPerThreadA, colsPerThreadB);

    MPI_Datatype recvType;
    MPI_Type_vector(rowsPerThreadA, colsPerThreadB, N3, MPI_DOUBLE, &recvType);
    MPI_Type_commit(&recvType);
 
    if(rank == 0) {
        for(int rank = 1; rank < size; ++rank){
            MPI_Recv(matrC + (rank / cols) * rowsPerThreadA * N3 + (rank % cols) * colsPerThreadB, 1, recvType, rank, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(int i = 0; i < rowsPerThreadA * colsPerThreadB; i++) {
            matrC[i / colsPerThreadB * N3 + i % colsPerThreadB] = partC[i];
        }
    }else{
        MPI_Send(partC, rowsPerThreadA * colsPerThreadB, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
    }
    double finish = MPI_Wtime();
    if(rank == 0){
        printf("time taken: %lf\n", finish - start);        
        printf("correct answer: %s\n", checkAnswer(matrC)? "true": "false");
    }
    
    if(rank == 0){
        delete[] matrA;
        delete[] matrB;
        delete[] matrC;
    }
    delete[] partA;
    delete[] partB;
    delete[] partC;
    
    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    MPI_Finalize();
    return 0;
}