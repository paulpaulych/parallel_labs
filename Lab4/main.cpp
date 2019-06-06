#include <iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<sys/time.h>
#include <mpi.h>

namespace{
    const double a = 200;
    const double Ni = 30;
    const double Nj = 30;
    const double Nk = 30;
    const double eps = 10e-8;
 
    const double Dx = 2.0;
    const double Dy = 2.0;
    const double Dz = 2.0;
    
    double hx;
    double hy;
    double hz;
    
    double powy;
    double powz;
    double powx;
    
    double denominator;
    double *func[2];
    double *slice[2];
    
    int sizeX, sizeY, sizeZ;
    int tmpf, finishFlag;   

    MPI_Request sendRequest[2] = {};
    MPI_Request recvRequest[2] = {};
}

double Phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}
 
double Rho(double x, double y, double z) {
    return 6.0 - a * Phi(x, y, z);
}
 
void initData(int *lines, int *offset, int rank) {
    for (int i = 0, start = offset[rank]; i < sizeX; i++, start++) {
        for (int j = 0; j < sizeY; j++) {
            for (int k = 0; k < sizeZ; k++) {
                if ((start != 0) && (j != 0) && (k != 0) && (start != Ni) && (j != Nj) && (k != Nk)) {
                    func[0][i*sizeY*sizeZ + j * sizeZ + k] = 0.0;
                    func[1][i*sizeY*sizeZ + j * sizeZ + k] = 0.0;
                }
                else {
                    func[0][i * sizeY * sizeZ + j * sizeZ + k] = Phi(start * hx, j * hy, k * hz);
                    func[1][i * sizeY * sizeZ + j * sizeZ + k] = Phi(start * hx, j * hy, k * hz);
                }
            }
        }
    }
}
 
void calcuateCenter(int *lines, int *offsets, int rank, int fooFlag) {
    for (int i = 1; i < sizeX - 1; i++) {
        for (int j = 1; j < sizeY - 1; j++) {
            for (int k = 1; k < sizeZ - 1; k++) {
                double fi, fj, fk;
                fi = (func[fooFlag][(i + 1)*sizeY*sizeZ + j * sizeZ + k] + func[fooFlag][(i - 1)*sizeY*sizeZ + j * sizeZ + k]) / powx;
                fj = (func[fooFlag][i*sizeY*sizeZ + (j + 1)* sizeZ + k] + func[fooFlag][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[fooFlag][i*sizeY*sizeZ + j * sizeZ + (k + 1)] + func[fooFlag][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
 
                func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank])*hx, j*hy, k*hz)) / denominator;
 
                if (fabs(func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank])*hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
        }
    }
}
 
void calculateEdges(int *lines, int *offsets, int rank, int size, int fooFlag) {
    for (int j = 1; j < sizeY - 1; j++) {
        for (int k = 1; k < sizeZ - 1; k++) {
            if (rank != 0) {
                int i = 0;
                double fi, fj, fk;
                fi = (func[fooFlag][(i + 1)*sizeY*sizeZ + j * sizeZ + k] + slice[0][j*sizeZ + k]) / powx;
                fj = (func[fooFlag][i*sizeY*sizeZ + (j + 1)*sizeZ + k] + func[fooFlag][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[fooFlag][i*sizeY*sizeZ + j * sizeZ + (k + 1)] + func[fooFlag][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
                func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank]) * hx, j*hy, k*hz)) / denominator;
                if (fabs(func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] - Phi(offsets[rank] * hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
 
            if (rank != size - 1) {
                double fi, fj, fk;
                int i = lines[rank] - 1;
                fi = (slice[1][j*sizeZ + k] + func[fooFlag][(i - 1)*sizeY*sizeZ + j * sizeZ + k]) / powx;
                fj = (func[fooFlag][i*sizeY*sizeZ + (j + 1)*sizeZ + k] + func[fooFlag][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[fooFlag][i* sizeY*sizeZ + j * sizeZ + (k + 1)] + func[fooFlag][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
                func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank]) * hx, j*hy, k*hz)) / denominator;
                if (fabs(func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank]) * hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
        }
    }
}
 
void sendSlice(int rank, int size, int fooFlag, int *lines) {
    if (rank != 0) {
        MPI_Isend(&(func[fooFlag][0]), sizeZ*sizeY, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &sendRequest[0]);
        MPI_Irecv(slice[0], sizeZ * sizeY, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &recvRequest[1]);
    }
    if (rank != size - 1) {
        MPI_Isend(&(func[fooFlag][(lines[rank] - 1) * sizeY * sizeZ]), sizeZ * sizeY, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &sendRequest[1]);
        MPI_Irecv(slice[1], sizeZ*sizeY, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recvRequest[0]);
    }
}
 
void waitForSlice(int rank, int size) {
    if (rank != 0) {
        MPI_Wait(&recvRequest[1], MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequest[0], MPI_STATUS_IGNORE);
    }
    if (rank != size - 1) {
        MPI_Wait(&recvRequest[0], MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequest[1], MPI_STATUS_IGNORE);
    }
}
 
void maxDiff(int *lines, int *offsets, int rank, int fooFlag) {
    double max = 0.0;
    double tmp;
    double tmpmax = 0.0;
    for (int i = 1; i < sizeX - 1; ++i) {
        for (int j = 1; j < sizeY - 1; ++j) {
            for (int k = 1; k < sizeZ - 1; ++k){
                if ((tmp = fabs(func[1-fooFlag][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank]) * hx, j*hy, k *hz))) > max){
                    max = tmp;
                }
            }        
        }
    }
    MPI_Allreduce(&max, &tmpmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) {
        max = tmpmax;
        std::cout << "Max difference = " << max << std::endl;
    }
}
 
int main(int argc, char **argv) {
    int size, rank;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    if(size > Ni){
        if(rank == 0){
            fprintf(stderr, "Too large number of processes\n");
        }
        return EXIT_FAILURE;
    }

    if((Ni < 2)||(Nj < 2)||(Nk < 2)){
        if(rank == 0){
            fprintf(stderr, "Wrong sizes\n");
        }
        return EXIT_FAILURE;    
    }

    int * lines = new int[size];
    int * offsets = new int[size];
    int currentLine = 0;
    int height = Ni + 1;
    int tmp = size - (height % size);
    for (int i = 0, currentLine = 0; i < size; i++) {
        lines[i] = (i < tmp) ? (height / size) : (height / size + 1);
        offsets[i] = currentLine;
        currentLine += lines[i];
    }
 
    sizeX = lines[rank];
    sizeY = Nj + 1;
    sizeZ = Nk + 1;

    try{
        slice[0] = new double[sizeZ*sizeY];
        slice[1] = new double[sizeZ*sizeY];
        func[0] = new double[sizeX*sizeY*sizeZ];
        func[1] = new double[sizeX*sizeY*sizeZ];
    }catch(std::exception &exc){
        std::cout << exc.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    hx = Dx / Ni;
    hy = Dy / Nj;
    hz = Dz / Nk;
 
    powx = hx * hx;
    powy = hy * hy;
    powz = hz * hz;
 
    denominator = 2 / powx + 2 / powy + 2 / powz + a;
 
    initData(lines, offsets, rank);
 
    int  fooFlag = 1;
    double startTime = MPI_Wtime();
    int iterations = 0;
    do {
        if (rank == 0) {
            iterations++;
        }
        finishFlag = 1;
        fooFlag = 1 - fooFlag;
    
        sendSlice(rank, size, fooFlag, lines);
        calcuateCenter(lines, offsets, rank, fooFlag);
        waitForSlice(rank, size);
        calculateEdges(lines, offsets, rank, size, fooFlag);
 
        MPI_Allreduce(&finishFlag, &tmpf, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        finishFlag = tmpf;

    } while (!finishFlag);
 
    maxDiff(lines, offsets, rank, fooFlag);
 
    if (rank == 0) {
        std::cout << "Time taken: " << MPI_Wtime() - startTime << std::endl;
        std::cout << "iterations: " << iterations << std::endl;
    }
 
    delete[] slice[0];
    delete[] slice[1];
    delete[] func[0];
    delete[] func[1];
    delete[] offsets;
    delete[] lines;
 
    MPI_Finalize();
 
    return EXIT_SUCCESS;
}