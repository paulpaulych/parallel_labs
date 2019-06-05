#include <iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<sys/time.h>
#include <mpi.h>

namespace{
    const double a = 200;
    const double Ni = 200;
    const double Nj = 200;
    const double Nk = 200;
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
    double *buffer[2];
    
    int sizeX, sizeY, sizeZ;
    int tmpf, finishFlag;   

    MPI_Request sendRequest[2] = {};
    MPI_Request recRequest[2] = {};
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
 
void calcuateCenter(int *lines, int *offsets, int rank, int L0, int L1) {
    for (int i = 1; i < sizeX - 1; i++) {
        for (int j = 1; j < sizeY - 1; j++) {
            for (int k = 1; k < sizeZ - 1; k++) {
                double fi, fj, fk;
                fi = (func[L0][(i + 1)*sizeY*sizeZ + j * sizeZ + k] + func[L0][(i - 1)*sizeY*sizeZ + j * sizeZ + k]) / powx;
                fj = (func[L0][i*sizeY*sizeZ + (j + 1)* sizeZ + k] + func[L0][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[L0][i*sizeY*sizeZ + j * sizeZ + (k + 1)] + func[L0][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
 
                func[L1][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank])*hx, j*hy, k*hz)) / denominator;
 
                if (fabs(func[L1][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank])*hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
        }
    }
}
 
void calculateEdges(int *lines, int *offsets, int rank, int size, int L0, int L1) {
    for (int j = 1; j < sizeY - 1; j++) {
        for (int k = 1; k < sizeZ - 1; k++) {
            if (rank != 0) {
                int i = 0;
                double fi, fj, fk;
                fi = (func[L0][(i + 1)*sizeY*sizeZ + j * sizeZ + k] + buffer[0][j*sizeZ + k]) / powx;
                fj = (func[L0][i*sizeY*sizeZ + (j + 1)*sizeZ + k] + func[L0][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[L0][i*sizeY*sizeZ + j * sizeZ + (k + 1)] + func[L0][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
                func[L1][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank]) * hx, j*hy, k*hz)) / denominator;
                if (fabs(func[L1][i*sizeY*sizeZ + j * sizeZ + k] - Phi(offsets[rank] * hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
 
            if (rank != size - 1) {
                double fi, fj, fk;
                int i = lines[rank] - 1;
                fi = (buffer[1][j*sizeZ + k] + func[L0][(i - 1)*sizeY*sizeZ + j * sizeZ + k]) / powx;
                fj = (func[L0][i*sizeY*sizeZ + (j + 1)*sizeZ + k] + func[L0][i*sizeY*sizeZ + (j - 1)*sizeZ + k]) / powy;
                fk = (func[L0][i* sizeY*sizeZ + j * sizeZ + (k + 1)] + func[L0][i*sizeY*sizeZ + j * sizeZ + (k - 1)]) / powz;
                func[L1][i*sizeY*sizeZ + j * sizeZ + k] = (fi + fj + fk - Rho((i + offsets[rank]) * hx, j*hy, k*hz)) / denominator;
                if (fabs(func[L1][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank]) * hx, j*hy, k*hz)) > eps){
                    finishFlag = 0;
                }
            }
        }
    }
}
 
void sendShadows(int rank, int size, int L0, int *lines) {
    if (rank != 0) {
        MPI_Isend(&(func[L0][0]), sizeZ*sizeY, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &sendRequest[0]);
        MPI_Irecv(buffer[0], sizeZ * sizeY, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &recRequest[1]);
    }
    if (rank != size - 1) {
        MPI_Isend(&(func[L0][(lines[rank] - 1) * sizeY * sizeZ]), sizeZ * sizeY, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &sendRequest[1]);
        MPI_Irecv(buffer[1], sizeZ*sizeY, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recRequest[0]);
    }
}
 
void waitForShadows(int rank, int size) {
    if (rank != 0) {
        MPI_Wait(&recRequest[1], MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequest[0], MPI_STATUS_IGNORE);
    }
    if (rank != size - 1) {
        MPI_Wait(&recRequest[0], MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequest[1], MPI_STATUS_IGNORE);
    }
}
 
void maxDiff(int *lines, int *offsets, int rank, int L0, int L1) {
    double max = 0.0;
    double tmp;
    double tmpmax = 0.0;
    for (int i = 1; i < sizeX - 1; ++i) {
        for (int j = 1; j < sizeY - 1; ++j) {
            for (int k = 1; k < sizeZ - 1; ++k){
                if ((tmp = fabs(func[L1][i*sizeY*sizeZ + j * sizeZ + k] - Phi((i + offsets[rank]) * hx, j*hy, k *hz))) > max){
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
        buffer[0] = new double[sizeZ*sizeY];
        buffer[1] = new double[sizeZ*sizeY];
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
 
    int  L0 = 1, L1 = 0;
    double startTime = MPI_Wtime();
    int iterations = 0;
    do {
        if (rank == 0) {
            iterations++;
        }
        finishFlag = 1;
        L0 = 1 - L0;
        L1 = 1 - L1;
    
        sendShadows(rank, size, L0, lines);
        calcuateCenter(lines, offsets, rank, L0, L1);
        waitForShadows(rank, size);
        calculateEdges(lines, offsets, rank, size, L0, L1);
 
        MPI_Allreduce(&finishFlag, &tmpf, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        finishFlag = tmpf;

    } while (!finishFlag);
 
    maxDiff(lines, offsets, rank, L0, L1);
 
    if (rank == 0) {
        std::cout << "Time taken: " << MPI_Wtime() - startTime << std::endl;
        std::cout << "iterations: " << iterations << std::endl;
    }
 
    delete[] buffer[0];
    delete[] buffer[1];
    delete[] func[0];
    delete[] func[1];
    delete[] offsets;
    delete[] lines;
 
    MPI_Finalize();
 
    return EXIT_SUCCESS;
}