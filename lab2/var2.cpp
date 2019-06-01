#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include <stdlib.h>
 
#define N 3264
 
const double EPSILON = 1e-9;
const double TAU = 1e-5;

double * initMatrix() {
    double * A = new double[N*N];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[N*i + j] = (i == j) ? 2.0 : 1.0;
        }
    }
    return A;
}
 
double* initVector(double initValue) {
    double* vector = new double[N];
 
    for (int i = 0; i < N; ++i) {
        vector[i] = initValue;
    }
 
    return vector;
}
 
double getSum(const double* vector) {
    static double result = 0;
    result  = 0;
 
    #pragma omp for reduction(+:result)  
    for (int i = 0; i < N; ++i) {
        result += vector[i] * vector[i];
    }
 
    return result;
}
 
void subtractVectors(double* leftVector, const double* rightVector) {
 
    #pragma omp for  
    for (int i = 0; i < N; ++i) {
        leftVector[i] -= rightVector[i];
    }
}
 
void multiply(const double* matrix, const double* vector, double* result) {
    double tmpSum = 0;
 
    #pragma omp for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            tmpSum += matrix[N * i + j] * vector[j];
        }
        result[i] = tmpSum;
        tmpSum = 0;
    }
}
 
void multiplyWithConstant(double* vector, const double value) {
 
    #pragma omp for
    for (int i = 0; i < N; ++i) {
        vector[i] = vector[i] *  value;
    }
}
 
void simpleIteration(double *A, double * b, double * x, int threadsNumber) {
    omp_set_num_threads(threadsNumber);
 
    double resultVector[N];
    double sum = 0;
    int flag = 0;
 
    double start = omp_get_wtime();
    int bNorm = sqrt(getSum(b));
    double  result = 0;
 
    #pragma omp parallel
    {
        do{
            multiply(A, x, resultVector);
            subtractVectors(resultVector, b);
            double value = 0;
           
            #pragma omp for reduction(+:result)
                for (int i = 0; i < N; ++i) {
                    result += resultVector[i] * resultVector[i];
                }
               
            sum = getSum(resultVector);
 
            #pragma omp single
            {
                flag = (sqrt(sum)/bNorm > EPSILON) ? 1 : 0;
            }
            multiplyWithConstant(resultVector, TAU);
            subtractVectors(x, resultVector);
        } while(flag);
    }
    double end = omp_get_wtime();
    printf("TIME TAKEN: %lf\n", end - start);
}
 
 
int main(int argc ,char* argv[]) {
 
    if(argc != 2){
        printf("wrong args number");
        return 1;
    }

    double * A = initMatrix();
    double * b = initVector(N + 1);
    double * x = initVector(0);

    simpleIteration(A, b, x, atoi(argv[1]));
 
    return 0;
}