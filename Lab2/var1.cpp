#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include <stdlib.h>
 
#define N 3264
 
const double EPSILON = 1e-9;
const double TAU = 1e-5;

double * initMatrix() {
    double * A = new double[N * N];
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[N*i + j] = (i == j) ? 2.0 : 1.0 ;
        }
    }

    return A;
}
 
double* initVector(double initValue) {
    double* vector = new double[N];
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        vector[i] = initValue;
    }
    
    return vector;
}
 
double getNorm(const double* vector) {
    double result = 0;
 
    #pragma omp parallel for reduction(+:result)  
    
    for (int i = 0; i < N; ++i) {
        result += vector[i] * vector[i];
    }
 
    return sqrt(result);
}
 
void subVectors(double* a, const double* b) {
 
    #pragma omp parallel for  
    for (int i = 0; i < N; ++i) {
        a[i] -= b[i];
    }

}
 
void multiply(const double* matrix, const double* vector, double* result) {
    double tmpSum;
 
    #pragma omp parallel for private(tmpSum)
    for (int i = 0; i < N; ++i) {
        tmpSum = 0;
        for (int j = 0; j < N; ++j) {
            tmpSum += matrix[N * i + j] * vector[j];
        }
        result[i] = tmpSum;
    }
   
}
 
void multiplyWithConstant(double* vector, const double value) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        vector[i] = vector[i] *  value;
    }
}
 
void simpleIteration(double * A, double * b, double * x, int count) {
    omp_set_num_threads(count);
 
    double result[N];
    double tmp = 0;
 
    double start = omp_get_wtime();
    int normB = getNorm(b);
 
    do{
        multiply(A, x, result);
        subVectors(result, b);
        tmp = getNorm(result) / normB;
        printf("%lf\n", tmp);
        multiplyWithConstant(result, TAU);
        subVectors(x, result);
    } while(tmp > EPSILON);
 
    double end = omp_get_wtime();
    printf("TIME TAKEN: %lf\n", end - start);
   
}

int main(int argc ,char* argv[]) {
 
    double * A = initMatrix();
    double * b = initVector(N + 1);
    double * x = initVector(0);
    simpleIteration(A, b, x, atoi(argv[1]));
    

    return 0;
}
