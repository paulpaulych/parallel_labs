#include <cmath>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <mpi/mpi.h>
#include <queue>
#include <random>
#include <unistd.h>
#include "tpool.h"
#include <thread>

int main(int argc, char ** argv) {
    ThreadPool tp(5);

    for(int i = 0; i < 5; i++){
        tp.push(Task{1, 0});
    }

    sleep(5);
    printf("waked up\n");
    printf("empty: %s\n", tp.empty() ? "true" : "false");

    for(int i = 0; i < 5; i++){
        tp.push(Task(1, 0));
    }
    tp.stop();

    printf("asd\n");

    return EXIT_SUCCESS;
}
