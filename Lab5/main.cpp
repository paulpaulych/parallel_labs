#include <cmath>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <mpi/mpi.h>
#include <queue>
#include <unistd.h>

using std::queue;

enum Option{ A, B, C };

const int ITERATIONS_CNT = 2;//M
const int TASKS_CNT = 20;//N
const int RHO = 2;
const int GEN_OPTION = Option::B;

const int REQUEST_TAG = 12;
const int ANSWER_TAG = 13;
const int DATA_TAG = 14;

pthread_t dataThreadID;
pthread_mutex_t taskMutex;
queue<int> tasks;
int rank = 0, size = 0, iteration = 0, task = 0;//number of current task
//long long int sleepTime = 0;

class ThreadPool{

private:
    queue<int> tasks;

public:


};

int newThread();
void * dataThread(void *);
void * taskThread(void *);
int calculateCount(int);
void generateTask(int);
void doTask(int&);
bool getTask(int);

int main() {
    int initThreadsFlag = 0;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &initThreadsFlag);

    if (initThreadsFlag != MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Gan't start MPI with THREAD_MULTIPLE \n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(nullptr));
    pthread_mutex_init(&taskMutex, nullptr);

    double start = MPI_Wtime();

    if (newThread() != EXIT_SUCCESS) {
        pthread_mutex_destroy(&taskMutex);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    double end = MPI_Wtime();
    printf("[MAIN thread from #%d process]: Time is %lf; \n", rank, end - start);
    //printf("[MAIN thread from #%d process]: Count operations is %lld \n", rank, sleepTime);

    pthread_mutex_destroy(&taskMutex);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

void * taskThread(void *args) {
    bool * neighbours = new bool[size]();

    //pthread_mutex_lock(&taskMutex);
    for (iteration = 0; iteration < ITERATIONS_CNT; ++iteration) {
        //pthread_mutex_unlock(&taskMutex);

        //printf("[Task thread from #%d process]: Start %d iteration \n" , rank, iteration);

        //pthread_mutex_lock(&taskMutex);
        int count = calculateCount(GEN_OPTION);
        if (count < 0) {
            printf("[Task thread from #%d process]: ERROR while task generating! \n", rank);
            iteration = ITERATIONS_CNT;
            continue;
        }
        generateTask(count);
        //pthread_mutex_unlock(&taskMutex);

        for (int i = 0; i < size; ++i) {
            neighbours[i] = true;
        }
        neighbours[rank] = false;

        //pthread_mutex_lock(&taskMutex);
        task = 0;
        for (int i = 1; i < size;) {
            pthread_mutex_lock(&taskMutex);
            while (!tasks.empty()) {
                int tmp = tasks.front();
                ++task;
                tasks.pop();

                pthread_mutex_unlock(&taskMutex);
                doTask(tmp);
                printf("[Task thread from #%d process]: Has done task #%d (%lu remained) ... \n",
                       rank, task, tasks.size());

                pthread_mutex_lock(&taskMutex);
            }
            pthread_mutex_unlock(&taskMutex);
            //FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            int neighbourRank = (rank + i) % size;
//            if (neighbours[neighbourRank] == false) {
//                ++i;
//                continue;
//            }

            neighbours[neighbourRank] = getTask(neighbourRank);
            i+=(neighbours[neighbourRank] == true ? 0 : 1);

//            pthread_mutex_lock(&taskMutex);
        }
//        pthread_mutex_unlock(&taskMutex);

        printf("[Task thread from #%d process]: End of main iteration cycle #%d. \n", rank, iteration);

        MPI_Barrier(MPI_COMM_WORLD);
        //pthread_mutex_lock(&taskMutex);
    }
    //pthread_mutex_unlock(&taskMutex);

    int req = 0;
    MPI_Send(&req, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);
    delete[] neighbours;
    return nullptr;
}

void * dataThread(void *args) {
    //pthread_mutex_lock(&taskMutex);
    while (true) {
        //pthread_mutex_unlock(&taskMutex);
        MPI_Status status;
        int receive;
        MPI_Recv(&receive, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD,&status);
        if (receive == 0) break;
        pthread_mutex_lock(&taskMutex);
        printf("[Data thread from #%d process]: Got request from %d process \n"
            ,rank, status.MPI_SOURCE);

        if (tasks.empty()) {
            int answer = 0;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            pthread_mutex_unlock(&taskMutex);
            printf("[Data thread from #%d process]: Haven't data to send to %d process \n"
                ,rank, status.MPI_SOURCE);
            //pthread_mutex_lock(&taskMutex);
            continue;
        }
        printf("[Data thread from #%d process]: Send data to %d process \n",rank, status.MPI_SOURCE);

        int answer = 1, send = tasks.front();
        MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
        MPI_Send(&send, 1, MPI_INT, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD);
        tasks.pop();
        pthread_mutex_unlock(&taskMutex);
    }
    //pthread_mutex_unlock(&taskMutex);
    printf("[Data thread from #%d process]: The end of work \n", rank);
    return nullptr;
}

int calculateCount(int option) {
    int ret = 0;
    switch (option) {
        case Option::A:
            ret = TASKS_CNT / size;
            break;
        case Option::B:
            ret = TASKS_CNT / (2 * size) *
                    std::min( (rank + iteration) % size, (2 * size - rank - iteration % size) % size );
            break;
        case Option::C:
            ret = (rank == iteration % size ? TASKS_CNT : 0);
            break;
        default:
            ret = -1;
            break;
    }
    return ret;
}

void generateTask(int count) {
    printf("[Task thread from #%d process]: Generation %d task\n"
            ,rank, count);
    //clearing task queue
    while (!tasks.empty()) tasks.pop();
    //if (count <= 0) return;

    int weight = RHO * (1 + abs(rank - (iteration % size)));

    for (int j = 0; j < count; ++j) {
        int task = rank + (1 + rand() % weight);
        tasks.push(task);
    }
}

void doTask(int & weight) {
    sleep(weight);
    return;
}

bool getTask(int neighbourRank) {
    printf("[Task thread from #%d process]: Trying to get data from %d process \n"
        ,rank, neighbourRank);

    int request = 1;
    MPI_Send(&request, 1, MPI_INT, neighbourRank, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&request, 1, MPI_INT, neighbourRank, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(request == 0) {
        printf("[Task thread from #%d process]: Have no task \n"
                ,rank);
        return false;
    }

    int task;
    MPI_Recv(&task, 1, MPI_INT, neighbourRank, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    pthread_mutex_lock(&taskMutex);
    tasks.push(task);
    pthread_mutex_unlock(&taskMutex);
    printf("[Task thread from #%d process]: New task(%d) got from %d process\n"
            ,rank, task, neighbourRank);

    return true;
}

int newThread() {
    pthread_attr_t attr;

    if (pthread_attr_init(&attr) != EXIT_SUCCESS) {
        fprintf(stderr, "Cannot initialize attributes");
        return EXIT_FAILURE;
    };

    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) != EXIT_SUCCESS) {
        fprintf(stderr, "Error in setting attributes");
        return EXIT_FAILURE;
    }

    if (pthread_create(&dataThreadID, &attr, dataThread, nullptr) != EXIT_SUCCESS) {
        fprintf(stderr, "Cannot create thread");
        return EXIT_FAILURE;
    }

    pthread_attr_destroy(&attr);

    taskThread(nullptr);

    if (pthread_join(dataThreadID, nullptr) != EXIT_SUCCESS) {
        fprintf(stderr, "Cannot join a thread");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}