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

enum Option{ A, B, C };

const int ITERATIONS_CNT = 2;//M
const int TASKS_CNT = 20;//N
const int GEN_OPTION = Option::B;

const int THREADS_PER_PROCESS = 1;

const int REQUEST_TAG = 12;
const int ANSWER_TAG = 13;
const int TASK_TAG = 14;

const bool FINISH = 0;
const bool NOT_FINISH = 1;
const int I_AM_DONE = 2;
const int TASK_REQUEST = 3;
const int TRUE_RESPONCE = 4;
const int FALSE_RESPONCE = 5;

int rank = 0, size = 0;

std::atomic_bool isWorking;
std::mutex workingM;
std::condition_variable woekingCV;
bool * finishedProcesses;

bool checkFinish(){
    for(int i = 0; i < size; i++){
        if(!finishedProcesses[i]){
            return false;
        }
    }
    return true;
}

void notifyAboutFinish(){
    for(int i = 1; i < size; i++){
        int message = FINISH;
        MPI_Send(&message, 1, MPI_INT, i, REQUEST_TAG, MPI_COMM_WORLD);
        isWorking = false;
    }
}

class TaskStealer{
private:
    ThreadPool * tp_;
    std::thread thread_;
public:
    TaskStealer(ThreadPool * tp){
        tp_ = tp;
        thread_ = std::thread(std::bind(&TaskStealer::run, this));
    }
    TaskStealer(){
        thread_ = std::thread(std::bind(&TaskStealer::run, this));
    }
    
    ~TaskStealer(){
        thread_.join();
    }

    void run(){
        while(isWorking){
            printf("[TaskStealer from #%d process]: Waiting...\n", rank);
            tp_->waitUntilEmpty();
            printf("[TaskStealer from #%d process]: Woke up.\n", rank);
            int i;
            for(i = (rank + 1)%size; i%size != rank; i = (i+1)%size){
                if(stealTask(i)){
                    break;
                }
            }
            if(i % size == rank){
                if(rank == 0){
                    finishedProcesses[0] = true;
                    if(checkFinish()){
                        printf("[MessageHandler from #%d process]: Notifiing.. \n", rank);
                        notifyAboutFinish();
                        printf("[MessageHandler from #%d process]: Notified.. \n", rank);
                        break;
                    }
                }else{
                    int message = I_AM_DONE;
                    printf("[TaskStealer from #%d process]: Sending I_AM_DONE\n", rank);
                    MPI_Send(&message, 1, MPI_INT, 0, REQUEST_TAG, MPI_COMM_WORLD);
                }
                printf("[TaskStealer from #%d process]: I am done\n", rank);
                return;
            }
        }
    }

private:
    bool stealTask(int neighbourRank){
        printf("[TaskStealer from #%d process]: Trying to get data from %d process \n", rank, neighbourRank);
        int message = TASK_REQUEST;
        MPI_Status status;
        MPI_Send(&message, 1, MPI_INT, neighbourRank, REQUEST_TAG, MPI_COMM_WORLD);
        printf("[TaskStealer from #%d process]: waiting for answer from process %d\n", rank, neighbourRank);
        int recieve = 0;
        MPI_Recv(&recieve, 1, MPI_INT, neighbourRank, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("[TaskStealer from #%d process]: Got answer %d\n", rank, neighbourRank);
        if(recieve == FALSE_RESPONCE) {
            printf("[TaskStealer from #%d process]: process %d has NOT tasks \n", rank, neighbourRank);
            return false;
        }
        if(recieve == TRUE_RESPONCE) {
            printf("[TaskStealer from #%d process]: process %d has tasks!!!\n", rank, neighbourRank);
            int weight;
            MPI_Recv(&weight, 1, MPI_INT, neighbourRank, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            tp_->push(Task(weight, rank));
            printf("[Task thread from #%d process]: New task got from %d process\n", rank, neighbourRank);    
            return true;
        }
        printf("[TaskStealer from #%d process]: got WRONG MESSAGE from %d process.\n", rank, neighbourRank);
        return true;
    }
};

class MessageHandler{
private:
    ThreadPool * tp_;
    std::thread thread_;
public:
    MessageHandler(ThreadPool * tp){
        tp_ = tp;
        thread_ = std::thread(std::bind(&MessageHandler::run, this));
    }

    ~MessageHandler(){
       thread_.join();
    }

    void run(){
        while (isWorking) {
            handleMessage();
            if(rank == 0){
                for(int i = 0; i < size; ++i){
                    printf("%d ", finishedProcesses[i]);
                }
                printf("\n");
            }
            if(isWorking &&(rank == 0) && checkFinish()){
                printf("[MessageHandler from #%d process]: Notifiing.. \n", rank);
                notifyAboutFinish();
                printf("[MessageHandler from #%d process]: Notified.. \n", rank);
                break;
            }
        }
        printf("[MessageHandler from #%d process]: The end of work \n", rank);
        tp_->stop();
    }

    void handleMessage(){
        MPI_Status status;
        int receive;
        //printf("[MessageHandler from #%d process]: Waiting...\n", rank);
        MPI_Recv(&receive, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        if((rank == 0) && (receive == I_AM_DONE)){
            printf("[MessageHandler from #%d process]: Got I_AM_DONE from %d process\n", rank, status.MPI_SOURCE);
            finishedProcesses[status.MPI_SOURCE] = true;
            return;
        }
        if(receive == FINISH){
            if(rank == size - 1){
                int message = FINISH;
                MPI_Send(&message, 1, MPI_INT, 0, REQUEST_TAG, MPI_COMM_WORLD);
            }
            printf("[MessageHandler from #%d process]: Got FINISH from %d process\n", rank, status.MPI_SOURCE);
            isWorking = false;
            return;
        }
        if(receive == TASK_REQUEST){
            if (tp_->empty()) {
                int answer = FALSE_RESPONCE;
                //printf("[MessageHandler from #%d process]: Haven't tasks to send to %d process \n", rank, status.MPI_SOURCE);
                MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                //printf("[MessageHandler from #%d process]: Sent answer to %d process...\n",rank, status.MPI_SOURCE);   
                return;
            }
            int answer = TRUE_RESPONCE;
            //printf("[MessageHandler from #%d process]: Sending answer to %d process...\n",rank, status.MPI_SOURCE);   
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            //printf("[MessageHandler from #%d process]: Sent answer to %d process...\n",rank, status.MPI_SOURCE);   
            printf("[MessageHandler thread from #%d process]: Sending task to %d process...\n",rank, status.MPI_SOURCE);
            int weight = tp_->popTask();
            MPI_Send(&weight, 1, MPI_INT, status.MPI_SOURCE, TASK_TAG, MPI_COMM_WORLD);
            return;
        }
        printf("[MessageHandler from #%d process]: Got WRONG MESSAGE from %d process...\n",rank, status.MPI_SOURCE);       
        return;
    }
};

int calculateTasksNum(int option);

void calculate(int tasksNum){
    ThreadPool threadPool(THREADS_PER_PROCESS); 
    isWorking = true;
    printf("rank : %d\n", rank); 
    
    for(int i = 0; i < tasksNum; ++i){
        std::random_device rd;
        int weight = (rank+1) * 3;
        threadPool.push(Task(weight, rank));
    }
    if(size > 1){
        TaskStealer taskStealer(&threadPool);
        MessageHandler messageHandler(&threadPool);
    }
    if(size == 1){
        threadPool.waitUntilEmpty();
        threadPool.stop();
    }
}

int main(int argc, char ** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided != MPI_THREAD_MULTIPLE){    
        printf("MPI_Init_thread() failed\n");
        return -1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(argc != 2){
        fprintf(stderr, "wrong args\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    int tasksNum = atoi(argv[1]);
    double start = MPI_Wtime();
    if(rank == 0){
        finishedProcesses = new bool[size]{false};
    }
    
    calculate(tasksNum);
    if(rank == 0){
        delete[] finishedProcesses;
    }
    
    double end = MPI_Wtime();
    if(rank == 0){
        printf("[MAIN thread from #%d process]: Time is %lf; \n", rank, end - start);
    }
    
    MPI_Finalize();

    return EXIT_SUCCESS;
}