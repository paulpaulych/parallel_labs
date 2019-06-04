#include <iostream>
#include <thread>
#include <queue>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <unistd.h>
#include <vector>
#include <functional>
#include <sys/times.h>

#define EMPTY_QUEUE (-1)

class Task{

private:
    int weight_ = 0;
    int rank_ = 0;

public: 

    Task(int weight__, int rank__){
        weight_ = weight__;
        rank_ = rank__;
    }
    
    int getWeight(){
        return weight_;
    }   

    void setRank(int rank__){
        rank_ = rank__;
    }

    double Random(){
        double x;
        x = (double)rand() / (RAND_MAX / 2) - 1;
        return x;
    }

    void doTask(){
        printf("[Task thread from #%d process]: task started \n", rank_);
        struct tms start, end;
        long clocks_per_sec = sysconf(_SC_CLK_TCK);
        long clocks;

        times(&start);
        long N = 3000000 * weight_;
        long i, m = 0;
        double pi, x, y;
        for (i = 0; i < N; i++){
            x = Random();
            y = Random();
            if (x * x + y * y <= 1){
                m++;
            }
        }

        pi = 4 * (double)m / N;
        times(&end);
        clocks = end.tms_utime - start.tms_utime; 
        printf("[Task thread from #%d process]: Task finished: %fl. Time taken: %lf sec.\n", rank_, pi,
            (double)clocks / clocks_per_sec );
    }
};

class ThreadPool {
private:
    int tasksDone = 0;

	std::condition_variable stealerCV_;
    std::mutex stealerM_;
    
	std::mutex m_;
	std::condition_variable cv_;
	std::queue<Task> q_;
	std::vector<std::thread> threads_;
	std::atomic<bool> isWorking_;
public:
	ThreadPool(unsigned int threads_count);
	
    ~ThreadPool();

 //   void setStealerCV(std::condition_variable * cv);

    void push(const Task& task);
	
    void waitUntilEmpty();

    bool empty();

    void stop();
	
	int popTask();

    void worker();
};