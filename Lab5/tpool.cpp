#include <iostream>
#include <thread>
#include <queue>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>
#include <functional>
#include "tpool.h"

ThreadPool::ThreadPool(unsigned int threads_count) {
    for (unsigned int i = 0; i < threads_count; ++i){
        threads_.push_back(std::thread(std::bind(&ThreadPool::worker, this)));
    }
    isWorking_ = true;
}

ThreadPool::~ThreadPool(){
    for (auto& t : threads_){
        t.join();
    }
    printf("~ThreadPool(). Total tasks done: %d\n", tasksDone);
}

void ThreadPool::push(const Task& task){
    std::lock_guard<std::mutex> lock(m_);
    q_.push(task);
    cv_.notify_one();
}

void ThreadPool::stop(){
    isWorking_ = false;
    cv_.notify_all();
}

bool ThreadPool::empty(){
    return q_.empty();
}

int ThreadPool::popTask(){
    std::lock_guard<std::mutex> lock(m_);
    if(!q_.empty()){
        int ret = q_.front().getWeight();
        q_.pop();
        return ret;
    }
    return EMPTY_QUEUE;
}

void ThreadPool::waitUntilEmpty(){
    std::unique_lock<std::mutex> lock(stealerM_);
    stealerCV_.wait(lock);
}

void ThreadPool::worker(){
    while (isWorking_){
        std::unique_lock<std::mutex> lock(m_);
        while(q_.empty() && isWorking_){
            cv_.wait(lock);
        }
        if(isWorking_){
            Task task = q_.front();
            q_.pop();
            lock.unlock();
            task.doTask();
            tasksDone++;
            if(empty()){
                stealerCV_.notify_one();
            }
        }
    }
    printf("[Worker]: Finished.\n");    
}