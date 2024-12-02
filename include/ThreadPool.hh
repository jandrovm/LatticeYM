#ifndef THREADPOOL_HH
#define THREADPOOL_HH

#include <thread>
#include <queue>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <vector>

class ThreadPool{
    public:
        ThreadPool(size_t nThreads);
        ~ThreadPool();
        void Enqueue(std::function<void()>);
        size_t GetThreadNumber(){ return threadNum; }
    private:
        bool stop = false;
        size_t threadNum;
        std::vector<std::thread> workerThreads;
        
        //tasks to be enqueued
        std::queue<std::function<void()> > tasks;
        
        //Synchronises access to shared data (whatever that means)
        std::mutex queueMutex;

        //Condition variable for locking the queue
        std::condition_variable cv;
};

#endif