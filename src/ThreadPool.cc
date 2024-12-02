#include "ThreadPool.hh"

using namespace std;

ThreadPool::ThreadPool(size_t nThreads = thread::hardware_concurrency()){
    threadNum = nThreads;
    for(size_t i = 0; i < nThreads; ++i){
        workerThreads.emplace_back([this]{
            
            while(true){
                function<void()> task;{
                    // The reason for putting the below code
                    // here is to unlock the queue before
                    // executing the task so that other
                    // threads can perform enqueue tasks

                    //Locks queue so that data can be shared safely
                    unique_lock<mutex> lock(queueMutex);

                    //Waits until there's a task to execute or the pool 
                    //is stopped
                    cv.wait(lock, [this]{
                        return !tasks.empty() || stop;
                    });

                    //Exits if the pool is stopped
                    if(stop && tasks.empty()){
                        return;
                    }
                    //Get next task from queue
                    task = move(tasks.front());
                    tasks.pop();
                }

                task();
            }
        });
    }
}
ThreadPool::~ThreadPool(){
    {
        //Queue bust be loked before changing stop flag
        unique_lock<mutex> lock(queueMutex);
        stop = true;
    }
    cv.notify_all();

    for(auto& thread : workerThreads){
        thread.join();
    }
}
void ThreadPool::Enqueue(function<void()> task){
    {
        unique_lock<mutex> lock(queueMutex);
        tasks.emplace(move(task));
    }
    cv.notify_one();
}