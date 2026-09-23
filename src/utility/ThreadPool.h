#pragma once

#include <thread>
#include <mutex>
#include <queue>
#include <vector>
#include <functional>
#include <condition_variable>

class ThreadPool {
public:
    ThreadPool(int max_threads);
    ~ThreadPool();

    void start();
    void QueueJob(std::function<void()> job, int id = 0);
    void Stop();
    void Pause();
    void Resume();
    bool Busy();
    bool CheckBusyThreads();
    
    // Wait until all queued jobs and active thread tasks are completed
    void waitFinished();

    int getMaxThreads() const { return m_MaxThreads; }

private:
    void ThreadLoop();

    bool Paused = false;
    int m_BusyThreads = 0;
    bool Terminate = false;

    std::mutex m_Queue;
    std::mutex m_Busy;
    std::condition_variable m_MutexCondition;
    std::condition_variable m_WaitCondition;

    std::vector<std::thread> m_threads;
    std::queue<std::pair<std::function<void()>, int>> jobs;

    int m_MaxThreads;
};