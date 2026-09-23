#include "ThreadPool.h"

ThreadPool::ThreadPool(int max_threads)
    : m_MaxThreads(max_threads), m_BusyThreads(0) {}

ThreadPool::~ThreadPool() {
    if (!Terminate) {
        Stop();
    }
}

void ThreadPool::start() {
    Terminate = false;
    for (int i = 0; i < m_MaxThreads; i++) {
        m_threads.emplace_back(std::thread(&ThreadPool::ThreadLoop, this));
    }
}

void ThreadPool::QueueJob(std::function<void()> job, int id) {
    {
        std::unique_lock<std::mutex> lock(m_Queue);
        jobs.push({ job, id });
    }
    m_MutexCondition.notify_one();
}

void ThreadPool::Stop() {
    {
        std::unique_lock<std::mutex> lock(m_Queue);
        Terminate = true;
    }
    m_MutexCondition.notify_all();
    for (std::thread& active_thread : m_threads) {
        if (active_thread.joinable()) {
            active_thread.join();
        }
    }
    m_threads.clear();
}

void ThreadPool::Pause() {
    Paused = true;
}

void ThreadPool::Resume() {
    Paused = false;
}

bool ThreadPool::Busy() {
    std::unique_lock<std::mutex> lock(m_Queue);
    return !jobs.empty() || (m_BusyThreads > 0);
}

bool ThreadPool::CheckBusyThreads() {
    std::unique_lock<std::mutex> lock(m_Busy);
    return m_BusyThreads > 0;
}

void ThreadPool::waitFinished() {
    std::unique_lock<std::mutex> lock(m_Busy);
    m_WaitCondition.wait(lock, [this]() {
        std::unique_lock<std::mutex> qLock(m_Queue);
        return jobs.empty() && (m_BusyThreads == 0);
    });
}

void ThreadPool::ThreadLoop() {
    while (true) {
        if (Paused) { 
            std::this_thread::yield();
            continue; 
        }

        std::pair<std::function<void()>, int> job;
        {
            std::unique_lock<std::mutex> lock(m_Queue);
            m_MutexCondition.wait(lock, [this] {
                return !jobs.empty() || Terminate;
            });

            if (Terminate && jobs.empty()) {
                return;
            }

            job = jobs.front();
            jobs.pop();
        }

        {
            std::unique_lock<std::mutex> lock(m_Busy);
            m_BusyThreads++;
        }

        job.first();

        {
            std::unique_lock<std::mutex> lock(m_Busy);
            m_BusyThreads--;
            if (m_BusyThreads == 0) {
                m_WaitCondition.notify_all();
            }
        }
    }
}