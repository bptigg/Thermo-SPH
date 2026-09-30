#include "ThreadPool.h"

ThreadPool::ThreadPool(int max_threads)
    : m_MaxThreads(max_threads), m_BusyThreads(0), Started(false) {}

ThreadPool::~ThreadPool() {
    if (!Terminate) {
        Stop();
    }
}

void ThreadPool::start() {
    Terminate = false;
    if(Started) { return; }
    m_threads.reserve(m_MaxThreads);
    for (int i = 0; i < m_MaxThreads; i++) {
        m_threads.emplace_back(std::thread(&ThreadPool::ThreadLoop, this));
    }
    Started = true;
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
        if (!Started) return;
        Terminate = true;
    }
    m_MutexCondition.notify_all();
    for (std::thread& active_thread : m_threads) {
        if (active_thread.joinable()) {
            active_thread.join();
        }
    }
    m_threads.clear();
    {
        std::unique_lock<std::mutex> lock(m_Queue);
        Started = false;
    }
}

void ThreadPool::Pause() {
    std::unique_lock<std::mutex> lock(m_Queue);
    Paused = true;
}

void ThreadPool::Resume() {
    {
        std::unique_lock<std::mutex> lock(m_Queue);
        Paused = false;
    }
    m_MutexCondition.notify_all();
}

bool ThreadPool::Busy() {
    std::unique_lock<std::mutex> lock(m_Queue);
    return !jobs.empty() || (m_BusyThreads > 0);
}

bool ThreadPool::CheckBusyThreads() {
    std::unique_lock<std::mutex> lock(m_Busy);
    return m_BusyThreads > 0;
}

void ThreadPool::parallel_for(size_t start, size_t end, const std::function<void(size_t, size_t)> &func, size_t chunkSize)
{
    if (start >= end) return;
    size_t total = end - start;
    size_t numChunks = (total + chunkSize - 1) / chunkSize;
    
    std::atomic<size_t> remainingChunks(numChunks);
    std::mutex completionMutex;
    std::condition_variable completionCv;

    for (size_t c = 0; c < numChunks; ++c) {
        size_t chunkStart = start + c * chunkSize;
        size_t chunkEnd = std::min(chunkStart + chunkSize, end);

        QueueJob([func, chunkStart, chunkEnd, &remainingChunks, &completionMutex, &completionCv]() {
            func(chunkStart, chunkEnd);

            if (--remainingChunks == 0) {
                std::lock_guard<std::mutex> lock(completionMutex);
                completionCv.notify_all();
            }
        });
    }
    std::unique_lock<std::mutex> lock(completionMutex);
    completionCv.wait(lock, [&]() { return remainingChunks.load() == 0; });
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