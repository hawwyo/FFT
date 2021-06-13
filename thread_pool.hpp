#ifndef __THREAD_POOL_HPP_INCLUDED__
#define __THREAD_POOL_HPP_INCLUDED__

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <complex>
#include <queue>
#include <atomic>

class thread_pool
{
private:

    struct my_func {
        std::function<void( std::vector< std::complex<double> >* , int , int , int , int , int , int, bool  )> f;
        std::vector< std::complex<double> >* a;
        int start_i, start_j, cnt, len, i_offset, j_ofset;
        bool inv;
    };
    
    std::vector< std::thread > all_threads;
    int _number_of_thread;

    void infinite_loop_function();
    void shutdown();
        
    std::mutex mutex;
    std::mutex threadpool_mutex;
    std::condition_variable condition;
    std::queue< my_func > job_queue;
    bool terminate_pool = 0;
    bool stopped = 0;

    std::atomic<int> cnt_of_unfinished;

public:
    thread_pool(int number_of_thread);
    ~thread_pool();

    void add_job(std::function<void( std::vector< std::complex<double> >* , int , int , int , int , int, int, bool )> job, std::vector< std::complex<double> >* a, int start_i, int start_j, int cnt, int len, int i_offset, int j_ofset, bool inv);
    const std::atomic<int> &cnt_of_jobs_left() {
        return cnt_of_unfinished;
    }
};

thread_pool::thread_pool( int number_of_threads ) : _number_of_thread(number_of_threads)
{
    all_threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        all_threads[i] = std::thread([this](){infinite_loop_function();});
    }
    cnt_of_unfinished = 0;
}

void thread_pool::add_job(std::function<void( std::vector< std::complex<double> >*, int , int , int , int , int , int, bool  )> job, std::vector< std::complex<double> >* a, int start_i, int start_j, int cnt, int len, int i_offset, int j_ofset, bool inv) {
    {
        std::unique_lock<std::mutex> lock(mutex);
        my_func f;
        f.f = job;
        f.a = a;
        f.start_i = start_i;
        f.start_j = start_j;
        f.len = len;
        f.cnt = cnt;
        f.i_offset = i_offset;
        f.j_ofset = j_ofset;
        f.inv = inv;
        // job_queue.push({ job, a, start_i, start_j, cnt, len, inv });
        job_queue.push(f);
    cnt_of_unfinished++;
    condition.notify_one();
    }
}

void thread_pool::infinite_loop_function() {

    while(true)
    {
        my_func job;
        {
            std::unique_lock<std::mutex> lock(mutex);

            condition.wait(lock, [&]{return !job_queue.empty() || terminate_pool;});

            if (job_queue.empty()) break;

            job = job_queue.front();
            job_queue.pop();
        }
        job.f( job.a, job.start_i, job.start_j, job.cnt, job.len, job.i_offset, job.j_ofset, job.inv ); 
        cnt_of_unfinished--;
    }

}

thread_pool::~thread_pool()
{
    shutdown();
}

void thread_pool::shutdown() {
    {
        std::unique_lock<std::mutex> lock(threadpool_mutex);
        terminate_pool = true;
    } // use this flag in condition.wait

    condition.notify_all(); // wake up all threads.

    // Join all threads.
    for(std::thread &every_thread : all_threads)
    {   every_thread.join();}

    all_threads.clear();  
    stopped = true; // use this flag in destructor, if not set, call shutdown()
}


#endif