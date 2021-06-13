#include <iostream>

#include "thread_pool.hpp"

using namespace std;

int main() {

    
    thread_pool *pool = new thread_pool(2);

    cout << pool->cnt_of_jobs_left() << endl;

    delete pool;

    


    return 0;
}