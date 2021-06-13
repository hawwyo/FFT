#include <bits/stdc++.h>
#include "tbb/tbb.h"

using namespace std;

int main()
{
    
    tbb::parallel_invoke(
        [](){ this_thread::sleep_for(10ms); cout << "TBB 1\n";},
        [](){cout << "TBB 2\n";}
    );
    
    return 0;
}
