#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <thread>

#include "thread_pool.hpp"
#include "time_measure.hpp"

#define cd complex<double>

using namespace std;

vector<int> bit_rev;
vector< vector< cd > > W;
vector< vector< cd > > W_1;
const double pi = acos(-1);
const int number_of_treads = 4;
const int log_cs = 25;

thread_pool *pool;

int log2ceil(int n) {
    return sizeof(n) * 8 - __builtin_clz(n) - 1;
}

void W_calc(int n) {
    int sz = log2ceil(n);

    W.resize(sz + 1);
    W_1.resize(sz + 1);
    W[sz].resize(n);
    W_1[sz].resize(n);
    W[sz][0] = 1;
    W_1[sz][0] = 1;

    double ang = 2 * pi / n;
    double ang_1 = -2 * pi / n;

    cd w(cos(ang), sin(ang));
    cd w_1(cos(ang_1), sin(ang_1));

    for (int j = 1; j < n; j++) {
        W[sz][j] = W[sz][j - 1] * w;
        W_1[sz][j] = W_1[sz][j - 1] * w_1;
    }

    for (int i = sz - 1; i > 0; i--) {
        W[i].resize( (1 << i) );
        W_1[i].resize( (1 << i) );

        for (int j = 0; j < W[i].size(); j++) {
            W[i][j] = W[i + 1][j << 1];
            W_1[i][j] = W_1[i + 1][j << 1];
        }    
    }
}

void bit_rev_calc(int n) {
    bit_rev.resize(n);
    int bit;
    for (int i = 1, j = 0; i < n; i++) {
        bit = (n >> 1);
        while ( j & bit ) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;

        bit_rev[i] = j;
    }
}

void fft_full_cache_per_thread(vector< cd >* a, int start_i, int start_j, int cnt, int len, int i_offset, int j_offset, bool inv) {
    
    // cout << this_thread::get_id() << endl;

    int half_len = len / 2;
    int log_len = log2ceil(len);
    
    for (int i = start_i; cnt; i += len) {
        for (int j = start_j; j < half_len && cnt; j += j_offset, cnt--) {
            start_j = 0;

            cd u = (*a)[i + j + i_offset];
            cd v = (*a)[i + j + half_len + i_offset] * ( inv ? W[log_len][j + i_offset] : W_1[log_len][j + i_offset] );

            (*a)[i + j + i_offset] = u + v;
            (*a)[i + j + half_len + i_offset] = u - v;
        }
    }
}

void fft_full_cache(vector< cd > *a, int left, int right, int len_scale, int i_offset, bool inv) {    
    int n = right - left + 1;    
    int start_len = (1 << len_scale);
    int cnt_of_operations = (n >> len_scale);
    int operations_per_thread = cnt_of_operations / number_of_treads;
    if (!operations_per_thread) {
        operations_per_thread = 1;
    }
    int j_offset = (1 << (len_scale - 1));
    
    // cout << operations_per_thread << endl;
    for (int len = start_len; len <= n; len <<= 1) {
        int half_len = (len >> 1);
        if ( half_len > operations_per_thread * j_offset) {
            for(int i = left; i < right; i += len){
                for(int j = 0; j < half_len; j += operations_per_thread * j_offset ){
                    // cout << i << " ! " << j << endl;
                    pool->add_job(fft_full_cache_per_thread, a, i, j, operations_per_thread, len, i_offset, j_offset, inv);
                }
            }
        }
        else {
            for (int i = left; i < right; i += (operations_per_thread << len_scale) ) {
                // cout << i << " !! " << 0 << endl;
                pool->add_job(fft_full_cache_per_thread, a, i, 0, operations_per_thread, len, i_offset, j_offset, inv);
            }
        }
        while ( pool->cnt_of_jobs_left() != 0 ) {}
    }
}

void fft_semi_cache(vector< cd > *a, int left, int right, int len_scale, int inv) {

    for (int i0 = 0; i0 < ( 1 << (len_scale - 1) ); i0++) {
        // cout << left + i0 << " " << right + i0 << " " << len_scale << endl;
        fft_full_cache(a, left, right, len_scale, i0, inv);

    }

}

void fft(vector< cd > &a, bool inv) {
    int n = a.size();

    for (int i = 0; i < n; i++) {
        if (i < bit_rev[i]) {
            swap( a[i], a[ bit_rev[i] ] );
        }
    }

    // fft_full_cache(&a, 0, n - 1, 1, inv);

    int log_n = log2ceil(n);

    int k0 = min(log_n, log_cs);
    while ( (log_n - k0) % log_cs ) k0++;
    int step_cs = (1 << k0);
    k0++;
    for (int i = 0; i < n; i += step_cs) {
        fft_full_cache(&a, i, i + step_cs - 1, 1, 0, inv);
    }

    // cout << k0 << " here" << endl;

    for (; k0 < log_n; k0 += log_cs) {

        step_cs = (1 << (k0 + log_cs - 1));

        for (int i = 0; i < n; i += step_cs) {
            // cout << i << " " << i + step_cs - 1 << endl;
            fft_semi_cache(&a, i, i + step_cs - 1, k0, inv);
        }

    }
    

    if (inv) {
        for (cd &item : a) {
            item /= n;
        }
    }
}

vector< int > mul(vector< int > &a, vector< int > &b) {

    pool = new thread_pool(number_of_treads);

    vector< complex<double> > ca(a.begin(), a.end());
    vector< complex<double> > cb(b.begin(), b.end());

    int n = 1;
    while (n < max( a.size(), b.size() )) {
        n *= 2;
    }
    n *= 2;

    bit_rev_calc(n);
    W_calc(n);

    ca.resize(n);
    cb.resize(n);

    time_measure tm1;

    fft(ca, 0);
    

    fft(cb, 0);

    vector< complex<double> > c(n);
    for (int i = 0; i < n; i++) {
        c[i] = ca[i] * cb[i];
    }

    fft(c, 1);
    // return {};

    vector< int > res(n);

    for (int i = 0; i < n; i++) {
        res[i] = int( round(c[i].real()) );
    }

    tm1.stop();
    cout << "! " << tm1.get() << endl; 

    delete pool;

    return res;

}



int main() {

    int n = 3;
    vector< int > a(n);
    vector< int > b(n);

    for (int i = 0; i < n; i++) {
        // a[i] = rand();
        // b[i] = rand();
        a[i] = 1;
        b[i] = 1;
    }

    vector< int > c;
    // time_measure tm;
    c = mul(a, b);
    // tm.stop();

    // cout << tm.get();

    // cout << endl;
    for (int i : c) {
        cout << i << " ";
    }
    cout << endl;

    return 0;
}
