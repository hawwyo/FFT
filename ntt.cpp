#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

#include "time_measure.hpp"

using namespace std;

// const int mod = 998244353; // 119 * (1 << 23) + 1
// const int root = 3;
// const int root_1 = 332748118;
// const int root_len = (1 << 23);
// const int inv_2 = 499122177;

const int mod = 7340033;
const int root = 5;
const int root_1 = 4404020;
const int root_len = 1<<20;

int power(int n, int a) {
    int res = 1;
    while (n) {
        if (n & 1) {
            res = (res * 1ll * a) % mod;
            n--;
        }
        else {
            a = (a * 1ll * a) % mod;
            n >>= 1;
        }
    }
    return res;
}

void fft( vector< int > &a, bool inv ) {

    int n = a.size();
    if (n == 1) {
        return;
    }

    vector< int > a1(n / 2), a2(n / 2);

    for (int i = 0; i < n; i++) {
        if ( (i & 1) == 0 ) {
            a1[i / 2] = a[i];
        }
        else {
            a2[i / 2] = a[i];
        }
    }

    fft(a1, inv);
    fft(a2, inv);

    int cur = 1;
    int w = (inv ? root_1 : root);

    for (int i = n; i < root_len; i *= 2) {
        w = (w * 1ll * w) % mod;
    }
    // int w = power( (mod - 1) / n, root );

    // if (inv) {
    //     w = power(mod - 2, w);
    // }

    for (int i = 0; i < n; i++) {
        if (i < n / 2) {
            a[i] = ( a1[i] * 1ll + cur * 1ll * a2[i] ) % mod;
        }
        else {
            a[i] = ( a1[i - n / 2] * 1ll + cur * 1ll * a2[i - n / 2] ) % mod;
        }

        if (inv) {
            a[i] = (a[i] * 1ll * 3670017) % mod;
            // a[i] = (a[i] * 1ll * power(mod - 2, 2)) % mod;
        }

        cur = (cur * 1ll * w) % mod;
    }

}

vector< int > mul(vector< int > &a, vector< int > &b) {
    vector< int > ca(a.begin(), a.end());
    vector< int > cb(b.begin(), b.end());

    int n = 1;
    while (n < max( a.size(), b.size() )) {
        n *= 2;
    }
    n *= 2;

    ca.resize(n);
    cb.resize(n);

    fft(ca, 0);
    fft(cb, 0);

    vector< int > c(n);
    for (int i = 0; i < n; i++) {
        c[i] = (ca[i] * 1ll * cb[i]) % mod;
    }

    fft(c, 1);

    vector< int > res(n);

    for (int i = 0; i < n; i++) {
        res[i] = c[i];
    }

    return res;

}


int main() {
    
    int n = 1000000;
    vector< int > a(n);
    vector< int > b(n);

    for (int i = 0; i < n; i++) {
        a[i] = rand();
        b[i] = rand();
        // a[i] = 1;
        // b[i] = 1;
    }

    time_measure tm;
    int sz = mul(a, b).size();
    tm.stop();

    cout << tm.get() << endl;

    // for (int i : mul(a, b)) {
    //     cout << i << " ";
    // }
    
}