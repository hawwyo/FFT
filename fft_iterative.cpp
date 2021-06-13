#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

#include "time_measure.hpp"
#define cd complex<double>

using namespace std;

vector<int> bit_rev;

const double pi = acos(-1);

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

void fft(vector< cd > &a, bool inv) {
    int n = a.size();

    for (int i = 0; i < n; i++) {
        if (i < bit_rev[i]) {
            swap( a[i], a[ bit_rev[i] ] );
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        int half_len = len / 2;

        double ang = 2 * pi / len;
        if (inv) ang = -ang;
        cd w( cos(ang), sin(ang) );

        for (int i = 0; i < n; i += len) {
            cd cur = 1;
            for (int j = 0; j < half_len; j++) {

                cd left = a[i + j], right = a[i + j + half_len];
                a[i + j           ] = left + cur * right;
                a[i + j + half_len] = left - cur * right;

                cur *= w;
            }
        }

    }

    if (inv) {
        for (cd &item : a) {
            item /= n;
        }
    }
}

vector< int > mul(vector< int > &a, vector< int > &b) {
    vector< complex<double> > ca(a.begin(), a.end());
    vector< complex<double> > cb(b.begin(), b.end());

    int n = 1;
    while (n < max( a.size(), b.size() )) {
        n *= 2;
    }
    n *= 2;

    bit_rev_calc(n);

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

    vector< int > res(n);

    for (int i = 0; i < n; i++) {
        res[i] = int( round(c[i].real()) );
    }

    tm1.stop();
    cout << "! " << tm1.get() << endl; 

    return res;

}



int main() {

    int n = 3;
    vector< int > a(n);
    vector< int > b(n);

    for (int i = 0; i < n; i++) {
        // a[i] = rand() % 1000;
        // b[i] = rand() % 1000;
        a[i] = 1;
        b[i] = 1;
    }



    vector< int > c;
    time_measure tm;
    c = mul(a, b);
    tm.stop();

    // ofstream in_file("input.txt");
    // ofstream out_file("output.txt");

    // for (int i = 0; i < n; i++) {
    //     in_file << a[i] << " ";
    // }
    // in_file << endl;
    // for (int i = 0; i < n; i++) {
    //     in_file << b[i] << " ";
    // }
    // for (int i = 0; i < c.size(); i++) {
    //     out_file << c[i] << " ";
    // }

    // cout << tm.get();

    cout << endl;
    for (int i : c) {
        cout << i << " ";
    }

    return 0;
}
