#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

#include "time_measure.hpp"

using namespace std;

const double pi = acos(-1);

void fft( vector< complex<double> > &a, bool inv ) {

    int n = a.size();
    if (n == 1) {
        return;
    }

    vector< complex< double > > a1(n / 2), a2(n / 2);

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

    double ang = 2 * pi / n * (inv ? -1 : 1);

    complex< double > w( cos(ang), sin(ang) ), cur(1);

    for (int i = 0; i < n; i++) {
        if (i < n / 2) {
            a[i] = a1[i] + cur * a2[i];
        }
        else {
            a[i] = a1[i - n / 2] + cur * a2[i - n / 2];
        }

        if (inv) {
            a[i] /= 2;
        }

        cur *= w;
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

    ca.resize(n);
    cb.resize(n);

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

    return res;

}


int main() {
    
    int n = 100000;
    vector< int > a(n);
    vector< int > b(n);

    for (int i = 0; i < n; i++) {
        a[i] = rand();
        b[i] = rand();
    }

    time_measure tm;
    int sz = mul(a, b).size();
    tm.stop();

    cout << sz << " " << tm.get();
    
}