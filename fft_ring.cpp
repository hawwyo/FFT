#include <bits/stdc++.h>

using namespace std;

const long long mod = 1484783617; // 177 * (2 ^ 23) + 1
const long long e = 5;

struct field {
    long long a, b;

    field operator + (const field &o) {
        return { (a + o.a) % mod, (b + o.b) % mod };
    }
    field operator - (const field &o) {
        return { ((a - o.a) % mod + mod) % mod, ((b - o.b) % mod + mod) % mod };
    }
    field operator * (const field &o) {
        return { ( (a * o.a) % mod + (b * o.b) % mod * e) % mod, ( (a * o.b) % mod + (b * o.a) % mod ) % mod };
    }
};


const field root = {0, 1};
const field root_1 = {0, 593913447};
const long long root_len = (1ll << 24);

field power(long long n, field a) {
    field res = {1, 0};
    while (n) {

        if (n & 1) {
            res = res * a;
            n--;
        }
        else {
            a = a * a;
            n >>= 1;
        }
    }
    return res;
}


vector<int> bit_rev;

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


void fft(vector< field > &a, bool inv) {
    int n = a.size();

    for (int i = 0; i < n; i++) {
        if (i < bit_rev[i]) {
            swap( a[i], a[ bit_rev[i] ] );
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        int half_len = len / 2;

        field w = {0, 0};
        if (!inv) {
            w = power( 2 * (mod - 1) / len, root );
        }
        else {
            w = power( 2 * (mod - 1) / len, root_1 );
        }

//        field tmp = power(len, w);
//        cout << len << " " << tmp.a << " " << tmp.b << "   w = " << w.a << " " << w.b << endl;

        for (int i = 0; i < n; i += len) {

            field cur = {1, 0};

            for (int j = 0; j < half_len; j++) {

                field left = a[i + j], right = a[i + j + half_len];
                a[i + j] = left + cur * right;
                a[i + j + half_len] = left - cur * right;

                cur = cur * w;

            }
        }

    }


//    field n_t = {n, 0};
//    field n_tmp = n_1 * n_t;

    if (inv) {
        field n_1 = {n, 0};
        n_1 = power( mod - 2, n_1 );
        for (field &item : a) {
            item = item * n_1;
        }
    }
}

void mul( vector< int > &a, vector< int > &b, vector< int > &res ) {
    long long n = 1;
    while (n < a.size() && n < b.size()) {
        n <<= 1;
    }
    n <<= 1;

    bit_rev_calc(n);

    vector< field > na(n, {0, 0});
    vector< field > nb(n, {0, 0});
    vector< field > nc(n, {0, 0});
//
    for (int i = 0; i < a.size(); i++) {
        na[i] = { a[i], 0 };
    }
    for (int i = 0; i < b.size(); i++) {
        nb[i] = { b[i], 0 };
    }
//
    fft(na, 0);

//    for (field t : na) {
//        cout << t.a << " " << t.b << endl;
//    }
//    return;

    fft(nb, 0);
//
    for (int i = 0; i < n; i++) {
        nc[i] = na[i] * nb[i];
    }
//
//    // TODO: inverse fft
    fft(nc, 1);
//
    res.resize(n);
    for (int i = 0; i < n; i++) {
        res[i] = nc[i].a;
    }
}

int main()
{


    int n = 5;
    vector< int > a(n);
    vector< int > b(n);
    for (int i = 0; i < n; i++) {
        a[i] = 1;
        b[i] = 1;
    }

    vector< int > res;
    mul(a, b, res);

    for (int i : res) {
        cout << i << " ";
    }


    return 0;
}
