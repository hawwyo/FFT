 #ifndef DICHO_HPP_INCLUDED
 #define DICHO_HPP_INCLUDED
 
 struct Product {
    template < class T > T operator()(const T& a, const T& b) {
        return a * b;
    }
 };

 template < typename T, typename T1 = int, typename Op = Product >
 T dicho(T a, int n, T1 e = 1, Op op = Op() ) {
    T v = e;
    for ( ; n; n >>= 1, a = op(a, a) ) {
        if ( n & 1 )
            v = op(v, a);
    }
    return v;
 }

 #endif
