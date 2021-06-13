 #ifndef PRIME_FIELD_INCLUDED
 #define PRIME_FIELD_INCLUDED

 #include <ostream>
 #include <mutex>
 #include <atomic>
 #include <cassert>

 #include "dicho.hpp"
 
 template < unsigned p >
 class PrimeField {
    unsigned value;
    template < typename T >
    static unsigned from_int(T v) {
        return v < 0 ? p - (-v) % p : v % p;
    }
    static PrimeField<p> get_primitive_element() {
        PrimeField<p> xi( 2 + rand() % (p-2) );
        while ( dicho(xi, (p-1)/2, 1) == 1 )
            xi = 2 + rand() % (p-2);
        return xi;
    }
 public:
    template < typename T >
    PrimeField(T v) : value ( from_int(v) ) {
    }
    PrimeField() : value (0) {
    }
    PrimeField<p>& operator+=(PrimeField<p> v) {
        static_assert( (p >> 31) == 0 );
        value += v.value;
        if ( value >= p )
            value -= p;
        return *this;
    }
    static size_t order() {
        return p;
    }
    PrimeField<p>& operator-=(PrimeField<p> v) {
        static_assert( (p >> 31) == 0 );
        if ( value < v.value )
            value += p;
        value -= v.value;
        return *this;
    }
    PrimeField<p>& operator*=(PrimeField<p> v) {
        value = static_cast<unsigned long long>(value) * v.value % p;
        return *this;
    }
    PrimeField<p>& operator/=(PrimeField<p> a) {
        return *this = dicho(a, p-2, *this);
    }
    PrimeField<p>& operator-() {
        value = p - value;
        return *this;
    }
    PrimeField<p>& operator++() {
        return *this += 1;
    }
    PrimeField<p>& operator--() {
        return *this -= 1;
    }
    PrimeField<p> operator++(int) {
        PrimeField<p> a(*this);
        *this += 1;
        return a;
    }
    PrimeField<p> operator--(int) {
        PrimeField<p> a(*this);
        *this -= 1;
        return a;
    }
    unsigned to_int() const {
        return value;
    }
    // operator unsigned() const {
    //     return value;
    // }
    // operator bool() const {
    //     return value;
    // }
    static PrimeField<p> primitive_element() {
        static PrimeField<p> xi = get_primitive_element();
        return xi;
    }
 };

 template < unsigned p, typename T >
 PrimeField<p> operator+ ( PrimeField<p> a, T b ) {
    return a += b;
 }

 template < unsigned p, typename T >
 PrimeField<p> operator- ( PrimeField<p> a, T b ) {
    return a -= b;
 }

 template < unsigned p, typename T >
 PrimeField<p> operator* ( PrimeField<p> a, T b ) {
    return a *= b;
 }

 template < unsigned p, typename T >
 PrimeField<p> operator/ ( PrimeField<p> a, T b ) {
    return a /= b;
 }

 template < unsigned p >
 bool operator== ( PrimeField<p> a, PrimeField<p> b ) {
    return a.to_int() == b.to_int();
 }

 template < unsigned p, typename T >
 bool operator== ( PrimeField<p> a, T b ) {
    return a == PrimeField<p>(b);
 }

 template < unsigned p, typename T >
 bool operator!= ( PrimeField<p> a, T b ) {
    return !(a == b);
 }

 template < unsigned p >
 std::ostream& operator<< ( std::ostream& out, PrimeField<p> a ) {
    return out << a.to_int();
 }

 template < unsigned p, int t, typename T  >
 void gen_roots_of_unity(T* &arr_w) {

    assert( (p-1) % (1 << t) == 0 );

    const unsigned N = 1 << t;

    const PrimeField<p> xi = PrimeField<p>::primitive_element();
    PrimeField<p> omega = dicho(xi, (p-1)/N);

    PrimeField<p>* w = new PrimeField<p>[N];
    w[0] = 1;
    for ( int i = 1; i < N; i++ )
        w[i] = w[i-1] * omega;

    arr_w[t] = w;
 }

 template < unsigned p, int t, typename T  >
 T& get_or_gen_roots_of_unity(T* &w) {
    static std::once_flag f;
    if ( w[t] )
        return w[t];
    std::call_once(f, gen_roots_of_unity<p, t, T>, w);
    return w[t];
 }

 template < unsigned p >
 PrimeField<p>* get_roots_of_unity(unsigned n) {

    static std::atomic<PrimeField<p>*> arr_w[27] = {0};
    std::atomic<PrimeField<p>*>* w = arr_w;

    switch ( n ) {
        case 1 <<  0: return get_or_gen_roots_of_unity<p,  0>(w);
        case 1 <<  1: return get_or_gen_roots_of_unity<p,  1>(w);
        case 1 <<  2: return get_or_gen_roots_of_unity<p,  2>(w);
        case 1 <<  3: return get_or_gen_roots_of_unity<p,  3>(w);
        case 1 <<  4: return get_or_gen_roots_of_unity<p,  4>(w);
        case 1 <<  5: return get_or_gen_roots_of_unity<p,  5>(w);
        case 1 <<  6: return get_or_gen_roots_of_unity<p,  6>(w);
        case 1 <<  7: return get_or_gen_roots_of_unity<p,  7>(w);
        case 1 <<  8: return get_or_gen_roots_of_unity<p,  8>(w);
        case 1 <<  9: return get_or_gen_roots_of_unity<p,  9>(w);
        case 1 << 10: return get_or_gen_roots_of_unity<p, 10>(w);
        case 1 << 11: return get_or_gen_roots_of_unity<p, 11>(w);
        case 1 << 12: return get_or_gen_roots_of_unity<p, 12>(w);
        case 1 << 13: return get_or_gen_roots_of_unity<p, 13>(w);
        case 1 << 14: return get_or_gen_roots_of_unity<p, 14>(w);
        case 1 << 15: return get_or_gen_roots_of_unity<p, 15>(w);
        case 1 << 16: return get_or_gen_roots_of_unity<p, 16>(w);
        case 1 << 17: return get_or_gen_roots_of_unity<p, 17>(w);
        case 1 << 18: return get_or_gen_roots_of_unity<p, 18>(w);
        case 1 << 19: return get_or_gen_roots_of_unity<p, 19>(w);
        case 1 << 20: return get_or_gen_roots_of_unity<p, 20>(w);
        case 1 << 21: return get_or_gen_roots_of_unity<p, 21>(w);
        case 1 << 22: return get_or_gen_roots_of_unity<p, 22>(w);
        case 1 << 23: return get_or_gen_roots_of_unity<p, 23>(w);
        case 1 << 24: return get_or_gen_roots_of_unity<p, 24>(w);
        case 1 << 25: return get_or_gen_roots_of_unity<p, 25>(w);
        case 1 << 26: return get_or_gen_roots_of_unity<p, 26>(w);
    }
    assert ( false );
    return 0;
 }

 #endif
