 #include <iostream>
 #include <vector>
 #include <cstdlib>
 #include <cassert>

 #include "prime-field.hpp"
 #include "time_measure.hpp"
 #include "fft_2.hpp"

 template < unsigned N >
 struct const_log2 {
	static_assert ( !(N & N-1) );
	static const unsigned value = 1 + const_log2< (N >> 1) >::value;
 };
 template < >
 struct const_log2<1> {
	static const unsigned value = 0; 
 };

 int main() {
	
	const unsigned p1 = 27 * (1 << 26) + 1;
	typedef PrimeField<p1> T;
	
	std::cout << " sizeof(T) = " << sizeof(T) << std::endl;
	
	const unsigned l = const_log2<512 / sizeof(T)>::value;
	const unsigned c = const_log2<32*8*1024 / const_pow2(l)>::value;
	const unsigned n = 24;
	
	std::cout << " l = " << l << " c = " << c << std::endl;
	
	std::vector<T> y(pow2(n)), z(pow2(n));
	for ( unsigned i = 0; i < pow2(n); ++i )
		y[i] = z[i] = T(rand());
	
	T omega = dicho(T::primitive_element(), (T::order() - 1) / pow2(n) );
	
	if ( y.size() <= 16 ) { for ( auto& t : y ) std::clog << ' ' << t; std::clog << " <<== y\n"; }
	if ( z.size() <= 16 ) { for ( auto& t : z ) std::clog << ' ' << t; std::clog << " <<== z\n"; }
	
	time_measure tm;

	fft_stage2_math(&z[0], omega, n, n);
	std::cout << " Hello!" << std::endl;
	fft_stage2_cache<l, c>(&y[0], omega, n);

	tm.stop();

	std::cout << tm.get() << std::endl;
	
	if ( y.size() <= 16 ) { for ( auto& t : y ) std::clog << ' ' << t; std::clog << " <<== y\n"; }
	if ( z.size() <= 16 ) { for ( auto& t : z ) std::clog << ' ' << t; std::clog << " <<== z\n"; }
	
	for ( unsigned i = 0; i < pow2(n); ++i )
		assert( y[i] == z[i] );
	
	return 0;
 }
