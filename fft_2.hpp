 #include <vector>
 
 unsigned pow2(unsigned v) {
	 return 1U << v;
 }
 constexpr unsigned const_pow2(unsigned n) {
	return 1U << n; 
 }
 
 template < typename T >
 void fft_stage2_math(T* y, const T& omega, const unsigned n, const unsigned t) {
	std::vector<T> w(n);
	w[n-1] = omega;
	for ( unsigned s = n; --s; )
		w[s-1] = w[s] * w[s];
	for ( unsigned s = 0; s < t; ++s ) {
		T wi(1);
		for ( unsigned i = 0; i < pow2(n); ++i, wi *= w[s] )
			if ( i & pow2(s) );
			else {
				auto j = i + pow2(s);
				auto yi = y[i] + wi * y[j];
				auto yj = y[i] - wi * y[j];
				y[i] = yi;
				y[j] = yj;
			}
	}
 }
 
 template < unsigned l, typename T >
 void fill_vector(T* v, T first, const T& step) {
	for ( unsigned i = 0; i < const_pow2(l); ++i, first *= step )
		v[i] = first;
 }
 
 template < unsigned l, typename T >
 void fill_const_vector(T* v, T first) {
	for ( unsigned i = 0; i < l; ++i )
		v[i] = first;
 }
 
 template < unsigned l, typename T >
 void mul_const(T* u, const T& v) {
	for ( unsigned i = 0; i < const_pow2(l); ++i )
		u[i] *= v;
 }
 
 template < unsigned l, typename T >
 void mul_pairwise(T* u, const T* v, const T* w) {
	for ( unsigned i = 0; i < l; ++i )
		u[i] = v[i] * w[i];
 }
 
 template < unsigned l, typename T >
 void add_pairwise(T* u, const T* v, const T* w) {
	for ( unsigned i = 0; i < const_pow2(l); ++i )
		u[i] = v[i] + w[i];
 }
 
 template < unsigned l, typename T >
 void sub_pairwise(T* u, const T* v, const T* w) {
	for ( unsigned i = 0; i < const_pow2(l); ++i )
		u[i] = v[i] - w[i];
 }
 
 template < unsigned l, unsigned c, typename T >
 void fft_stage2_cache(T* y, const T& omega, const unsigned n) {
	const auto n0 = n - (n - l) / c * c;
	fft_stage2_math(y, omega, n, n0);
	std::vector<T> w(n);
	w[n-1] = omega;
	for ( unsigned s = n; --s; )
		w[s-1] = w[s] * w[s];
	for ( unsigned t = n0; t < n; t += c )
		for ( unsigned jc = 0; jc < pow2(n); jc += pow2(t+c) )
			for ( unsigned il = 0; il < pow2(t); il += const_pow2(l) )
				for ( unsigned s = t; s < t + c; ++s ) {
					T vec_w[const_pow2(l)];
					fill_vector<l>(vec_w, dicho(w[s], il), w[s]);
					for ( unsigned js = 0; js < pow2(t+c); js += pow2(s+1) ) {
						T vec_v[const_pow2(l)];
						std::copy(&vec_w[0], &vec_w[const_pow2(l)], &vec_v[0]);
						for ( unsigned it = 0; it < pow2(s); it += pow2(t), mul_const<l>(vec_v, w[s-t]) ) {
							unsigned i = it + js + il + jc;
							unsigned j = i + pow2(s);
							T tmp[const_pow2(l)];
							mul_pairwise<const_pow2(l)>(tmp, vec_v, y+j);
							sub_pairwise<l>(y+j, y+i, tmp);
							add_pairwise<l>(y+i, y+i, tmp);
						}
					}
				}
 }
