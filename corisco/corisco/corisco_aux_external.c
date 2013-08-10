// Copyright 2011 Nicolau Leal Werneck, Anna Helena Reali Costa and
// Universidade de São Paulo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include<stdio.h>

// file corisco_aux_external.c
//
// Contains function SSE_rsqrt to calculate reciprocal of square root
// using SSE helper instruction plus a single Newton-Raphson iteration.
//
// Totally ripped off from example posted by dude at stackoverflow
// To use with Cython, compile with gcc -fPIC -c rsqrt.c -o rsqrt.o
// stolen and ressocialized by nwerneck@gmail.com - 16/09/2010

#include<emmintrin.h>
#include<xmmintrin.h>

#include<math.h>
#include<tgmath.h>


/** These macros use inline asm to apply the rsqrt isntruction to calculate the reverse square root ofa  value in either a variable that will be the output, or an expression result, etc (the "2" version).*/
#define NIC_SSE_RSQRT(A)   asm ("rsqrtps %[a], %[a]" : [a] "=x"((A)) : "x"((A)) );
#define NIC_SSE_RSQRT2(B,A)   asm ("rsqrtps %0, %0" : "=x"((B)) : "x"((A)) );

#define VEC_ASSIGN(var, val) (var).f[0]=(val);(var).f[1]=(val);(var).f[2]=(val);(var).f[3]=(val);

typedef float v4sf __attribute__ ((  vector_size (16) )); // vector of four single floats
  
union f4vector 
{
  v4sf v;
  float f[4];
};




/* Uses rsqrt instruction intrinsic to calculate reciprocal of square root */
inline void SSE_rsqrt( float *pOut, float *pIn )
{
    _mm_store_ss( pOut, _mm_rsqrt_ss( _mm_load_ss( pIn ) ) );
}

/* Gets rsqrt estimate then perform as ingle step of the
   Newton-Raphson iteration. */
inline void SSE_rsqrt_NR( float *pOut, float *pIn )
{
    _mm_store_ss( pOut, _mm_rsqrt_ss( _mm_load_ss( pIn ) ) );
    // compiles to movss, sqrtss, movss
    *pOut *= ((3.0f - *pOut * *pOut * *pIn) * 0.5f);
}

/* Normalizes a vector using approximate square root calculation. */
inline void external_normalize(double* x, double* y)
{
  float fni = x[0]*x[0]+y[0]*y[0];
  float fn;
  SSE_rsqrt_NR(&fn, &fni);
  x[0] = fn * x[0];
  y[0] = fn * y[0];
}





/** SSE Vector operations */
inline void myvec_sumacc(float*a, float*b) {
  (*(union f4vector*)a).v +=   (*(union f4vector*)b).v;
}
inline void myvec_mulacc(float*a, float*b) {
  (*(union f4vector*)a).v *=   (*(union f4vector*)b).v;
}
inline void myvec_pos_lim(float*a) {
  static float z[4] = {0.0f,0.0f,0.0f,0.0f};
  (*(__m128*)a) = _mm_max_ps((*(__m128*)a) , (*(__m128*)z) );
}
inline void myvec_abs(float*a) {
  (*(__m128*)a) = _mm_max_ps((*(__m128*)a) , -(*(__m128*)a) );
}
inline void myvec_copy(float*a, float*b) {
  (*(__m128*)a) = (*(__m128*)b);
}
inline void myvec_rsqrt(float*a) {
  (*(__m128*)a) = _mm_rsqrt_ps((*(__m128*)a));
}
inline void myvec_rcp(float*a) {
  (*(__m128*)a) = _mm_rcp_ps((*(__m128*)a));
}

/** Some floating-pount ioperations */
inline void myabs(float* x){
  *x = (*x>0)?*x:-*x;
}
inline void mypos_lim(float* x){
  *x = (*x>0)?*x:0;
}
/** */





/** Reads clock time (it's a x86_64 thing) */
static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
