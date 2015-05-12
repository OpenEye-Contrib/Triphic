//
// file popcount_ssse2.cc
// David Cosgrove
// AstraZeneca
// 1st December 2014
//
// This file contains function popcount_ssse3, which counts the number of set
// bits in an array of unsigned ints.  It is taken from Andrew Dalke's popcnt
// benchmarking program:
// For more details see
//   http://dalkescientific.com/writings/diary/archive/2008/07/03/hakmem_and_other_popcounts.html
//   http://dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html
//
// This code is from Supplemental Information to "Anatomy of High-Performance 2D
// Similarity Calculations", Haque IS, Pande VS, and Walters WP,
// Journal of Chemical Information and Modeling 2011
// (dx.doi.org/10.1021/ci200235e|J.Chem.Inf.Model., 2011, 51, 23452351)
// It is under a 3-clause BSD license:
// The I below is Andrew Dalke!

/****** The following code comes with the 3-clause BSD license  ******/
/* It is CPU specific because it uses SSE2 instructions. */
/* It should work under Visual Studio; what's the right SSE2 macro check
       for that compiler?  */

/* Imran notes:
    Note that you do need to compile the SSE popcounts with
    -fno-strict-aliasing, or they'll return incorrect results (the
    epilogues play a little fast and loose with pointer
    casting). Compiling for 64-bit will probably help the performance
    too, since the inner loops are unrolled and were optimized for the
    extra 8 registers. */
/* I found that -fno-strict-aliasing does not affect this benchmark */


/* Source: Supplemental Information to "Anatomy of High-Performance 2D
   Similarity Calculations", Haque IS, Pande VS, and Walters WP,
   Journal of Chemical Information and Modeling 2011 */

/*  Written by Imran S. Haque (ihaque@cs.stanford.edu)

 Copyright (c) 2011 Stanford University.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

     * Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
     * Neither the name of Stanford University nor the names of its contributors
       may be used to endorse or promote products derived from this software without
       specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

// This code is used as a backup when popcount_ssse3 can't be used, due to
// lack of compiler support etc.  At present, the AP cluster is on RH5 whose
// system compiler, gcc 4.1 does seem to support ssse3.

#if defined(__GNUC__) && defined(__SSE2__)
#include <emmintrin.h>

// ****************************************************************************
static inline __m128i _mm_srli1_epi8(__m128i v) {
    // Avg returns (a_i + b_i + 1) >> 1; so make b_i = -1 to compute (a >> 1)
    const __m128i neg1 = _mm_set1_epi8(0xFF);
    return _mm_avg_epu8(v,neg1);
}

// ****************************************************************************
static inline __m128i _mm_srli2_epi8(__m128i v) {
    // Faster than mask-and-shift approach used in srli4
    return _mm_srli1_epi8(_mm_srli1_epi8(v));
}

// ****************************************************************************
static inline __m128i _mm_srli4_epi8(__m128i v) {
    // Using ANDNOT lets us share the magic number with what's used for the popcount itself
    const __m128i mask = _mm_set1_epi32(0x0F0F0F0F);
    return _mm_srli_epi16(_mm_andnot_si128(mask,v),4);
}

// ****************************************************************************
static inline __m128i popcount_sse2_8_helper_1(unsigned* buf, int N) {
    __m128i* vbuf = (__m128i*)buf;
    __m128i total,count0,count1,count2,count3;
    int i;
    __m128i B0 = _mm_set1_epi32(0x55555555);
    __m128i B1 = _mm_set1_epi32(0x33333333);
    __m128i B2 = _mm_set1_epi32(0x0F0F0F0F);
    total = _mm_setzero_si128();
    for (i = 0; i < N; i+=4) {
      __m128i v0 = _mm_load_si128(vbuf+i+0);
      __m128i v1 = _mm_load_si128(vbuf+i+1);
      __m128i v2 = _mm_load_si128(vbuf+i+2);
      __m128i v3 = _mm_load_si128(vbuf+i+3);

      count0 = _mm_sub_epi8(v0, _mm_and_si128(B0, _mm_srli1_epi8(v0)));
      count1 = _mm_sub_epi8(v1, _mm_and_si128(B0, _mm_srli1_epi8(v1)));
      count2 = _mm_sub_epi8(v2, _mm_and_si128(B0, _mm_srli1_epi8(v2)));
      count3 = _mm_sub_epi8(v3, _mm_and_si128(B0, _mm_srli1_epi8(v3)));

      count0 = _mm_add_epi8(_mm_and_si128(B1,count0),
          _mm_and_si128(B1,_mm_srli2_epi8(count0)));
      count1 = _mm_add_epi8(_mm_and_si128(B1,count1),
          _mm_and_si128(B1,_mm_srli2_epi8(count1)));
      count2 = _mm_add_epi8(_mm_and_si128(B1,count2),
          _mm_and_si128(B1,_mm_srli2_epi8(count2)));
      count3 = _mm_add_epi8(_mm_and_si128(B1,count3),
          _mm_and_si128(B1,_mm_srli2_epi8(count3)));

      count0 = _mm_and_si128(B2,_mm_add_epi8(count0,
               _mm_srli4_epi8(count0)));
      count1 = _mm_and_si128(B2,_mm_add_epi8(count1,
               _mm_srli4_epi8(count1)));
      count2 = _mm_and_si128(B2,_mm_add_epi8(count2,
               _mm_srli4_epi8(count2)));
      count3 = _mm_and_si128(B2,_mm_add_epi8(count3,
               _mm_srli4_epi8(count3)));

      total = _mm_add_epi8(total,_mm_add_epi8(_mm_add_epi8(count0,count1),
                _mm_add_epi8(count2,count3)));
    }

    // Reduce 16*8b->{-,-,-,16b,-,-,-,16b}
    const __m128i ZERO = _mm_setzero_si128();
    return _mm_sad_epu8(total,ZERO);
}

// ****************************************************************************
static inline int popcount_sse2_8bit(unsigned* buf,int n) {
  int N = n/4;
  __m128i count32 = _mm_setzero_si128();
  // 2^5 loop iters might overflow 8-bit counter, so
  // cap it at 2^4 iters per chunk
  const int inner_maxits = 16;
  while (N > inner_maxits) {
    count32 = _mm_add_epi32(count32,popcount_sse2_8_helper_1(buf,inner_maxits));
    buf += inner_maxits*4;
    N -= inner_maxits;
  }
  if (N > 0) count32 = _mm_add_epi32(count32,popcount_sse2_8_helper_1(buf,N));

  // Layout coming from PSADBW accumulation is 2*{0,32}: 0 S1 0 S0
  int count;
  _mm_store_ss((float*)&count,(__m128)(_mm_add_epi32(count32,_mm_shuffle_epi32(count32,_MM_SHUFFLE(2,2,2,2)))));
  return count;
}

#endif
