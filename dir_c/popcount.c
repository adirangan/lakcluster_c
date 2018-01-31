static inline unsigned int popcount_uchar_array(unsigned char *wp,unsigned long wl)
{
  unsigned char *we=&(wp[wl]);
  unsigned int p=0;
  while (wp<we){ p += popcount_uchar[*wp++];}
  return p;
}

typedef union {
  __m128i vi;
  unsigned long u8[2];
  unsigned int u4[4];
} __uni16;

// SSE2 implementations of Lauradoux/Walisch popcount, combined with xor to
// handle Hamming distance, and masking to handle missingness.
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX / BITCT) is a multiple of 3, and no greater than 30.

static inline long long int popcount(__m128i** mem1p, __m128i** maskp, __m128i** maskp_end) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      count2 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      half1 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow
    // 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    d += (unsigned int)(acc.u8[0] + acc.u8[1]);
  } while ((*maskp) < (*maskp_end));
  return d;
}

static inline long long int popcount_and(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end) {
  int verbose=0;
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);  
    if (verbose){ printf(" %% do start *mem1p %p *mem2p %p *maskp %p maskp_tmp %p *maskp_end %p\n",(*mem1p),(*mem2p),(*maskp),maskp_tmp,(*maskp_end));}
    acc.vi = _mm_setzero_si128();
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      if (verbose){ printf(" %% %% *maskp %p\n",(*maskp));}
      count2 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      if (verbose){ printf(" %% %% *maskp %p\n",(*maskp));}
      half1 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      if (verbose){ printf(" %% %% *maskp %p\n",(*maskp));}
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
      if (verbose){ printf(" %% %% *maskp %p\n",(*maskp));}
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    d += (unsigned int)(acc.u8[0] + acc.u8[1]);
    if (verbose){ printf(" %% do final *mem1p %p *mem2p %p *maskp %p maskp_tmp %p *maskp_end %p\n",(*mem1p),(*mem2p),(*maskp),maskp_tmp,(*maskp_end));}
  } while ((*maskp) < (*maskp_end));
  return d;
}

static inline long long int popcount_xor(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      count2 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half1 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow
    // 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    d += (unsigned int)(acc.u8[0] + acc.u8[1]);
  } while ((*maskp) < (*maskp_end));
  return d;
}

static inline long long int popcount_notxorxor(__m128i** mem1p, __m128i** memSp,__m128i** mem2p, __m128i** maskp, __m128i** maskp_end) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_andnot_si128(_mm_xor_si128(_mm_xor_si128(*((*mem1p)++),*((*memSp)++)), *((*mem2p)++)), *((*maskp)++));
      count2 = _mm_andnot_si128(_mm_xor_si128(_mm_xor_si128(*((*mem1p)++),*((*memSp)++)), *((*mem2p)++)), *((*maskp)++));
      half1  = _mm_andnot_si128(_mm_xor_si128(_mm_xor_si128(*((*mem1p)++),*((*memSp)++)), *((*mem2p)++)), *((*maskp)++));
      half2  = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1  = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow
    // 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    d += (unsigned int)(acc.u8[0] + acc.u8[1]);
  } while ((*maskp) < (*maskp_end));
  return d;
}

static inline double popcount_lf(__m128i** mem1p, __m128i** maskp, __m128i** maskp_end, double **dinp) {
  /* Here we calculate popcount(*mem1p), with interval (j+0)*POPLENGTH-->(j+1)*POPLENGTH-1 weighted by (*dinp)[j]. ;
     We use mask from *maskp to *maskp_end. ;
     Note that, as stated above, the vector *dinp is given only on a coarse-grid with grid-step POPLENGTH. ;
  */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  double dout=0; 
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      count2 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      half1 = _mm_and_si128(*((*mem1p)++), *((*maskp)++));
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow
    // 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    dout += (double)(*(*dinp)++) * (double)((unsigned int)(acc.u8[0] + acc.u8[1]));
  } while ((*maskp) < (*maskp_end));
  return dout;
}

static inline double popcount_xor_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end, double **dinp) {
  /* Here we calculate popcount(xor(*mem1p,*mem2p)), with interval (j+0)*POPLENGTH-->(j+1)*POPLENGTH-1 weighted by (*dinp)[j]. ;
     We use mask from *maskp to *maskp_end. ;
     Note that, as stated above, the vector *dinp is given only on a coarse-grid with grid-step POPLENGTH. ;
  */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  double dout=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      count2 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half1 = _mm_and_si128(_mm_xor_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow
    // 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    dout += (double)(*(*dinp)++) * (double)((unsigned int)(acc.u8[0] + acc.u8[1]));
  } while ((*maskp) < (*maskp_end));
  return dout;
}

static inline double popcount_and_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end, double **dinp) {
  /* Here we calculate popcount(and(*mem1p,*mem2p)), with interval (j+0)*POPLENGTH-->(j+1)*POPLENGTH-1 weighted by (*dinp)[j]. ;
     We use mask from *maskp to *maskp_end. ;
     Note that, as stated above, the vector *dinp is given only on a coarse-grid with grid-step POPLENGTH. ;
  */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* maskp_tmp = NULL;
  double dout=0;
  do{
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);  
    acc.vi = _mm_setzero_si128();
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      count2 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half1 = _mm_and_si128(_mm_and_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1 and
      // count2 store a partial bitcount covering themselves AND another bit from
      // elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
    // are guaranteed to be <= 120, thus adding two together does not overflow 255.
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    dout += (double)(*(*dinp)++) * (double)((unsigned int)(acc.u8[0] + acc.u8[1]));
  } while ((*maskp) < (*maskp_end));
  return dout;
}

static inline double popcount_pm0(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end) {
  /* here we calculate transpose(*mem1p) * (*mem2p) ; */
  /* We use mask from *maskp to *maskp_end. ; */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1_pos, count2_pos, half1_pos, half2_pos;
  __uni16 acc_pos;
  __m128i count1_neg, count2_neg, half1_neg, half2_neg;
  __uni16 acc_neg;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  double dout=0;
  do {
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc_pos.vi = _mm_setzero_si128(); acc_neg.vi = _mm_setzero_si128();
    do {
      count1_pos = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      count1_neg = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      count2_pos = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      count2_neg = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half1_pos  = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      half1_neg  = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half2_pos = _mm_and_si128(_mm_srli_epi64(half1_pos, 1), m1);
      half1_pos = _mm_and_si128(half1_pos, m1);
      count1_pos = _mm_sub_epi64(count1_pos, _mm_and_si128(_mm_srli_epi64(count1_pos, 1), m1));
      count2_pos = _mm_sub_epi64(count2_pos, _mm_and_si128(_mm_srli_epi64(count2_pos, 1), m1));
      count1_pos = _mm_add_epi64(count1_pos, half1_pos);
      count2_pos = _mm_add_epi64(count2_pos, half2_pos);
      count1_pos = _mm_add_epi64(_mm_and_si128(count1_pos, m2), _mm_and_si128(_mm_srli_epi64(count1_pos, 2), m2));
      count1_pos = _mm_add_epi64(count1_pos, _mm_add_epi64(_mm_and_si128(count2_pos, m2), _mm_and_si128(_mm_srli_epi64(count2_pos, 2), m2)));
      acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_add_epi64(_mm_and_si128(count1_pos, m4), _mm_and_si128(_mm_srli_epi64(count1_pos, 4), m4)));
      half2_neg = _mm_and_si128(_mm_srli_epi64(half1_neg, 1), m1);
      half1_neg = _mm_and_si128(half1_neg, m1);
      count1_neg = _mm_sub_epi64(count1_neg, _mm_and_si128(_mm_srli_epi64(count1_neg, 1), m1));
      count2_neg = _mm_sub_epi64(count2_neg, _mm_and_si128(_mm_srli_epi64(count2_neg, 1), m1));
      count1_neg = _mm_add_epi64(count1_neg, half1_neg);
      count2_neg = _mm_add_epi64(count2_neg, half2_neg);
      count1_neg = _mm_add_epi64(_mm_and_si128(count1_neg, m2), _mm_and_si128(_mm_srli_epi64(count1_neg, 2), m2));
      count1_neg = _mm_add_epi64(count1_neg, _mm_add_epi64(_mm_and_si128(count2_neg, m2), _mm_and_si128(_mm_srli_epi64(count2_neg, 2), m2)));
      acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_add_epi64(_mm_and_si128(count1_neg, m4), _mm_and_si128(_mm_srli_epi64(count1_neg, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc_pos.vi = _mm_add_epi64(_mm_and_si128(acc_pos.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_pos.vi, 8), m8));
    acc_neg.vi = _mm_add_epi64(_mm_and_si128(acc_neg.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_neg.vi, 8), m8));
#else
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 8)), m8);
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 8)), m8);
#endif
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 16)), m16);
    acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 32));
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 16)), m16);
    acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 32));
    d = (int)(acc_pos.u8[0] + acc_pos.u8[1] - acc_neg.u8[0] - acc_neg.u8[1]);
    dout += (double)(d);
  } while ((*maskp) < (*maskp_end));
  return dout;
}

static inline long long int popcount_pmpm0(__m128i** mem0p,__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end) {
  /* here we calculate transpose(*mem0p .* *mem1p) * (*mem2p), using mask *maskp-->*maksp_end */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1_pos, count2_pos, half1_pos, half2_pos;
  __uni16 acc_pos;
  __m128i count1_neg, count2_neg, half1_neg, half2_neg;
  __uni16 acc_neg;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  do {
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc_pos.vi = _mm_setzero_si128(); acc_neg.vi = _mm_setzero_si128();
    do {
      count1_pos = _mm_and_si128(_mm_andnot_si128(_mm_xor_si128(*((*mem0p)  ), *((*mem1p)  )), *((*mem2p)  )), *((*maskp)  ));
      count1_neg = _mm_and_si128(   _mm_and_si128(_mm_xor_si128(*((*mem0p)++), *((*mem1p)++)), *((*mem2p)++)), *((*maskp)++));
      count2_pos = _mm_and_si128(_mm_andnot_si128(_mm_xor_si128(*((*mem0p)  ), *((*mem1p)  )), *((*mem2p)  )), *((*maskp)  ));
      count2_neg = _mm_and_si128(   _mm_and_si128(_mm_xor_si128(*((*mem0p)++), *((*mem1p)++)), *((*mem2p)++)), *((*maskp)++));
      half1_pos  = _mm_and_si128(_mm_andnot_si128(_mm_xor_si128(*((*mem0p)  ), *((*mem1p)  )), *((*mem2p)  )), *((*maskp)  ));
      half1_neg  = _mm_and_si128(   _mm_and_si128(_mm_xor_si128(*((*mem0p)++), *((*mem1p)++)), *((*mem2p)++)), *((*maskp)++));
      half2_pos = _mm_and_si128(_mm_srli_epi64(half1_pos, 1), m1);
      half1_pos = _mm_and_si128(half1_pos, m1);
      count1_pos = _mm_sub_epi64(count1_pos, _mm_and_si128(_mm_srli_epi64(count1_pos, 1), m1));
      count2_pos = _mm_sub_epi64(count2_pos, _mm_and_si128(_mm_srli_epi64(count2_pos, 1), m1));
      count1_pos = _mm_add_epi64(count1_pos, half1_pos);
      count2_pos = _mm_add_epi64(count2_pos, half2_pos);
      count1_pos = _mm_add_epi64(_mm_and_si128(count1_pos, m2), _mm_and_si128(_mm_srli_epi64(count1_pos, 2), m2));
      count1_pos = _mm_add_epi64(count1_pos, _mm_add_epi64(_mm_and_si128(count2_pos, m2), _mm_and_si128(_mm_srli_epi64(count2_pos, 2), m2)));
      acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_add_epi64(_mm_and_si128(count1_pos, m4), _mm_and_si128(_mm_srli_epi64(count1_pos, 4), m4)));
      half2_neg = _mm_and_si128(_mm_srli_epi64(half1_neg, 1), m1);
      half1_neg = _mm_and_si128(half1_neg, m1);
      count1_neg = _mm_sub_epi64(count1_neg, _mm_and_si128(_mm_srli_epi64(count1_neg, 1), m1));
      count2_neg = _mm_sub_epi64(count2_neg, _mm_and_si128(_mm_srli_epi64(count2_neg, 1), m1));
      count1_neg = _mm_add_epi64(count1_neg, half1_neg);
      count2_neg = _mm_add_epi64(count2_neg, half2_neg);
      count1_neg = _mm_add_epi64(_mm_and_si128(count1_neg, m2), _mm_and_si128(_mm_srli_epi64(count1_neg, 2), m2));
      count1_neg = _mm_add_epi64(count1_neg, _mm_add_epi64(_mm_and_si128(count2_neg, m2), _mm_and_si128(_mm_srli_epi64(count2_neg, 2), m2)));
      acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_add_epi64(_mm_and_si128(count1_neg, m4), _mm_and_si128(_mm_srli_epi64(count1_neg, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc_pos.vi = _mm_add_epi64(_mm_and_si128(acc_pos.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_pos.vi, 8), m8));
    acc_neg.vi = _mm_add_epi64(_mm_and_si128(acc_neg.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_neg.vi, 8), m8));
#else
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 8)), m8);
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 8)), m8);
#endif
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 16)), m16);
    acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 32));
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 16)), m16);
    acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 32));
    d += (int)(unsigned int)(acc_pos.u8[0] + acc_pos.u8[1]) - (int)(unsigned int)(acc_neg.u8[0] + acc_neg.u8[1]);
  } while ((*maskp) < (*maskp_end));
  return d;
}

static inline double popcount_pm0_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end,double **dinp) {
  /* here we calculate transpose(*mem1p) * (*mem2p), with interval (j+0)*POPLENGTH-->(j+1)*POPLENGTH-1 weighted by (*dinp)[j]. ; */
  /*      We use mask from *maskp to *maskp_end. ; */
  /*      Note that, as stated above, the vector *dinp is given only on a coarse-grid with grid-step POPLENGTH. ; */
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1_pos, count2_pos, half1_pos, half2_pos;
  __uni16 acc_pos;
  __m128i count1_neg, count2_neg, half1_neg, half2_neg;
  __uni16 acc_neg;
  __m128i* maskp_tmp = NULL;
  long long int d=0;
  double dout=0;
  do {
    maskp_tmp = &((*maskp)[MULTIPLEX_2DIST / 128]);
    acc_pos.vi = _mm_setzero_si128(); acc_neg.vi = _mm_setzero_si128();
    do {
      count1_pos = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      count1_neg = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      count2_pos = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      count2_neg = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half1_pos  = _mm_and_si128(   _mm_and_si128(*((*mem1p)  ), *((*mem2p)  )), *((*maskp)  ));
      half1_neg  = _mm_and_si128(_mm_andnot_si128(*((*mem1p)++), *((*mem2p)++)), *((*maskp)++));
      half2_pos = _mm_and_si128(_mm_srli_epi64(half1_pos, 1), m1);
      half1_pos = _mm_and_si128(half1_pos, m1);
      count1_pos = _mm_sub_epi64(count1_pos, _mm_and_si128(_mm_srli_epi64(count1_pos, 1), m1));
      count2_pos = _mm_sub_epi64(count2_pos, _mm_and_si128(_mm_srli_epi64(count2_pos, 1), m1));
      count1_pos = _mm_add_epi64(count1_pos, half1_pos);
      count2_pos = _mm_add_epi64(count2_pos, half2_pos);
      count1_pos = _mm_add_epi64(_mm_and_si128(count1_pos, m2), _mm_and_si128(_mm_srli_epi64(count1_pos, 2), m2));
      count1_pos = _mm_add_epi64(count1_pos, _mm_add_epi64(_mm_and_si128(count2_pos, m2), _mm_and_si128(_mm_srli_epi64(count2_pos, 2), m2)));
      acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_add_epi64(_mm_and_si128(count1_pos, m4), _mm_and_si128(_mm_srli_epi64(count1_pos, 4), m4)));
      half2_neg = _mm_and_si128(_mm_srli_epi64(half1_neg, 1), m1);
      half1_neg = _mm_and_si128(half1_neg, m1);
      count1_neg = _mm_sub_epi64(count1_neg, _mm_and_si128(_mm_srli_epi64(count1_neg, 1), m1));
      count2_neg = _mm_sub_epi64(count2_neg, _mm_and_si128(_mm_srli_epi64(count2_neg, 1), m1));
      count1_neg = _mm_add_epi64(count1_neg, half1_neg);
      count2_neg = _mm_add_epi64(count2_neg, half2_neg);
      count1_neg = _mm_add_epi64(_mm_and_si128(count1_neg, m2), _mm_and_si128(_mm_srli_epi64(count1_neg, 2), m2));
      count1_neg = _mm_add_epi64(count1_neg, _mm_add_epi64(_mm_and_si128(count2_neg, m2), _mm_and_si128(_mm_srli_epi64(count2_neg, 2), m2)));
      acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_add_epi64(_mm_and_si128(count1_neg, m4), _mm_and_si128(_mm_srli_epi64(count1_neg, 4), m4)));
    } while ((*maskp) < (maskp_tmp));
#if MULTIPLEX_DIST > 960
    acc_pos.vi = _mm_add_epi64(_mm_and_si128(acc_pos.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_pos.vi, 8), m8));
    acc_neg.vi = _mm_add_epi64(_mm_and_si128(acc_neg.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_neg.vi, 8), m8));
#else
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 8)), m8);
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 8)), m8);
#endif
    acc_pos.vi = _mm_and_si128(_mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 16)), m16);
    acc_pos.vi = _mm_add_epi64(acc_pos.vi, _mm_srli_epi64(acc_pos.vi, 32));
    acc_neg.vi = _mm_and_si128(_mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 16)), m16);
    acc_neg.vi = _mm_add_epi64(acc_neg.vi, _mm_srli_epi64(acc_neg.vi, 32));
    d = (int)(acc_pos.u8[0] + acc_pos.u8[1] - acc_neg.u8[0] - acc_neg.u8[1]);
    dout += (double)(*(*dinp)++) * (double)(d);
  } while ((*maskp) < (*maskp_end));
  return dout;
}
