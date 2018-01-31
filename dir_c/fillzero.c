inline void fill_uchar_zero(unsigned char* iarr, size_t size) { unsigned char*iptr = &(iarr[size]); while (iarr<iptr){ *iarr++=0;}}
inline void fill_uchar_ones(unsigned char* iarr, size_t size) { unsigned char*iptr = &(iarr[size]); while (iarr<iptr){ *iarr++=255;}}
inline void fill_long_zero(long* larr, size_t size) { long* lptr = &(larr[size]); while (larr < lptr) { *larr++ = 0;} }
