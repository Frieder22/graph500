#if !defined(BITMAP)
#define BITMAP

#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#define ulong_bits 64
#define ulong_mask &63
#define ulong_divide >>6

static inline void Bitmap_Set(unsigned long long *array, uint32_t vertex){
    array[vertex ulong_divide] |= (1ULL << (vertex ulong_mask));
}

static inline bool Bitmap_Test(unsigned long long *array, uint32_t vertex){
    return (array[vertex ulong_divide] & (1ULL << (vertex ulong_mask))) != 0;
}

static inline void Bitmap_Clean(unsigned long long *array, size_t size){
    memset(array, 0, size * sizeof(unsigned long long));
}


#endif // BITMAP