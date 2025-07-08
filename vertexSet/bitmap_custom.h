#if !defined(BITMAP)
#define BITMAP

#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#define ulong_bits 64
#define ulong_mask &63
#define ulong_divide >>6

/**
 * Sets a value in the bitmap. If the value is already set, the 
 * value keeps set.
 * @param array array, that represents the bitmap.
 * @param vertex value, that should be set
 */
static inline void Bitmap_Set(unsigned long long* const array, const uint32_t vertex){
    array[vertex ulong_divide] |= (1ULL << (vertex ulong_mask));
}

/**
 * Sets the value, to a word shifted by 'shift'. The shift is done in
 * mod size_bitarray (circular). To get the correct value later, 
 * Bitmap_Test_Shifted() musst be used with the same shift input.
 * @param array array, that represents the bitmap.
 * @param vertex value, that should be set
 * @param shift size of shift
 * @param size_bitarray number of words the bitarray has (elements in array)
 */
static inline void Bitmap_Set_Shifted(unsigned long long* const array, const uint32_t vertex, const int shift, const int size_bitarray){
    array[((vertex ulong_divide) - shift + size_bitarray) % size_bitarray] |= (1ULL << (vertex ulong_mask));
}

/**
 * Tests, if a value is set in bitmap.
 * @param array array, that represents the bitmap.
 * @param vertex value to test.
 */
static inline bool Bitmap_Test( const unsigned long long* const array, const uint32_t vertex){
    return (array[vertex ulong_divide] & (1ULL << (vertex ulong_mask))) != 0;
}

/**
 * Tets a value, that was set with Bitmap_Set_Shifted(). The param shift
 * musst be the same here.
 * @param array array, that represents the bitmap.
 * @param vertex value, that should be checked
 * @param shift size of shift
 * @param size_bitarray number of words the bitarray has (elements in array)
 */
static inline bool Bitmap_Test_Shifted(const unsigned long long* const array, const uint32_t vertex, const int shift, const int size_bitarray){
    return (array[((vertex ulong_divide) - shift + size_bitarray) % size_bitarray] & (1ULL << (vertex ulong_mask))) != 0;
}

/**
 * Resets the bitmap (deletes all set elements).
 * @param array array, that represents the bitmap.
 * @param size_bitarray number of words the bitarray has (elements in array)
 */
static inline void Bitmap_Clean(unsigned long long* const array, const size_t size_bitarray){
    memset(array, 0, size_bitarray * sizeof(unsigned long long));
}


#endif // BITMAP