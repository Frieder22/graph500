#define ulong_bits_CUSTOM 64
#define ulong_mask_CUSTOM &63
#define ulong_shift_CUSTOM >>6
#define SET_VISITED_CUSTOM(array, v) do {array[VERTEX_LOCAL((v)) ulong_shift] |= (1UL << (VERTEX_LOCAL((v)) ulong_mask));} while (0)
#define SET_VISITEDLOC_CUSTOM(array, v) do {array[(v) ulong_shift] |= (1ULL << ((v) ulong_mask));} while (0)
#define TEST_VISITED_CUSTOM(array, v) ((array[VERTEX_LOCAL((v)) ulong_shift] & (1UL << (VERTEX_LOCAL((v)) ulong_mask))) != 0)
#define TEST_VISITEDLOC_CUSTOM(array, v) ((array[(v) ulong_shift] & (1ULL << ((v) ulong_mask))) != 0)
#define CLEAN_VISITED_CUSTOM(array, size)  memset(array,0,size*sizeof(unsigned long));