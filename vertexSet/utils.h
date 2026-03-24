#if !defined(UTILS)
#define UTILS
static inline int Vertexset_Log2_floor(int x){
    int i;
    for (i = 0; i < 32; i++) {
        x = x >> 1;
        if (x <= 0){
            break;
        }
    }
    return i;
};

static inline int Vertexset_reverseBits(int num, int N_bits){
    int reverse_num = 0;
    for (int i = 0; i < N_bits; i++){
        if ((num & (1<<i))) {
            reverse_num |= 1 << (N_bits - 1 - i);
        }
    }
    return reverse_num;
    
}

static inline int Vertexset_nearestLog2(int p){
    int count = 0;
    while (p != 1) {
        p >>= 1;
        count++;
    }
    return (1<<count);
}
#endif // UTILS
