/*
Internal_mxArray.h
Matlab version: 2010B
*/

typedef struct {
    void *reserved;
    int reserved1[2];
    void *reserved2;
    size_t number_of_dims;
    unsigned int reserved3;
    struct {
        unsigned int flag0 : 1;
        unsigned int flag1 : 1;
        unsigned int flag2 : 1;
        unsigned int flag3 : 1;
        unsigned int flag4 : 1;
        unsigned int flag5 : 1;
        unsigned int flag6 : 1;
        unsigned int flag7 : 1;
        unsigned int flag7a: 1;
        unsigned int flag8 : 1;
        unsigned int flag9 : 1;
        unsigned int flag10 : 1;
        unsigned int flag11 : 4;
        unsigned int flag12 : 8;
        unsigned int flag13 : 8;
    } flags;
    size_t reserved4[2];
    union {
        struct {
            void *pdata;
            void *pimag_data;
            void *reserved5;
            size_t reserved6[3];
        } number_array;
    } data;
} Internal_mxArray; 