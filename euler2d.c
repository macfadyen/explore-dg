#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_COMMAND_LEN 1024
#define MAX_DG_ORDER 11
#define NUM_FIELDS 3
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))
#define TRY(n)                                                                 \
    do {                                                                       \
        int res = n;                                                           \
        if (res)                                                               \
            return res;                                                        \
    } while (0)

#define A_GAUSS_GRID   0
#define A_LOBATTO_GRID 1
#define A_WGTS         2
#define A_WGTS_CACHE   3
#define A_WGTS_DELTA   4
#define A_PRIM         5
#define A_CONS         6
#define A_CONS_DELTA   7
#define A_CONS_SRF_I   8
#define A_CONS_SRF_J   9
#define A_CONS_SRF_K   10
#define A_FLUX_I       11
#define A_FLUX_J       12
#define A_FLUX_K       13
#define A_FLUX_GOD_I   14
#define A_FLUX_GOD_J   15
#define A_FLUX_GOD_K   16
#define A_TRZN         17
#define A_NUM_ARRAYS   18
static double* global_array[A_NUM_ARRAYS];


// user-configurable global variables
static int num_zones_i = 10;
static int num_zones_j = 10;
static int num_zones_k = 10;
static int dg_order = 3;

struct Array7
{
    double *ptr;
    int strides[7];
};

#define GET7(a, i, j, k, r, s, t, q) \
 a.ptr[0 \
 + (i) * a.strides[0] \
 + (j) * a.strides[1] \
 + (k) * a.strides[2] \
 + (r) * a.strides[3] \
 + (s) * a.strides[4] \
 + (t) * a.strides[5] \
 + (q) * a.strides[6]]


void array_shape(int array, int* shape)
{
    switch (array) {

    // grid points
    case A_GAUSS_GRID:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_order;
        shape[4] = dg_order;
        shape[5] = dg_order;
        shape[6] = 3; // (x, y, z) coordinates of Gauss point
        return;
    case A_LOBATTO_GRID:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_order + 1;
        shape[4] = dg_order + 1;
        shape[5] = dg_order + 1;
        shape[6] = 3; // (x, y, z) coordinates of Lobatto points
        return;

    // DG weights
    case A_WGTS:
    case A_WGTS_CACHE:
    case A_WGTS_DELTA:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = NUM_FIELDS;
        shape[4] = dg_order;
        shape[5] = dg_order;
        shape[6] = dg_order;
        return;

    // hydrodynamics fields
    case A_PRIM:
    case A_CONS:
    case A_CONS_DELTA:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_order;
        shape[4] = dg_order;
        shape[5] = dg_order;
        shape[6] = NUM_FIELDS;
        return;

    // surface conserved states
    case A_CONS_SRF_I:
        shape[0] = num_zones_i + 1;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = 2;
        shape[4] = dg_order;
        shape[5] = dg_order;
        shape[6] = NUM_FIELDS;
        return;
    case A_CONS_SRF_J:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j + 1;
        shape[2] = num_zones_k;
        shape[3] = dg_order;
        shape[4] = 2;
        shape[5] = dg_order;
        shape[6] = NUM_FIELDS;
        return;
    case A_CONS_SRF_K:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k + 1;
        shape[3] = dg_order;
        shape[4] = dg_order;
        shape[5] = 2;
        shape[6] = NUM_FIELDS;
        return;

    // fluxes
    case A_FLUX_GOD_I:
    case A_FLUX_I:
        shape[0] = num_zones_i + 1;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = 1;
        shape[4] = dg_order;
        shape[5] = dg_order;
        shape[6] = NUM_FIELDS;
        return;
    case A_FLUX_GOD_J:
    case A_FLUX_J:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j + 1;
        shape[2] = num_zones_k;
        shape[3] = dg_order;
        shape[4] = 1;
        shape[5] = dg_order;
        shape[6] = NUM_FIELDS;
        return;
    case A_FLUX_GOD_K:
    case A_FLUX_K:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k + 1;
        shape[3] = dg_order;
        shape[4] = dg_order;
        shape[5] = 1;
        shape[6] = NUM_FIELDS;
        return;
    }

    assert(0);
}

struct Array7 array7_make(int array)
{
    struct Array7 a;
    int shape[7];
    array_shape(array, shape);
    a.ptr = global_array[array];
    a.strides[6] = 1;
    a.strides[5] = a.strides[6] * shape[6];
    a.strides[4] = a.strides[5] * shape[5];
    a.strides[3] = a.strides[4] * shape[4];
    a.strides[2] = a.strides[3] * shape[3];
    a.strides[1] = a.strides[2] * shape[2];
    return a;
}

double gauss_quadrature_node(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return +0.0000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return -0.5773502691896;
        case 1: return +0.5773502691896;
        }
        break;
    case 3:
        switch (index) {
        case 0: return -0.7745966692415;
        case 1: return +0.0000000000000;
        case 2: return +0.7745966692415;
        }
        break;
    case 4:
        switch (index) {
        case 0: return -0.8611363115941;
        case 1: return -0.3399810435849;
        case 2: return +0.3399810435849;
        case 3: return +0.8611363115941;
        }
        break;
    case 5:
        switch (index) {
        case 0: return -0.9061798459387;
        case 1: return -0.5384693101057;
        case 2: return +0.0000000000000;
        case 3: return +0.5384693101057;
        case 4: return +0.9061798459387;
        }
        break;
    }
    return 0.0;
}

double gauss_quadrature_weight(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return +2.0000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return +1.0000000000000;
        case 1: return +1.0000000000000;
        }
        break;
    case 3:
        switch (index) {
        case 0: return +0.5555555555556;
        case 1: return +0.8888888888889;
        case 2: return +0.5555555555556;
        }
        break;
    case 4:
        switch (index) {
        case 0: return +0.3478548451375;
        case 1: return +0.6521451548625;
        case 2: return +0.6521451548625;
        case 3: return +0.3478548451375;
        }
        break;
    case 5:
        switch (index) {
        case 0: return +0.2369268850562;
        case 1: return +0.4786286704994;
        case 2: return +0.5688888888889;
        case 3: return +0.4786286704994;
        case 4: return +0.2369268850562;
        }
        break;
    }
    return 0.0;
}

double lobatto_node(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return -1.0000000000000;
        case 1: return +1.0000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return -1.0000000000000;
        case 1: return +0.0000000000000;
        case 2: return +1.0000000000000;
        }
        break;
    case 3:
        switch (index) {
        case 0: return -1.0000000000000;
        case 1: return -0.4472135955000;
        case 2: return +0.4472135955000;
        case 3: return +1.0000000000000;
        }
        break;
    case 4:
        switch (index) {
        case 0: return -1.0000000000000;
        case 1: return -0.6546536707080;
        case 2: return +0.0000000000000;
        case 3: return +0.6546536707080;
        case 4: return +1.0000000000000;
        }
        break;
    case 5:
        switch (index) {
        case 0: return -1.0000000000000;
        case 1: return -0.7650553239295;
        case 2: return -0.2852315164806;
        case 3: return +0.2852315164806;
        case 4: return +0.7650553239295;
        case 5: return +1.0000000000000;
        }
        break;
    }
    return 0.0;
}

double lobatto_weight(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return +1.0000000000000;
        case 1: return +1.0000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return +0.3333333333333;
        case 1: return +1.3333333333333;
        case 2: return +0.3333333333333;
        }
        break;
    case 3:
        switch (index) {
        case 0: return +0.1666666666667;
        case 1: return +0.8333333333333;
        case 2: return +0.8333333333333;
        case 3: return +0.1666666666667;
        }
        break;
    case 4:
        switch (index) {
        case 0: return +0.1000000000000;
        case 1: return +0.5444444444444;
        case 2: return +0.7111111111111;
        case 3: return +0.5444444444444;
        case 4: return +0.1000000000000;
        }
        break;
    case 5:
        switch (index) {
        case 0: return +0.0666666666667;
        case 1: return +0.3784749562978;
        case 2: return +0.5548583770355;
        case 3: return +0.5548583770355;
        case 4: return +0.3784749562978;
        case 5: return +0.0666666666667;
        }
        break;
    }
    return 0.0;
}

int main()
{
    return 0;
}
