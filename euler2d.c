#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>




// ============================================================================
// Utility macros
// ============================================================================
#define MAX_COMMAND_LEN 1024
#define MAX_DG_ORDER 11
#define NUM_FIELDS 5
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




// ============================================================================
// Runtime configuration
// ============================================================================
static int num_zones_i = 10;
static int num_zones_j = 10;
static int num_zones_k = 10;
static int dg_order = 3;




// ============================================================================
// 7 dimensional array data structure and API.
// ============================================================================
struct Array7
{
    double *ptr;
    int which;
    int shape[7];
    int strides[7];
};

// Array codes
// --------------------------------------------------------
#define A_GAUSS_GRID   1
#define A_LOBATTO_GRID 2
#define A_WGTS         3
#define A_WGTS_CACHE   4
#define A_WGTS_DELTA   5
#define A_PRIM         6
#define A_CONS         7
#define A_CONS_DELTA   8
#define A_CONS_SRF_I   9
#define A_CONS_SRF_J   10
#define A_CONS_SRF_K   11
#define A_FLUX_I       12
#define A_FLUX_J       13
#define A_FLUX_K       14
#define A_FLUX_GOD_I   15
#define A_FLUX_GOD_J   16
#define A_FLUX_GOD_K   17
#define A_TRZN         18
#define A_NUM_ARRAYS   19
static double* global_array[A_NUM_ARRAYS]; // C: elements are initialized to NULL

// Function prototypes
// --------------------------------------------------------
struct Array7 array7_make(int array);
const char* array_name(int array);
void array_shape(int array, int* shape);
size_t array_len(int array);
struct Array7 array_require_alloc(int array);
int array_clear(int array);

// Function implementations
// --------------------------------------------------------
struct Array7 array7_make(int array)
{
    assert(global_array[array]);
    struct Array7 a;
    array_shape(array, a.shape);
    a.which = array;
    a.ptr = global_array[array];
    a.strides[6] = 1;
    a.strides[5] = a.strides[6] * a.shape[6];
    a.strides[4] = a.strides[5] * a.shape[5];
    a.strides[3] = a.strides[4] * a.shape[4];
    a.strides[2] = a.strides[3] * a.shape[3];
    a.strides[1] = a.strides[2] * a.shape[2];
    return a;
}

const char* array_name(int array)
{
    switch (array) {
    case A_GAUSS_GRID: return "GRID";
    case A_LOBATTO_GRID: return "GRID";
    case A_WGTS: return "WGTS";
    case A_WGTS_CACHE: return "CACHE";
    case A_WGTS_DELTA: return "DELTA";
    case A_PRIM: return "PRIM";
    case A_CONS: return "CONS";
    case A_CONS_DELTA: return "DELTA";
    case A_CONS_SRF_I: return "SRF_I";
    case A_CONS_SRF_J: return "SRF_J";
    case A_CONS_SRF_K: return "SRF_K";
    case A_FLUX_I: return "FLUX_I";
    case A_FLUX_J: return "FLUX_J";
    case A_FLUX_K: return "FLUX_K";
    case A_FLUX_GOD_I: return "GOD_I";
    case A_FLUX_GOD_J: return "GOD_J";
    case A_FLUX_GOD_K: return "GOD_K";
    case A_TRZN: return "TRZN";
    }
    assert(0);
}

void array_shape(int array, int* shape)
{
    switch (array) {

    // grid points
    // ----------------------------------------------------
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
    // ----------------------------------------------------
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
    // ----------------------------------------------------
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
    // ----------------------------------------------------
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
    // ----------------------------------------------------
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

size_t array_len(int array)
{
    int n[7];
    array_shape(array, n);
    return n[0] * n[1] * n[2] * n[3] * n[4] * n[5] * n[6];
}

struct Array7 array_require_alloc(int array)
{
    double** ptr = &global_array[array];

    if (*ptr == NULL) {
        size_t elem = array_len(array);
        size_t size = elem * sizeof(double);
        *ptr = malloc(size);

        for (size_t a = 0; a < elem; ++a) {
            (*ptr)[elem] = 0.0;
        }
    }
    return array7_make(array);
}

int array_clear(int array)
{
    free(global_array[array]);
    global_array[array] = NULL;
    return 0;
}

#define GET7(a, i, j, k, r, s, t, q) \
 a.ptr[(( \
    (i < 0 || i >= a.shape[0]) || \
    (j < 0 || j >= a.shape[1]) || \
    (k < 0 || k >= a.shape[2]) || \
    (r < 0 || r >= a.shape[3]) || \
    (s < 0 || s >= a.shape[4]) || \
    (t < 0 || t >= a.shape[5]) || \
    (q < 0 || q >= a.shape[6])) && bounds_error_set(a.which, i, j, k, r, s, t, q)) \
 + (i) * a.strides[0] \
 + (j) * a.strides[1] \
 + (k) * a.strides[2] \
 + (r) * a.strides[3] \
 + (s) * a.strides[4] \
 + (t) * a.strides[5] \
 + (q) * a.strides[6]]




// ============================================================================
// Error reporting for out-of-bounds index.
// ============================================================================
static int BOUNDS_ERROR_ARRAY = 0;
static int BOUNDS_ERROR_INDEX[7];

static int bounds_error_set(int array, int i, int j, int k, int r, int s, int t, int q)
{
    BOUNDS_ERROR_ARRAY = array;
    BOUNDS_ERROR_INDEX[0] = i;
    BOUNDS_ERROR_INDEX[1] = j;
    BOUNDS_ERROR_INDEX[2] = k;
    BOUNDS_ERROR_INDEX[3] = r;
    BOUNDS_ERROR_INDEX[4] = s;
    BOUNDS_ERROR_INDEX[5] = t;
    BOUNDS_ERROR_INDEX[6] = q;
    return 0;
}

static int bounds_error_report()
{
    if (BOUNDS_ERROR_ARRAY == 0)
    {
        return 0;
    }
    else
    {
        printf("out-of-bounds at index (%d %d %d %d %d %d %d) on array %s\n",
            BOUNDS_ERROR_INDEX[0],
            BOUNDS_ERROR_INDEX[1],
            BOUNDS_ERROR_INDEX[2],
            BOUNDS_ERROR_INDEX[3],
            BOUNDS_ERROR_INDEX[4],
            BOUNDS_ERROR_INDEX[5],
            BOUNDS_ERROR_INDEX[6],
            array_name(BOUNDS_ERROR_ARRAY));
        return 1;
    }
}




// ============================================================================
// 7 dimensional array data structure and functions to work with it.
// ============================================================================
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
    struct Array7 a = array_require_alloc(A_CONS);
    double x = GET7(a, 0, 0, 0, 0, 0, 0, 5);

    if (bounds_error_report())
    {
        exit(1);
    }
    return 0;
}
