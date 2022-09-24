#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>




// --------------------------------------------------------
// Program control
// --------------------------------------------------------
int sim_reset();
int sim_set_terminal(const char *terminal_str);
int sim_set_binary_output(const char *terminal_str);




// --------------------------------------------------------
// Array API
// --------------------------------------------------------

// Return an array struct for the given array number.
//
// The array is returned even if its memory is not allocated, or it would have a
// zero-shape.
//
struct Array array_make(int array_num);

// Return the name of the array for the given array number.
//
const char *array_name(int array_num);

// Return the 7-component shape of the array with the given array number.
//
void array_shape(int array_num, int shape[7]);

// Return the number of elements in the array with the given array number.
//
// Multiply this number by sizeof(double) to get the allocation size for the
// corresponding array.
//
size_t array_len(int array_num);

// Require that the given array is allocated, and return its struct instance.
//
struct Array array_require_alloc(int array_num);

// Deallocate an array.
//
int array_clear(int array_num);

// Write an array to the current binary output stream.
//
int array_write(int array_num);

// Indicate an array was indexed out of bounds.
//
// This function writes a message to stderr and exits the program
//
int array_bounds_error_exit(int array, int i, int j, int k, int r, int s, int t,
                            int q);




// =============================================================================
// Utility macros
// =============================================================================
#define MAX_COMMAND_LEN 1024
#define MAX_DG_ORDER    11
#define NUM_FIELDS      5
#define MIN2(a, b)      ((a) < (b) ? (a) : (b))
#define MAX2(a, b)      ((a) > (b) ? (a) : (b))
#define MIN3(a, b, c)   MIN2(a, MIN2(b, c))
#define MAX3(a, b, c)   MAX2(a, MAX2(b, c))
#define TRY(n)                                                                 \
    do {                                                                       \
        int res = n;                                                           \
        if (res)                                                               \
            return res;                                                        \
    } while (0)

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

// All the allocations in use throughout the program.
//
// Note: since this is a static array, the pointers are initialized to NULL per
// C standard.
static double *global_array[A_NUM_ARRAYS];




// =============================================================================
// Runtime configuration
// =============================================================================
static FILE *terminal;
static FILE *binary_output;
static int dg_order;
static int num_zones_i;
static int num_zones_j;
static int num_zones_k;
static double domain_x0;
static double domain_x1;
static double domain_y0;
static double domain_y1;
static double domain_z0;
static double domain_z1;




// =============================================================================
// User interaction
// =============================================================================
int sim_reset()
{
    sim_set_terminal("stdout");
    sim_set_binary_output(NULL);

    dg_order = 3;
    num_zones_i = 10;
    num_zones_j = 10;
    num_zones_k = 1;
    domain_x0 = 0.0;
    domain_x1 = 1.0;
    domain_y0 = 0.0;
    domain_y1 = 1.0;
    domain_z0 = 0.0;
    domain_z1 = 1.0;

    for (int n = 0; n < A_NUM_ARRAYS; ++n) {
        array_clear(n);
    }
    return 0;
}

int sim_set_terminal(const char *terminal_str)
{
    FILE *new_terminal = NULL;

    if (terminal != NULL && terminal != stdout) {
        fclose(terminal);
    }
    if (strcmp(terminal_str, "stdout") == 0) {
        terminal = stdout;
    } else if ((new_terminal = fopen(terminal_str, "w")) == NULL) {
        fprintf(stderr, "[error] unable to open terminal file %s\n",
                terminal_str);
        return 1;
    }
    terminal = new_terminal;
    return 0;
}

int sim_set_binary_output(const char *file_str)
{
    FILE *new_binary_output = NULL;

    if (binary_output != NULL) {
        fclose(binary_output);
    } else if (file_str &&
               (new_binary_output = fopen(file_str, "wb")) == NULL) {
        fprintf(stderr, "[error] unable to open binary output file %s\n",
                file_str);
        return 1;
    }
    binary_output = new_binary_output;
    return 0;
}




// =============================================================================
// 7 dimensional array data structure and API.
// =============================================================================
struct Array {
    double *ptr;
    int array_num;
    int shape[7];
    int strides[7];
};

// Function implementations
// --------------------------------------------------------
struct Array array_make(int array_num)
{
    struct Array a;
    array_shape(array_num, a.shape);
    a.array_num = array_num;
    a.ptr = global_array[array_num];
    a.strides[6] = 1;
    a.strides[5] = a.strides[6] * a.shape[6];
    a.strides[4] = a.strides[5] * a.shape[5];
    a.strides[3] = a.strides[4] * a.shape[4];
    a.strides[2] = a.strides[3] * a.shape[3];
    a.strides[1] = a.strides[2] * a.shape[2];
    a.strides[0] = a.strides[1] * a.shape[1];
    return a;
}

const char *array_name(int array_num)
{
    switch (array_num) {
    case A_GAUSS_GRID: return "GAUSS_GRID";
    case A_LOBATTO_GRID: return "LOBATTO_GRID";
    case A_WGTS: return "WGTS";
    case A_WGTS_CACHE: return "WGTS_CACHE";
    case A_WGTS_DELTA: return "WGTS_DELTA";
    case A_PRIM: return "PRIM";
    case A_CONS: return "CONS";
    case A_CONS_DELTA: return "CONS_DELTA";
    case A_CONS_SRF_I: return "CONS_SRF_I";
    case A_CONS_SRF_J: return "CONS_SRF_J";
    case A_CONS_SRF_K: return "CONS_SRF_K";
    case A_FLUX_I: return "FLUX_I";
    case A_FLUX_J: return "FLUX_J";
    case A_FLUX_K: return "FLUX_K";
    case A_FLUX_GOD_I: return "FLUX_GOD_I";
    case A_FLUX_GOD_J: return "FLUX_GOD_J";
    case A_FLUX_GOD_K: return "FLUX_GOD_K";
    case A_TRZN: return "TRZN";
    default: return NULL;
    }
}

void array_shape(int array, int shape[7])
{
    int dg_r = num_zones_i > 1 ? dg_order : 1;
    int dg_s = num_zones_j > 1 ? dg_order : 1;
    int dg_t = num_zones_k > 1 ? dg_order : 1;

    switch (array) {

    // grid points
    // ----------------------------------------------------
    case A_GAUSS_GRID:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_r;
        shape[4] = dg_s;
        shape[5] = dg_t;
        shape[6] = 3; // (x, y, z) coordinates of Gauss point
        return;
    case A_LOBATTO_GRID:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_r + 1;
        shape[4] = dg_s + 1;
        shape[5] = dg_t + 1;
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
        shape[3] = dg_r;
        shape[4] = dg_s;
        shape[5] = dg_t;
        return;

    // hydrodynamics fields
    // ----------------------------------------------------
    case A_PRIM:
    case A_CONS:
    case A_CONS_DELTA:
        shape[0] = num_zones_i;
        shape[1] = num_zones_j;
        shape[2] = num_zones_k;
        shape[3] = dg_r;
        shape[4] = dg_s;
        shape[5] = dg_t;
        shape[6] = NUM_FIELDS;
        return;

    // surface conserved states
    // ----------------------------------------------------
    case A_CONS_SRF_I:
        if (num_zones_i > 1) {
            shape[0] = num_zones_i + 1;
            shape[1] = num_zones_j;
            shape[2] = num_zones_k;
            shape[3] = 2;
            shape[4] = dg_s;
            shape[5] = dg_t;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;
    case A_CONS_SRF_J:
        if (num_zones_j > 1) {
            shape[0] = num_zones_i;
            shape[1] = num_zones_j + 1;
            shape[2] = num_zones_k;
            shape[3] = dg_r;
            shape[4] = 2;
            shape[5] = dg_t;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;
    case A_CONS_SRF_K:
        if (num_zones_k > 1) {
            shape[0] = num_zones_i;
            shape[1] = num_zones_j;
            shape[2] = num_zones_k + 1;
            shape[3] = dg_r;
            shape[4] = dg_s;
            shape[5] = 2;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;

    // fluxes
    // ----------------------------------------------------
    case A_FLUX_GOD_I:
    case A_FLUX_I:
        if (num_zones_i > 1) {
            shape[0] = num_zones_i + 1;
            shape[1] = num_zones_j;
            shape[2] = num_zones_k;
            shape[3] = 1;
            shape[4] = dg_s;
            shape[5] = dg_t;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;
    case A_FLUX_GOD_J:
    case A_FLUX_J:
        if (num_zones_j > 1) {
            shape[0] = num_zones_i;
            shape[1] = num_zones_j + 1;
            shape[2] = num_zones_k;
            shape[3] = dg_r;
            shape[4] = 1;
            shape[5] = dg_t;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;
    case A_FLUX_GOD_K:
    case A_FLUX_K:
        if (num_zones_k > 1) {
            shape[0] = num_zones_i;
            shape[1] = num_zones_j;
            shape[2] = num_zones_k + 1;
            shape[3] = dg_r;
            shape[4] = dg_s;
            shape[5] = 1;
            shape[6] = NUM_FIELDS;
            return;
        }
        break;
    }

    shape[0] = 0;
    shape[1] = 0;
    shape[2] = 0;
    shape[3] = 0;
    shape[4] = 0;
    shape[5] = 0;
    shape[6] = 0;
}

size_t array_len(int array_num)
{
    int n[7];
    array_shape(array_num, n);
    return n[0] * n[1] * n[2] * n[3] * n[4] * n[5] * n[6];
}

struct Array array_require_alloc(int array_num)
{
    double *ptr = global_array[array_num];

    if (ptr == NULL) {
        size_t elem = array_len(array_num);
        size_t size = elem * sizeof(double);
        ptr = malloc(size);

        for (size_t a = 0; a < elem; ++a) {
            ptr[elem] = 0.0;
        }
        global_array[array_num] = ptr;
    }
    return array_make(array_num);
}

int array_clear(int array_num)
{
    free(global_array[array_num]);
    global_array[array_num] = NULL;
    return 0;
}

int array_write(int array_num)
{
    struct Array a = array_make(array_num);

    if (a.ptr == NULL || binary_output == NULL) {
        return 1;
    }

    fwrite(a.ptr, sizeof(double), array_len(array_num), binary_output);
    fflush(binary_output);
    return 0;
}

int array_bounds_error_exit(int array, int i, int j, int k, int r, int s, int t,
                            int q)
{
    fprintf(stderr,
            "out-of-bounds at index (%d %d %d %d %d %d %d) on array%s\n", i, j,
            k, r, s, t, q, array_name(array));
    exit(1);
    return 0;
}

#define GET7(a, i, j, k, r, s, t, q)                                           \
    a.ptr[(((i < 0 || i >= a.shape[0]) || (j < 0 || j >= a.shape[1]) ||        \
            (k < 0 || k >= a.shape[2]) || (r < 0 || r >= a.shape[3]) ||        \
            (s < 0 || s >= a.shape[4]) || (t < 0 || t >= a.shape[5]) ||        \
            (q < 0 || q >= a.shape[6])) &&                                     \
           array_bounds_error_exit(a.array_num, i, j, k, r, s, t, q)) +        \
          (i)*a.strides[0] + (j)*a.strides[1] + (k)*a.strides[2] +             \
          (r)*a.strides[3] + (s)*a.strides[4] + (t)*a.strides[5] +             \
          (q)*a.strides[6]]

#define GETP6(a, i, j, k, r, s, t) (&GET7(a, i, j, k, r, s, t, 0))

#define FOR_EACH_IJKRST_(ni, nj, nk, nr, ns, nt)                               \
    for (int i = 0; i < ni; ++i)                                               \
        for (int j = 0; j < nj; ++j)                                           \
            for (int k = 0; k < nk; ++k)                                       \
                for (int r = 0; r < nr; ++r)                                   \
                    for (int s = 0; s < ns; ++s)                               \
                        for (int t = 0; t < nt; ++t)

#define FOR_EACH_IJKRST(a)                                                     \
    FOR_EACH_IJKRST_(a.shape[0], a.shape[1], a.shape[2], a.shape[3],           \
                     a.shape[4], a.shape[5])




// =============================================================================
// Functions to generate stencil data
// =============================================================================
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




// =============================================================================
// Grid functions
// =============================================================================
int grid_gauss_init()
{
    printf("[grid_gauss_init]\n");

    struct Array grid = array_require_alloc(A_GAUSS_GRID);
    int ng_i = num_zones_i > 1 ? 1 : 0;
    int ng_j = num_zones_j > 1 ? 1 : 0;
    int ng_k = num_zones_k > 1 ? 1 : 0;
    int dg_r = num_zones_i > 1 ? dg_order : 1;
    int dg_s = num_zones_j > 1 ? dg_order : 1;
    int dg_t = num_zones_k > 1 ? dg_order : 1;
    double dx = (domain_x1 - domain_x0) / (num_zones_i - 2 * ng_i);
    double dy = (domain_y1 - domain_y0) / (num_zones_j - 2 * ng_j);
    double dz = (domain_z1 - domain_z0) / (num_zones_k - 2 * ng_k);

    FOR_EACH_IJKRST (grid) {
        double xsi_x = gauss_quadrature_node(dg_r, r);
        double xsi_y = gauss_quadrature_node(dg_s, s);
        double xsi_z = gauss_quadrature_node(dg_t, t);
        GET7(grid, i, j, k, r, s, t, 0) =
            domain_x0 + dx * (i - ng_i + 0.5 * (1.0 + xsi_x));
        GET7(grid, i, j, k, r, s, t, 1) =
            domain_y0 + dy * (j - ng_j + 0.5 * (1.0 + xsi_y));
        GET7(grid, i, j, k, r, s, t, 2) =
            domain_z0 + dz * (k - ng_k + 0.5 * (1.0 + xsi_z));
    }
    return 0;
}

int grid_lobatto_init()
{
    printf("[grid_lobatto_init]\n");

    struct Array grid = array_require_alloc(A_LOBATTO_GRID);
    int ng_i = num_zones_i > 1 ? 1 : 0;
    int ng_j = num_zones_j > 1 ? 1 : 0;
    int ng_k = num_zones_k > 1 ? 1 : 0;
    int dg_r = num_zones_i > 1 ? dg_order : 1;
    int dg_s = num_zones_j > 1 ? dg_order : 1;
    int dg_t = num_zones_k > 1 ? dg_order : 1;
    double dx = (domain_x1 - domain_x0) / (num_zones_i - 2 * ng_i);
    double dy = (domain_y1 - domain_y0) / (num_zones_j - 2 * ng_j);
    double dz = (domain_z1 - domain_z0) / (num_zones_k - 2 * ng_k);

    FOR_EACH_IJKRST (grid) {
        double xsi_x = lobatto_node(dg_r, r);
        double xsi_y = lobatto_node(dg_s, s);
        double xsi_z = lobatto_node(dg_t, t);
        GET7(grid, i, j, k, r, s, t, 0) =
            domain_x0 + dx * (i - ng_i + 0.5 * (1.0 + xsi_x));
        GET7(grid, i, j, k, r, s, t, 1) =
            domain_y0 + dy * (j - ng_j + 0.5 * (1.0 + xsi_y));
        GET7(grid, i, j, k, r, s, t, 2) =
            domain_z0 + dz * (k - ng_k + 0.5 * (1.0 + xsi_z));
    }
    return 0;
}




// =============================================================================
// Initial conditions
// =============================================================================
int prim_init()
{
    if (!global_array[A_GAUSS_GRID] || !global_array[A_PRIM]) {
        return 1;
    }

    struct Array grid = array_make(A_GAUSS_GRID);
    struct Array prim = array_make(A_PRIM);

    FOR_EACH_IJKRST (grid) {
        double *x = GETP6(grid, i, j, k, r, s, t);
        double *p = GETP6(prim, i, j, k, r, s, t);

        if (x[0] < 0.5) {
            p[0] = 1.0;
            p[1] = 0.0;
            p[2] = 0.0;
            p[3] = 0.0;
            p[4] = 1.0;
        } else {
            p[0] = 0.1;
            p[1] = 0.0;
            p[2] = 0.0;
            p[3] = 0.0;
            p[4] = 0.125;
        }
    }
    return 0;
}




// =============================================================================
// Main
// =============================================================================
int main()
{
    sim_reset();

    grid_gauss_init();
    grid_lobatto_init();
    prim_init();

    sim_set_binary_output("gauss.bin");
    array_write(A_GAUSS_GRID);

    sim_set_binary_output("lobatto.bin");
    array_write(A_LOBATTO_GRID);

    sim_set_binary_output("prim.bin");
    array_write(A_PRIM);
    sim_reset();

    printf("[done]\n");
    return 0;
}
