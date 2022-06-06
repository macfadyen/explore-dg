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

static FILE* terminal = NULL;
static int num_zones = 20;
static int order = 3;
static int rk_order = 1;
static double domain_x0 = 0.0;
static double domain_x1 = 1.0;
static double time_phys = 0.0;
static double time_step = 0.0;
static double rk_parameter = 0.0;
static double tci_threshold = 0.001;
static double cfl_parameter = 0.1;
static double time_final = 0.1;
static double adiabatic_gamma = 5.0 / 3.0;

struct timespec timer_start()
{
    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    return start_time;
}

long timer_end(struct timespec start_time)
{
    struct timespec end_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    long nanos = (end_time.tv_sec - start_time.tv_sec) * 1000000000L
        + (end_time.tv_nsec - start_time.tv_nsec);
    return nanos;
}

long factorial(long n)
{
    if (n <= 1) {
        return 1;
    }
    return n * factorial(n - 1);
}

long choose(long n, long k)
{
    return factorial(n) / factorial(k) / factorial(n - k);
}

double legendre_polynomial(int n, double x)
{
    double p = 0.0;

    for (int k = 0; k <= n; ++k) {
        p += choose(n, k) * choose(n + k, k) * pow(0.5 * (x - 1), k);
    }
    return p * sqrt(2 * n + 1);
}

double legendre_polynomial_derivative(int n, double x)
{
    double p = 0.0;

    for (int k = 0; k <= n; ++k) {
        p += choose(n, k) * choose(n + k, k) * 0.5 * k
            * pow(0.5 * (x - 1), k - 1);
    }
    return p * sqrt(2 * n + 1);
}

double gauss_quadrature_node(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return +0.00000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return -0.57735026918963;
        case 1: return +0.57735026918963;
        }
        break;
    case 3:
        switch (index) {
        case 0: return -0.77459666924148;
        case 1: return +0.00000000000000;
        case 2: return +0.77459666924148;
        }
        break;
    case 4:
        switch (index) {
        case 0: return -0.86113631159405;
        case 1: return -0.33998104358486;
        case 2: return +0.33998104358486;
        case 3: return +0.86113631159405;
        }
        break;
    case 5:
        switch (index) {
        case 0: return -0.90617984593866;
        case 1: return -0.53846931010568;
        case 2: return +0.00000000000000;
        case 3: return +0.53846931010568;
        case 4: return +0.90617984593866;
        }
        break;
    case 6:
        switch (index) {
        case 0: return -0.93246951420315;
        case 1: return -0.66120938646626;
        case 2: return -0.23861918608320;
        case 3: return +0.23861918608320;
        case 4: return +0.66120938646626;
        case 5: return +0.93246951420315;
        }
        break;
    case 7:
        switch (index) {
        case 0: return -0.94910791234276;
        case 1: return -0.74153118559939;
        case 2: return -0.40584515137740;
        case 3: return +0.00000000000000;
        case 4: return +0.40584515137740;
        case 5: return +0.74153118559939;
        case 6: return +0.94910791234276;
        }
        break;
    case 8:
        switch (index) {
        case 0: return -0.96028985649754;
        case 1: return -0.79666647741363;
        case 2: return -0.52553240991633;
        case 3: return -0.18343464249565;
        case 4: return +0.18343464249565;
        case 5: return +0.52553240991633;
        case 6: return +0.79666647741363;
        case 7: return +0.96028985649754;
        }
        break;
    case 9:
        switch (index) {
        case 0: return -0.96816023950763;
        case 1: return -0.83603110732664;
        case 2: return -0.61337143270059;
        case 3: return -0.32425342340381;
        case 4: return +0.00000000000000;
        case 5: return +0.32425342340381;
        case 6: return +0.61337143270059;
        case 7: return +0.83603110732664;
        case 8: return +0.96816023950763;
        }
        break;
    case 10:
        switch (index) {
        case 0: return -0.97390652851717;
        case 1: return -0.86506336668898;
        case 2: return -0.67940956829902;
        case 3: return -0.43339539412925;
        case 4: return -0.14887433898163;
        case 5: return +0.14887433898163;
        case 6: return +0.43339539412925;
        case 7: return +0.67940956829902;
        case 8: return +0.86506336668898;
        case 9: return +0.97390652851717;
        }
        break;
    case 11:
        switch (index) {
        case 0: return -0.97822865814606;
        case 1: return -0.88706259976810;
        case 2: return -0.73015200557405;
        case 3: return -0.51909612920681;
        case 4: return -0.26954315595234;
        case 5: return +0.00000000000000;
        case 6: return +0.26954315595234;
        case 7: return +0.51909612920681;
        case 8: return +0.73015200557405;
        case 9: return +0.88706259976810;
        case 10: return +0.97822865814606;
        }
        break;
    }
    assert(0);
}

double gauss_quadrature_weight(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0: return +2.00000000000000;
        }
        break;
    case 2:
        switch (index) {
        case 0: return +1.00000000000000;
        case 1: return +1.00000000000000;
        }
        break;
    case 3:
        switch (index) {
        case 0: return +0.55555555555556;
        case 1: return +0.88888888888889;
        case 2: return +0.55555555555556;
        }
        break;
    case 4:
        switch (index) {
        case 0: return +0.34785484513745;
        case 1: return +0.65214515486255;
        case 2: return +0.65214515486255;
        case 3: return +0.34785484513745;
        }
        break;
    case 5:
        switch (index) {
        case 0: return +0.23692688505619;
        case 1: return +0.47862867049937;
        case 2: return +0.56888888888889;
        case 3: return +0.47862867049937;
        case 4: return +0.23692688505619;
        }
        break;
    case 6:
        switch (index) {
        case 0: return +0.17132449237917;
        case 1: return +0.36076157304814;
        case 2: return +0.46791393457269;
        case 3: return +0.46791393457269;
        case 4: return +0.36076157304814;
        case 5: return +0.17132449237917;
        }
        break;
    case 7:
        switch (index) {
        case 0: return +0.12948496616887;
        case 1: return +0.27970539148928;
        case 2: return +0.38183005050512;
        case 3: return +0.41795918367347;
        case 4: return +0.38183005050512;
        case 5: return +0.27970539148928;
        case 6: return +0.12948496616887;
        }
        break;
    case 8:
        switch (index) {
        case 0: return +0.10122853629038;
        case 1: return +0.22238103445337;
        case 2: return +0.31370664587789;
        case 3: return +0.36268378337836;
        case 4: return +0.36268378337836;
        case 5: return +0.31370664587789;
        case 6: return +0.22238103445337;
        case 7: return +0.10122853629038;
        }
        break;
    case 9:
        switch (index) {
        case 0: return +0.08127438836157;
        case 1: return +0.18064816069486;
        case 2: return +0.26061069640294;
        case 3: return +0.31234707704000;
        case 4: return +0.33023935500126;
        case 5: return +0.31234707704000;
        case 6: return +0.26061069640294;
        case 7: return +0.18064816069486;
        case 8: return +0.08127438836157;
        }
        break;
    case 10:
        switch (index) {
        case 0: return +0.06667134430869;
        case 1: return +0.14945134915058;
        case 2: return +0.21908636251598;
        case 3: return +0.26926671931000;
        case 4: return +0.29552422471475;
        case 5: return +0.29552422471475;
        case 6: return +0.26926671931000;
        case 7: return +0.21908636251598;
        case 8: return +0.14945134915058;
        case 9: return +0.06667134430869;
        }
        break;
    case 11:
        switch (index) {
        case 0: return +0.05566856711617;
        case 1: return +0.12558036946490;
        case 2: return +0.18629021092773;
        case 3: return +0.23319376459199;
        case 4: return +0.26280454451025;
        case 5: return +0.27292508677790;
        case 6: return +0.26280454451025;
        case 7: return +0.23319376459199;
        case 8: return +0.18629021092773;
        case 9: return +0.12558036946490;
        case 10: return +0.05566856711617;
        }
        break;
    }
    assert(0);
}

double rk_parameter_for_step(int s, int rk_order)
{
    switch (rk_order) {
    case 1:
        switch (s) {
        case 0: return 0.0;
        }
        break;
    case 2:
        switch (s) {
        case 0: return 0.0;
        case 1: return 0.5;
        }
        break;
    case 3:
        switch (s) {
        case 0: return 0.0;
        case 1: return 3.0 / 4.0;
        case 2: return 1.0 / 3.0;
        }
        break;
    }
    assert(0);
}

void hydro_prim_to_cons(double* prim, double* cons)
{
    double rho = prim[0];
    double vel = prim[1];
    double pre = prim[2];
    cons[0] = rho;
    cons[1] = rho * vel;
    cons[2] = rho * vel * vel * 0.5 + pre / (adiabatic_gamma - 1.0);
}

void hydro_cons_to_prim(double* cons, double* prim)
{
    double d = cons[0];
    double s = cons[1];
    double e = cons[2];
    double pre = (e - 0.5 * s * s / d) * (adiabatic_gamma - 1.0);

    if (d <= 0.0) {
        fprintf(stderr, "[error] negative density\n");
        exit(1);
    }
    if (pre <= 0.0) {
        fprintf(stderr, "[error] negative pressure\n");
        exit(1);
    }

    prim[0] = d;
    prim[1] = s / d;
    prim[2] = pre;
}

void hydro_flux(double* prim, double* flux)
{
    double rho = prim[0];
    double vel = prim[1];
    double pre = prim[2];
    double e = rho * vel * vel * 0.5 + pre / (adiabatic_gamma - 1.0);
    flux[0] = rho * vel;
    flux[1] = rho * vel * vel + pre;
    flux[2] = (e + pre) * vel;
}

void hydro_flux2(double* prim, double* cons, double* flux)
{
    double vel = prim[1];
    double pre = prim[2];
    flux[0] = cons[0] * vel;
    flux[1] = cons[1] * vel + pre;
    flux[2] = (cons[2] + pre) * vel;
}

double hydro_sound_speed(double* prim)
{
    double rho = prim[0];
    double pre = prim[2];
    return sqrt(adiabatic_gamma * pre / rho);
}

void hydro_riemann_hlle(
    double* pl, double* pr, double* ul, double* ur, double* fhat)
{
    double fl[NUM_FIELDS];
    double fr[NUM_FIELDS];
    double vel_l = pl[1];
    double vel_r = pr[1];
    double cs_l = hydro_sound_speed(pl);
    double cs_r = hydro_sound_speed(pr);
    double lam_pl = vel_l + cs_l;
    double lam_ml = vel_l - cs_l;
    double lam_pr = vel_r + cs_r;
    double lam_mr = vel_r - cs_r;
    double ap = max3(0.0, lam_pl, lam_pr);
    double am = min3(0.0, lam_ml, lam_mr);

    hydro_flux2(pr, ur, fr);
    hydro_flux2(pl, ul, fl);

    for (int q = 0; q < NUM_FIELDS; ++q) {
        fhat[q]
            = (ap * fl[q] - am * fr[q] + ap * am * (ur[q] - ul[q])) / (ap - am);
    }
}



#define A_GRID 0
#define A_WGTS 1
#define A_WGTS_CACHE 2
#define A_WGTS_DELTA 3
#define A_PRIM 4
#define A_CONS 5
#define A_CONS_DELTA 6
#define A_CONS_SRF 7
#define A_FLUX 8
#define A_FLUX_GOD 9
#define A_TRZN 10
#define ARRAY_COUNT 12

#define CURRENT 0
#define INVALID 1

static double* global_array[ARRAY_COUNT];
static int global_array_status[ARRAY_COUNT];

const char* array_name(int array)
{
    switch (array) {
    case A_GRID: return "grid";
    case A_WGTS: return "wgts";
    case A_WGTS_DELTA: return "wgts_delta";
    case A_PRIM: return "prim";
    case A_CONS: return "cons";
    case A_CONS_DELTA: return "cons_delta";
    case A_CONS_SRF: return "cons_srf";
    case A_FLUX: return "flux";
    case A_FLUX_GOD: return "flux_god";
    case A_TRZN: return "trzn";
    }
    assert(0);
}

void array_shape(int array, int* shape)
{
    int ni = num_zones;
    int nq = NUM_FIELDS;
    int nr = order;
    int nl = order;

    switch (array) {
    case A_GRID:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = 1;
        return;
    case A_WGTS:
    case A_WGTS_CACHE:
    case A_WGTS_DELTA:
        shape[0] = ni;
        shape[1] = nq;
        shape[2] = nl;
        return;
    case A_PRIM:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case A_CONS:
    case A_CONS_DELTA:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case A_CONS_SRF:
        shape[0] = ni;
        shape[1] = 2;
        shape[2] = nq;
        return;
    case A_FLUX:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case A_FLUX_GOD:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case A_TRZN:
        shape[0] = ni;
        shape[1] = 1;
        shape[2] = 1;
        return;
    }
    assert(0);
}

void array_stride(int array, int* stride)
{
    int shape[3];
    array_shape(array, shape);
    stride[2] = 1;
    stride[1] = stride[2] * shape[2];
    stride[0] = stride[1] * shape[1];
}

double* array_ptr_stride(int array, int* stride)
{
    array_stride(array, stride);
    return global_array[array];
}

size_t array_len(int array)
{
    int shape[3];
    array_shape(array, shape);
    return shape[0] * shape[1] * shape[2];
}

void array_alloc_if_needed(int array)
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
}

int array_clear(int array)
{
    free(global_array[array]);
    global_array[array] = NULL;
    global_array_status[array] = INVALID;
    return 0;
}

int array_require_current(int array)
{
    if (global_array_status[array] != CURRENT) {
        fprintf(stderr, "[error] %s is not current\n", array_name(array));
        return 1;
    }
    return 0;
}

int array_set_current(int array)
{
    global_array_status[array] = CURRENT;
    return 0;
}

int array_invalidate(int array)
{
    global_array_status[array] = INVALID;
    return 0;
}

int array_print(int array)
{
    if (array_require_current(array)) {
        return 1;
    }

    double* data = global_array[array];
    int n[3];
    int s[3];
    array_shape(array, n);
    array_stride(array, s);

    for (int i = 0; i < n[0]; ++i) {
        for (int r = 0; r < n[1]; ++r) {
            for (int q = 0; q < n[2]; ++q) {
                fprintf(
                    terminal, "%+.8f ", data[i * s[0] + r * s[1] + q * s[2]]);
            }
            fprintf(terminal, "\n");
        }
    }
    return 0;
}

int project(int w_array, int u_array, int use_tci)
{
    if (array_require_current(u_array)) {
        return 1;
    }
    if (use_tci && array_require_current(A_TRZN)) {
        return 1;
    }
    array_alloc_if_needed(w_array);

    int num_poly = order;
    int num_points = order;
    int num_fields = NUM_FIELDS;
    int us[3];
    int ws[3];
    double* t = global_array[A_TRZN];
    double* u = array_ptr_stride(u_array, us);
    double* w = array_ptr_stride(w_array, ws);

    double phi[MAX_DG_ORDER * MAX_DG_ORDER];
    double xsi[MAX_DG_ORDER];
    double wgt[MAX_DG_ORDER];

    for (int r = 0; r < num_points; ++r) {
        xsi[r] = gauss_quadrature_node(order, r);
        wgt[r] = gauss_quadrature_weight(order, r);

        for (int l = 0; l < num_poly; ++l) {
            phi[r * order + l] = legendre_polynomial(l, xsi[r]);
        }
    }

    for (int i = 0; i < num_zones; ++i) {
        if (use_tci && t[i] <= tci_threshold) {
            continue;
        }
        for (int q = 0; q < num_fields; ++q) {
            for (int l = 0; l < num_poly; ++l) {
                double wiql = 0.0;

                for (int r = 0; r < num_points; ++r) {
                    double wr = wgt[r];
                    double prl = phi[r * order + l];
                    double* uirq = &u[i * us[0] + r * us[1] + q * us[2]];
                    wiql += *uirq * prl * wr * 0.5;
                }
                w[i * ws[0] + q * ws[1] + l * ws[2]] = wiql;
            }
        }
    }
    return array_set_current(w_array);
}

int grid_init()
{
    array_alloc_if_needed(A_GRID);

    double* grid = global_array[A_GRID];
    double dx = (domain_x1 - domain_x0) / (num_zones - 2);
    double x0 = domain_x0 - dx;

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {
            double xsi = gauss_quadrature_node(order, r);
            grid[i * order + r] = x0 + dx * (i + 0.5 * (1.0 + xsi));
        }
    }
    return array_set_current(A_GRID);
}

int prim_init(void (*prim_func)(double x, double*))
{
    if (array_require_current(A_GRID)) {
        return 1;
    }
    array_alloc_if_needed(A_PRIM);

    int num_points = order;
    int xs[3];
    int ps[3];
    double* x = array_ptr_stride(A_GRID, xs);
    double* p = array_ptr_stride(A_PRIM, ps);

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < num_points; ++r) {
            double* xir = &x[i * xs[0] + r * xs[1]];
            double* pir = &p[i * ps[0] + r * ps[1]];
            prim_func(*xir, pir);
        }
    }
    return array_set_current(A_PRIM);
}

void prim_func_sod(double x, double* prim)
{
    if (x < 0.5) {
        prim[0] = 1.0;
        prim[1] = 0.0;
        prim[2] = 1.0;
    } else {
        prim[0] = 0.125;
        prim[1] = 0.0;
        prim[2] = 0.1;
    }
}

int prim_init_sod()
{
    return prim_init(prim_func_sod);
}

void prim_func_dwave(double x, double* prim)
{
    prim[0] = 1.0 + 0.5 * sin(2 * M_PI * x);
    prim[1] = 1.0;
    prim[2] = 1.0;
}

int prim_init_dwave()
{
    return prim_init(prim_func_dwave);
}

// Riemann Tests from Table 4.1 of Toro (2009)
void prim_func_test1(double x, double* prim)
{
    // run until t=0.25
    if (x < 0.5) {
        prim[0] = 1.0;
        prim[1] = 0.0;
        prim[2] = 1.0;
    } else {
        prim[0] = 0.125;
        prim[1] = 0.0;
        prim[2] = 0.1;
    }
}

int prim_init_test1()
{
    return prim_init(prim_func_test1);
}

void prim_func_test2(double x, double* prim)
{
    // run until t=0.15
    if (x < 0.5) {
        prim[0] =  1.0;
        prim[1] = -2.0;
        prim[2] =  0.4;
    } else {
        prim[0] =  1.0;
        prim[1] =  2.0;
        prim[2] =  0.4;
    }
}

int prim_init_test2()
{
    return prim_init(prim_func_test2);
}

void prim_func_test3(double x, double* prim)
{
    // run until t=0.012
    if (x < 0.5) {
        prim[0] =  1.0;
        prim[1] =  0.0;
        prim[2] =  1000.0;
    } else {
        prim[0] =  1.0;
        prim[1] =  0.0;
        prim[2] =  0.01;
    }
}

int prim_init_test3()
{
    return prim_init(prim_func_test3);
}

void prim_func_test4(double x, double* prim)
{
    // run until t=0.035
    if (x < 0.5) {
        prim[0] =  1.0;
        prim[1] =  0.0;
        prim[2] =  0.01;
    } else {
        prim[0] =  1.0;
        prim[1] =  0.0;
        prim[2] =  100.0;
    }
}

int prim_init_test4()
{
    return prim_init(prim_func_test4);
}

void prim_func_test5(double x, double* prim)
{
    // run until t=0.035
    if (x < 0.5) {
        prim[0] =  5.99924;
        prim[1] =  19.5975;
        prim[2] =  460.894;
    } else {
        prim[0] =  5.99242;
        prim[1] = -6.19633;
        prim[2] =  46.0950;
    }
}

int prim_init_test5()
{
    return prim_init(prim_func_test5);
}

void prim_func_test_noh(double x, double* prim)
{
    // run until t=1.0 USE Gamma = 5/3
    if (x < 0.5) {
        prim[0] = 1.0;
        prim[1] = 1.0;
        prim[2] = 1.0e-6;
    } else {
        prim[0] = 1.0;
        prim[1] = -1.0;
        prim[2] = 1.0e-6;
    }
}

int prim_init_test_noh()
{
    return prim_init(prim_func_test_noh);
}

int prim_from_cons()
{
    if (array_require_current(A_CONS)) {
        return 1;
    }
    array_alloc_if_needed(A_PRIM);

    int num_points = order;
    int us[3];
    int ps[3];
    double* u = array_ptr_stride(A_CONS, us);
    double* p = array_ptr_stride(A_PRIM, ps);

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < num_points; ++r) {
            double* uirq = &u[i * us[0] + r * us[1]];
            double* pirq = &p[i * ps[0] + r * ps[1]];
            hydro_cons_to_prim(uirq, pirq);
        }
    }
    return array_set_current(A_PRIM);
}

int cons_from_prim()
{
    if (array_require_current(A_PRIM)) {
        return 1;
    }
    array_alloc_if_needed(A_CONS);

    int num_points = order;
    int ps[3];
    int us[3];
    double* p = array_ptr_stride(A_PRIM, ps);
    double* u = array_ptr_stride(A_CONS, us);

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < num_points; ++r) {
            double* pirq = &p[i * ps[0] + r * ps[1]];
            double* uirq = &u[i * us[0] + r * us[1]];
            hydro_prim_to_cons(pirq, uirq);
        }
    }
    return array_set_current(A_CONS);
}

int cons_from_wgts()
{
    if (array_require_current(A_WGTS)) {
        return 1;
    }
    array_alloc_if_needed(A_CONS);

    int num_poly = order;
    int num_points = order;
    int num_fields = NUM_FIELDS;
    int us[3];
    int ws[3];
    double* u = array_ptr_stride(A_CONS, us);
    double* w = array_ptr_stride(A_WGTS, ws);

    double phi[MAX_DG_ORDER * MAX_DG_ORDER];

    for (int r = 0; r < num_points; ++r) {
        for (int l = 0; l < num_poly; ++l) {
            phi[r * order + l]
                = legendre_polynomial(l, gauss_quadrature_node(order, r));
        }
    }

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < num_points; ++r) {
            for (int q = 0; q < num_fields; ++q) {
                double uirq = 0.0;

                for (int l = 0; l < num_poly; ++l) {
                    double prl = phi[r * order + l];
                    double* wiql = &w[i * ws[0] + q * ws[1] + l * ws[2]];
                    uirq += *wiql * prl;
                }
                u[i * us[0] + r * us[1] + q * us[2]] = uirq;
            }
        }
    }
    return array_set_current(A_CONS);
}

int cons_srf_from_wgts()
{
    if (array_require_current(A_WGTS)) {
        return 1;
    }
    array_alloc_if_needed(A_CONS_SRF);

    int num_poly = order;
    int num_fields = NUM_FIELDS;
    int us[3];
    int ws[3];
    double* u = array_ptr_stride(A_CONS_SRF, us);
    double* w = array_ptr_stride(A_WGTS, ws);

    double phi_srf[2 * MAX_DG_ORDER];

    for (int l = 0; l < num_poly; ++l) {
        phi_srf[0 * order + l] = legendre_polynomial(l, -1.0);
        phi_srf[1 * order + l] = legendre_polynomial(l, +1.0);
    }

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < 2; ++r) {
            for (int q = 0; q < num_fields; ++q) {
                double uirq = 0.0;

                for (int l = 0; l < num_poly; ++l) {
                    double prl = phi_srf[r * order + l];
                    double* wiql = &w[i * ws[0] + q * ws[1] + l * ws[2]];
                    uirq += *wiql * prl;
                }
                u[i * us[0] + r * us[1] + q * us[2]] = uirq;
            }
        }
    }
    return array_set_current(A_CONS_SRF);
}

int cons_delta_from_flux_god()
{
    if (array_require_current(A_FLUX_GOD) || array_require_current(A_GRID)) {
        return 1;
    }
    array_alloc_if_needed(A_CONS_DELTA);

    int num_points = order;
    int num_fields = NUM_FIELDS;
    double dt = time_step;

    double* f = global_array[A_FLUX_GOD];
    double* x = global_array[A_GRID];
    double* du = global_array[A_CONS_DELTA];

    for (int r = 1; r < num_zones * num_points - 1; ++r) {
        double* fimh = &f[(r + 0) * num_fields];
        double* fiph = &f[(r + 1) * num_fields];
        double ximh = 0.5 * (x[r - 1] + x[r + 0]);
        double xiph = 0.5 * (x[r + 0] + x[r + 1]);

        for (int q = 0; q < num_fields; ++q) {
            du[r * num_fields + q] = -(fiph[q] - fimh[q]) * dt / (xiph - ximh);
        }
    }
    return array_set_current(A_CONS_DELTA);
}

int cons_apply_bc()
{
    if (array_require_current(A_CONS)) {
        return 1;
    }
    int num_points = order;
    int num_fields = NUM_FIELDS;
    int us[3];
    double* u = array_ptr_stride(A_WGTS, us);

    for (int r = 0; r < num_points; ++r) {
        for (int q = 0; q < num_fields; ++q) {
            int a = r * us[1] + q * us[2];
            u[0 * us[0] + a] = u[(num_zones - 2) * us[0] + a];
            u[(num_zones - 1) * us[0] + a] = u[1 * us[0] + a];
        }
    }
    return 0;
}

int wgts_delta_from_dg()
{
    if (array_require_current(A_FLUX) || array_require_current(A_FLUX_GOD)
        || array_require_current(A_TRZN)) {
        return 1;
    }
    array_alloc_if_needed(A_WGTS_DELTA);

    int num_poly = order;
    int num_points = order;
    int num_fields = NUM_FIELDS;
    int gs[3];
    int fs[3];
    int ws[3];
    double* t = global_array[A_TRZN];
    double* g = array_ptr_stride(A_FLUX_GOD, gs);
    double* f = array_ptr_stride(A_FLUX, fs);
    double* dw = array_ptr_stride(A_WGTS_DELTA, ws);
    double dx = (domain_x1 - domain_x0) / (num_zones - 2);
    double dt = time_step;

    double phi_srf[2 * MAX_DG_ORDER];
    double dph[MAX_DG_ORDER * MAX_DG_ORDER];
    double xsi[MAX_DG_ORDER];
    double wgt[MAX_DG_ORDER];

    for (int l = 0; l < num_poly; ++l) {
        phi_srf[0 * order + l] = legendre_polynomial(l, -1.0);
        phi_srf[1 * order + l] = legendre_polynomial(l, +1.0);
    }

    for (int r = 0; r < num_points; ++r) {
        xsi[r] = gauss_quadrature_node(order, r);
        wgt[r] = gauss_quadrature_weight(order, r);

        for (int l = 0; l < num_poly; ++l) {
            dph[r * order + l] = legendre_polynomial_derivative(l, xsi[r]);
        }
    }

    for (int i = 1; i < num_zones - 1; ++i) {
        if (t[i] > tci_threshold) {
            continue;
        }
        for (int q = 0; q < num_fields; ++q) {
            double* fimh = &g[(i + 0) * gs[0]];
            double* fiph = &g[(i + 1) * gs[0]];

            for (int l = 0; l < num_poly; ++l) {
                double dwiql = 0.0;

                for (int r = 0; r < num_points; ++r) {
                    double firq = f[i * fs[0] + r * fs[1] + q * fs[2]];
                    double dprl = dph[r * order + l];
                    dwiql += firq * dprl * wgt[r];
                }
                dwiql += fimh[q] * phi_srf[0 * order + l];
                dwiql -= fiph[q] * phi_srf[1 * order + l];

                dw[i * ws[0] + q * ws[1] + l * ws[2]] = dwiql * dt / dx;
            }
        }
    }
    return array_set_current(A_WGTS_DELTA);
}

int wgts_cache_from_wgts()
{
    if (array_require_current(A_WGTS)) {
        return 1;
    }
    array_alloc_if_needed(A_WGTS_CACHE);

    double* w = global_array[A_WGTS];
    double* w0 = global_array[A_WGTS_CACHE];
    memcpy(w0, w, array_len(A_WGTS) * sizeof(double));

    return array_set_current(A_WGTS_CACHE);
}

int wgts_add_delta()
{
    if (array_require_current(A_WGTS) || array_require_current(A_WGTS_CACHE)
        || array_require_current(A_WGTS_DELTA)) {
        return 1;
    }

    size_t elem = array_len(A_WGTS);
    double* w = global_array[A_WGTS];
    double* w0 = global_array[A_WGTS_CACHE];
    double* dw = global_array[A_WGTS_DELTA];
    double rk = rk_parameter;

    for (size_t a = 0; a < elem; ++a) {
        w[a] = w0[a] * rk + (w[a] + dw[a]) * (1.0 - rk);
    }

    array_invalidate(A_CONS);
    array_invalidate(A_PRIM);
    array_invalidate(A_FLUX);
    array_invalidate(A_CONS_SRF);
    array_invalidate(A_FLUX_GOD);
    return 0;
}

int wgts_from_cons()
{
    return project(A_WGTS, A_CONS, 0);
}

int wgts_apply_bc()
{
    if (array_require_current(A_WGTS)) {
        return 1;
    }

    int num_poly = order;
    int num_fields = NUM_FIELDS;
    int ws[3];
    double* w = array_ptr_stride(A_WGTS, ws);

    for (int q = 0; q < num_fields; ++q) {
        for (int l = 0; l < num_poly; ++l) {
            int a = q * ws[1] + l * ws[2];
            w[0 * ws[0] + a] = w[(num_zones - 2) * ws[0] + a];
            w[(num_zones - 1) * ws[0] + a] = w[1 * ws[0] + a];
        }
    }
    return 0;
}

int wgts_delta_from_cons_delta()
{
    return project(A_WGTS_DELTA, A_CONS_DELTA, 1);
}

int flux_from_prim()
{
    if (array_require_current(A_PRIM)) {
        return 1;
    }
    array_alloc_if_needed(A_FLUX);

    int num_points = order;
    int ps[3];
    int fs[3];
    double* p = array_ptr_stride(A_PRIM, ps);
    double* f = array_ptr_stride(A_FLUX, fs);

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < num_points; ++r) {
            double* pir = &p[i * ps[0] + r * ps[1]];
            double* fir = &f[i * ps[0] + r * ps[1]];
            hydro_flux(pir, fir);
        }
    }
    return array_set_current(A_FLUX);
}

int flux_god_compute_fv()
{
    if (array_require_current(A_PRIM) || array_require_current(A_CONS)) {
        return 1;
    }
    array_alloc_if_needed(A_FLUX_GOD);

    int num_points = order;
    int num_fields = NUM_FIELDS;

    double* u = global_array[A_CONS];
    double* p = global_array[A_PRIM];
    double* f = global_array[A_FLUX_GOD];

    for (int r = 0; r < num_zones * num_points - 1; ++r) {
        double* pl = &p[(r + 0) * num_fields];
        double* pr = &p[(r + 1) * num_fields];
        double* ul = &u[(r + 0) * num_fields];
        double* ur = &u[(r + 1) * num_fields];
        double* fhat = &f[(r + 1) * num_fields];
        hydro_riemann_hlle(pl, pr, ul, ur, fhat);
    }
    return array_set_current(A_FLUX_GOD);
}

int flux_god_compute_dg()
{
    if (array_require_current(A_CONS_SRF) || array_require_current(A_TRZN)) {
        return 1;
    }
    array_alloc_if_needed(A_FLUX_GOD);

    int us[3];
    int fs[3];
    double* t = global_array[A_TRZN];
    double* u = array_ptr_stride(A_CONS_SRF, us);
    double* f = array_ptr_stride(A_FLUX_GOD, fs);
    double pl[NUM_FIELDS];
    double pr[NUM_FIELDS];

    for (int i = 0; i < num_zones - 1; ++i) {
        if (t[i] > tci_threshold || t[i + 1] > tci_threshold) {
            continue;
        }
        double* ul = &u[(i + 0) * us[0] + 1 * us[1]];
        double* ur = &u[(i + 1) * us[0] + 0 * us[1]];
        double* fhat = &f[(i + 1) * fs[0]];

        hydro_cons_to_prim(ul, pl);
        hydro_cons_to_prim(ur, pr);
        hydro_riemann_hlle(pl, pr, ul, ur, fhat);
    }
    return array_set_current(A_FLUX_GOD);
}

int trzn_compute()
{
    if (array_require_current(A_WGTS)) {
        return 1;
    }
    array_alloc_if_needed(A_TRZN);

    int k = order - 1;
    int indicator_field = 0; // density
    int ts[3];
    int ws[3];
    double* t = array_ptr_stride(A_TRZN, ts);
    double* w = array_ptr_stride(A_WGTS, ws);

    for (int i = 0; i < num_zones; ++i) {
        double w0 = w[i * ws[0] + indicator_field * ws[1] + 0 * ws[2]];
        double wk = w[i * ws[0] + indicator_field * ws[1] + k * ws[2]];
        t[i * ts[0]] = fabs(wk / w0);
    }
    return array_set_current(A_TRZN);
}

int stencil_print()
{
    printf("stencil data at order %d (xsi, wgt, phi, dph):\n", order);

    for (int n = 0; n < order; ++n) {
        printf("p%d\n", n);
        for (int r = 0; r < order; ++r) {
            double x = gauss_quadrature_node(order, r);
            double w = gauss_quadrature_weight(order, r);
            double y = legendre_polynomial(n, x);
            double z = legendre_polynomial_derivative(n, x);
            printf("%+.8f %+.8f %+.8f %+.8f\n", x, w, y, z);
        }
    }
    return 0;
}

int run()
{
    int iteration = 0;
    double dx = (domain_x1 - domain_x0) / (num_zones - 2);
    double wavespeed = 2.0;
    time_step = dx / wavespeed * cfl_parameter;

    TRY(grid_init());
    TRY(prim_init_test5());
    TRY(cons_from_prim());
    TRY(wgts_from_cons());

    while (time_phys < time_final) {
        struct timespec start = timer_start();

        // 1. compute troubled zone indicator
        // 2. compute low-order Godunov flux on internal and zone interfaces
        // 3. replace with high-order Godunov flux on zone interfaces, except
        //    on faces adjacent to a troubled zone
        // 4. compute low-order time derivative Lx of conserved fields in
        //    troubled zones
        // 5. compute weights L corresponding to Lx
        // 6. compute high-order time derivative M of weights
        // 7. add either L dt or M dt to weights, depending on whether it's a
        //    troubled zone

        TRY(trzn_compute());
        TRY(wgts_cache_from_wgts());

        for (int s = 0; s < rk_order; ++s) {
            rk_parameter = rk_parameter_for_step(s, rk_order);
            TRY(cons_srf_from_wgts());
            TRY(flux_god_compute_fv());
            TRY(flux_god_compute_dg());
            TRY(cons_delta_from_flux_god());
            TRY(wgts_delta_from_cons_delta());
            TRY(flux_from_prim());
            TRY(wgts_delta_from_dg());
            TRY(wgts_add_delta());
            TRY(cons_from_wgts());
            TRY(prim_from_cons());
        }
        array_invalidate(A_TRZN);
        array_invalidate(A_WGTS_CACHE);

        double seconds = timer_end(start) * 1e-9;
        time_phys += time_step;
        iteration += 1;
        printf("[%04d] t = %.3f Mzps=%.3f\n", iteration, time_phys,
            num_zones / seconds * 1e-6);
    }
    return 0;
}

int set_order(int new_order)
{
    if (new_order > MAX_DG_ORDER) {
        fprintf(stderr, "[error] maximum order is %d\n", MAX_DG_ORDER);
        return 1;
    }
    for (int i = 0; i < ARRAY_COUNT; ++i) {
        if (global_array[i] != NULL) {
            fprintf(stderr,
                "[error] cannot set order when %s is "
                "allocated\n",
                array_name(i));
            return 1;
        }
    }
    order = new_order;
    return 0;
}

int set_rk_order(int new_rk_order)
{
    if (new_rk_order < 1 || new_rk_order > 3) {
        fprintf(stderr, "[error] rk_order must be between 1 and 3, got %d\n",
            new_rk_order);
        return 1;
    }
    rk_order = new_rk_order;
    return 0;
}

int set_num_zones(int new_num_zones)
{
    for (int i = 0; i < ARRAY_COUNT; ++i) {
        if (global_array[i] != NULL) {
            fprintf(stderr,
                "[error] cannot set num_zones when %s is "
                "allocated\n",
                array_name(i));
            return 1;
        }
    }
    num_zones = new_num_zones + 2; // tack on guard zones here
    return 0;
}

int set_tci_threshold(double tci)
{
    tci_threshold = tci;
    return 0;
}

int set_cfl_parameter(double cfl)
{
    cfl_parameter = cfl;
    return 0;
}

int set_time_final(double tmax)
{
    time_final = tmax;
    return 0;
}

int set_adiabatic_gamma(double gamma)
{
    adiabatic_gamma = gamma;
    return 0;
}

int set_terminal(const char* terminal_str)
{
    FILE* new_terminal = NULL;

    if (terminal != stdout) {
        fclose(terminal);
    }
    if (strcmp(terminal_str, "stdout") == 0) {
        terminal = stdout;
    } else if ((new_terminal = fopen(terminal_str, "w")) == NULL) {
        fprintf(stderr, "[error] unable to open file %s\n", terminal_str);
        return 1;
    }
    terminal = new_terminal;
    return 0;
}

int load_commands_from_file(const char* filename)
{
    int load_command(const char* cmd);

    FILE* input = NULL;
    char buffer[MAX_COMMAND_LEN];

    if (strcmp(filename, "stdin") == 0) {
        input = stdin;
    } else if ((input = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "[error] unable to open file %s\n", filename);
        return 1;
    }

    while (fgets(buffer, MAX_COMMAND_LEN, input)) {
        buffer[strcspn(buffer, "\n")] = 0;
        if (strcmp("done", buffer) == 0) {
            return 0;
        }
        if (load_command(buffer) && input != stdin) {
            return 1;
        }
    }
    return 0;
}

int load_commands_from_array(int argc, const char** argv)
{
    int load_command(const char* cmd);

    for (int n = 1; n < argc; ++n) {
        if (load_command(argv[n])) {
            return 1;
        }
    }
    return 0;
}

int load_command(const char* cmd)
{
    if (cmd[0] == '#' || strlen(cmd) == 0)
        return 0;
    if (strcmp(cmd, "stencil:print") == 0)
        return stencil_print();

    if (strcmp(cmd, "grid:print") == 0)
        return array_print(A_GRID);
    if (strcmp(cmd, "grid:clear") == 0)
        return array_clear(A_GRID);
    if (strcmp(cmd, "grid:init") == 0)
        return grid_init();

    if (strcmp(cmd, "prim:print") == 0)
        return array_print(A_PRIM);
    if (strcmp(cmd, "prim:clear") == 0)
        return array_clear(A_PRIM);
    if (strcmp(cmd, "prim:init_sod") == 0)
        return prim_init_sod();
    if (strcmp(cmd, "prim:init_dwave") == 0)
        return prim_init_dwave();
    if (strcmp(cmd, "prim:from_cons") == 0)
        return prim_from_cons();

    if (strcmp(cmd, "cons:print") == 0)
        return array_print(A_CONS);
    if (strcmp(cmd, "cons:clear") == 0)
        return array_clear(A_CONS);
    if (strcmp(cmd, "cons:from_prim") == 0)
        return cons_from_prim();
    if (strcmp(cmd, "cons:from_wgts") == 0)
        return cons_from_wgts();

    if (strcmp(cmd, "cons_srf:print") == 0)
        return array_print(A_CONS_SRF);
    if (strcmp(cmd, "cons_srf:clear") == 0)
        return array_clear(A_CONS_SRF);
    if (strcmp(cmd, "cons_srf:from_wgts") == 0)
        return cons_srf_from_wgts();

    if (strcmp(cmd, "flux:print") == 0)
        return array_print(A_FLUX);
    if (strcmp(cmd, "flux:from_prim") == 0)
        return flux_from_prim();

    if (strcmp(cmd, "wgts:print") == 0)
        return array_print(A_WGTS);
    if (strcmp(cmd, "wgts:clear") == 0)
        return array_clear(A_WGTS);
    if (strcmp(cmd, "wgts:from_cons") == 0)
        return wgts_from_cons();
    if (strcmp(cmd, "wgts:add_delta") == 0)
        return wgts_add_delta();

    if (strcmp(cmd, "cons_delta:from_flux_god") == 0)
        return cons_delta_from_flux_god();
    if (strcmp(cmd, "wgts_delta:from_cons_delta") == 0)
        return wgts_delta_from_cons_delta();

    if (strcmp(cmd, "flux_god:print") == 0)
        return array_print(A_FLUX_GOD);
    if (strcmp(cmd, "flux_god:clear") == 0)
        return array_clear(A_FLUX_GOD);
    if (strcmp(cmd, "flux_god:compute_fv") == 0)
        return flux_god_compute_fv();
    if (strcmp(cmd, "flux_god:compute_dg") == 0)
        return flux_god_compute_dg();

    if (strcmp(cmd, "trzn:print") == 0)
        return array_print(A_TRZN);
    if (strcmp(cmd, "trzn:clear") == 0)
        return array_clear(A_TRZN);
    if (strcmp(cmd, "trzn:compute") == 0)
        return trzn_compute();

    if (strcmp(cmd, "run") == 0)
        return run();

    if (strncmp(cmd, "order=", 6) == 0)
        return set_order(atoi(cmd + 6));
    if (strncmp(cmd, "rk=", 3) == 0)
        return set_rk_order(atoi(cmd + 3));

    if (strncmp(cmd, "num_zones=", 10) == 0)
        return set_num_zones(atoi(cmd + 10));
    if (strncmp(cmd, "tci=", 4) == 0)
        return set_tci_threshold(atof(cmd + 4));
    if (strncmp(cmd, "cfl=", 4) == 0)
        return set_cfl_parameter(atof(cmd + 4));
    if (strncmp(cmd, "tmax=", 5) == 0)
        return set_time_final(atof(cmd + 5));
    if (strncmp(cmd, "gamma=", 6) == 0)
        return set_adiabatic_gamma(atof(cmd + 6));
    if (strncmp(cmd, "terminal=", 9) == 0)
        return set_terminal(cmd + 9);
    if (strncmp(cmd, "load:", 5) == 0)
        return load_commands_from_file(cmd + 5);
    if (strcmp(cmd, "done") == 0)
        return 1;
    if (strcmp(cmd, "quit") == 0)
        return 1;

    fprintf(stderr, "[error] unrecognized: %s\n", cmd);
    return 1;
}

int main(int argc, const char** argv)
{
    for (int i = 0; i < ARRAY_COUNT; ++i) {
        global_array[i] = NULL;
        global_array_status[i] = INVALID;
    }

    terminal = stdout;
    load_commands_from_array(argc, argv);

    for (int i = 0; i < ARRAY_COUNT; ++i) {
        array_clear(i);
    }
    return 0;
}
