#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_COMMAND_LEN 1024
#define MAX_DG_ORDER 11
#define NUM_FIELDS 3
#define ADIABATIC_GAMMA (5.0 / 3.0)
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))

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
    return p;
}

double legendre_polynomial_derivative(int n, double x)
{
    double p = 0.0;

    for (int k = 0; k <= n; ++k) {
        p += choose(n, k) * choose(n + k, k) * 0.5 * k
            * pow(0.5 * (x - 1), k - 1);
    }
    return p;
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

void hydro_prim_to_cons(double* prim, double* cons)
{
    double rho = prim[0];
    double vel = prim[1];
    double pre = prim[2];
    cons[0] = rho;
    cons[1] = rho * vel;
    cons[2] = rho * vel * vel * 0.5 + pre / (ADIABATIC_GAMMA - 1.0);
}

void hydro_cons_to_prim(double* cons, double* prim)
{
    double d = cons[0];
    double s = cons[1];
    double e = cons[2];
    double pre = (e - 0.5 * s * s / d) * (ADIABATIC_GAMMA - 1.0);

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
    double e = rho * vel * vel * 0.5 + pre / (ADIABATIC_GAMMA - 1.0);
    flux[0] = rho * vel;
    flux[1] = rho * vel * vel + pre;
    flux[2] = (e + pre) * vel;
}

double hydro_sound_speed(double* prim)
{
    double rho = prim[0];
    double pre = prim[2];
    return sqrt(ADIABATIC_GAMMA * pre / rho);
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

    hydro_flux(pr, fr);
    hydro_flux(pl, fl);

    for (int q = 0; q < NUM_FIELDS; ++q) {
        fhat[q]
            = (ap * fl[q] - am * fr[q] + ap * am * (ur[q] - ul[q])) / (ap - am);
    }
}

static FILE* terminal = NULL;
static int num_zones = 20;
static int order = 3;
static double domain_x0 = 0.0;
static double domain_x1 = 1.0;

static double time = 0.0;
static double timestep = 0.0;
static double* grid = NULL;
static double* weights = NULL;
static double* gauss_prim = NULL;
static double* gauss_cons = NULL;
static double* gauss_flux = NULL;
static double* surface_cons = NULL;
static double* godunov_flux = NULL;

#define ARRAY_GRID 31
#define ARRAY_WEIGHTS 32
#define ARRAY_GAUSS_PRIM 33
#define ARRAY_GAUSS_CONS 34
#define ARRAY_GAUSS_FLUX 35
#define ARRAY_SURFACE_CONS 36
#define ARRAY_RIEMANN_FLUX 37
#define ARRAY_GODUNOV_FLUX 38

double** array_pointer(int array)
{
    switch (array) {
    case ARRAY_GRID: return &grid;
    case ARRAY_WEIGHTS: return &weights;
    case ARRAY_GAUSS_PRIM: return &gauss_prim;
    case ARRAY_GAUSS_CONS: return &gauss_cons;
    case ARRAY_GAUSS_FLUX: return &gauss_flux;
    case ARRAY_SURFACE_CONS: return &surface_cons;
    case ARRAY_GODUNOV_FLUX: return &godunov_flux;
    }
    assert(0);
}

const char* array_name(int array)
{
    switch (array) {
    case ARRAY_GRID: return "grid";
    case ARRAY_WEIGHTS: return "weights";
    case ARRAY_GAUSS_PRIM: return "gauss_prim";
    case ARRAY_GAUSS_CONS: return "gauss_cons";
    case ARRAY_GAUSS_FLUX: return "gauss_flux";
    case ARRAY_SURFACE_CONS: return "surface_cons";
    case ARRAY_GODUNOV_FLUX: return "godunov_flux";
    }
    assert(0);
}

void array_shape(int array, int* shape)
{
    // grid: ni x nr x 1
    // weights: ni x nq x nl
    // gauss_prim: ni x nr x nq
    // gauss_cons: ni x nr x nq
    // gauss_flux: ni x nr x nq
    // surface_cons: ni x 2  x nq
    // godunov_flux: ni x nr x nq

    int ni = num_zones;
    int nq = NUM_FIELDS;
    int nr = order;
    int nl = order;

    switch (array) {
    case ARRAY_GRID:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = 1;
        return;
    case ARRAY_WEIGHTS:
        shape[0] = ni;
        shape[1] = nq;
        shape[2] = nl;
        return;
    case ARRAY_GAUSS_PRIM:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case ARRAY_GAUSS_CONS:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case ARRAY_GAUSS_FLUX:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
        return;
    case ARRAY_SURFACE_CONS:
        shape[0] = ni;
        shape[1] = 2;
        shape[2] = nq;
        return;
    case ARRAY_GODUNOV_FLUX:
        shape[0] = ni;
        shape[1] = nr;
        shape[2] = nq;
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

size_t array_len(int array)
{
    int shape[3];
    array_shape(array, shape);
    return shape[0] * shape[1] * shape[2];
}

void array_alloc_if_needed(int array)
{
    size_t elem = array_len(array);
    size_t size = elem * sizeof(double);
    double** ptr = array_pointer(array);

    if (*ptr == NULL) {
        *ptr = malloc(size);

        for (size_t a = 0; a < elem; ++a) {
            (*ptr)[elem] = 0.0;
        }
    }
}

int array_print(int array)
{
    int n[3];
    int s[3];

    array_shape(array, n);
    array_stride(array, s);
    double* data = *array_pointer(array);

    if (data == NULL) {
        fprintf(stderr, "[error] %s is not initialized\n", array_name(array));
        return 1;
    }
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

int grid_print()
{
    return array_print(ARRAY_GRID);
}

int grid_clear()
{
    free(grid);
    grid = NULL;
    return 0;
}

int grid_init()
{
    array_alloc_if_needed(ARRAY_GRID);

    double x0 = domain_x0;
    double x1 = domain_x1;
    double dx = (x1 - x0) / num_zones;

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {
            double xsi = gauss_quadrature_node(order, r);
            grid[i * order + r] = x0 + dx * (i + 0.5 * (1.0 + xsi));
        }
    }
    return 0;
}

int prim_print()
{
    return array_print(ARRAY_GAUSS_PRIM);
}

int prim_clear()
{
    free(gauss_prim);
    gauss_prim = NULL;
    return 0;
}

int prim_init(void (*prim_func)(double x, double*))
{
    if (grid == NULL) {
        fprintf(stderr, "[error] grid is not initialized\n");
        return 1;
    }
    int xn[3];
    int pn[3];
    int xs[3];
    int ps[3];

    array_shape(ARRAY_GRID, xn);
    array_shape(ARRAY_GAUSS_PRIM, pn);
    array_stride(ARRAY_GRID, xs);
    array_stride(ARRAY_GAUSS_PRIM, ps);
    array_alloc_if_needed(ARRAY_GAUSS_PRIM);

    double* x = grid;
    double* p = gauss_prim;

    for (int i = 0; i < pn[0]; ++i) {
        for (int r = 0; r < pn[1]; ++r) {
            double* xir = &x[i * xs[0] + r * xs[1]];
            double* pir = &p[i * ps[0] + r * ps[1]];
            prim_func(*xir, pir);
        }
    }
    return 0;
}

void prim_func_sod(double x, double* prim)
{
    if (x < 0.5) {
        prim[0] = 1.0;
        prim[1] = 0.0;
        prim[2] = 1.0;
    } else {
        prim[0] = 0.1;
        prim[1] = 0.0;
        prim[2] = 0.125;
    }
}

void prim_func_dwave(double x, double* prim)
{
    prim[0] = 1.0 + 0.5 * sin(2 * M_PI * x);
    prim[1] = 0.0;
    prim[2] = 1.0;
}

int prim_init_sod()
{
    return prim_init(prim_func_sod);
}

int prim_init_dwave()
{
    return prim_init(prim_func_dwave);
}

int prim_from_cons()
{
    if (gauss_cons == NULL) {
        fprintf(stderr, "[error] cons must be initialized\n");
        return 1;
    }
    int n[3];
    int s[3];

    array_shape(ARRAY_GAUSS_PRIM, n);
    array_stride(ARRAY_GAUSS_PRIM, s);
    array_alloc_if_needed(ARRAY_GAUSS_PRIM);

    for (int i = 0; i < n[0]; ++i) {
        for (int r = 0; r < n[1]; ++r) {
            double* uirq = &gauss_cons[i * s[0] + r * s[1]];
            double* pirq = &gauss_prim[i * s[0] + r * s[1]];
            hydro_cons_to_prim(uirq, pirq);
        }
    }
    return 0;
}

int cons_print()
{
    return array_print(ARRAY_GAUSS_CONS);
}

int cons_clear()
{
    free(gauss_cons);
    gauss_cons = NULL;
    return 0;
}

int cons_from_prim()
{
    if (gauss_prim == NULL) {
        fprintf(stderr, "[error] prim must be initialized\n");
        return 1;
    }
    int n[3];
    int s[3];

    array_shape(ARRAY_GAUSS_CONS, n);
    array_stride(ARRAY_GAUSS_CONS, s);
    array_alloc_if_needed(ARRAY_GAUSS_CONS);

    for (int i = 0; i < n[0]; ++i) {
        for (int r = 0; r < n[1]; ++r) {
            double* pirq = &gauss_prim[i * s[0] + r * s[1]];
            double* uirq = &gauss_cons[i * s[0] + r * s[1]];
            hydro_prim_to_cons(pirq, uirq);
        }
    }
    return 0;
}

int cons_from_wgts()
{
    if (weights == NULL) {
        fprintf(stderr, "[error] weights must be initialized\n");
        return 1;
    }
    int un[3];
    int wn[3];
    int us[3];
    int ws[3];

    array_shape(ARRAY_GAUSS_CONS, un);
    array_shape(ARRAY_WEIGHTS, wn);
    array_stride(ARRAY_GAUSS_CONS, us);
    array_stride(ARRAY_WEIGHTS, ws);
    array_alloc_if_needed(ARRAY_GAUSS_CONS);

    double* u = gauss_cons;
    double* w = weights;

    for (int i = 0; i < un[0]; ++i) {
        for (int r = 0; r < un[1]; ++r) {
            for (int q = 0; q < un[2]; ++q) {
                double xr = gauss_quadrature_node(order, r);
                double wr = gauss_quadrature_weight(order, r);
                double uirq = 0.0;

                for (int l = 0; l < wn[2]; ++l) {
                    double plr = legendre_polynomial(l, xr);
                    double* wiql = &w[i * ws[0] + q * ws[1] + l * ws[2]];
                    uirq += *wiql * plr;
                }
                u[i * us[0] + r * us[1] + q * us[2]] = uirq;
            }
        }
    }
    return 0;
}

int cons_add_gflux()
{
    if (gauss_cons == NULL) {
        fprintf(stderr, "[error] cons must be initialized\n");
        return 1;
    }
    if (godunov_flux == NULL) {
        fprintf(stderr, "[error] gflux must be initialized\n");
        return 1;
    }
    if (grid == NULL) {
        fprintf(stderr, "[error] grid must be initialized\n");
        return 1;
    }

    int n[3];
    int s[3];
    int num_points = order;
    int num_fields = NUM_FIELDS;
    double dt = timestep;

    for (int r = 1; r < num_zones * num_points - 1; ++r) {
        double* fimh = &godunov_flux[(r + 0) * num_fields];
        double* fiph = &godunov_flux[(r + 1) * num_fields];

        double ximh = 0.5 * (grid[r - 1] + grid[r + 0]);
        double xiph = 0.5 * (grid[r + 0] + grid[r + 1]);

        for (int q = 0; q < num_fields; ++q) {
            gauss_cons[r * num_fields + q]
                -= (fiph[q] - fimh[q]) * dt / (xiph - ximh);
        }
    }
    return 0;
}

int wgts_print()
{
    return array_print(ARRAY_WEIGHTS);
}

int wgts_clear()
{
    free(weights);
    weights = NULL;
    return 0;
}

int wgts_from_cons()
{
    if (gauss_cons == NULL) {
        fprintf(stderr, "[error] cons must be initialized\n");
        return 1;
    }
    int un[3];
    int wn[3];
    int us[3];
    int ws[3];

    array_shape(ARRAY_GAUSS_CONS, un);
    array_shape(ARRAY_WEIGHTS, wn);
    array_stride(ARRAY_GAUSS_CONS, us);
    array_stride(ARRAY_WEIGHTS, ws);
    array_alloc_if_needed(ARRAY_WEIGHTS);

    double* u = gauss_cons;
    double* w = weights;

    for (int i = 0; i < un[0]; ++i) {
        for (int q = 0; q < un[2]; ++q) {
            for (int l = 0; l < wn[2]; ++l) {
                double wiql = 0.0;

                for (int r = 0; r < un[1]; ++r) {
                    double xr = gauss_quadrature_node(order, r);
                    double wr = gauss_quadrature_weight(order, r);
                    double plr = legendre_polynomial(l, xr);
                    double* uirq = &u[i * us[0] + r * us[1] + q * us[2]];
                    wiql += *uirq * plr * wr * 0.5;
                }
                w[i * ws[0] + q * ws[1] + l * ws[2]] = wiql;
            }
        }
    }
    return 0;
}

int gflux_print()
{
    return array_print(ARRAY_GODUNOV_FLUX);
}

int gflux_compute()
{
    if (gauss_prim == NULL) {
        fprintf(stderr, "[error] prim must be initialized\n");
        return 1;
    }
    if (gauss_cons == NULL) {
        fprintf(stderr, "[error] cons must be initialized\n");
        return 1;
    }

    int n[3];
    int s[3];
    int num_points = order;
    int num_fields = NUM_FIELDS;
    array_alloc_if_needed(ARRAY_GODUNOV_FLUX);

    for (int r = 0; r < num_zones * num_points - 1; ++r) {
        double* pl = &gauss_prim[(r + 0) * num_fields];
        double* pr = &gauss_prim[(r + 1) * num_fields];
        double* ul = &gauss_cons[(r + 0) * num_fields];
        double* ur = &gauss_cons[(r + 1) * num_fields];
        double* fhat = &godunov_flux[(r + 1) * num_fields];
        hydro_riemann_hlle(pl, pr, ul, ur, fhat);
    }
    return 0;
}

int stencil_print()
{
    printf("stencil data at order %d (xsi, wgt, phi, dph):\n", order);

    for (int n = 0; n < order; ++n) {
        printf("p%d\n", n);
        for (int r = 0; r < order; ++r) {
            double x = gauss_quadrature_node(order, r);
            double w = gauss_quadrature_weight(order, r);
            double y = legendre_polynomial(n, x) * sqrt(2 * n + 1);
            double z = legendre_polynomial_derivative(n, x) * sqrt(2 * n + 1);
            printf("%+.8f %+.8f %+.8f %+.8f\n", x, w, y, z);
        }
    }
    return 0;
}

int set_order(int new_order)
{
    if (new_order > MAX_DG_ORDER) {
        fprintf(stderr, "[error] maximum order is %d\n", MAX_DG_ORDER);
        return 1;
    }
    if (grid != NULL) {
        fprintf(stderr, "[error] cannot set order when grid is allocated\n");
        return 1;
    }
    order = new_order;
    return 0;
}

int set_num_zones(int new_num_zones)
{
    if (grid || gauss_prim || gauss_cons || weights) {
        fprintf(stderr,
            "[error] cannot set num_zones when grid|prim|cons|wgts are "
            "allocated\n");
        return 1;
    }
    num_zones = new_num_zones;
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
        return grid_print();
    if (strcmp(cmd, "grid:clear") == 0)
        return grid_clear();
    if (strcmp(cmd, "grid:init") == 0)
        return grid_init();
    if (strcmp(cmd, "prim:print") == 0)
        return prim_print();
    if (strcmp(cmd, "prim:clear") == 0)
        return prim_clear();
    if (strcmp(cmd, "prim:init_sod") == 0)
        return prim_init_sod();
    if (strcmp(cmd, "prim:init_dwave") == 0)
        return prim_init_dwave();
    if (strcmp(cmd, "prim:init_dwave") == 0)
        return prim_init_dwave();
    if (strcmp(cmd, "prim:from_cons") == 0)
        return prim_from_cons();
    if (strcmp(cmd, "cons:clear") == 0)
        return cons_clear();
    if (strcmp(cmd, "cons:from_prim") == 0)
        return cons_from_prim();
    if (strcmp(cmd, "cons:from_wgts") == 0)
        return cons_from_wgts();
    if (strcmp(cmd, "cons:add_gflux") == 0)
        return cons_add_gflux();
    if (strcmp(cmd, "wgts:print") == 0)
        return wgts_print();
    if (strcmp(cmd, "wgts:clear") == 0)
        return wgts_clear();
    if (strcmp(cmd, "wgts:from_cons") == 0)
        return wgts_from_cons();
    if (strcmp(cmd, "gflux:print") == 0)
        return gflux_print();
    if (strcmp(cmd, "gflux:compute") == 0)
        return gflux_compute();

    if (strcmp(cmd, "toddler:simulate") == 0)
    {
        int iter = 0;
        double x0 = domain_x0;
        double x1 = domain_x1;
        double dx = (x1 - x0) / num_zones;
        timestep = dx / 1.0 * 0.05;

        load_command("grid:init");
        load_command("prim:init_sod");
        load_command("cons:from_prim");

        while (time < 0.1) {
            load_command("gflux:compute");
            load_command("cons:add_gflux");
            load_command("prim:from_cons");
            time += timestep;
            iter += 1;
            printf("[%04d] t = %.3f\n", iter, time);
        }
        return 0;
    }

    if (strncmp(cmd, "order=", 6) == 0)
        return set_order(atoi(cmd + 6));
    if (strncmp(cmd, "num_zones=", 10) == 0)
        return set_num_zones(atoi(cmd + 10));
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
    terminal = stdout;
    load_commands_from_array(argc, argv);
    grid_clear();
    prim_clear();
    cons_clear();
    wgts_clear();
    return 0;
}
