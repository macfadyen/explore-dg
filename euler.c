#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_DG_ORDER 11
#define NUM_FIELDS 3
#define ADIABATIC_GAMMA (5.0 / 3.0)

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

double gauss_quadrature_weight(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0:
            return +2.00000000000000;
        }
    case 2:
        switch (index) {
        case 0:
            return +1.00000000000000;
        case 1:
            return +1.00000000000000;
        }
    case 3:
        switch (index) {
        case 0:
            return +0.55555555555556;
        case 1:
            return +0.88888888888889;
        case 2:
            return +0.55555555555556;
        }
    case 4:
        switch (index) {
        case 0:
            return +0.34785484513745;
        case 1:
            return +0.65214515486255;
        case 2:
            return +0.65214515486255;
        case 3:
            return +0.34785484513745;
        }
    case 5:
        switch (index) {
        case 0:
            return +0.23692688505619;
        case 1:
            return +0.47862867049937;
        case 2:
            return +0.56888888888889;
        case 3:
            return +0.47862867049937;
        case 4:
            return +0.23692688505619;
        }
    case 6:
        switch (index) {
        case 0:
            return +0.17132449237917;
        case 1:
            return +0.36076157304814;
        case 2:
            return +0.46791393457269;
        case 3:
            return +0.46791393457269;
        case 4:
            return +0.36076157304814;
        case 5:
            return +0.17132449237917;
        }
    case 7:
        switch (index) {
        case 0:
            return +0.12948496616887;
        case 1:
            return +0.27970539148928;
        case 2:
            return +0.38183005050512;
        case 3:
            return +0.41795918367347;
        case 4:
            return +0.38183005050512;
        case 5:
            return +0.27970539148928;
        case 6:
            return +0.12948496616887;
        }
    case 8:
        switch (index) {
        case 0:
            return +0.10122853629038;
        case 1:
            return +0.22238103445337;
        case 2:
            return +0.31370664587789;
        case 3:
            return +0.36268378337836;
        case 4:
            return +0.36268378337836;
        case 5:
            return +0.31370664587789;
        case 6:
            return +0.22238103445337;
        case 7:
            return +0.10122853629038;
        }
    case 9:
        switch (index) {
        case 0:
            return +0.08127438836157;
        case 1:
            return +0.18064816069486;
        case 2:
            return +0.26061069640294;
        case 3:
            return +0.31234707704000;
        case 4:
            return +0.33023935500126;
        case 5:
            return +0.31234707704000;
        case 6:
            return +0.26061069640294;
        case 7:
            return +0.18064816069486;
        case 8:
            return +0.08127438836157;
        }
    case 10:
        switch (index) {
        case 0:
            return +0.06667134430869;
        case 1:
            return +0.14945134915058;
        case 2:
            return +0.21908636251598;
        case 3:
            return +0.26926671931000;
        case 4:
            return +0.29552422471475;
        case 5:
            return +0.29552422471475;
        case 6:
            return +0.26926671931000;
        case 7:
            return +0.21908636251598;
        case 8:
            return +0.14945134915058;
        case 9:
            return +0.06667134430869;
        }
    case 11:
        switch (index) {
        case 0:
            return +0.05566856711617;
        case 1:
            return +0.12558036946490;
        case 2:
            return +0.18629021092773;
        case 3:
            return +0.23319376459199;
        case 4:
            return +0.26280454451025;
        case 5:
            return +0.27292508677790;
        case 6:
            return +0.26280454451025;
        case 7:
            return +0.23319376459199;
        case 8:
            return +0.18629021092773;
        case 9:
            return +0.12558036946490;
        case 10:
            return +0.05566856711617;
        }
    }
    return 0.0;
}

double gauss_quadrature_node(int order, int index)
{
    switch (order) {
    case 1:
        switch (index) {
        case 0:
            return +0.00000000000000;
        }
    case 2:
        switch (index) {
        case 0:
            return -0.57735026918963;
        case 1:
            return +0.57735026918963;
        }
    case 3:
        switch (index) {
        case 0:
            return -0.77459666924148;
        case 1:
            return +0.00000000000000;
        case 2:
            return +0.77459666924148;
        }
    case 4:
        switch (index) {
        case 0:
            return -0.86113631159405;
        case 1:
            return -0.33998104358486;
        case 2:
            return +0.33998104358486;
        case 3:
            return +0.86113631159405;
        }
    case 5:
        switch (index) {
        case 0:
            return -0.90617984593866;
        case 1:
            return -0.53846931010568;
        case 2:
            return +0.00000000000000;
        case 3:
            return +0.53846931010568;
        case 4:
            return +0.90617984593866;
        }
    case 6:
        switch (index) {
        case 0:
            return -0.93246951420315;
        case 1:
            return -0.66120938646626;
        case 2:
            return -0.23861918608320;
        case 3:
            return +0.23861918608320;
        case 4:
            return +0.66120938646626;
        case 5:
            return +0.93246951420315;
        }
    case 7:
        switch (index) {
        case 0:
            return -0.94910791234276;
        case 1:
            return -0.74153118559939;
        case 2:
            return -0.40584515137740;
        case 3:
            return +0.00000000000000;
        case 4:
            return +0.40584515137740;
        case 5:
            return +0.74153118559939;
        case 6:
            return +0.94910791234276;
        }
    case 8:
        switch (index) {
        case 0:
            return -0.96028985649754;
        case 1:
            return -0.79666647741363;
        case 2:
            return -0.52553240991633;
        case 3:
            return -0.18343464249565;
        case 4:
            return +0.18343464249565;
        case 5:
            return +0.52553240991633;
        case 6:
            return +0.79666647741363;
        case 7:
            return +0.96028985649754;
        }
    case 9:
        switch (index) {
        case 0:
            return -0.96816023950763;
        case 1:
            return -0.83603110732664;
        case 2:
            return -0.61337143270059;
        case 3:
            return -0.32425342340381;
        case 4:
            return +0.00000000000000;
        case 5:
            return +0.32425342340381;
        case 6:
            return +0.61337143270059;
        case 7:
            return +0.83603110732664;
        case 8:
            return +0.96816023950763;
        }
    case 10:
        switch (index) {
        case 0:
            return -0.97390652851717;
        case 1:
            return -0.86506336668898;
        case 2:
            return -0.67940956829902;
        case 3:
            return -0.43339539412925;
        case 4:
            return -0.14887433898163;
        case 5:
            return +0.14887433898163;
        case 6:
            return +0.43339539412925;
        case 7:
            return +0.67940956829902;
        case 8:
            return +0.86506336668898;
        case 9:
            return +0.97390652851717;
        }
    case 11:
        switch (index) {
        case 0:
            return -0.97822865814606;
        case 1:
            return -0.88706259976810;
        case 2:
            return -0.73015200557405;
        case 3:
            return -0.51909612920681;
        case 4:
            return -0.26954315595234;
        case 5:
            return +0.00000000000000;
        case 6:
            return +0.26954315595234;
        case 7:
            return +0.51909612920681;
        case 8:
            return +0.73015200557405;
        case 9:
            return +0.88706259976810;
        case 10:
            return +0.97822865814606;
        }
    }
    return 0.0;
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

static int num_zones = 20;
static int order = 3;
static double domain_x0 = 0.0;
static double domain_x1 = 1.0;

static double* grid = NULL;
static double* solution_primitive = NULL;
static double* solution_conserved = NULL;
static double* solution_weights = NULL;

int grid_print()
{
    if (grid == NULL) {
        printf("[error] grid is not initialized\n");
        return 1;
    }
    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {
            printf("%.6f ", grid[i * order + r]);
        }
        printf("\n");
    }
    return 0;
}

int grid_clear()
{
    free(grid);
    grid = NULL;
    return 0;
}

int grid_init()
{
    if (grid != NULL) {
        printf("[error] grid is already initialized\n");
        return 1;
    }

    double x0 = domain_x0;
    double x1 = domain_x1;
    double dx = (x1 - x0) / num_zones;
    grid = malloc(num_zones * order * sizeof(double));

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
    int num_points = order;
    int num_fields = 3;
    int sr = 1;
    int sq = sr * num_points;
    int si = sq * num_fields;

    double* p = solution_primitive;

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {
            for (int q = 0; q < num_fields; ++q) {
                printf("%+.8f ", p[i * si + r * sr + q * sq]);
            }
            printf("\n");
        }
    }
    return 0;
}

int prim_clear()
{
    free(solution_primitive);
    solution_primitive = NULL;
    return 0;
}

int prim_init_sod()
{
    if (grid == NULL) {
        printf("[error] grid is not initialized\n");
        return 1;
    }

    int num_points = order;
    int num_fields = NUM_FIELDS;
    int x_sr = 1;
    int x_si = x_sr * num_points;
    int p_sr = 1;
    int p_sq = p_sr * num_points;
    int p_si = p_sq * num_fields;

    if (solution_primitive == NULL) {
        solution_primitive
            = malloc(num_zones * num_fields * num_points * sizeof(double));
    }

    double* x = grid;
    double* p = solution_primitive;

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {
            double xir = x[i * x_si + r * x_sr];
            double prim[NUM_FIELDS];

            if (xir < 0.5) {
                prim[0] = 1.0;
                prim[1] = 0.0;
                prim[2] = 1.0;
            } else {
                prim[0] = 0.1;
                prim[1] = 0.0;
                prim[2] = 0.125;
            }
            for (int q = 0; q < num_fields; ++q) {
                p[i * p_si + r * p_sr + q * p_sq] = prim[q];
            }
        }
    }
    return 0;
}

int cons_clear()
{
    free(solution_conserved);
    solution_conserved = NULL;
    return 0;
}

int cons_from_prim()
{
    if (solution_primitive == NULL) {
        printf("[error] prim must be initialized\n");
        return 1;
    }

    int num_points = order;
    int num_fields = NUM_FIELDS;
    int sr = 1;
    int sq = sr * num_points;
    int si = sq * num_fields;

    if (solution_conserved == NULL) {
        solution_conserved
            = malloc(num_zones * num_fields * num_points * sizeof(double));
    }

    double* p = solution_primitive;
    double* u = solution_conserved;

    double prim[NUM_FIELDS];
    double cons[NUM_FIELDS];

    for (int i = 0; i < num_zones; ++i) {
        for (int r = 0; r < order; ++r) {

            for (int q = 0; q < num_fields; ++q) {
                prim[q] = p[i * si + r * sr + q * sq];
            }
            hydro_prim_to_cons(prim, cons);

            for (int q = 0; q < num_fields; ++q) {
                u[q] = p[i * si + r * sr + q * sq] = cons[q];
            }
        }
    }
    return 0;
}

int wgts_clear()
{
    free(solution_weights);
    solution_weights = NULL;
    return 0;
}

int wgts_from_cons()
{
    if (solution_conserved == NULL) {
        printf("[error] cons must be initialized\n");
        return 1;
    }

    int num_poly = order;
    int num_points = order;
    int num_fields = 3;

    int u_sr = 1;
    int u_sq = u_sr * num_points;
    int u_si = u_sq * num_fields;

    int w_sl = 1;
    int w_sq = w_sl * num_poly;
    int w_si = w_sq * num_fields;

    if (solution_weights == NULL) {
        solution_weights
            = malloc(num_zones * num_fields * num_points * sizeof(double));
    }
    for (int a = 0; a < num_zones * num_fields * num_points; ++a) {
        solution_weights[a] = 0.0;
    }

    double* u = solution_conserved;
    double* w = solution_weights;

    for (int i = 0; i < num_zones; ++i) {
        for (int q = 0; q < num_fields; ++q) {
            for (int l = 0; l < num_poly; ++l) {
                for (int r = 0; r < num_points; ++r) {
                    double uiqr = u[i * u_si + q * u_sq + r * u_sr];
                    double xr = gauss_quadrature_node(order, r);
                    double wr = gauss_quadrature_weight(order, r);
                    double plr = legendre_polynomial(l, xr);
                    w[i * w_si + q * w_sq + l * w_sl] += uiqr * plr * wr;
                }
            }
        }
    }
    return 0;
}

int stencil_print()
{
    printf("stencil data at order %d (xsi, wgt, phi):\n", order);

    for (int n = 0; n < order; ++n) {
        printf("p%d\n", n);
        for (int r = 0; r < order; ++r) {
            double x = gauss_quadrature_node(order, r);
            double w = gauss_quadrature_weight(order, r);
            double y = legendre_polynomial(n, x) * sqrt(2 * n + 1);
            printf("%+.8f %+.8f %+.8f\n", x, w, y);
        }
    }
    return 0;
}

int set_order(int new_order)
{
    if (new_order > MAX_DG_ORDER) {
        printf("[error] maximum order is %d\n", MAX_DG_ORDER);
        return 1;
    }
    if (grid != NULL) {
        printf("[error] cannot set order when grid is allocated\n");
        return 1;
    }
    order = new_order;
    printf("[set] order=%d\n", order);
    return 0;
}

int set_num_zones(int new_num_zones)
{
    if (grid != NULL) {
        printf("[error] cannot set num_zones when grid is allocated\n");
        return 1;
    }
    num_zones = new_num_zones;
    printf("[set] num_zones=%d\n", num_zones);
    return 0;
}

int load_command(const char* cmd)
{
    if (strcmp(cmd, "stencil:print") == 0) {
        return stencil_print();
    }
    if (strcmp(cmd, "grid:print") == 0) {
        return grid_print();
    }
    if (strcmp(cmd, "grid:clear") == 0) {
        return grid_clear();
    }
    if (strcmp(cmd, "grid:init") == 0) {
        return grid_init();
    }
    if (strcmp(cmd, "prim:print") == 0) {
        return prim_print();
    }
    if (strcmp(cmd, "prim:clear") == 0) {
        return prim_clear();
    }
    if (strcmp(cmd, "prim:init_sod") == 0) {
        return prim_init_sod();
    }
    if (strcmp(cmd, "cons:clear") == 0) {
        return cons_clear();
    }
    if (strcmp(cmd, "cons:from_prim") == 0) {
        return cons_from_prim();
    }
    if (strcmp(cmd, "wgts:clear") == 0) {
        return wgts_clear();
    }
    if (strcmp(cmd, "wgts:from_cons") == 0) {
        return wgts_from_cons();
    }
    if (strncmp(cmd, "order=", 6) == 0) {
        return set_order(atoi(cmd + 6));
    }
    if (strncmp(cmd, "num_zones=", 10) == 0) {
        return set_num_zones(atoi(cmd + 10));
    }
    printf("[error] unrecognized: %s\n", cmd);
    return 1;
}

int main(int argc, const char** argv)
{
    for (int n = 1; n < argc; ++n) {
        if (load_command(argv[n])) {
            break;
        }
    }
    grid_clear();
    prim_clear();
    cons_clear();
    wgts_clear();
    return 0;
}
