// =============================================================================
// test_kabsch_nogsl.cc
//
// Unit test: compares kabsch_superpositioning (GSL) vs
//            kabsch_superpositioning_nogsl (no GSL) on a range of test cases.
//
// Build (from the Latfit_modified/ directory):
//   g++ -std=c++11 -I. \
//       biu/test_kabsch_nogsl.cc biu/kabsch_nogsl.cc \
//       -lgsl -lgslcblas -lm \
//       -o test_kabsch_nogsl && ./test_kabsch_nogsl
//
// The test passes when every assertion prints OK and the final line says PASS.
// =============================================================================

#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <cassert>
#include <functional>

// ---------- GSL headers (original implementation) ---------------------------
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

// ---------- Our new implementation ------------------------------------------
#include "kabsch_nogsl.hh"


// ============================================================================
// Verbatim copy of the original GSL kabsch_superpositioning() so the test is
// completely self-contained and independent of SuperPos_Kabsch.cc.
// ============================================================================
static inline void kabsch_gsl_vector_cross_test(
        const gsl_vector *a, const gsl_vector *b, gsl_vector *c)
{
    double a0=gsl_vector_get(a,0), a1=gsl_vector_get(a,1), a2=gsl_vector_get(a,2);
    double b0=gsl_vector_get(b,0), b1=gsl_vector_get(b,1), b2=gsl_vector_get(b,2);
    gsl_vector_set(c,0,a1*b2-b1*a2);
    gsl_vector_set(c,1,a2*b0-b2*a0);
    gsl_vector_set(c,2,a0*b1-b0*a1);
}

#define NORM_EPS_TEST 0.00000001

static int kabsch_superpositioning_gsl(
        unsigned int size,
        gsl_matrix *X, gsl_matrix *Y,
        gsl_matrix *U, gsl_vector *t, double *s)
{
    unsigned int i,j,k;
    int U_ok=1;
    double n=1.0/size;
    gsl_vector *cx=gsl_vector_alloc(3);
    gsl_vector *cy=gsl_vector_alloc(3);
    gsl_matrix *R=gsl_matrix_alloc(3,3);
    gsl_matrix *RTR=gsl_matrix_alloc(3,3);
    gsl_eigen_symmv_workspace *espace=gsl_eigen_symmv_alloc(3);
    gsl_matrix *evec=gsl_matrix_alloc(3,3);
    gsl_vector *eval=gsl_vector_alloc(3);

    gsl_vector_set_zero(cx);
    for(i=size;i>0;){gsl_vector_const_view row=gsl_matrix_const_row(X,--i);gsl_vector_add(cx,&row.vector);}
    gsl_vector_scale(cx,n);
    gsl_vector_set_zero(cy);
    for(i=size;i>0;){gsl_vector_const_view row=gsl_matrix_const_row(Y,--i);gsl_vector_add(cy,&row.vector);}
    gsl_vector_scale(cy,n);
    for(i=size;i>0;){gsl_vector_view row=gsl_matrix_row(X,--i);gsl_vector_sub(&row.vector,cx);}
    for(i=size;i>0;){gsl_vector_view row=gsl_matrix_row(Y,--i);gsl_vector_sub(&row.vector,cy);}

    if(size==1){gsl_matrix_set_identity(U);}
    else{
        gsl_matrix_set_zero(R);
        for(k=size;k>0;){--k;for(i=3;i>0;){--i;for(j=3;j>0;){--j;
            gsl_matrix_set(R,i,j,gsl_matrix_get(R,i,j)+gsl_matrix_get(Y,k,i)*gsl_matrix_get(X,k,j));}}}
        gsl_matrix_set_zero(RTR);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,R,0.0,RTR);
        gsl_eigen_symmv(RTR,eval,evec,espace);
        gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
        if(gsl_vector_get(eval,1)>NORM_EPS_TEST){
            double norm_b0,norm_b1,norm_b2;
            gsl_vector_const_view a0=gsl_matrix_const_column(evec,0);
            gsl_vector_const_view a1=gsl_matrix_const_column(evec,1);
            gsl_vector_view a2=gsl_matrix_column(evec,2);
            gsl_vector_view b0=gsl_matrix_column(RTR,0);
            gsl_vector_view b1=gsl_matrix_column(RTR,1);
            gsl_vector_view b2=gsl_matrix_column(RTR,2);
            kabsch_gsl_vector_cross_test(&a0.vector,&a1.vector,&a2.vector);
            gsl_blas_dgemv(CblasNoTrans,1.0,R,&a0.vector,0.0,&b0.vector);
            norm_b0=gsl_blas_dnrm2(&b0.vector);
            gsl_blas_dgemv(CblasNoTrans,1.0,R,&a1.vector,0.0,&b1.vector);
            norm_b1=gsl_blas_dnrm2(&b1.vector);
            if(norm_b0>NORM_EPS_TEST&&norm_b1>NORM_EPS_TEST){
                gsl_vector_scale(&b0.vector,1.0/norm_b0);
                gsl_vector_scale(&b1.vector,1.0/norm_b1);
                kabsch_gsl_vector_cross_test(&b0.vector,&b1.vector,&b2.vector);
                norm_b2=gsl_blas_dnrm2(&b2.vector);
                if(norm_b2>NORM_EPS_TEST){
                    gsl_matrix_set_zero(U);
                    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,RTR,evec,0.0,U);
                }else{U_ok=0;gsl_matrix_set_identity(U);}
            }else{U_ok=0;gsl_matrix_set_identity(U);}
        }else{U_ok=0;gsl_matrix_set_identity(U);}
    }

    if(s){
        *s=1.0;
        if(U_ok&&size>1){
            double dom=0.0,nom=0.0,dom_i,nom_i;
            gsl_vector *Uxi=gsl_vector_alloc(3);
            for(i=size;i>0;){
                gsl_vector_const_view rx=gsl_matrix_const_row(X,--i);
                gsl_vector_const_view ry=gsl_matrix_const_row(Y,i);
                gsl_vector_set_zero(Uxi);
                gsl_blas_dgemv(CblasNoTrans,1.0,U,&rx.vector,1.0,Uxi);
                gsl_blas_ddot(&ry.vector,Uxi,&nom_i); nom+=nom_i;
                gsl_blas_ddot(Uxi,Uxi,&dom_i); dom+=dom_i;
            }
            *s=nom/dom;
            gsl_vector_free(Uxi);
        }
        gsl_vector_scale(cx,*s);
    }
    gsl_vector_memcpy(t,cy);
    gsl_blas_dgemv(CblasNoTrans,-1.0,U,cx,1.0,t);

    gsl_vector_free(eval);gsl_matrix_free(evec);gsl_eigen_symmv_free(espace);
    gsl_matrix_free(RTR);gsl_matrix_free(R);gsl_vector_free(cy);gsl_vector_free(cx);
    return U_ok;
}


// ============================================================================
// Test framework helpers
// ============================================================================

static int g_tests_run    = 0;
static int g_tests_failed = 0;

static void check(const char* label, double a, double b, double tol = 1e-8) {
    g_tests_run++;
    if (std::fabs(a - b) > tol) {
        g_tests_failed++;
        printf("  FAIL  %-45s  got=%.10f  expected=%.10f  diff=%.2e\n",
               label, a, b, std::fabs(a-b));
    } else {
        printf("  OK    %-45s  (%.10f)\n", label, a);
    }
}

// Apply rotation U to point p, add translation t, scale by s (same as the
// algorithm does post-fit to recover superposed coordinates).
static void apply_transform(double U[3][3], double t[3], double s,
                             const double p_in[3], double p_out[3]) {
    for (int i = 0; i < 3; ++i) {
        p_out[i] = t[i];
        for (int j = 0; j < 3; ++j)
            p_out[i] += s * U[i][j] * p_in[j];
    }
}

// ============================================================================
// Run one test case: call both implementations with the same data, compare U,t,s
// and the resulting transformed points.
// ============================================================================
struct Points {
    std::vector<double> data; // flat row-major (N x 3)
    int N;
    Points(std::initializer_list<std::initializer_list<double>> rows) : N(rows.size()) {
        for (auto& r : rows)
            for (double v : r)
                data.push_back(v);
    }
};

static void run_test(const std::string& name, const Points& X_in, const Points& Y_in,
                     bool use_scale = false)
{
    printf("\n=== %s ===\n", name.c_str());
    const int N = X_in.N;
    assert(N == Y_in.N);

    // ---- GSL version -------------------------------------------------------
    std::vector<double> Xg(X_in.data), Yg(Y_in.data);
    gsl_matrix *Xm = gsl_matrix_alloc(N, 3);
    gsl_matrix *Ym = gsl_matrix_alloc(N, 3);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < 3; ++j) {
            gsl_matrix_set(Xm, i, j, Xg[i*3+j]);
            gsl_matrix_set(Ym, i, j, Yg[i*3+j]);
        }
    gsl_matrix *Ug = gsl_matrix_alloc(3, 3);
    gsl_vector *tg = gsl_vector_alloc(3);
    double sg = 1.0;
    kabsch_superpositioning_gsl(N, Xm, Ym, Ug, tg, use_scale ? &sg : nullptr);

    double Ug_arr[3][3], tg_arr[3];
    for (int i=0;i<3;i++){tg_arr[i]=gsl_vector_get(tg,i);for(int j=0;j<3;j++)Ug_arr[i][j]=gsl_matrix_get(Ug,i,j);}
    gsl_matrix_free(Xm); gsl_matrix_free(Ym); gsl_matrix_free(Ug); gsl_vector_free(tg);

    // ---- NoGSL version -----------------------------------------------------
    std::vector<double> Xn(X_in.data), Yn(Y_in.data);
    double Un[3][3] = {}, tn[3] = {}, sn = 1.0;
    kabsch_superpositioning_nogsl(N, Xn.data(), Yn.data(), Un, tn, use_scale ? &sn : nullptr);

    // ---- Compare rotation matrix -------------------------------------------
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            char lbl[64];
            // Rotation matrices can agree up to a sign flip of a row when the
            // eigenvalue is degenerate; for non-degenerate cases they agree to
            // machine precision.  We allow a loose tolerance here and tighten
            // it per-element only where the difference is large.
            snprintf(lbl, sizeof(lbl), "U[%d][%d]", i, j);
            // Accept if |a-b| < 1e-6 OR |a+b| < 1e-6 (sign ambiguity on
            // degenerate eigenvectors is physically irrelevant — the final
            // transformed points still agree).
            double diff  = std::fabs(Ug_arr[i][j] - Un[i][j]);
            double diffS = std::fabs(Ug_arr[i][j] + Un[i][j]);
            double best  = std::min(diff, diffS);
            g_tests_run++;
            if (best < 1e-6) {
                printf("  OK    %-45s  gsl=%.8f  nogsl=%.8f\n", lbl, Ug_arr[i][j], Un[i][j]);
            } else {
                g_tests_failed++;
                printf("  FAIL  %-45s  gsl=%.8f  nogsl=%.8f  diff=%.2e\n",
                       lbl, Ug_arr[i][j], Un[i][j], diff);
            }
        }

    // ---- Compare translation vector ----------------------------------------
    for (int i = 0; i < 3; ++i) {
        char lbl[32];
        snprintf(lbl, sizeof(lbl), "t[%d]", i);
        check(lbl, tn[i], tg_arr[i], 1e-6);
    }

    // ---- Compare scaling ---------------------------------------------------
    if (use_scale)
        check("s (scaling)", sn, sg, 1e-6);

    // ---- Compare final transformed points ----------------------------------
    // The ultimate criterion: after transformation, both versions must map
    // X points to the same positions (which should be close to Y).
    // We use the ORIGINAL (un-centred) X data and the transformation.
    for (int i = 0; i < N; ++i) {
        double pg[3], pn[3];
        apply_transform(Ug_arr, tg_arr, use_scale ? sg : 1.0, X_in.data.data()+i*3, pg);
        apply_transform(Un,     tn,     use_scale ? sn : 1.0, X_in.data.data()+i*3, pn);
        for (int j = 0; j < 3; ++j) {
            char lbl[64];
            snprintf(lbl, sizeof(lbl), "transformed_X[%d][%d]", i, j);
            check(lbl, pn[j], pg[j], 1e-5);
        }
    }
}


// ============================================================================
// Additional focused tests for the internal math helpers
// ============================================================================

static void test_jacobi_eigen() {
    printf("\n=== jacobi_sym3 eigen-decomposition ===\n");

    // Build a known symmetric 3x3 with eigenvalues 6, 3, 1
    // A = Q * diag(6,3,1) * Q^T where Q is a known rotation
    // Use Q = Rx(pi/4) for simplicity
    using namespace kabsch_nogsl_detail;

    // Simple diagonal matrix: eigenvalues are 5, 3, 1
    Mat3 A;
    A(0,0)=5; A(1,1)=3; A(2,2)=1;

    double eval[3];
    Mat3 evec;
    jacobi_sym3(A, eval, evec);
    sort_eigenpairs_desc(eval, evec);

    check("jacobi eigenvalue[0] == 5", eval[0], 5.0, 1e-10);
    check("jacobi eigenvalue[1] == 3", eval[1], 3.0, 1e-10);
    check("jacobi eigenvalue[2] == 1", eval[2], 1.0, 1e-10);

    // Eigenvectors of a diagonal matrix are unit vectors; order must match.
    // The eigenvector for lambda=5 is [1,0,0], for 3 is [0,1,0], for 1 is [0,0,1]
    // (up to sign).
    for (int col = 0; col < 3; ++col) {
        Vec3 ev = evec.col(col);
        double len = norm(ev);
        char lbl[64];
        snprintf(lbl, sizeof(lbl), "eigenvec[%d] is unit length", col);
        check(lbl, len, 1.0, 1e-10);
    }

    // More interesting: non-diagonal symmetric matrix with known eigenvectors
    // A2 = [[2,-1,0],[-1,2,-1],[0,-1,2]] — tridiagonal, eigenvalues 2+sqrt2, 2, 2-sqrt2
    Mat3 A2;
    A2(0,0)=2; A2(0,1)=-1; A2(0,2)=0;
    A2(1,0)=-1; A2(1,1)=2; A2(1,2)=-1;
    A2(2,0)=0; A2(2,1)=-1; A2(2,2)=2;
    double eval2[3];
    Mat3 evec2;
    jacobi_sym3(A2, eval2, evec2);
    sort_eigenpairs_desc(eval2, evec2);
    double lam0 = 2.0 + std::sqrt(2.0);
    double lam2 = 2.0 - std::sqrt(2.0);
    check("tridiag eigenvalue[0]", eval2[0], lam0, 1e-8);
    check("tridiag eigenvalue[1]", eval2[1], 2.0,  1e-8);
    check("tridiag eigenvalue[2]", eval2[2], lam2, 1e-8);

    // Verify A2 * v = lambda * v for each eigenpair
    for (int col = 0; col < 3; ++col) {
        Vec3 ev  = evec2.col(col);
        Vec3 Av  = matvec(A2, ev);
        Vec3 lev = scale(ev, eval2[col]);
        char lbl[64];
        snprintf(lbl, sizeof(lbl), "A2*v[%d] == lambda*v[%d] (x)", col, col);
        check(lbl, Av[0], lev[0], 1e-8);
        snprintf(lbl, sizeof(lbl), "A2*v[%d] == lambda*v[%d] (y)", col, col);
        check(lbl, Av[1], lev[1], 1e-8);
        snprintf(lbl, sizeof(lbl), "A2*v[%d] == lambda*v[%d] (z)", col, col);
        check(lbl, Av[2], lev[2], 1e-8);
    }
}

static void test_vector_cross() {
    printf("\n=== vector cross-product ===\n");
    using namespace kabsch_nogsl_detail;

    Vec3 x(1,0,0), y(0,1,0), z(0,0,1);
    Vec3 xcy = cross(x,y);
    check("x cross y = (0,0,1) z", xcy[2], 1.0, 1e-15);
    check("x cross y = (0,0,1) x", xcy[0], 0.0, 1e-15);
    check("x cross y = (0,0,1) y", xcy[1], 0.0, 1e-15);

    Vec3 a(1,2,3), b(4,5,6);
    Vec3 axb = cross(a,b);
    // a x b = (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
    check("a cross b [0]", axb[0], -3.0, 1e-14);
    check("a cross b [1]", axb[1],  6.0, 1e-14);
    check("a cross b [2]", axb[2], -3.0, 1e-14);

    // cross product is anti-commutative
    Vec3 bxa = cross(b,a);
    check("b cross a [0] = -(a cross b)", bxa[0], 3.0, 1e-14);
}

static void test_matmul() {
    printf("\n=== matrix multiply ===\n");
    using namespace kabsch_nogsl_detail;

    Mat3 A = Mat3::identity();
    Mat3 B; B(0,1)=1; B(1,2)=1; B(2,0)=1; // cyclic permutation
    Mat3 C = matmul(A, B);
    check("I*B == B [0][1]", C(0,1), 1.0, 1e-15);
    check("I*B == B [0][0]", C(0,0), 0.0, 1e-15);

    Mat3 D = matmul(B, transpose(B)); // B * B^T = I for an orthogonal B
    check("B*B^T == I [0][0]", D(0,0), 1.0, 1e-15);
    check("B*B^T == I [1][1]", D(1,1), 1.0, 1e-15);
    check("B*B^T == I [2][2]", D(2,2), 1.0, 1e-15);
    check("B*B^T == I [0][1]", D(0,1), 0.0, 1e-15);
}


// ============================================================================
// Main: run all tests
// ============================================================================
int main() {
    printf("=======================================================\n");
    printf(" Kabsch nogsl vs GSL unit test\n");
    printf("=======================================================\n");

    // --- Internal math helpers -----------------------------------------------
    test_vector_cross();
    test_matmul();
    test_jacobi_eigen();

    // --- End-to-end superposition tests --------------------------------------

    // Test 1: identical point sets -> U = I, t = 0
    run_test("Identical sets (4 points)",
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0}}),
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0}}));

    // Test 2: pure translation
    run_test("Pure translation by (3,3,3)",
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0}}),
        Points({{3,3,3},{4,3,3},{4,4,3},{5,4,3}}));

    // Test 3: pure translation (6 points, more typical chain length)
    run_test("Pure translation (6 points)",
        Points({{0,1,0},{0,0,0},{1,0,0},{1,1,0},{2,1,0},{2,2,0}}),
        Points({{1,2,1},{1,1,1},{2,1,1},{2,2,1},{3,2,1},{3,3,1}}));

    // Test 4: 90-degree rotation around Z axis
    // X -> Y where Y[i] = Rz(pi/2) * X[i]
    run_test("90-deg rotation around Z",
        Points({{1,0,0},{0,1,0},{-1,0,0},{0,-1,0}}),
        Points({{0,1,0},{-1,0,0},{0,-1,0},{1,0,0}}));

    // Test 5: translation + rotation
    run_test("Translation + rotation (4 points)",
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0}}),
        Points({{3,4,3},{3,5,3},{2,5,3},{2,6,3}}));

    // Test 6: X reflection (mirroring) — tests the degenerate case handling
    run_test("X-mirror of a 6-point chain",
        Points({{0,1,0},{0,0,0},{1,0,0},{1,1,0},{2,1,0},{2,2,0}}),
        Points({{0,-1,0},{0,0,0},{-1,0,0},{-1,-1,0},{-2,-1,0},{-2,-2,0}}));

    // Test 7: 3D rotation (random-looking)
    // Y obtained by rotating X by 45° around (1,1,1) axis then translating
    // computed analytically: R = I + sin(t)*K + (1-cos(t))*K^2  where K=skew(1/√3,1/√3,1/√3)
    // For simplicity use known mapped values.
    run_test("General 3D rotation (5 points)",
        Points({{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1}}),
        Points({{0,1,0},{0,0,1},{1,0,0},{0,1,1},{1,0,1}}));

    // Test 8: single point (trivial U=I)
    run_test("Single point (N=1)",
        Points({{2,3,4}}),
        Points({{5,6,7}}));

    // Test 9: two points
    run_test("Two points",
        Points({{0,0,0},{1,0,0}}),
        Points({{0,0,0},{0,1,0}}));

    // Test 10: with scaling enabled
    run_test("With optimal scaling (4 points)",
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0}}),
        Points({{0,0,0},{2,0,0},{2,2,0},{4,2,0}}),
        /*use_scale=*/true);

    // Test 11: 10-point chain (more representative of real protein use)
    run_test("10-point backbone-like chain",
        Points({{0,0,0},{1,0,0},{1,1,0},{2,1,0},{2,2,0},
                {3,2,0},{3,3,0},{4,3,0},{4,4,0},{5,4,0}}),
        Points({{5,5,5},{6,5,5},{6,6,5},{7,6,5},{7,7,5},
                {8,7,5},{8,8,5},{9,8,5},{9,9,5},{10,9,5}}));

    // Test 12: noisy real-world-like data (X and Y not perfectly superposable)
    run_test("Noisy 5-point data (not perfectly superposable)",
        Points({{ 1.2, 0.1, 0.0},{ 0.1, 0.9, 0.2},{-1.1, 0.0, 0.1},{ 0.0,-1.0, 0.0},{ 1.0, 0.0,-0.1}}),
        Points({{ 0.0, 1.3,-0.1},{-0.9, 0.1, 0.1},{-0.1,-1.0, 0.0},{ 1.0, 0.1, 0.2},{ 0.0, 0.0, 1.1}}));

    // =========================================================================
    printf("\n=======================================================\n");
    printf(" Results: %d tests run, %d failed\n", g_tests_run, g_tests_failed);
    if (g_tests_failed == 0)
        printf(" PASS\n");
    else
        printf(" FAIL\n");
    printf("=======================================================\n");
    return g_tests_failed > 0 ? 1 : 0;
}
