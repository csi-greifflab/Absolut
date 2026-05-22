// =============================================================================
// kabsch_nogsl.cc
//
// GSL-free implementation of the Kabsch superposition algorithm.
// See kabsch_nogsl.hh for full documentation.
//
// Only standard headers are required:
//   <cmath>   — sqrt, fabs, fmod
//   <cstring> — memcpy (not actually used; left for clarity)
// =============================================================================

#include "kabsch_nogsl.hh"

#include <cmath>
#include <algorithm>   // std::swap

using namespace kabsch_nogsl_detail;

#define KABSCH_NOGSL_NORM_EPS 0.00000001

// Convenience macros for row-major flat array indexing
#define X_(i,j)  X[(i)*3+(j)]
#define Y_(i,j)  Y[(i)*3+(j)]

int kabsch_superpositioning_nogsl(
        unsigned int    size,
        double*         X,
        double*         Y,
        double          U[3][3],
        double          t[3],
        double*         s
) {
    int U_ok = 1;

    // ------------------------------------------------------------------
    // 1. Compute centroids
    // ------------------------------------------------------------------
    Vec3 cx, cy;
    for (unsigned int i = 0; i < size; ++i) {
        cx[0] += X_(i,0);  cx[1] += X_(i,1);  cx[2] += X_(i,2);
        cy[0] += Y_(i,0);  cy[1] += Y_(i,1);  cy[2] += Y_(i,2);
    }
    double n = 1.0 / size;
    cx[0]*=n;  cx[1]*=n;  cx[2]*=n;
    cy[0]*=n;  cy[1]*=n;  cy[2]*=n;

    // ------------------------------------------------------------------
    // 2. Center both point sets
    // ------------------------------------------------------------------
    for (unsigned int i = 0; i < size; ++i) {
        X_(i,0) -= cx[0];  X_(i,1) -= cx[1];  X_(i,2) -= cx[2];
        Y_(i,0) -= cy[0];  Y_(i,1) -= cy[1];  Y_(i,2) -= cy[2];
    }

    // ------------------------------------------------------------------
    // 3. Build R = Y^T * X  (Kabsch cross-covariance matrix)
    //    R[i][j] = sum_k  Y[k][i] * X[k][j]
    // ------------------------------------------------------------------
    Mat3 R;
    for (unsigned int k = 0; k < size; ++k)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                R(i,j) += Y_(k,i) * X_(k,j);

    // ------------------------------------------------------------------
    // 4. RTR = R^T * R  (symmetric, used for eigenproblem)
    // ------------------------------------------------------------------
    Mat3 RT  = transpose(R);
    Mat3 RTR = matmul(RT, R);

    // ------------------------------------------------------------------
    // 5. Eigen-decompose RTR via Jacobi iteration
    //    eigenvalues -> eval[3],  eigenvectors (columns) -> evec
    // ------------------------------------------------------------------
    double eval[3];
    Mat3   evec;

    if (size == 1) {
        // Trivial rotation
        U[0][0]=U[1][1]=U[2][2]=1.0;
        U[0][1]=U[0][2]=U[1][0]=U[1][2]=U[2][0]=U[2][1]=0.0;
    } else {
        jacobi_sym3(RTR, eval, evec);

        // ------------------------------------------------------------------
        // 6. Sort eigenpairs descending (same as GSL_EIGEN_SORT_VAL_DESC)
        // ------------------------------------------------------------------
        sort_eigenpairs_desc(eval, evec);

        if (eval[1] > KABSCH_NOGSL_NORM_EPS) {
            // ---------------------------------------------------------------
            // 7. Build Kabsch's ak and bk frames
            //    ak = columns of evec (already computed)
            //    bk = R * ak  (normalised)
            //    a2 = a0 x a1   (complete right-handed frame)
            // ---------------------------------------------------------------
            Vec3 a0 = evec.col(0);
            Vec3 a1 = evec.col(1);
            Vec3 a2 = cross(a0, a1);       // a2 = a0 x a1
            evec.setCol(2, a2);            // overwrite 3rd column

            Vec3 b0 = matvec(R, a0);
            Vec3 b1 = matvec(R, a1);
            double norm_b0 = norm(b0);
            double norm_b1 = norm(b1);

            if (norm_b0 > KABSCH_NOGSL_NORM_EPS && norm_b1 > KABSCH_NOGSL_NORM_EPS) {
                b0 = scale(b0, 1.0/norm_b0);
                b1 = scale(b1, 1.0/norm_b1);
                Vec3 b2 = cross(b0, b1);
                double norm_b2 = norm(b2);

                if (norm_b2 > KABSCH_NOGSL_NORM_EPS) {
                    // ----------------------------------------------------------
                    // 8. U = B * A^T
                    //    B = [b0 | b1 | b2]  (columns), A = evec (columns)
                    // ----------------------------------------------------------
                    Mat3 B, A;
                    B.setCol(0, b0);  B.setCol(1, b1);  B.setCol(2, b2);
                    A = evec;         // columns are a0, a1, a2

                    Mat3 AT = transpose(A);
                    Mat3 Umat = matmul(B, AT);

                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            U[i][j] = Umat(i,j);
                } else {
                    U_ok = 0;
                    U[0][0]=U[1][1]=U[2][2]=1.0;
                    U[0][1]=U[0][2]=U[1][0]=U[1][2]=U[2][0]=U[2][1]=0.0;
                }
            } else {
                U_ok = 0;
                U[0][0]=U[1][1]=U[2][2]=1.0;
                U[0][1]=U[0][2]=U[1][0]=U[1][2]=U[2][0]=U[2][1]=0.0;
            }
        } else {
            U_ok = 0;
            U[0][0]=U[1][1]=U[2][2]=1.0;
            U[0][1]=U[0][2]=U[1][0]=U[1][2]=U[2][0]=U[2][1]=0.0;
        }
    }

    // ------------------------------------------------------------------
    // 9. Optionally compute optimal isotropic scaling
    //    s = <Y, U*X> / <U*X, U*X>   (after centring)
    // ------------------------------------------------------------------
    if (s) {
        *s = 1.0;
        if (U_ok && size > 1) {
            double nom = 0.0, dom = 0.0;
            for (unsigned int i = 0; i < size; ++i) {
                Vec3 xi(X_(i,0), X_(i,1), X_(i,2));
                Vec3 yi(Y_(i,0), Y_(i,1), Y_(i,2));
                // Uxi = U * xi
                Vec3 Uxi( U[0][0]*xi[0]+U[0][1]*xi[1]+U[0][2]*xi[2],
                           U[1][0]*xi[0]+U[1][1]*xi[1]+U[1][2]*xi[2],
                           U[2][0]*xi[0]+U[2][1]*xi[1]+U[2][2]*xi[2] );
                nom += dot(yi, Uxi);
                dom += dot(Uxi, Uxi);
            }
            if (dom > KABSCH_NOGSL_NORM_EPS) *s = nom / dom;
        }
        // Scale cx by s (mirrors GSL path used before computing t)
        cx[0] *= *s;  cx[1] *= *s;  cx[2] *= *s;
    }

    // ------------------------------------------------------------------
    // 10. Translation  t = cy - s * U * cx
    // ------------------------------------------------------------------
    Vec3 Ucx( U[0][0]*cx[0]+U[0][1]*cx[1]+U[0][2]*cx[2],
               U[1][0]*cx[0]+U[1][1]*cx[1]+U[1][2]*cx[2],
               U[2][0]*cx[0]+U[2][1]*cx[1]+U[2][2]*cx[2] );
    t[0] = cy[0] - Ucx[0];
    t[1] = cy[1] - Ucx[1];
    t[2] = cy[2] - Ucx[2];

    return U_ok;
}

#undef X_
#undef Y_
