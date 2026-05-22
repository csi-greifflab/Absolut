#ifndef BIU_KABSCH_NOGSL_HH_
#define BIU_KABSCH_NOGSL_HH_

// =============================================================================
// kabsch_nogsl.hh / kabsch_nogsl.cc
//
// Drop-in replacement for the GSL-based Kabsch superposition used in
// SuperPos_Kabsch.cc.  Implements the same interface as the original
// kabsch_superpositioning() but with zero external dependencies:
//   - No GSL headers, no GSL link flag.
//   - Uses only <cmath>, <algorithm> (C++11).
//

#include <cmath>
#include <algorithm>   // std::swap
// Algorithm
// ---------
// Given N point pairs (X[i], Y[i]):
//   1. Compute centroids cx, cy.
//   2. Center both sets.
//   3. Build the 3x3 cross-covariance matrix  R = Y^T * X.
//   4. Compute RTR = R^T * R  (symmetric, 3x3).
//   5. Eigen-decompose RTR via Jacobi iteration -> eigenvalues (eval),
//      eigenvectors (evec, stored as columns).
//   6. Sort eigenpairs by eigenvalue descending.
//   7. Build bk = R * ak / ||R * ak||  (Kabsch's orthonormal frame).
//   8. U = B * A^T.
//   9. Optionally compute optimal isotropic scaling s.
//  10. Translation  t = cy - s * U * cx.
//
// This reproduces exactly the same numerical steps as the original GSL version.
// =============================================================================

// Mat3 and Vec3: lightweight fixed-size helpers used only inside kabsch_nogsl.
// They are deliberately kept in a private namespace to avoid clashing with
// biu::Matrix or any other project type.

namespace kabsch_nogsl_detail {

    // -------------------------------------------------------------------------
    // Vec3 — a 3-element double vector
    // -------------------------------------------------------------------------
    struct Vec3 {
        double v[3];
        Vec3()  { v[0]=v[1]=v[2]=0.0; }
        Vec3(double a, double b, double c) { v[0]=a; v[1]=b; v[2]=c; }
        double& operator[](int i)       { return v[i]; }
        double  operator[](int i) const { return v[i]; }
    };

    inline Vec3 operator+(const Vec3& a, const Vec3& b) {
        return Vec3(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
    }
    inline Vec3 operator-(const Vec3& a, const Vec3& b) {
        return Vec3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
    }
    inline Vec3 operator*(double s, const Vec3& a) {
        return Vec3(s*a[0], s*a[1], s*a[2]);
    }
    inline double dot(const Vec3& a, const Vec3& b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    inline double norm(const Vec3& a) {
        return std::sqrt(dot(a,a));
    }
    inline Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3( a[1]*b[2] - b[1]*a[2],
                     a[2]*b[0] - b[2]*a[0],
                     a[0]*b[1] - b[0]*a[1] );
    }
    inline Vec3 scale(const Vec3& a, double s) {
        return Vec3(s*a[0], s*a[1], s*a[2]);
    }

    // -------------------------------------------------------------------------
    // Mat3 — a 3x3 double matrix stored row-major: m[row][col]
    // -------------------------------------------------------------------------
    struct Mat3 {
        double m[3][3];
        Mat3() { for(int i=0;i<3;i++) for(int j=0;j<3;j++) m[i][j]=0.0; }

        // Access
        double& operator()(int r, int c)       { return m[r][c]; }
        double  operator()(int r, int c) const { return m[r][c]; }

        // Column vector
        Vec3 col(int c) const { return Vec3(m[0][c], m[1][c], m[2][c]); }

        // Set column
        void setCol(int c, const Vec3& v) {
            m[0][c]=v[0]; m[1][c]=v[1]; m[2][c]=v[2];
        }

        // Identity
        static Mat3 identity() {
            Mat3 I;
            I(0,0)=I(1,1)=I(2,2)=1.0;
            return I;
        }
    };

    // Matrix multiply C = A * B
    inline Mat3 matmul(const Mat3& A, const Mat3& B) {
        Mat3 C;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                    C(i,j) += A(i,k)*B(k,j);
        return C;
    }

    // Matrix-vector multiply: y = M * x
    inline Vec3 matvec(const Mat3& M, const Vec3& x) {
        return Vec3( M(0,0)*x[0]+M(0,1)*x[1]+M(0,2)*x[2],
                     M(1,0)*x[0]+M(1,1)*x[1]+M(1,2)*x[2],
                     M(2,0)*x[0]+M(2,1)*x[1]+M(2,2)*x[2] );
    }

    // Transpose
    inline Mat3 transpose(const Mat3& A) {
        Mat3 T;
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) T(i,j)=A(j,i);
        return T;
    }

    // -------------------------------------------------------------------------
    // jacobi_sym3 — classical Jacobi iteration for a 3x3 symmetric matrix.
    //
    // On entry:  A is symmetric.
    // On exit:   eval[k] = eigenvalue k
    //            evec    = matrix whose COLUMNS are the corresponding
    //                      orthonormal eigenvectors  (same convention as GSL)
    //            eigenvalues are NOT sorted here; caller must sort.
    //
    // Convergence is guaranteed for 3x3; typically < 10 sweeps.
    // -------------------------------------------------------------------------
    inline void jacobi_sym3(Mat3 A, double eval[3], Mat3& evec) {
        // Start with identity for eigenvectors
        evec = Mat3::identity();

        const int MAX_ITER = 100;
        const double EPS   = 1e-12;

        for (int iter = 0; iter < MAX_ITER; ++iter) {
            // Find largest off-diagonal element
            double maxOff = 0.0;
            int p = 0, q = 1;
            for (int r = 0; r < 3; ++r) {
                for (int c = r+1; c < 3; ++c) {
                    if (std::fabs(A(r,c)) > maxOff) {
                        maxOff = std::fabs(A(r,c));
                        p = r; q = c;
                    }
                }
            }
            if (maxOff < EPS) break; // Converged

            // Compute the Jacobi rotation angle
            double theta = 0.5*(A(q,q) - A(p,p)) / A(p,q);
            double t;
            if (theta >= 0.0)
                t =  1.0 / (theta + std::sqrt(1.0 + theta*theta));
            else
                t = -1.0 / (-theta + std::sqrt(1.0 + theta*theta));

            double cosT = 1.0 / std::sqrt(1.0 + t*t);
            double sinT = t * cosT;

            // Update A in place: A' = G^T * A * G  (Jacobi rotation G)
            // This is the standard Golub & Van Loan update.
            double App = A(p,p) - t*A(p,q);
            double Aqq = A(q,q) + t*A(p,q);
            A(p,p) = App;
            A(q,q) = Aqq;
            A(p,q) = A(q,p) = 0.0;

            for (int r = 0; r < 3; ++r) {
                if (r == p || r == q) continue;
                double Arp = A(r,p);
                double Arq = A(r,q);
                A(r,p) = A(p,r) = cosT*Arp - sinT*Arq;
                A(r,q) = A(q,r) = sinT*Arp + cosT*Arq;
            }

            // Accumulate rotation in evec (columns = eigenvectors)
            for (int r = 0; r < 3; ++r) {
                double evp = evec(r,p);
                double evq = evec(r,q);
                evec(r,p) = cosT*evp - sinT*evq;
                evec(r,q) = sinT*evp + cosT*evq;
            }
        }

        // Extract diagonal as eigenvalues
        for (int i = 0; i < 3; ++i) eval[i] = A(i,i);
    }

    // -------------------------------------------------------------------------
    // sort_eigenpairs_desc — sort eigenvalues (and matching eigenvector columns)
    // in descending order.  Simple selection sort on 3 elements.
    // -------------------------------------------------------------------------
    inline void sort_eigenpairs_desc(double eval[3], Mat3& evec) {
        for (int i = 0; i < 2; ++i) {
            int maxIdx = i;
            for (int j = i+1; j < 3; ++j)
                if (eval[j] > eval[maxIdx]) maxIdx = j;
            if (maxIdx != i) {
                std::swap(eval[i], eval[maxIdx]);
                // Swap columns i and maxIdx in evec
                for (int r = 0; r < 3; ++r)
                    std::swap(evec(r,i), evec(r,maxIdx));
            }
        }
    }

} // namespace kabsch_nogsl_detail


// =============================================================================
// Public interface — mirrors the original C function in SuperPos_Kabsch.cc
//
// Parameters
// ----------
//   size   : number of point pairs
//   X      : [size x 3] array (row-major) — points to be moved; MODIFIED in
//             place (centred) by this call, matching GSL behaviour
//   Y      : [size x 3] array (row-major) — target points;       MODIFIED
//   U      : [3 x 3]  OUTPUT rotation matrix (row-major)
//   t      : [3]      OUTPUT translation vector
//   s      : pointer to OUTPUT scaling (may be nullptr to skip scaling)
//
// Returns 1 if rotation is valid, 0 if degenerate (U set to identity).
// =============================================================================
int kabsch_superpositioning_nogsl(
        unsigned int    size,
        double*         X,   // size x 3, row-major
        double*         Y,   // size x 3, row-major
        double          U[3][3],
        double          t[3],
        double*         s
);

#endif // BIU_KABSCH_NOGSL_HH_
