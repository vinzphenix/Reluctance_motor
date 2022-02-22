#include "motor.h"
#include <math.h>
#include <time.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// - TOOL FUNCTIONS - //////////////////////////////////////////////////////

double femMin(double *x, int n) {
    double myMin = x[0];
    int i;
    for (i = 1; i < n; i++)
        myMin = fmin(myMin, x[i]);
    return myMin;
}

double femMax(double *x, int n) {
    double myMax = x[0];
    int i;
    for (i = 1; i < n; i++)
        myMax = fmax(myMax, x[i]);
    return myMax;
}

void femError(char *text, int line, char *file) {
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femErrorScan(int test, int line, char *file) {
    if (test >= 0)
        return;

    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file) {
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - GLOBAL VARIABLES - /////////////////////////////////////////////////////

int reMeshInit = 0;
int solveInit = 0;

motorReMesh *globalReMesh; // Remaillage
double *globalX, *globalY;
double *theGlobalCoord;

int *number;    // numerotation des noeuds
int *dirichlet; // conditions frontières

double *ddH; // d^2(H)/dB^2  [size = 43]
double deltaTheta = 0.;
int currentCoil = 1;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - INTEGRATION RULE - /////////////////////////////////////////////////////

static const double _gaussTriXsi[3] = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTriEta[3] = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTriWeight[3] = {0.166666666666667, 0.166666666666667, 0.166666666666667};
static const int _nInteg = 3;

/*static const double _gaussTriXsi[1]     = { 0.333333333333333};
static const double _gaussTriEta[1]     = { 0.333333333333333};
static const double _gaussTriWeight[1]  = { 0.500000000000000};
static const int _nInteg = 1;*/

static const double _dphidxsi[3] = {-1.0, 1.0, 0.0};
static const double _dphideta[3] = {-1.0, 0.0, 1.0};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////// - COO MATRIX - ////////////////////////////////////////////////////////

// COO functions with *** where already implemented in a template of the course LINMA1170

int compareijkSLOW(const void *_a, const void *_b) { // ***
    ijk *a = (ijk *)_a;
    ijk *b = (ijk *)_b;
    if (a->i < b->i)
        return -1;
    if (a->i > b->i)
        return 1;
    if (a->j < b->j)
        return -1;
    if (a->j > b->j)
        return 1;
    return 0;
}

int compareijk(const void *_a, const void *_b) {
    if (((ijk *)_a)->i < ((ijk *)_b)->i)
        return -1;
    else if (((ijk *)_a)->i > ((ijk *)_b)->i)
        return 1;
    else if (((ijk *)_a)->j < ((ijk *)_b)->j)
        return -1;
    else if (((ijk *)_a)->j > ((ijk *)_b)->j)
        return 1;
    return 0;
}

COO *allocateCOO(size_t initial) { // ***
    COO *coo = (COO *)malloc(sizeof(COO));
    coo->compressed = 0;
    coo->nnz = 0;
    coo->alloc = initial;
    coo->data = (ijk *)malloc(initial * sizeof(ijk));
    return coo;
}

void freeCOO(COO *myCOO) { // ***
    free(myCOO->data);
    free(myCOO);
}

// add element to the matrix and increase size if needed
void addToCoo(COO *coo, size_t I, size_t J, double K) { // ***
    // if max size is reached: we double the size of the memory allocated
    if (coo->nnz == coo->alloc) {
        ijk *temp = (ijk *)malloc(2 * coo->alloc * sizeof(ijk));
        for (size_t i = 0; i < coo->alloc; i++) {
            temp[i].i = coo->data[i].i;
            temp[i].j = coo->data[i].j;
            temp[i].k = coo->data[i].k;
        }
        free(coo->data);
        coo->data = temp;
        coo->alloc *= 2;
    }
    coo->compressed = 0;
    coo->data[coo->nnz].i = I;
    coo->data[coo->nnz].j = J;
    coo->data[coo->nnz++].k = K;
}

// sort lexicographically then compress
void sortAndCompress(COO *coo) {
    if (coo->compressed == 1)
        return;

    qsort(coo->data, coo->nnz, sizeof(ijk), compareijk);

    size_t i, n, c = 0;
    ijk *data = coo->data;
    for (i = 1, n = 0; i < coo->nnz; i++) {
        if (compareijk(data + c, data + i)) { // new element
            c = ++n;
            data[n].i = data[i].i;
            data[n].j = data[i].j;
            data[n].k = data[i].k;
        } else
            data[n].k += data[i].k;
    }
    // update status
    coo->nnz = n + 1;
    coo->compressed = 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////// - CSR MATRIX - ////////////////////////////////////////////////////////

// convert sparse COO matrix to sparse CSR matrix
CSR *convertToCSR(COO *A, size_t size) {
    // A    : matrix in COO form
    // size : size of the square matrix
    // return : the matrix in CSR form

    ijk *data = A->data;
    size_t nnz = A->nnz;

    CSR *myCSR = (CSR *)malloc(sizeof(CSR));
    myCSR->nnz = nnz;
    myCSR->size = size;

    myCSR->IA = (size_t *)malloc((2 * size + 1 + nnz) * sizeof(size_t)); // size+1
    myCSR->JA = myCSR->IA + size + 1;                                    // nnz
    myCSR->DIAG = myCSR->JA + nnz;                                       // size

    myCSR->VA = (double *)malloc(nnz * sizeof(double));

    size_t i, q = 0;
    // for (i = 0 ; i < size + 1 ; i++)
    //     myCSR->IA[i] = 0;
    memset(myCSR->IA, 0, (size + 1) * sizeof(double));

    for (i = 0; i < nnz; i++) {
        if (data[i].i == data[i].j)
            myCSR->DIAG[q++] = i;
        myCSR->IA[data[i].i + 1]++;
        myCSR->JA[i] = data[i].j;
        myCSR->VA[i] = data[i].k;
    }

    for (i = 0; i < size; i++)
        myCSR->IA[i + 1] += myCSR->IA[i];

    return myCSR;
}

// duplicate the matrix to get the preconditionner matrix
CSR *duplicateCSR(CSR *M) {
    // M    : matrix in CSR form
    // return : matrix in CSR form

    CSR *myCSR = (CSR *)malloc(sizeof(CSR));
    myCSR->nnz = M->nnz;
    myCSR->size = M->size;
    myCSR->IA = (size_t *)malloc((M->size + 1 + M->nnz) * sizeof(size_t)); // size+1
    myCSR->JA = myCSR->IA + M->size + 1;                                   // nnz
    myCSR->VA = (double *)malloc(M->nnz * sizeof(double));

    memcpy(myCSR->IA, M->IA, (M->size + 1) * sizeof(size_t));
    memcpy(myCSR->JA, M->JA, (M->nnz) * sizeof(size_t));
    memcpy(myCSR->VA, M->VA, (M->nnz) * sizeof(double));
    return myCSR;
}

void freeCSR(CSR *myCSR) {
    free(myCSR->IA);
    free(myCSR->VA);
    free(myCSR);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - Matrix Operations - ////////////////////////////////////////////////////

// compute ILU
void ILU(CSR *LU) {
    // R. C . Mittal & A.H Al-Kurdi
    // An efficient method for constructing an ILU preconditioner
    // for solving large sparse nonsymmetric linear systems by the GMRES method

    size_t i, j, k, v, w, n = LU->size;
    size_t *DIAG = LU->DIAG;
    size_t *IA = LU->IA;
    size_t *JA = LU->JA;
    double *VA = LU->VA;

    double *point = (double *)calloc(n, sizeof(double)); // working array

    for (i = 1; i < n; i++) {
        for (v = IA[i] + 1; v < IA[i + 1]; v++)
            point[JA[v]] = v;

        for (v = IA[i]; v < DIAG[i]; v++) {
            j = JA[v];
            VA[v] /= VA[DIAG[j]];

            for (w = DIAG[j] + 1; w < IA[j + 1]; w++) {
                k = point[JA[w]];
                if (k > 0) // there is an element "above"
                    VA[k] -= VA[v] * VA[w];
            } // end for w
        }     // end for v

        for (v = IA[i] + 1; v < IA[i + 1]; v++)
            point[JA[v]] = 0;
    } // end for i

    free(point);
}

// solve LU x = b
void systemSolve(CSR *LU, double *x, double *b) {

    size_t i, j, n = LU->size;
    size_t *DIAG = LU->DIAG;
    size_t *IA = LU->IA;
    size_t *JA = LU->JA;
    double *VA = LU->VA;
    double sum;

    // Forward substitution
    x[0] = b[0];
    for (i = 1; i < n; i++) {
        sum = 0.;
        for (j = IA[i]; j < DIAG[i]; j++)
            sum += VA[j] * x[JA[j]];
        x[i] = b[i] - sum;
    }

    // Back substitution
    x[n - 1] /= VA[DIAG[n - 1]];
    for (i = n - 1; i-- > 0;) {
        // start at n-2 and finish at 0 included (complicated with size_t !)
        sum = 0.;
        for (j = DIAG[i] + 1; j < IA[i + 1]; j++)
            sum += VA[j] * x[JA[j]];
        x[i] = (x[i] - sum) / VA[DIAG[i]];
    }
}

// compute y = A * x
void matrixVectorProduct(CSR *A, double *x, double *y, size_t size) {
    size_t *IA = A->IA;
    size_t *JA = A->JA;
    double *VA = A->VA;
    size_t end;

    for (size_t i = 0; i < size; i++) {
        y[i] = 0.;
        end = IA[i + 1];
        for (size_t j = IA[i]; j < end; j++) {
            y[i] += VA[j] * x[JA[j]];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - ITERATIVE SOLVER - /////////////////////////////////////////////////////

// create the solver
sparseSolver *sparseSolverCreate(size_t size) {
    sparseSolver *mySolver = (sparseSolver *)malloc(sizeof(sparseSolver));
    mySolver->size = size;

    mySolver->x = (double *)malloc(2 * size * sizeof(double));
    mySolver->b = mySolver->x + size;

    mySolver->R = (double *)malloc(4 * size * sizeof(double));
    mySolver->Z = mySolver->R + size;
    mySolver->D = mySolver->R + 2 * size;
    mySolver->S = mySolver->R + 3 * size;

    return mySolver;
}

// set x and b of the system A x = b
// !!! NUMBER NEEDS to be a permutation
void setSystemVectors(motor *theMotor, sparseSolver *mySolver) {

    for (int i = 0; i < theMotor->size; i++) {
        mySolver->x[number[i]] = theMotor->a[i];
        mySolver->b[i] = 0.;
    }
}

void freeSparseSolver(sparseSolver *mySolver) {
    free(mySolver->R);
    free(mySolver->x);
    // free(mySolver->dirichlet); now in GLOBAL
    // free(mySolver->number); now in GLOBAL
    // freeCOO(mySolver->A); already DONE in convertToCSR
    freeCSR(mySolver->B);
    freeCSR(mySolver->M);
    free(mySolver);
}

// initialize the CG before iterating
void initCG(sparseSolver *mySolver) {

    int size = mySolver->size;
    sortAndCompress(mySolver->A);

    // saveMatrix(mySolver->A, size);  // to see the sparsity pattern

    mySolver->M = convertToCSR(mySolver->A, size); // transform to CSR
    freeCOO(mySolver->A);                          // don't need it anymore
    mySolver->B = duplicateCSR(mySolver->M);       // duplicate to perform preconditionner

    mySolver->rr = 0.;
    mySolver->norm_b = 0.;

    matrixVectorProduct(mySolver->B, mySolver->x, mySolver->S, size);

    for (int i = 0; i < size; i++) {
        mySolver->R[i] = mySolver->b[i] - mySolver->S[i];
        mySolver->norm_b += (mySolver->b[i]) * (mySolver->b[i]);
    }
    mySolver->norm_b = sqrt(mySolver->norm_b);

    ILU(mySolver->M);
    systemSolve(mySolver->M, mySolver->Z, mySolver->R);

    for (int i = 0; i < size; i++) {
        mySolver->D[i] = mySolver->Z[i];
        mySolver->rr += (mySolver->R[i]) * (mySolver->Z[i]);
    }
}

// Apply the CG algorithm
int PCG(sparseSolver *mySolver, size_t itermax, double precision) {

    CSR *A = mySolver->B;
    CSR *M = mySolver->M;
    double *R = mySolver->R; // residue
    double *Z = mySolver->Z; // result of M Z = R
    double *D = mySolver->D; // direction
    double *S = mySolver->S; // result of A * D
    double *x = mySolver->x; // solution of sytem

    double alpha, beta, rrNew;
    double rr = mySolver->rr;
    double norm_b = mySolver->norm_b;
    double norm_r = norm_b;

    size_t size = mySolver->size;
    size_t nnz = A->nnz;
    size_t iter, i;

    for (iter = 1; (iter < itermax) && (precision < norm_r / norm_b); iter++) {

        alpha = 0.;
        matrixVectorProduct(A, D, S, size);
        for (i = 0; i < size; i++)
            alpha += S[i] * D[i];
        alpha = rr / alpha;

        for (i = 0; i < size; i++) {
            x[i] += alpha * D[i];
            R[i] -= alpha * S[i];
        }

        systemSolve(M, Z, R);

        rrNew = 0.;
        norm_r = 0.;
        for (i = 0; i < size; i++) {
            rrNew += R[i] * Z[i];
            norm_r += R[i] * R[i];
        }
        norm_r = sqrt(norm_r);

        beta = rrNew / rr;
        rr = rrNew;
        for (i = 0; i < size; i++)
            D[i] = Z[i] + beta * D[i];
    }

    // printf("CONVERGED ? size = %ld \t nbIterations = %ld \t error = %le\n", size, iter, norm_r/norm_b);
    return iter;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - MATRIX ASSEMBLE - //////////////////////////////////////////////////////

// get all needed informations about nodes of the current element
void getLocalInfo(const motor *theMotor, const int *number, const int *dirichlet,
                  const int iElem, int *map, int *ctr, double *x, double *y, double *u) {
    motorMesh *theMesh = theMotor->mesh;
    int nLoc = theMesh->nLocalNode;

    for (int j = 0; j < nLoc; j++) {
        map[j] = theMesh->elem[iElem * nLoc + j];
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
        ctr[j] = dirichlet[map[j]];
        u[j] = theMotor->a[map[j]];
        map[j] = number[map[j]];
    }
}

// apply a dirichlet boundary condition a = 0
void applyConstraint(double *Aloc, double *Bloc, int myNode, int nLoc) {
    for (int i = 0; i < nLoc; i++) {
        Aloc[i * nLoc + myNode] = 0.;
        Aloc[myNode * nLoc + i] = 0.;
    }
    Aloc[myNode * nLoc + myNode] = 1.;
    Bloc[myNode] = 0.;
}

// insert local system in the full system
void insertLocalToGlobal(COO *A, double *b, double *Aloc, double *Bloc, int *map, int nLoc) {
    int i, j;
    for (i = 0; i < nLoc; i++) {
        b[map[i]] += Bloc[i];
        for (j = 0; j < nLoc; j++)
            addToCoo(A, map[i], map[j], Aloc[i * nLoc + j]);
    }
}

// Compute the matrix of the linear system
void assembleMotorMatrix(motor *theMotor, sparseSolver *theSolver) {
    motorMesh *theMesh = theMotor->mesh;

    int nLoc = theMesh->nLocalNode;
    int iElem, iInteg, i, j, map[3], ctr[3];
    double Xloc[3], Yloc[3], Uloc[3], phi[3], dphidx[3], dphidy[3];

    double xsi, eta, weight, dxdxsi, dxdeta, dydxsi, dydeta, jac, invJ;
    double js, invMu, dInvMudb2;

    double *Aloc = (double *)malloc(sizeof(double) * nLoc * nLoc);
    double *Bloc = (double *)malloc(sizeof(double) * nLoc);

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        invMu = 1. / theMotor->mu[theMesh->domain[iElem]];
        js = theMotor->js[theMesh->domain[iElem]];

        memset(Aloc, 0, nLoc * nLoc * sizeof(double));
        memset(Bloc, 0, nLoc * sizeof(double));
        getLocalInfo(theMotor, number, dirichlet, iElem, map, ctr, Xloc, Yloc, Uloc);

        for (iInteg = 0; iInteg < _nInteg; iInteg++) {
            xsi = _gaussTriXsi[iInteg];
            eta = _gaussTriEta[iInteg];
            weight = _gaussTriWeight[iInteg];
            phi[0] = 1. - xsi - eta;
            phi[1] = xsi;
            phi[2] = eta;

            dxdxsi = 0.;
            dxdeta = 0.;
            dydxsi = 0.;
            dydeta = 0.;

            for (i = 0; i < nLoc; i++) {
                dxdxsi += Xloc[i] * _dphidxsi[i];
                dxdeta += Xloc[i] * _dphideta[i];
                dydxsi += Yloc[i] * _dphidxsi[i];
                dydeta += Yloc[i] * _dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            invJ = 1. / jac;
            for (i = 0; i < nLoc; i++) {
                dphidx[i] = (_dphidxsi[i] * dydeta - _dphideta[i] * dydxsi) * invJ;
                dphidy[i] = (_dphideta[i] * dxdxsi - _dphidxsi[i] * dxdeta) * invJ;
            }

            for (i = 0; i < nLoc; i++) {
                for (j = 0; j < nLoc; j++)
                    Aloc[i * nLoc + j] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight * invMu;
                Bloc[i] += phi[i] * jac * js * weight;
            }
        } // end of iInteg

        for (i = 0; i < nLoc; i++)
            if (ctr[i] == 1)
                applyConstraint(Aloc, Bloc, i, nLoc);

        insertLocalToGlobal(theSolver->A, theSolver->b, Aloc, Bloc, map, nLoc);

    } // end of iElem

    free(Aloc);
    free(Bloc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////// - NON LINEAR - ////////////////////////////////////////////////////////

double *getCubicSplines(const double *B, const double *H, int n) {
    // Solve tridiagonal system (abc) x = d, where x is H''(B)
    int i;
    double a[n], b[n], c[n], d[n], h1, h2, w;
    double *x = (double *)calloc(n, sizeof(double));

    a[0] = 0.;
    a[n - 1] = 0.;
    b[0] = 1.;
    b[n - 1] = 1.;
    c[0] = 0.;
    c[n - 1] = 0.;
    d[0] = 0.;
    d[n - 1] = 0.;

    // Construct Tridiagonal matrix
    for (i = 1; i < n - 1; i++) {
        h1 = B[i] - B[i - 1];
        h2 = B[i + 1] - B[i];
        a[i] = h1 / 6.;
        b[i] = 2. * (h1 + h2) / 6.;
        c[i] = h2 / 6.;
        d[i] = (H[i + 1] - H[i]) / h2 - (H[i] - H[i - 1]) / h1;
    }

    // "gaussian elimination"
    for (i = 1; i < n; i++) {
        w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }

    // back substitution
    x[n - 1] = d[n - 1] / b[n - 1];
    for (i = n - 2; i > 0; i--) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }

    return x;
}

void getMu(const double *B, const double *H, int n, double *ddH, double b2, double db2, double *mu) {
    // at end of function :
    //   mu[0] = 1 / mu
    //   mu[1] = d/db^2 (1/mu)
    double w1, w2, Db;
    double b_mid = sqrt(b2);
    double b_prev = sqrt(b2 - db2);
    double b_next = sqrt(b2 + db2);

    for (int i = 0; i < n - 1; i++) {
        if (B[i] <= b_mid && b_mid < B[i + 1]) {
            Db = B[i + 1] - B[i];
            w1 = B[i + 1] - b_mid;
            w2 = b_mid - B[i];
            mu[0] = ddH[i] / (6. * Db) * w1 * w1 * w1 + ddH[i + 1] / (6. * Db) * w2 * w2 * w2 + (H[i] / Db - ddH[i] * Db / 6.) * w1 + (H[i + 1] / Db - ddH[i + 1] * Db / 6.) * w2; // here we get H(b)
            mu[0] /= b_mid;                                                                                                                                                        // we divide by b to get 1/mu
        }
        // Compute a centered finite difference for derviative approximation
        if (B[i] <= b_prev && b_prev < B[i + 1]) {
            Db = B[i + 1] - B[i];
            w1 = B[i + 1] - b_prev;
            w2 = b_prev - B[i];
            mu[1] = ddH[i] / (6. * Db) * w1 * w1 * w1 + ddH[i + 1] / (6. * Db) * w2 * w2 * w2 + (H[i] / Db - ddH[i] * Db / 6.) * w1 + (H[i + 1] / Db - ddH[i + 1] * Db / 6.) * w2;
            mu[1] /= b_prev;
        }

        if (B[i] <= b_next && b_next < B[i + 1]) {
            Db = B[i + 1] - B[i];
            w1 = B[i + 1] - b_next;
            w2 = b_next - B[i];
            mu[2] = ddH[i] / (6. * Db) * w1 * w1 * w1 + ddH[i + 1] / (6. * Db) * w2 * w2 * w2 + (H[i] / Db - ddH[i] * Db / 6.) * w1 + (H[i + 1] / Db - ddH[i + 1] * Db / 6.) * w2;
            mu[2] /= b_next;
            break;
        }
    }
    mu[1] = (mu[2] - mu[1]) / (2 * db2);
}

// very similar to assembleMotorMatrix
void assembleJacobian(motor *theMotor, sparseSolver *theSolver) {
    motorMesh *theMesh = theMotor->mesh;

    int nLoc = theMesh->nLocalNode;
    int iElem, iInteg, i, j, flag, map[3], ctr[3];
    double Xloc[3], Yloc[3], Uloc[3], phi[3], dphidx[3], dphidy[3];

    double xsi, eta, weight, dxdxsi, dxdeta, dydxsi, dydeta, jac, invJ;
    double js, invMu, dInvMudb2;

    double dadx, dady, b2; // used for nonlinear
    double term1, term2, product, db2 = 1e-8;
    const double *B = theMotor->hystereticCurveB;
    const double *H = theMotor->hystereticCurveH;
    int N = theMotor->nHystereticCurve;

    double *mu = (double *)malloc(3 * sizeof(double));

    double *Aloc = (double *)malloc(sizeof(double) * nLoc * nLoc);
    double *Bloc = (double *)malloc(sizeof(double) * nLoc);

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {

        invMu = 1. / theMotor->mu[theMesh->domain[iElem]];
        js = theMotor->js[theMesh->domain[iElem]];
        flag = (theMesh->domain[iElem] == 0 || theMesh->domain[iElem] == 8);

        for (i = 0; i < nLoc * nLoc; i++)
            Aloc[i] = 0.;
        for (i = 0; i < nLoc; i++)
            Bloc[i] = 0.;
        memset(Aloc, 0, nLoc * nLoc * sizeof(double));
        memset(Bloc, 0, nLoc * sizeof(double));
        getLocalInfo(theMotor, number, dirichlet, iElem, map, ctr, Xloc, Yloc, Uloc);

        for (iInteg = 0; iInteg < _nInteg; iInteg++) {
            xsi = _gaussTriXsi[iInteg];
            eta = _gaussTriEta[iInteg];
            weight = _gaussTriWeight[iInteg];
            phi[0] = 1. - xsi - eta;
            phi[1] = xsi;
            phi[2] = eta;

            dxdxsi = 0.;
            dxdeta = 0.;
            dydxsi = 0.;
            dydeta = 0.;
            dadx = 0.;
            dady = 0.;

            for (i = 0; i < nLoc; i++) {
                dxdxsi += Xloc[i] * _dphidxsi[i];
                dxdeta += Xloc[i] * _dphideta[i];
                dydxsi += Yloc[i] * _dphidxsi[i];
                dydeta += Yloc[i] * _dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            invJ = 1. / jac;
            for (i = 0; i < nLoc; i++) {
                dphidx[i] = (_dphidxsi[i] * dydeta - _dphideta[i] * dydxsi) * invJ;
                dphidy[i] = (_dphideta[i] * dxdxsi - _dphidxsi[i] * dxdeta) * invJ;
                dadx += Uloc[i] * dphidx[i];
                dady += Uloc[i] * dphidy[i];
            }

            if (flag) { // if nonlinear solver and mu variable
                getMu(B, H, N, ddH, dadx * dadx + dady * dady, db2, mu);
                for (i = 0; i < nLoc; i++) {
                    product = 0.;
                    for (j = 0; j < nLoc; j++) {
                        term1 = (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]);
                        term2 = (dphidx[i] * dadx + dphidy[i] * dady) * (dphidx[j] * dadx + dphidy[j] * dady);
                        product += term1 * Uloc[j];
                        Aloc[i * nLoc + j] += (mu[0] * term1 + 2. * mu[1] * term2) * jac * weight;
                    }
                    Bloc[i] += (phi[i] * js - mu[0] * product) * jac * weight;
                }
            } else { // if nonlinear solver and mu constant
                for (i = 0; i < nLoc; i++) {
                    product = 0.;
                    for (j = 0; j < nLoc; j++) {
                        term1 = (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]);
                        product += term1 * Uloc[j];
                        Aloc[i * nLoc + j] += invMu * term1 * jac * weight;
                    }
                    Bloc[i] += (phi[i] * js - invMu * product) * jac * weight;
                }
            }
        } // end of iInteg

        for (i = 0; i < nLoc; i++)
            if (ctr[i] == 1)
                applyConstraint(Aloc, Bloc, i, nLoc);

        insertLocalToGlobal(theSolver->A, theSolver->b, Aloc, Bloc, map, nLoc);

    } // end of iElem

    free(Aloc);
    free(Bloc);
    free(mu);
}

int newtonRaphson(motor *theMotor, sparseSolver *mySolver, size_t iterMaxNewton, double precisionNewton,
                  size_t iterMaxPCG, double precisionPCG, int flagSolveLinear) {

    size_t iterNewton, i, size = mySolver->size;
    double norm2_delta_x = 100.;
    int res;

    for (iterNewton = 0; (iterNewton < iterMaxNewton) && (precisionNewton < sqrt(norm2_delta_x)); iterNewton++) {

        // free matrices from previous iteration : maybe slower :-(
        if (0 < iterNewton || flagSolveLinear == 1) {
            freeCSR(mySolver->B); // reset CSR matrix (classic one)
            freeCSR(mySolver->M); // reset CSR matrix (matrix of the preconditionner)
        }
        mySolver->A = allocateCOO(3 * size); // already freed since there is a free(mySolver->A) in PCG of linear solver

        // reset x and b
        memset(mySolver->x, 0, size * sizeof(double));
        memset(mySolver->b, 0, size * sizeof(double));

        // assemble jacobian matrix : J Dx = b - Ax
        assembleJacobian(theMotor, mySolver);

        // execute PCG algorithm
        initCG(mySolver);
        res = PCG(mySolver, iterMaxPCG, precisionPCG);

        // increment the solution with the Dx obtained
        norm2_delta_x = 0.;
        for (i = 0; i < size; i++) {
            norm2_delta_x += mySolver->x[i] * mySolver->x[i];
            theMotor->a[i] += mySolver->x[number[i]];
        }

        // printf("\t nb iter = %ld\t norm_delta_x = %le\n", iterNewton, sqrt(norm2_delta_x));
    }

    // printf("nbIter = %ld norm DX = %le\n", iterNewton, sqrt(norm2_delta_x));
    return iterNewton;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////// - UPDATE MESH - ////////////////////////////////////////////////////////

void freeReMesh(motorReMesh *theReMesh) {
    free(theReMesh->rotorNodes);
    free(theReMesh->statorNodes);
    free(theReMesh);
}

int compareByNumber(const void *n1, const void *n2) {
    return *(int *)n1 - *(int *)n2;
}

int compareByAngle(const void *n1, const void *n2) {
    double angle1 = atan2(globalY[*(int *)n1], globalX[*(int *)n1]);
    double angle2 = atan2(globalY[*(int *)n2], globalX[*(int *)n2]);
    if (angle1 > angle2)
        return 1;
    else if (angle1 < angle2)
        return -1;
    return 0;
}

// split nodes of domain GAP : moving ones and static ones
// sort them by angle
void setupReMesh(motor *theMotor, motorReMesh *theReMesh) {

    // reach the Gap domain
    motorMesh *theMesh = theMotor->mesh;
    int shift, i, n, N, m, M, iDomain = 11;
    for (i = 0, shift = 0; i < iDomain; i++)
        shift += theMesh->nElemDomain[i];
    int *gapElem = &theMesh->elem[3 * shift];
    int nGapElem = theMesh->nElemDomain[iDomain];
    int *rotorNodes = (int *)malloc(sizeof(int) * 3 * nGapElem);
    int *statorNodes = (int *)malloc(sizeof(int) * 3 * nGapElem);

    // split nodes of the Gap
    for (i = 0, N = 0, M = 0; i < 3 * nGapElem; i++) {
        if (theMotor->movingNodes[gapElem[i]])
            rotorNodes[N++] = gapElem[i];
        else
            statorNodes[M++] = gapElem[i];
    }

    rotorNodes = realloc(rotorNodes, sizeof(int) * N);
    statorNodes = realloc(statorNodes, sizeof(int) * M);

    qsort(rotorNodes, N, sizeof(int), compareByNumber);
    qsort(statorNodes, M, sizeof(int), compareByNumber);

    // delete duplicates Rotor
    int currentNode = -1;
    for (i = 0, n = 0; i < N; i++) {
        if (currentNode == rotorNodes[i])
            continue;
        currentNode = rotorNodes[i];
        rotorNodes[n++] = currentNode;
    }
    // delete duplicates Stator
    currentNode = -1;
    for (i = 0, m = 0; i < M; i++) {
        if (currentNode == statorNodes[i])
            continue;
        currentNode = statorNodes[i];
        statorNodes[m++] = currentNode;
    }

    rotorNodes = realloc(rotorNodes, sizeof(int) * n);
    statorNodes = realloc(statorNodes, sizeof(int) * m);

    // sort by angle
    globalX = theMesh->X;
    globalY = theMesh->Y;
    qsort(rotorNodes, n, sizeof(int), compareByAngle);
    qsort(statorNodes, m, sizeof(int), compareByAngle);

    // fill the structure
    theReMesh->nStatorNodes = m;
    theReMesh->statorNodes = statorNodes;
    theReMesh->nRotorNodes = n;
    theReMesh->rotorNodes = rotorNodes;
    theReMesh->nGapElem = nGapElem;
    theReMesh->gapElem = gapElem;
}

double elementQuality(double *x, double *y) {
    double jac = (x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]);
    double per = hypot(x[1] - x[0], y[1] - y[0]) + hypot(x[2] - x[1], y[2] - y[1]) + hypot(x[0] - x[2], y[0] - y[2]);

    // perimeter squared is way better for small mesh since it favors much more small triangles
    return jac / (per * per);
}

void checkAndSwap(motor *theMotor, motorReMesh *theReMesh, double delta) {
    motorMesh *theMesh = theMotor->mesh;

    int *rotorNodes = theReMesh->rotorNodes;
    int *statorNodes = theReMesh->statorNodes;
    int n = theReMesh->nRotorNodes;
    int m = theReMesh->nStatorNodes;

    // Sort rotorNodes to get smaller angle in first position
    globalX = theMesh->X;
    globalY = theMesh->Y;
    qsort(rotorNodes, theReMesh->nRotorNodes, sizeof(int), compareByAngle);

    int i, idxR, idxS, consecutive = 0;
    double X[3], Y[3], qualIn, qualOut;
    X[0] = theMesh->X[rotorNodes[0]];
    Y[0] = theMesh->Y[rotorNodes[0]];
    X[1] = theMesh->X[statorNodes[0]];
    Y[1] = theMesh->Y[statorNodes[0]];

    for (i = 0, idxR = 0, idxS = 0; i < theReMesh->nGapElem; i++) {
        // very usefull that number of triangles remains constant
        X[2] = theMesh->X[statorNodes[(idxS + 1) % m]];
        Y[2] = theMesh->Y[statorNodes[(idxS + 1) % m]];
        qualOut = elementQuality(X, Y);

        X[2] = theMesh->X[rotorNodes[(idxR + 1) % n]];
        Y[2] = theMesh->Y[rotorNodes[(idxR + 1) % n]];
        qualIn = elementQuality(X, Y);

        if (qualIn < qualOut) {
            theReMesh->gapElem[3 * i] = rotorNodes[idxR];
            theReMesh->gapElem[3 * i + 1] = statorNodes[idxS];
            idxS = (idxS + 1) % m;
            theReMesh->gapElem[3 * i + 2] = statorNodes[idxS];
            X[1] = theMesh->X[statorNodes[idxS]];
            Y[1] = theMesh->Y[statorNodes[idxS]];
        } else {
            theReMesh->gapElem[3 * i] = rotorNodes[idxR];
            theReMesh->gapElem[3 * i + 1] = statorNodes[idxS];
            idxR = (idxR + 1) % n;
            theReMesh->gapElem[3 * i + 2] = rotorNodes[idxR];
            X[0] = theMesh->X[rotorNodes[idxR]];
            Y[0] = theMesh->Y[rotorNodes[idxR]];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////// - EDGES - ///////////////////////////////////////////////////////////

int motorEdgesCompare(const void *e0, const void *e1) {
    int diff = MIN(((motorEdge *)e0)->node[0], ((motorEdge *)e0)->node[1]) - MIN(((motorEdge *)e1)->node[0], ((motorEdge *)e1)->node[1]);
    if (diff < 0)
        return 1;
    if (diff > 0)
        return -1;

    diff = MAX(((motorEdge *)e0)->node[0], ((motorEdge *)e0)->node[1]) - MAX(((motorEdge *)e1)->node[0], ((motorEdge *)e1)->node[1]);
    if (diff < 0)
        return 1;
    if (diff > 0)
        return -1;
    return 0;
}

motorEdges *motorEdgesCreate(motorMesh *theMesh) {
    motorEdges *theEdges = (motorEdges *)malloc(sizeof(motorEdges));
    int nLoc = theMesh->nLocalNode;
    int i, j, n = theMesh->nElem * nLoc;
    motorEdge *edges = (motorEdge *)malloc(n * sizeof(motorEdge));
    theEdges->mesh = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;

    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i * nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc];
        }
    }

    qsort(theEdges->edges, theEdges->nEdge, sizeof(motorEdge), motorEdgesCompare);

    int index = 0;
    int nBoundary = 0;
    for (i = 0; i < theEdges->nEdge; i++, index++) {
        if (i == theEdges->nEdge - 1 || motorEdgesCompare(&edges[i], &edges[i + 1]) != 0) {
            edges[index] = edges[i];
            nBoundary++;
        } else {
            edges[index] = edges[i];
            edges[index].elem[1] = edges[i + 1].elem[0];
            i = i + 1;
        }
    }

    theEdges->edges = realloc(edges, index * sizeof(motorEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void motorEdgesFree(motorEdges *theEdges) {
    free(theEdges->edges);
    free(theEdges);
}

int compareByCoord(const void *nodeOne, const void *nodeTwo) {
    double diff = theGlobalCoord[*(int *)nodeOne] - theGlobalCoord[*(int *)nodeTwo];
    if (diff < 0)
        return 1;
    if (diff > 0)
        return -1;
    return 0;
}

void motorMeshRenumber(motorMesh *theMesh, int *number) {
    int *inverse = (int *)malloc(sizeof(int) * theMesh->nNode);
    for (int i = 0; i < theMesh->nNode; i++)
        inverse[i] = i;

    if (1)
        theGlobalCoord = theMesh->X;
    if (0)
        theGlobalCoord = theMesh->Y;

    qsort(inverse, theMesh->nNode, sizeof(int), compareByCoord);
    for (int i = 0; i < theMesh->nNode; i++)
        number[inverse[i]] = i;
    // for (int i = 0 ; i < theMesh->nNode ; i++) number[i] = i; // EASY TRICK FOR NO RENUMBER
    free(inverse);
}

void findBoundary(motorMesh *theMesh) {

    for (int i = 0; i < theMesh->nNode; i++)
        dirichlet[i] = 0;

    motorEdges *theEdges = motorEdgesCreate(theMesh);
    for (int i = 0; i < theEdges->nEdge; i++) {
        if (theEdges->edges[i].elem[1] < 0) {
            dirichlet[theEdges->edges[i].node[0]] = 1;
            dirichlet[theEdges->edges[i].node[1]] = 1;
        }
    }
    motorEdgesFree(theEdges);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// - TORQUE FUNCTIONS - /////////////////////////////////////////////////////

void getInfoTorque(const motorMesh *theMesh, double *a, const int iElem, int *map, double *x, double *y, double *u) {

    int nLoc = theMesh->nLocalNode;
    for (int j = 0; j < nLoc; j++) {
        map[j] = theMesh->elem[iElem * nLoc + j];
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
        u[j] = a[map[j]];
    }
}

double computeIntegral(motorMesh *theMesh, double *a, motor *theMotor) {

    int node, shift, i, n, N, iDomain = 10;
    int iElem, iInteg, map[3];
    for (i = 0, shift = 0; i < iDomain; i++)
        shift += theMesh->nElemDomain[i];
    int nElem = theMesh->nElemDomain[iDomain];
    int nLoc = theMesh->nLocalNode;

    double Xloc[3], Yloc[3], Uloc[3], phi[3], dphidx[3], dphidy[3];
    double xsi, eta, weight, dxdxsi, dxdeta, dydxsi, dydeta, jac, invJ;
    double x, y, r, dadx, dady, dadr, dadt;

    double I = 0.;

    for (iElem = shift; iElem < shift + nElem; iElem++) {
        getInfoTorque(theMesh, a, iElem, map, Xloc, Yloc, Uloc);

        for (iInteg = 0; iInteg < _nInteg; iInteg++) {

            xsi = _gaussTriXsi[iInteg];
            eta = _gaussTriEta[iInteg];
            weight = _gaussTriWeight[iInteg];

            x = 0.;
            y = 0.;
            dxdxsi = 0.;
            dxdeta = 0.;
            dydxsi = 0.;
            dydeta = 0.;
            dadx = 0.;
            dady = 0.;

            phi[0] = 1. - xsi - eta;
            phi[1] = xsi;
            phi[2] = eta;

            for (i = 0; i < nLoc; i++) {
                x += Xloc[i] * phi[i];
                y += Yloc[i] * phi[i];
                dxdxsi += Xloc[i] * _dphidxsi[i];
                dxdeta += Xloc[i] * _dphideta[i];
                dydxsi += Yloc[i] * _dphidxsi[i];
                dydeta += Yloc[i] * _dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            invJ = 1. / jac;

            for (i = 0; i < nLoc; i++) {
                dphidx[i] = (_dphidxsi[i] * dydeta - _dphideta[i] * dydxsi) * invJ;
                dphidy[i] = (_dphideta[i] * dxdxsi - _dphidxsi[i] * dxdeta) * invJ;
                dadx += Uloc[i] * dphidx[i];
                dady += Uloc[i] * dphidy[i];
            }

            r = hypot(x, y);
            dadr = dadx * x / r + dady * y / r;
            dadt = dadx * (-y) + dady * x;
            I += dadr * dadt * jac * weight;
        }
    }
    return I;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// - SAVE FUNCTIONS - //////////////////////////////////////////////////////

void saveMu(motor *theMotor) {
    const double *B = theMotor->hystereticCurveB;
    const double *H = theMotor->hystereticCurveH;
    int n = theMotor->nHystereticCurve;

    double *ddH = getCubicSplines(B, H, n);
    double u, h, b2, mu[3];

    FILE *ptr = fopen("../data/mu.txt", "w");
    double db2 = 1e-8;
    for (int i = 1; i < 1000; i++) {
        b2 = pow(2. * i / 1000., 2);

        getMu(theMotor->hystereticCurveB, theMotor->hystereticCurveH,
              theMotor->nHystereticCurve, ddH, b2, db2, mu);
        fprintf(ptr, "%le %le %le\n", b2, mu[0], mu[1]);
    }
    fclose(ptr);
}

int firstWrite = 0;
double globalAngle = 0.;
//  globalAngle = theMotor->theta;  // in computeMagneticPotential : just needed for python plots
void saveMatrix(COO *A, size_t n) {
    FILE *ptr;
    if (firstWrite == 0) {
        ptr = fopen("cooMatrixZERO.txt", "w");
        firstWrite = 1;
    } else
        ptr = fopen("cooMatrixZERO.txt", "a");

    fprintf(ptr, "%ld %ld %lf\n", n, A->nnz, globalAngle);

    for (size_t i = 0; i < A->nnz; i++) {
        fprintf(ptr, "%ld %ld %le\n", A->data[i].i, A->data[i].j, A->data[i].k);
    }
    fclose(ptr);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// - MAIN FUNCTIONS - //////////////////////////////////////////////////////

void motorAdaptMesh(motor *theMotor, double delta) {
    motorMesh *theMesh = theMotor->mesh;
    double x, y;
    theMotor->theta += delta;
    deltaTheta = delta;

    for (int i = 0; i < theMesh->nNode; ++i) {
        if (theMotor->movingNodes[i] == 1) {
            x = theMesh->X[i] * cos(delta) - theMesh->Y[i] * sin(delta);
            y = theMesh->X[i] * sin(delta) + theMesh->Y[i] * cos(delta);
            theMesh->X[i] = x;
            theMesh->Y[i] = y;
        }
    }

    motorReMesh *theReMesh;
    if (reMeshInit == 0) {
        theReMesh = malloc(sizeof(motorReMesh));
        setupReMesh(theMotor, theReMesh);

        double rIn = 0.;
        double rOut = 0.;
        FILE *ptr = fopen("radiusOut.txt", "w");
        int node;
        for (int i = 0; i < theReMesh->nStatorNodes; i++) {
            node = theReMesh->statorNodes[i];
            rOut += hypot(theMesh->X[node], theMesh->Y[node]);
            fprintf(ptr, "%le\n", hypot(theMesh->X[node], theMesh->Y[node]));
        }
        fclose(ptr);
        ptr = fopen("radiusIn.txt", "w");
        for (int i = 0; i < theReMesh->nRotorNodes; i++) {
            node = theReMesh->rotorNodes[i];
            rIn += hypot(theMesh->X[node], theMesh->Y[node]);
            fprintf(ptr, "%le\n", hypot(theMesh->X[node], theMesh->Y[node]));
        }
        fclose(ptr);
        printf("r out = %lf\n", rOut / theReMesh->nStatorNodes);
        printf("r in  = %lf\n", rIn / theReMesh->nRotorNodes);
        printf("Delta r = %.6le\n", rOut / theReMesh->nStatorNodes - rIn / theReMesh->nRotorNodes);
        printf("Delta r = %.6le\n", 1e-3 / 3.);

        globalReMesh = theReMesh;
        reMeshInit = 1;
    }

    theReMesh = globalReMesh;
    checkAndSwap(theMotor, theReMesh, delta);
}

double motorComputeCouple(motor *theMotor) {
    double d = 0.33e-3; // a bit ugly but it seems to be the exact value
    double L = theMotor->L;
    double mu0 = 4e-7 * M_PI;

    double I = computeIntegral(theMotor->mesh, theMotor->a, theMotor);

    return -L * I / (mu0 * d);
}

void motorComputeCurrent(motor *theMotor) {

    int i;
    double js; // = 8.8464*1e5;
    double phi0 = M_PI / 180. * 18.33;
    double phi = M_PI / 6.;
    double th = fmod(fmod(theMotor->theta + phi0, M_PI / 2.) + M_PI / 2., M_PI / 2.);

    for (i = 1; i < 7; i++) {
        if (theMotor->js[i] > 1) {
            js = theMotor->js[i];
            currentCoil = i;
        }
        theMotor->js[i] = 0.;
    }

    if ((th < phi)) {
        theMotor->js[1] = js;
        theMotor->js[2] = -js;
    } else if ((phi < th) && (th < 2. * phi)) {
        theMotor->js[5] = js;
        theMotor->js[6] = -js;
    } else if ((2 * phi < th) && (th < 3. * phi)) {
        theMotor->js[3] = js;
        theMotor->js[4] = -js;
    }
    return;
}

void motorComputeMagneticPotential(motor *theMotor) {

    int res, flagSolveLinear;

    if (solveInit == 0) {
        dirichlet = (int *)malloc(theMotor->mesh->nNode * sizeof(int));
        number = (int *)malloc(theMotor->mesh->nNode * sizeof(int));
        findBoundary(theMotor->mesh);

        solveInit = 1;
        flagSolveLinear = 1;

        ddH = getCubicSplines(theMotor->hystereticCurveB, theMotor->hystereticCurveH, theMotor->nHystereticCurve);
    } else {
        flagSolveLinear = (theMotor->js[currentCoil] < 1.) || (10. < deltaTheta);
    }

    if (theMotor->nonLinearFlag) {
        theMotor->mu[0] = theMotor->hystereticCurveB[1] / theMotor->hystereticCurveH[1];
        theMotor->mu[8] = theMotor->hystereticCurveB[1] / theMotor->hystereticCurveH[1];
    }

    // Create the solver
    sparseSolver *theSolver = sparseSolverCreate(theMotor->mesh->nNode);

    // renumbering to increase efficency
    motorMeshRenumber(theMotor->mesh, number);

    if (theMotor->nonLinearFlag == 0 || flagSolveLinear) {
        // set first get of PCG as the last solution
        setSystemVectors(theMotor, theSolver);

        // Assemble matrix of linear sytem
        theSolver->A = allocateCOO(3 * theMotor->size);
        assembleMotorMatrix(theMotor, theSolver);

        // Solve the system with PCG
        initCG(theSolver);
        res = PCG(theSolver, 3000, 1e-8);

        // Update solution
        for (int i = 0; i < theMotor->size; i++)
            theMotor->a[i] = theSolver->x[number[i]];
    }

    // if necessary apply Newton's method to solve nonlinear system
    if (theMotor->nonLinearFlag)
        res = newtonRaphson(theMotor, theSolver, 20, 1e-8, 3000, 1e-10, flagSolveLinear);

    // free the solver for the next iteration
    freeSparseSolver(theSolver);

    return;
}

void motorFree(motor *theMotor) {
    free(theMotor->mu);
    free(theMotor->js);
    free(theMotor->a);
    free(theMotor->movingNodes);
    free(theMotor);

    if (globalReMesh)
        freeReMesh(globalReMesh);
    if (number)
        free(number);
    if (dirichlet)
        free(dirichlet);
    if (ddH)
        free(ddH);

    // usefull when we change the mesh during the animation
    solveInit = 0;
    reMeshInit = 0;
}
