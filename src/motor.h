#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define Error(a) femError(a, __LINE__, __FILE__)
#define Warning(a) femWarning(a, __LINE__, __FILE__)
#define FALSE 0
#define TRUE 1

#define MAXSTRING 50
typedef char string[MAXSTRING];

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
    int *number;
} femMesh;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
    int nDomain;
    int *nElemDomain;
    string *nameDomain;
    int *domain;
} motorMesh;

typedef struct {
    int size;
    motorMesh *mesh;
    double *a;
    double time;
    double theta;
    double omega;
    int *movingNodes;
    double inertia;
    double L;
    double *js;
    double *mu;
    int nonLinearFlag;
    int nHystereticCurve;
    const double *hystereticCurveH;
    const double *hystereticCurveB;
} motor;

typedef struct {
    motor *motor;
    femMesh **domain;
    const char **filename;
    char **messages;
    char easterEgg[6];
    int completed;
    int reallyCompleted;
    int nbMessages;
    char action;
    int mode;
    int subMode;
    int nbSubModes;
    int informations;
    int fieldLines;
    int meshChange;
    int iteration;
    int running;
    int flag;
    int manCurrent;
    int displayMesh;
    double discreteTime;
    double duration;
    double startTime;
    double time;
    double stopTime;
    double fullTime;
    double timeStep;
    double fastFactor;
    double funStartTime;
    double delayTime;
    double C;
} animation;

typedef struct {
    int nStatorNodes; // number moving nodes of the Gap
    int *statorNodes; // moving nodes of the Gap

    int nRotorNodes; // number moving nodes of the Gap
    int *rotorNodes; // moving nodes of the Gap

    int nGapElem;
    int *gapElem; // pointer to the elements on the Gap
} motorReMesh;

typedef struct {
    int elem[2];
    int node[2];
} motorEdge;

typedef struct {
    motorMesh *mesh;
    motorEdge *edges; // List of edges
    int nEdge;        // Number of edges
    int nBoundary;    // Number of edges on the boundary
} motorEdges;

typedef struct ijk { // Contains one triplet (index, index, value)
    size_t i, j;     // Indices
    double k;        // Value
} ijk;

typedef struct COO { // Sparse matrix "coordinate list"
    size_t nnz;      // Number of non zeros in the matrix
    size_t alloc;    // Number of entries that have been allocated
    int compressed;  // Status : sorted/compressed or not
    ijk *data;       // Content of the matrix [alloc]
} COO;

typedef struct CSR { // Sparse matrix "compressed sparse row"
    size_t nnz;      // Number of non zeros in the matrix
    size_t size;     // Number of rows/cols of the matrix
    size_t *DIAG;    // index of DIAG elements [size]
    size_t *IA;      // row index [size + 1]
    size_t *JA;      // col index [nnz]
    double *VA;      // values [nnz]
} CSR;

typedef struct {
    size_t size;   // Size of the system = number of nodes
    COO *A;        // Matrix of the linear system in COO format
    CSR *B;        // Matrix of the linear system in CSR format
    CSR *M;        // Matrix obtained after ILU0
    double *x;     // Solution of the system [size]
    double *b;     // RHS vector [size]
    double *R;     // Residue [size]
    double *Z;     // Residue Tilde [size]
    double *D;     // Direction [size]
    double *S;     // Result of Matrix-Vector product [size]
    double rr;     // R dot R
    double norm_b; // Norm 2 of b
    // int *number;        // Permutation list for minimal bandwidth [size]
    // int *dirichlet;     // Bitmap of nodes on the boundary [size]
} sparseSolver;

animation *animationCreate(const char **filename, int nbMessages, int nDomain);
void initAnim(animation *myAnim, int indexMesh, double duration, double dt, double omega0);
void freeAnimation(animation *myAnim);

void motorMeshWrite(const motorMesh *theMotorMesh, const char *filename);
motorMesh *motorMeshRead(const char *filename);
void motorMeshFree(motorMesh *myMotorMesh);
femMesh *motorDomainCreate(const motorMesh *theMotorMesh, int iDomain);

motor *motorCreate(motorMesh *theMotorMesh);
void motorPrintInfos(const motor *theMotor);
void motorComputeMagneticPotential(motor *myMotor);
double motorComputeCouple(motor *myMotor);
void motorAdaptMesh(motor *myMotor, double delta);
void motorComputeCurrent(motor *myMotor);
void motorFree(motor *theMotor);

int newtonRaphson(motor *theMotor, sparseSolver *mySolver, size_t iterMaxNewton, double precisionNewton,
                  size_t iterMaxPCG, double precisionPCG, int flagSolveLinear);
void assembleMotorMatrix(motor *theMotor, sparseSolver *theSolver);
int PCG(sparseSolver *mySolver, size_t itermax, double precision);

motorEdges *motorEdgesCreate(motorMesh *theMesh);
void motorEdgesFree(motorEdges *theEdges);
int motorEdgesCompare(const void *edgeOne, const void *edgeTwo);

void findNodes(motor *theMotor, motorReMesh *theReMesh);
void checkAndSwap(motor *theMotor, motorReMesh *theReMesh, double delta);
void freeReMesh(motorReMesh *theReMesh);

double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femWarning(char *text, int line, char *file);
