#include "motor.h"
#include <math.h>

int getFromCoord(motor *theMotor, double x, double y, double *sol) {
    motorMesh *theMesh = theMotor->mesh;
    double den, xsi, eta, xLoc[4], yLoc[4], uLoc[4];
    xLoc[0] = x;
    yLoc[0] = y;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {

        xLoc[1] = theMesh->X[theMesh->elem[3 * iElem + 0]];
        yLoc[1] = theMesh->Y[theMesh->elem[3 * iElem + 0]];
        xLoc[2] = theMesh->X[theMesh->elem[3 * iElem + 1]];
        yLoc[2] = theMesh->Y[theMesh->elem[3 * iElem + 1]];
        xLoc[3] = theMesh->X[theMesh->elem[3 * iElem + 2]];
        yLoc[3] = theMesh->Y[theMesh->elem[3 * iElem + 2]];

        den = (xLoc[2] - xLoc[1]) * (yLoc[3] - yLoc[1]) - (yLoc[2] - yLoc[1]) * (xLoc[3] - xLoc[1]);
        xsi = (yLoc[3] - yLoc[1]) * (x - xLoc[1]) + (xLoc[1] - xLoc[3]) * (y - yLoc[1]);
        eta = (yLoc[1] - yLoc[2]) * (x - xLoc[1]) + (xLoc[2] - xLoc[1]) * (y - yLoc[1]);

        if ((xsi >= 0.) && (eta >= 0.) && (xsi + eta <= den)) {
            xsi /= den;
            eta /= den;
            uLoc[1] = theMotor->a[theMesh->elem[3 * iElem + 0]];
            uLoc[2] = theMotor->a[theMesh->elem[3 * iElem + 1]];
            uLoc[3] = theMotor->a[theMesh->elem[3 * iElem + 2]];

            sol[0] = (1 - xsi - eta) * uLoc[1] + xsi * uLoc[2] + eta * uLoc[3];
            // sol[1] = (-uLoc[1] + uLoc[2]) * (yLoc[3] - yLoc[1]) / den + (-uLoc[1] + uLoc[3]) * (yLoc[1] - yLoc[2]) / den;
            // sol[2] = (-uLoc[1] + uLoc[2]) * (xLoc[1] - xLoc[3]) / den + (-uLoc[1] + uLoc[3]) * (xLoc[2] - xLoc[1]) / den;
            return 0;
            // return (1 - xsi - eta) * uLoc[1] + xsi * uLoc[2] + eta * uLoc[3];
        }
    }
    printf("we have an error houston : x = %lf \t y = %lf\n", x, y);
    return -1;
}

double *getFieldNewMesh(motor *theMotor, motor *oldMotor) {
    int res;
    double *aOld = (double *)malloc(theMotor->size * sizeof(double));

    double *X_new = theMotor->mesh->X;
    double *Y_new = theMotor->mesh->Y;

    for (int iNode = 0; iNode < theMotor->size; iNode++) {
        res = getFromCoord(oldMotor, X_new[iNode], Y_new[iNode], &aOld[iNode]);
        if (res == -1) {
            aOld[iNode] = theMotor->a[iNode];
        }
    }
    return aOld;
}

void computeError(motor *theMotor, double *aOLD, double *values) {

    motorMesh *theMesh = theMotor->mesh;

    double xLoc[3], yLoc[3], uLoc[3], uLocOLD[3], phi[3];
    double x, y, u, uOLD, dudx, dudy, dudxOLD, dudyOLD;
    double xsi, eta, w, J, invJ, dxdxsi, dxdeta, dydxsi, dydeta;
    int iElem, iInteg, node, k, nLoc = 3;

    static const double _dphidxsi[3] = {-1.0, 1.0, 0.0};
    static const double _dphideta[3] = {-1.0, 0.0, 1.0};
    static const double _gaussTriXsi[3] = {0.166666666666667, 0.666666666666667, 0.166666666666667};
    static const double _gaussTriEta[3] = {0.166666666666667, 0.166666666666667, 0.666666666666667};
    static const double _gaussTriWeight[3] = {0.166666666666667, 0.166666666666667, 0.166666666666667};
    static const int _nInteg = 3;

    double errorSoluceL2 = 0.0;
    double errorSoluceH1 = 0.0;

    double lc = 0.;
    double minlc = 1.;
    double maxlc = 0.;
    double sum = 0.;
    double sum2 = 0.;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {

        for (k = 0; k < nLoc; k++) {
            xLoc[k] = theMesh->X[theMesh->elem[3 * iElem + k]];
            yLoc[k] = theMesh->Y[theMesh->elem[3 * iElem + k]];
            uLoc[k] = theMotor->a[theMesh->elem[3 * iElem + k]];
            uLocOLD[k] = aOLD[theMesh->elem[3 * iElem + k]];
        }

        for (iInteg = 0; iInteg < _nInteg; iInteg++) {
            xsi = _gaussTriXsi[iInteg];
            eta = _gaussTriEta[iInteg];
            w = _gaussTriWeight[iInteg];

            x = 0.;
            y = 0.;
            u = 0.;
            uOLD = 0.;
            dudx = 0.;
            dudy = 0.;
            dudxOLD = 0.;
            dudyOLD = 0.;
            dxdxsi = 0.;
            dxdeta = 0.;
            dydxsi = 0.;
            dydeta = 0.;

            phi[0] = 1. - xsi - eta;
            phi[1] = xsi;
            phi[2] = eta;

            for (k = 0; k < nLoc; k++) {
                dxdxsi += xLoc[k] * _dphidxsi[k];
                dxdeta += xLoc[k] * _dphideta[k];
                dydxsi += yLoc[k] * _dphidxsi[k];
                dydeta += yLoc[k] * _dphideta[k];

                x += xLoc[k] * phi[k];
                y += yLoc[k] * phi[k];
                u += uLoc[k] * phi[k];
                uOLD += uLocOLD[k] * phi[k];
            }

            J = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            invJ = 1. / J;

            for (k = 0; k < nLoc; k++) {
                dudx += uLoc[k] * (dydeta * _dphidxsi[k] - dydxsi * _dphideta[k]) * invJ;
                dudy += uLoc[k] * (dxdxsi * _dphideta[k] - dxdeta * _dphidxsi[k]) * invJ;
                dudxOLD += uLocOLD[k] * (dydeta * _dphidxsi[k] - dydxsi * _dphideta[k]) * invJ;
                dudyOLD += uLocOLD[k] * (dxdxsi * _dphideta[k] - dxdeta * _dphidxsi[k]) * invJ;
            }

            errorSoluceL2 += (u - uOLD) * (u - uOLD) * w * J;
            errorSoluceH1 += ((dudx - dudxOLD) * (dudx - dudxOLD) + (dudy - dudyOLD) * (dudy - dudyOLD)) * w * J;
        }
        lc = sqrt(J / 2.);

        minlc = fmin(lc, minlc);
        maxlc = fmax(lc, maxlc);
        sum += lc;
        sum2 += lc * lc;
    }

    errorSoluceH1 = sqrt(errorSoluceL2 + errorSoluceH1);
    errorSoluceL2 = sqrt(errorSoluceL2);
    printf("\tH1 norm = %le\n", errorSoluceH1);
    printf("\tL2 norm = %le\n", errorSoluceL2);
    printf("min = %lf\t mean = %lf\t rms = %lf\n", minlc, sum / theMesh->nElem, sqrt(sum2 / theMesh->nElem));

    values[0] = theMesh->nElem;
    values[1] = theMesh->nNode;
    values[2] = errorSoluceL2;
    values[3] = errorSoluceH1;
    values[4] = minlc;
    values[5] = sum / theMesh->nElem;
    values[6] = sqrt(sum2 / theMesh->nElem);
    values[7] = maxlc;
}

int mainConvergence(int argc, char *argv[]) {

    char *filename[5];
    filename[0] = "../data/motor400.txt";
    filename[1] = "../data/motor838.txt";
    filename[2] = "../data/motor1667.txt";
    filename[3] = "../data/motor4424.txt";
    filename[4] = "../data/motor14608.txt";

    double values[8];

    FILE *ptr = fopen("convergence.txt", "a");

    // for (int i = 0 ; i < 5 ; i++) {

    printf("%d \t %d\n", atoi(argv[1]), atoi(argv[2]));

    motorMesh *motorMesh1 = motorMeshRead(filename[atoi(argv[1])]);
    motor *motor1 = motorCreate(motorMesh1);
    motorComputeMagneticPotential(motor1);

    motorMesh *motorMesh2 = motorMeshRead(filename[atoi(argv[2])]);
    motor *motor2 = motorCreate(motorMesh2);
    motorComputeMagneticPotential(motor2);

    // double *aOLD = getFieldNewMesh(motor2, motor1);
    double *aOLD = (double *)calloc(motor2->size, sizeof(double));
    computeError(motor2, aOLD, values);

    for (int i = 0; i < 8; i++)
        fprintf(ptr, "%le ", values[i]);
    fprintf(ptr, "\n");

    motorFree(motor1);
    motorFree(motor2);

    // printf("FILE N %d : %d\t %lf\n", i, theMotorMesh->nNode, 0.10/sqrt(theMotorMesh->nNode));
    // }

    fclose(ptr);

    return 0;
}
