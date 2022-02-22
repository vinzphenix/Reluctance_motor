/*
 * 
 * ==========================================================================
 * 
 *  Simulation par éléments finis d'un motor à reluctance variable
 *  Projet du cours LEPL1110 : année 20-21
 *
 *  Vincent Legat
 *  Michel Henry
 *  François Henrotte
 *  Benjamin Legat
 *  
 * ==========================================================================
 */

# define mainAnim main
//# define mainPerf main
//# define mainTest1 main

#include "glfem.h"
#include <time.h>

const char *filename[5] = {"../data/motor400.txt",
                           "../data/motor838.txt",
                           "../data/motor1667.txt",
                           "../data/motor4424.txt",
                           "../data/motor14608.txt"};

int mainAnim(void) {
    
    clock_t tic;
    
    int indexFile = 2;
    double duration = 1e-10 + 3.;
    double dt = 0.003;
    double omega0 = 200. * M_PI/30.;
    double C;
    
    animation *myAnim = animationCreate(filename, 6, 12); //nbMessages, nDomain
    initAnim(myAnim, indexFile, duration, dt, omega0); // indexMesh, duration, dt

    motor *theMotor = myAnim->motor;
    motorMesh *theMotorMesh = theMotor->mesh;

    double js   = 8.8464*1e5 * 1.;
    theMotor->js[1] = js;
    theMotor->js[2] = -js;
    theMotor->nonLinearFlag = 0;

    const char theHelpMessage[] = {
    "   [esc] : Exit\n"
    "\n"
    " <- / -> : Navigate between views \n"
    " up/down : Navigate between modes \n"
    "  space  : Play / Pause animation \n"
    "\n"
    "    C    : Display / hide credits \n"
    "    D    : Display / hide mesh \n"
    "    F    : Display / hide field lines \n"
    "    H    : Display / hide keyboard shortcuts \n"
    "\n"
    "    I    : More / less information \n"
    "    J    : first tap: manual coil selection\n"
    "           other tap: change coil and try make it turn yourself\n"
    "    K    : Show domains with colors \n"
    "    L    : Enable / disable nonlinear solver\n"
    "    M    : Enter / Leave mesh selection\n"
    "    R    : Reset animation, zoom and translations \n"
    "    S    : Show magnetic potential \n"
    };
    printf("\n%s\n", theHelpMessage);
    glfemWindowCreate("EPL1110 : Switched Reluctance Motor",900,700, theMotorMesh->nNode,theMotorMesh->X,theMotorMesh->Y);
    glfemWindowSetHelpMessage(theHelpMessage);
    
    do {
        glfemWindowUpdate();
        char action = glfemGetAction();
        myAnim->action = action;
    
        
        glfemSetColorLine(GLFEM_BLACK);
        glfemSetColor(GLFEM_BACKGROUND);
        glfemSetLineWidth(0.0001);
        
        if (myAnim->meshChange == 1) {
            if ((48 < action) && (action < 54) && (action - '0' - 1 != indexFile)) {
                indexFile = action - '0' - 1;
                motorFree(theMotor);
                motorMeshFree(theMotorMesh);
                initAnim(myAnim, indexFile, duration, dt, omega0);
                theMotor = myAnim->motor;
                theMotorMesh = theMotor->mesh;
                myAnim->meshChange == 0;
                myAnim->flag = 0;
                if (indexFile == 4) {
                    myAnim->funStartTime = glfwGetTime();
                    myAnim->flag = 2;
                }
            }
            if (action != 'C' && 54 <= action && action != 'I'
                && action != 'd' && action != 'u') {
                myAnim->flag = 1;
                myAnim->funStartTime = glfwGetTime();
            }
        }
        
        //   Calcul du potentiel magnétique A
        //   Calcul du couple
        //   Calcul de omega par l'équation de Newton
        //   Rotation du rotor et remaillage de la bande glissante
        //   Mise a jour des courants dans les inducteurs en fonction de l'angle
        
        updateParameters(myAnim);
        theMotor->time = myAnim->time;
        updateDisplay(myAnim);
        
        if (myAnim->discreteTime < myAnim->time) {
            
            tic = clock();
            motorComputeMagneticPotential(theMotor);
            //printf("\tTime for solve : %.5f [s]\n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            
            C = motorComputeCouple(theMotor);
            theMotor->omega += C * myAnim->timeStep / theMotor->inertia;
            
            motorAdaptMesh(theMotor, theMotor->omega * myAnim->timeStep);
            
            if (myAnim->manCurrent == 0)
                motorComputeCurrent(theMotor);
            
            myAnim->C = C;
            myAnim->iteration += 1;
            myAnim->discreteTime += myAnim->timeStep;
            myAnim->fullTime += myAnim->timeStep;
            //printf("Iteration = %4d \t t = %.2f \t\t \u03B8 = %.2f deg \t C = %lf\n\n",myAnim->iteration, myAnim->fullTime, fmod(180./M_PI * theMotor->theta, 360.), C);
        }
        
        updateMessages(myAnim);

        
    } while(!glfemWindowShouldClose());
    
    freeAnimation(myAnim);
    motorFree(theMotor);
    motorMeshFree(theMotorMesh);
    
    return 0;
}


// =========== Quelques fonctions fournies gracieusement par l'équipe didactique ====== //
//
// Attention : il n'est pas permis de les modifier :-)
// Attention : la structure de données du problème est figée et vous ne pouvez pas
//             la modifier : c'est cette structure que le correcteur automatique
//             utilisera pour tester votre programme
//
// =====================================================================================//

static const double _hystereticCurveH[43] = { 0, 10, 20, 30, 40, 50, 
      60, 70, 80, 90, 100, 125, 150, 175, 200, 250,
      300, 400, 500, 600,  700, 800, 900, 1000, 1250, 1500, 2000, 2500, 5000,
      7500,  10000, 15000, 20000, 59000, 174000, 514000, 1520000, 4470000,
      13200000, 38900000, 115000000, 339000000, 1000000000 }; 
static const double _hystereticCurveB[43] =  { 0.0,           
      0.194880829963, 0.377143018857, 0.537767739762, 0.672888260835,
      0.783043000477, 0.871342430831,0.941778611986, 0.998183303557, 1.04378111223,
      1.08110469369, 1.14963767549, 1.19607212343, 1.22964695907, 1.25515221835,
      1.29162498935, 1.31678879432, 1.35015120537, 1.37220092877, 1.38859114656,
      1.4017440574, 1.41287024565, 1.42264180514, 1.43146158921, 1.45082466146,
      1.46784549989, 1.49819370601, 1.52578650709, 1.64314027719, 1.73458485332,
      1.8039068939,1.89568786291, 1.95213815187, 2.1390774927, 2.45827909293,
      3.32303272825, 5.85485500678, 13.2701832298, 35.2114648741, 99.8027446541,
      291.062951228, 854.036370229, 2515.3105707 };   

  

motor *motorCreate(motorMesh *theMesh) {
    motor *theMotor = malloc(sizeof(motor));
    theMotor->mesh = theMesh;
    theMotor->size = theMesh->nNode;
    
    //
    //  Identification des noeuds mobiles 
    //
    
    theMotor->movingNodes = malloc(sizeof(int)*theMesh->nNode);
    for (int i=0; i < theMotor->size; i++) 
        theMotor->movingNodes[i] = 0;
    for (int i=0; i < theMesh->nElem; i++) {
        int domain = theMesh->domain[i];
        if (domain == 8 || domain == 9 || domain == 10 ) {
            int *elem = &(theMesh->elem[i*3]);
            for (int j=0; j < 3; j++) {
                theMotor->movingNodes[elem[j]] = 1;
            }
        }
    }
    //
    //  Initialisation des inconnues à la valeur X pour voir le maillage tourner
    //   
    theMotor->a = malloc(sizeof(double)*theMotor->size);
    for (int i = 0 ; i < theMotor->size ; i++){
            double th = atan2(theMotor->mesh->Y[i], theMotor->mesh->X[i]);
            double r = hypot(theMotor->mesh->Y[i], theMotor->mesh->X[i]);
            if (theMotor->movingNodes[i]) theMotor->a[i] = sin(3*th);
            else theMotor->a[i] = -sin(3*th);
        }    
    theMotor->theta = 0;
    theMotor->omega = 0;
    theMotor->time  = 0;
    
    //
    //  Parametres materiels
    //    
    
    double mu_0 = 4*3.141592653589793*1e-7; //kg m /(A**2 s**2)
    double mu_r = 1e3;                      
    double js   = 8.8464*1e5;               // A / m**2
    theMotor->inertia = 5*1e-4;             // kg m**2
    theMotor->L       = 0.06;               // m
    theMotor->js = malloc(sizeof(double)*theMesh->nDomain);
    theMotor->mu = malloc(sizeof(double)*theMesh->nDomain);  
    for(int i = 0; i < theMesh->nDomain; i++) {
        theMotor->js[i] = 0.0;
        theMotor->mu[i] = mu_0; }
    theMotor->mu[0] *= mu_r;
    theMotor->mu[8] *= mu_r;
    theMotor->js[1] = js;
    theMotor->js[2] = -js;
    theMotor->js[3] = 0.;
    theMotor->js[4] = 0.;
    theMotor->js[5] = 0.;
    theMotor->js[6] = 0.;
    
    //
    //  Bonus : magnetostatique non-lineaire
    //
    
    theMotor->nonLinearFlag = 0;
    theMotor->nHystereticCurve = 43;
    theMotor->hystereticCurveH = _hystereticCurveH;
    theMotor->hystereticCurveB = _hystereticCurveB;

    return theMotor;
}

void motorPrintInfos(const motor *theMotor) {
    int  size = theMotor->size;
    motorMesh *theMesh = theMotor->mesh;
    printf(" \n");
    printf(" ====== Switched Reluctance Motor Simulation ============================\n");
    for(int i = 0; i < theMesh->nDomain; i++) 
        printf("    Domain %2d : %-16s : nElem = %6d, mu = %10.3e, js = %10.3e \n",i,theMesh->nameDomain[i],
                    theMesh->nElemDomain[i],theMotor->mu[i],theMotor->js[i]);
    printf("                                 : mu permeability [kg m s2/A2] - js current density [A/m2]\n");
    printf("    Number of elements           : %d\n",theMesh->nElem);   
    printf("    Number of nodes              : %d\n",theMotor->size);  
    printf("    Flag for non linearities     : %d\n",theMotor->nonLinearFlag); 
    printf("    Rotor inertia                : %13.7e [kg m2]\n",theMotor->inertia); 
    printf("    Motor axial length           : %13.7e [m]\n",theMotor->L); 
    printf("    Time                         : %13.7e [s]\n",theMotor->time); 
    printf("    Angular position             : %13.7e [rad]\n",theMotor->theta); 
    printf("    Angular velocity             : %13.7e [rad/s]\n",theMotor->omega); 
    printf("=========================================================================\n");
}


motorMesh *motorMeshRead(const char *filename) {
    motorMesh *theMesh = malloc(sizeof(motorMesh));
    theMesh->nLocalNode = 3;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");   
    int i,j,trash,*elem;
     
    ErrorScan(fscanf(file, "Number of nodes %d\n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);     
    for (i = 0; i < theMesh->nNode; i++) 
         ErrorScan(fscanf(file,"%d : %le %le\n",&trash,&theMesh->X[i],&theMesh->Y[i])); 
    
    ErrorScan(fscanf(file, "Number of sub-domains %d \n", &theMesh->nDomain)); 
    theMesh->nameDomain = malloc(sizeof(string)*theMesh->nDomain);
    theMesh->nElemDomain = malloc(sizeof(int)*theMesh->nDomain);
    for(i = 0; i < theMesh->nDomain; i++) 
     	  ErrorScan(fscanf(file, "%6d : %s : %d\n",&trash,theMesh->nameDomain[i],&theMesh->nElemDomain[i]));
    	   
    ErrorScan(fscanf(file, "Number of triangles %d  \n", &theMesh->nElem));   
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    theMesh->domain = malloc(sizeof(int)*theMesh->nElem);
    for (i = 0; i < theMesh->nElem; i++) {
    	  elem = &(theMesh->elem[i*3]);
    	  ErrorScan(fscanf(file,"%d : %d %d %d %d\n",&trash,&elem[0],&elem[1],&elem[2],&theMesh->domain[i]));
    }
    
    fclose(file);
    return theMesh;
}


void motorMeshWrite(const motorMesh *theMesh, const char *filename) {
    int i,j,*elem;
    
    FILE* file = fopen(filename,"w");    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; i++) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]);
    }
    
    fprintf(file, "Number of sub-domains %d \n", theMesh->nDomain);  
    for (i = 0; i < theMesh->nDomain; i++) { 
     	  fprintf(file, "%6d : %-16s : %d \n",i,theMesh->nameDomain[i],theMesh->nElemDomain[i]);
    }
    
    fprintf(file, "Number of triangles %d \n", theMesh->nElem);
    for (i = 0; i < theMesh->nElem; i++) {
        elem = &(theMesh->elem[i*3]);
        fprintf(file,"%6d : %6d %6d %6d %6d \n",i,elem[0],elem[1],elem[2],theMesh->domain[i]);
    }
    
    fclose(file);
}

void motorMeshFree(motorMesh *theMesh) {
    free(theMesh->elem);
    free(theMesh->nElemDomain);
    free(theMesh->nameDomain);
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->domain);
    free(theMesh);
}

femMesh *motorDomainCreate(const motorMesh *theMotorMesh, int iDomain) {
    femMesh *theMesh = (femMesh *) malloc(sizeof(femMesh)); 
    theMesh->nLocalNode = 3;
    theMesh->nNode = theMotorMesh->nNode; 
    theMesh->X = theMotorMesh->X;
    theMesh->Y = theMotorMesh->Y;
    
    int shift = 0;
    for (int i=0; i < iDomain; i++)
        shift += theMotorMesh->nElemDomain[i];
    theMesh->elem = &theMotorMesh->elem[3*shift];
    theMesh->nElem = theMotorMesh->nElemDomain[iDomain];
    return theMesh;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////// - ANIMATION FUNCTIONS - ///////////////////////////////////////////////////

animation *animationCreate(const char **filename, int nbMessages, int nDomain) {
    animation *myAnim = malloc(sizeof(animation));
    
    char **theMessages = malloc(nbMessages * sizeof(char*));
    for (int i = 0 ; i < nbMessages ; i++)
        theMessages[i] = malloc(256 * sizeof(char));
    
    femMesh **thedomain = (femMesh **) malloc(nDomain * sizeof(femMesh *));
    
    myAnim->filename = filename;
    myAnim->nbMessages = nbMessages;
    myAnim->messages = theMessages;
    myAnim->domain = thedomain;
    
    sprintf(myAnim->easterEgg, "000000");
    myAnim->completed = 0;
    myAnim->reallyCompleted = 0;
    
    return myAnim;
}

void freeAnimation(animation *myAnim) {
    for (int i = 0 ; i < myAnim->nbMessages ; i++)
        free(myAnim->messages[i]);
    for (int i = 0 ; i < myAnim->motor->mesh->nDomain ; i++)
        free(myAnim->domain[i]);
        //myAnim->domain[i] = NULL;
    free(myAnim->messages);
    free(myAnim->domain);
    free(myAnim);
}

void initAnim(animation *myAnim, int indexMesh, double duration, double dt, double omega0) {
    motorMesh *theMotorMesh = motorMeshRead(myAnim->filename[indexMesh]);
    motor *theMotor = motorCreate(theMotorMesh);
    theMotor->omega = omega0;
    //motorPrintInfos(theMotor);
    
    for (int i = 0 ; i < theMotorMesh->nDomain ; i++)
        myAnim->domain[i] = motorDomainCreate(theMotorMesh, i);
    
    myAnim->motor = theMotor;
    myAnim->mode = 1;
    myAnim->subMode = 0;
    myAnim->nbSubModes = 2;
    myAnim->informations = 1;
    myAnim->meshChange = 0;
    myAnim->iteration = 0;
    myAnim->displayMesh = 0;
    myAnim->manCurrent = 0;
    myAnim->discreteTime = 0.;
    myAnim->startTime = 0.;
    myAnim->duration = duration - dt;
    myAnim->fullTime = 0.;
    myAnim->timeStep = dt;
    myAnim->stopTime = 0.;
    myAnim->fastFactor = 1.;
    myAnim->fieldLines = 1;
    myAnim->C = 0.;
    theMotor->omega=omega0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////// - ANALYSIS FUNCTIONS - ////////////////////////////////////////////////////

int mainPerf(int argc, char *argv[]) {
    int fileNb = 4;
    int iterNb = 1;
    
    clock_t tic, tac;
    //FILE *ptr = fopen("renumberZERO.txt", "a");
    //FILE *ptr = fopen("renumberONCE.txt", "a");
    //FILE *ptr = fopen("renumberFULL.txt", "a");
    //FILE *ptr = fopen("../data/nonLinear.txt", "a");
    

    motorMesh *theMotorMesh = motorMeshRead(filename[fileNb]);
    motor *theMotor = motorCreate(theMotorMesh);
    
    double js   = 8.8464 * 1.e5 * 1.;
    theMotor->js[1] = js;
    theMotor->js[2] = -js;
    int bobine;
    
    theMotor->nonLinearFlag = 0;
    
    for (int i = 0 ; i < iterNb ; i++) {
        
        for (int j = 1 ; j < 7 ; j++)
            if (theMotor->js[j] > 1.)
                bobine = j;
        
        //tic = clock();
        motorComputeMagneticPotential(theMotor);
        //tac = clock();
        //fprintf(ptr, "%d %le %le\n", theMotor->size, theMotor->theta, (tac - tic) / (1. * CLOCKS_PER_SEC));
        //printf("Time for SOLVE : %.3f [s]\n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
        
        //tic = clock();
        double C = motorComputeCouple(theMotor);
        //fprintf(ptr, "%d %d %d %le %le %le\n", theMotor->size, theMotor->nonLinearFlag, bobine, js, theMotor->theta, C);
        printf("angle = %5.2lf \t C = %lf\n", 180./M_PI * theMotor->theta, C);
        //printf("Time for COUPLE : %.3f [s]\t C = %le\n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC, C);
        
        //tic = clock();
        motorAdaptMesh(theMotor, -0.5*M_PI/180.);
        //printf("Time for ADAPT : %.3f [s]\n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
        
        motorComputeCurrent(theMotor);
        
        printf("---------------------------------------------------------------------\n");
    }
    
    //fclose(ptr);
    motorFree(theMotor);
    motorMeshFree(theMotorMesh);
    
    return 0;
}


int mainTest1(void) {

    int nb = 3;
    int coil = 3;
    int flag = 1;

    motorMesh *theMotorMesh = motorMeshRead(filename[nb]);
    motor *theMotor = motorCreate(theMotorMesh);
    theMotor->nonLinearFlag = flag;
    
    double C, js   = 8.8464*1e5;
    
    for (int i = 1 ; i < 7 ; i++) theMotor->js[i] = 0.;
    if (coil == 1) {
        theMotor->js[1] = js;
        theMotor->js[2] = -js;
    } else if (coil == 2) {
        theMotor->js[3] = js;
        theMotor->js[4] = -js;
    } else {
        theMotor->js[5] = js;
        theMotor->js[6] = -js;
    }
    
    FILE *ptr = fopen("../data/torque.txt", "a");
    
    for (int i = 0 ; i < 0 ; i++) {
        motorComputeMagneticPotential(theMotor);
        C = motorComputeCouple(theMotor);
        
        fprintf(ptr, "%d %d %le %le\n", coil, flag, -180./M_PI * theMotor->theta, C);
        printf("iter %3d \t angle %lf \t C = %lf\n", i, -180./M_PI * theMotor->theta, C);
        
        motorAdaptMesh(theMotor, -0.5 * M_PI/180.);
    }
    fclose(ptr);
    
    return 0;
}

