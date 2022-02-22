/*
 *  glfem.c - BOV version
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */
 
 
#include "glfem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))

double   zoom_init;
double   translate_init[2];
double   w_init;
double   h_init;
double   min_colormap;
double   max_colormap;
float    current_text_color[4] = {1.0,0.0,0.0,1.0};
float    current_color[4]      = {0.0,0.0,0.0,1.0};
float    current_line_color[4] = {0.0,0.0,0.0,1.0};
float 	 current_line_width = 0.005;

typedef bov_window_t glfemWindow;
glfemWindow* theCurrentWindow = NULL;
float w;
float h;

float LGRED_TEXT[4]   = {1.0,0.0,0.0,0.65};
float RED_TEXT[4]   = {1.0,0.0,0.0,1.0};
float BLACK_TEXT[4] = {0.0,0.0,0.0,1.0};
float TRSP_TEXT[4] = {1.0,1.0,1.0,0.0};
float SHADE_TEXT[4] = {.85,.85,.85,.8};
int PIXEL = 2;

void glfemSetColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_color[i] = color[i]; }
}

void glfemSetTextColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_text_color[i] = color[i]; }
}

void glfemSetLineWidth(float width) 
{
    current_line_width = width;
}

void glfemSetColorLine(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_line_color[i] = color[i]; }
}

int glfemGetKey(char theKey)
{
    return (glfwGetKey(theCurrentWindow->self,theKey) == GLFW_PRESS) ;
}


static char theCurrentAction = '0';

char glfemGetAction(void) {
    if (theCurrentAction != '0')  {
          char theLastCurrentAction = theCurrentAction;
          theCurrentAction = '0';
          return theLastCurrentAction;}
    else  return theCurrentAction;
}

//
// Attention : GLFW ne detecte pas le type du clavier et donc : il ne faut utiliser
// que les lettres communes de AZERTY et QWERTY :-)
//
// Il n’y a pas tant de lettres communes (20 lettres). 
// Dans l’ordre d’apparition : e, r, t, u, i, o, p, s, d, f, g, h, j, k, l, x, c, v, b, n
//

static void glfemKeyCallback(GLFWwindow* self,
            int key, int scancode, int action,int mods) 
{
    bov_window_t* window = (bov_window_t*) glfwGetWindowUserPointer(self);
    if (action==GLFW_PRESS || action==GLFW_REPEAT) {
        switch (key) {
          case GLFW_KEY_ESCAPE :
            glfwSetWindowShouldClose(self,GL_TRUE);
            break;
          case GLFW_KEY_H :
            if(window->help_needed==0) window->help_needed = 1;
            else                       window->help_needed = 0;
            break;
          case GLFW_KEY_R :
            window->param.zoom = zoom_init;
            window->param.translate[0] = translate_init[0];
            window->param.translate[1] = translate_init[1];
            break;}}
    
    if      (key == GLFW_KEY_C && action == GLFW_PRESS)  {theCurrentAction = 'C';}
    else if (key == GLFW_KEY_E && action == GLFW_PRESS)  {theCurrentAction = 'E';}
    else if (key == GLFW_KEY_D && action == GLFW_PRESS)  {theCurrentAction = 'D';}
    else if (key == GLFW_KEY_F && action == GLFW_PRESS)  {theCurrentAction = 'F';}
    else if (key == GLFW_KEY_H && action == GLFW_PRESS)  {theCurrentAction = 'H';}
    else if (key == GLFW_KEY_I && action == GLFW_PRESS)  {theCurrentAction = 'I';}
    else if (key == GLFW_KEY_J && action == GLFW_PRESS)  {theCurrentAction = 'J';}
    else if (key == GLFW_KEY_K && action == GLFW_PRESS)  {theCurrentAction = 'K';}
    else if (key == GLFW_KEY_L && action == GLFW_PRESS)  {theCurrentAction = 'L';}
    else if (key == GLFW_KEY_M && action == GLFW_PRESS)  {theCurrentAction = 'M';}
    else if (key == GLFW_KEY_R && action == GLFW_PRESS)  {theCurrentAction = 'R';}
    else if (key == GLFW_KEY_S && action == GLFW_PRESS)  {theCurrentAction = 'S';}

    else if (key == GLFW_KEY_SLASH && action == GLFW_PRESS)  {theCurrentAction = 'A';}
    else if (key == GLFW_KEY_EQUAL && action == GLFW_PRESS)  {theCurrentAction = 'Z';}
    
    else if (key == GLFW_KEY_1 && action == GLFW_PRESS)  {theCurrentAction = '1';}
    else if (key == GLFW_KEY_2 && action == GLFW_PRESS)  {theCurrentAction = '2';}
    else if (key == GLFW_KEY_3 && action == GLFW_PRESS)  {theCurrentAction = '3';}
    else if (key == GLFW_KEY_4 && action == GLFW_PRESS)  {theCurrentAction = '4';}
    else if (key == GLFW_KEY_5 && action == GLFW_PRESS)  {theCurrentAction = '5';}
    else if (key == GLFW_KEY_6 && action == GLFW_PRESS)  {theCurrentAction = '6';}
    else if (key == GLFW_KEY_7 && action == GLFW_PRESS)  {theCurrentAction = '7';}
    else if (key == GLFW_KEY_8 && action == GLFW_PRESS)  {theCurrentAction = '8';}
    else if (key == GLFW_KEY_9 && action == GLFW_PRESS)  {theCurrentAction = '9';}
    
    else if (key == GLFW_KEY_SEMICOLON && action == GLFW_PRESS)  {theCurrentAction = 'M';}
    else if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)  {theCurrentAction = 'P';}
    else if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS)  {theCurrentAction = 'r';}
    else if (key == GLFW_KEY_LEFT && action == GLFW_PRESS)  {theCurrentAction = 'l';}
    else if (key == GLFW_KEY_UP && action == GLFW_PRESS)  {theCurrentAction = 'u';}
    else if (key == GLFW_KEY_DOWN && action == GLFW_PRESS)  {theCurrentAction = 'd';}
    else if (key==GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(self,GL_TRUE);
    
}

void glfemWindowCreate(const char *windowName,int w, int h,int n,double *x, double *y) {
    bov_window_t *window = bov_window_new(w,h, windowName);
    theCurrentWindow = window;
    bov_window_set_color(window, (GLfloat[4]){0.9, 0.9, 0.8, 0.0});
    
    //
    // Defining the current viewport and stores it as the reference
    //
  
    /*double minX  = femMin(x,n)* -0.15;
    double maxX  = femMax(x,n)*0.45;
    double minY  = femMin(y,n)*-0.3;
    double maxY  = femMax(y,n)*0.4;*/
    double minX  = femMin(x,n);
    double maxX  = femMax(x,n);
    double minY  = femMin(y,n);
    double maxY  = femMax(y,n);
    double sizeX = (maxX-minX)/1.45;
    double meanX = (maxX+minX)/2.0; 
    double sizeY = (maxY-minY)/1.45;
    double meanY = (maxY+minY)/2.0;
    
    double ratio = (GLfloat) h / (GLfloat) w;
    double size = fmax(sizeX,sizeY);
    double left,right,top,bottom;
    if (ratio > 1.0) {
        left = meanX - size;
        right = meanX + size;
        bottom = meanY - size*ratio;
        top = meanY + size*ratio;  }   
    else {
        left = meanX - size/ratio;
        right = meanX + size/ratio;
        bottom = meanY - size;
        top = meanY + size;  }   
    if ((fabs(top-bottom)) >= (fabs(left-right))) {
        window->param.zoom = 1/(0.45*(fabs(top-bottom))); }   
    else {
        window->param.zoom = 1/(0.38*(fabs(left-right))); }   
  
    window->param.translate[0] = -(left + right) / 2.0;
    window->param.translate[1] = -(bottom + top) / 2.0;
  
    zoom_init = window->param.zoom;
    translate_init[0] = window->param.translate[0];
    translate_init[1] = window->param.translate[1];
    w_init = window->size[0];
    h_init = window->size[1];
    w = w_init;
    h = h_init;
    //
    // Default call back and help message
    //
  
    glfwSetKeyCallback(window->self, glfemKeyCallback);
    glfemWindowSetHelpMessage((const char[]) {
    "   [esc]   Exit\n"
    "    R      Reset zoom and translation\n"
    "    H      Display/hide keyboard shortcuts\n"});
}

void glfemWindowSetHelpMessage(const char *message) {    
        
    bov_text_delete(theCurrentWindow->help);
    theCurrentWindow->help = bov_text_new((const GLubyte*)(message),GL_STATIC_DRAW);

    bov_text_set_space_type(theCurrentWindow->help, PIXEL_SPACE);
    bov_text_set_fontsize(theCurrentWindow->help, 20.0f); 
    bov_text_set_boldness(theCurrentWindow->help, 0.1f);
    bov_text_set_outline_width(theCurrentWindow->help, 0.5f);
    bov_text_set_color(theCurrentWindow->help,GLFEM_BLACK);
}

void glfemWindowResetSize(){
		theCurrentWindow->param.zoom = zoom_init;
		theCurrentWindow->param.translate[0] = translate_init[0];
		theCurrentWindow->param.translate[1] = translate_init[1];
}

void glfemWindowUpdate()
{
    float w_ = theCurrentWindow->param.res[0];
    float w2 = theCurrentWindow->size[0];
	float ratio = w_/w2;
	
	w = theCurrentWindow->size[0];
	h = theCurrentWindow->size[1];
	  
    bov_text_set_fontsize(theCurrentWindow->help, ratio*20.0f);
    bov_text_set_pos(theCurrentWindow->help, (GLfloat[2]){0., h*0.55} );
    bov_window_update(theCurrentWindow);
}

void glfemWindowUpdateAndWait(){
    bov_text_set_pos(theCurrentWindow->help, 
            (GLfloat[2]){-20.0f, theCurrentWindow->size[1] - 30.0f} );
    bov_window_update_and_wait_events(theCurrentWindow);
    
}

void glfemWindowFree() 
{
    bov_window_delete(theCurrentWindow);
}

int glfemWindowShouldClose() 
{
    return bov_window_should_close(theCurrentWindow);
}

void glfemDrawMessage(char *message, double pos[2]){

    float w = theCurrentWindow->param.res[0];
    float w2 = theCurrentWindow->size[0];
	float ratio = w/w2;
    float ratioH = h/h_init;
    float ratioW = w/w_init;
  
    bov_text_t* text = bov_text_new((const GLubyte *)message, GL_STATIC_DRAW);
    text->param =  (bov_text_param_t) {
      .fillColor = {current_text_color[0],current_text_color[1],current_text_color[2],current_text_color[3]},
      .outlineColor = {1.0f ,1.0f, 1.0f, 2.0f},
      .pos = {0.0f, 0.0f},
      .shift = {0.0f, 0.0f},
      .fontSize = ratio*20.0f,
      .boldness = 0.0f,
      .outlineWidth = 0.0f,
      .spaceType = PIXEL_SPACE};
    bov_text_set_pos(text, (GLfloat[2]) {pos[0]*ratioW, pos[1]*ratioH});  
    bov_text_draw(theCurrentWindow, text);
    bov_text_delete(text);
}

void glfemDrawMsgParams(char *message, double pos[2], float color[4], double size, int type){
    float ratioH = h/h_init;
    float ratioW = w/w_init;
    
    bov_text_t* text = bov_text_new((const GLubyte *)message, GL_STATIC_DRAW);
    bov_text_set_space_type(text, type);
    bov_text_set_fontsize(text, size);
    bov_text_set_outline_width(text, 0.f);
    bov_text_set_color(text, color);
    bov_text_set_pos(text, (GLfloat[2]) {pos[0]*ratioW, pos[1]*ratioH});  
    bov_text_draw(theCurrentWindow, text);
    bov_text_delete(text);
}


void glfemDrawNodes(double *x, double *y, int n) 
{
    GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
    for(int i = 0; i < n; ++i){
      coord[i][0] = x[i];
      coord[i][1] = y[i];}
    bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
    bov_points_set_color(points,current_color);
    bov_points_set_width(points, current_line_width/zoom_init);
    bov_points_draw(theCurrentWindow,points, 0, BOV_TILL_END);
    bov_points_delete(points);
    free(coord);
}

void glfemDrawElement(double *x, double *y, int n)
{
    GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
    for(int i = 0; i < n; ++i){
      coord[i][0] = x[i];
      coord[i][1] = y[i];}
    bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
    bov_points_set_color(points,current_color);
    bov_points_set_width(points,current_line_width/zoom_init);
    bov_points_set_outline_width(points, current_line_width/zoom_init);
    bov_points_set_outline_color(points, current_line_color);

    bov_line_loop_draw(theCurrentWindow, points, 0, BOV_TILL_END);
    bov_points_delete(points);
    free(coord);
}

void glfemPlotMesh(femMesh *theMesh)
{
    int i,j,*nodes;
    int nLocalNode = theMesh->nLocalNode;

    GLfloat (*data)[2] = malloc(sizeof(data[0])*nLocalNode*theMesh->nElem);
    
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*nLocalNode]);
        for (j=0; j < nLocalNode; ++j) {
            data[i*nLocalNode+j][0] = theMesh->X[nodes[j]];
            data[i*nLocalNode+j][1] = theMesh->Y[nodes[j]]; }}
      
    bov_points_t* points = bov_points_new(data,nLocalNode*theMesh->nElem,GL_STATIC_DRAW);
  //   bov_points_set_color(points, GLFEM_BACKGROUND);
    bov_points_set_color(points, current_color);
    bov_points_set_width(points,current_line_width/zoom_init);
    bov_points_set_outline_width(points, 5*current_line_width/zoom_init);
    bov_points_set_outline_color(points, current_line_color);
    
    bov_triangles_draw(theCurrentWindow, points, 0, BOV_TILL_END);
    
    bov_points_delete(points);
    free(data);
}

void glfemSetScale(femMesh *theMesh, double *u)
{
    max_colormap = femMax(u,theMesh->nNode);
    min_colormap = femMin(u,theMesh->nNode);
}

double glfemScale(double minimum, double maximum, double value)
{
    if (value < minimum)        return 0;
    if (minimum == maximum)     return minimum;
    return (value - minimum) / fabs(maximum - minimum);
}

void glfemPlotSolution(femMesh* theMesh, double *u){
    int i,j,*nodes;
    int nLocalNode = theMesh->nLocalNode;
   
    for(int i = 0; i < theMesh->nElem; ++i){
        nodes = &(theMesh->elem[i*nLocalNode]);
        GLfloat (*data)[3] = malloc(sizeof(data[0])*3);
        for (j=0; j < 3; ++j) {
           data[j][0] = theMesh->X[nodes[j]];
           data[j][1] = theMesh->Y[nodes[j]]; 
           data[j][2] = glfemScale(min_colormap,max_colormap,u[nodes[j]]);} 

        bov_points_t* points = bov_points_new_with_value(data, 3, GL_STATIC_DRAW);   
        
        bov_points_set_color(points, GLFEM_BLACK);
        bov_points_set_width(points,current_line_width/zoom_init);
        bov_points_set_outline_color(points, GLFEM_BLACK);
        bov_points_set_outline_width(points, 5*current_line_width/zoom_init);
        bov_triangles_draw(theCurrentWindow, points, 0, BOV_TILL_END);
        bov_points_delete(points);
        free(data);}
}

void getInfo(const motorMesh *theMesh, double *soluce, const int iElem, int *map, double *x, double *y, double *u) {

    int nLocal = theMesh->nLocalNode;    
    for (int j = 0 ; j < nLocal ; ++j) {
        map[j] = theMesh->elem[iElem*nLocal + j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]];
        u[j]   = soluce[map[j]];
    }
}

double *getGradientPLOT(motor *theMotor){

    motorMesh *theMesh = theMotor->mesh;
    double *soluce = theMotor->a;
    double *dudx = calloc(theMotor->size, sizeof(double));
    double *dudy = calloc(theMotor->size, sizeof(double));
    int *nbElem = calloc(theMotor->size, sizeof(int));
    
    int nLocalNode = theMesh->nLocalNode;
    int iElem, i, j, map[3];
    double Xloc[3], Yloc[3], Uloc[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3];
    
    double xsi, eta, dxdxsi, dxdeta, dydxsi, dydeta, jac, invJ;
    double dadx, dady, gradLoc;
    
    xsi = 1./3.;
    eta = 1./3.;
    dphidxsi[0] = -1.0; dphidxsi[1] = 1.0; dphidxsi[2] = 0.0;
    dphideta[0] = -1.0; dphideta[1] = 0.0; dphideta[2] = 1.0;
    
    for (iElem = 0 ; iElem < theMesh->nElem ; iElem++) {
    
        getInfo(theMesh, soluce, iElem, map, Xloc, Yloc, Uloc);
        
        dxdxsi = 0 ; dxdeta = 0; 
        dydxsi = 0 ; dydeta = 0;
        
        for (i = 0 ; i < nLocalNode ; i++) {
            dxdxsi += Xloc[i]*dphidxsi[i]; dxdeta += Xloc[i]*dphideta[i];   
            dydxsi += Yloc[i]*dphidxsi[i]; dydeta += Yloc[i]*dphideta[i];}
            
        jac = dxdxsi * dydeta - dxdeta * dydxsi;
        invJ = 1. / jac;
        for (i = 0 ; i < nLocalNode ; i++) {
            dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) * invJ;
            dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) * invJ;
            dudx[map[i]] += Uloc[i] * dphidx[i];
            dudy[map[i]] += Uloc[i] * dphidy[i];
            nbElem[map[i]] ++;
        }
    }
    
    for (i = 0 ; i < theMotor->size ; i++){
        dudx[i] = sqrt(dudx[i]*dudx[i] + dudy[i]*dudy[i]) / (1.*nbElem[i]);
    }
    
    
            
    //if (sqrt(dudx[map[i]]*dudx[map[i]] + dudy[map[i]]*dudy[map[i]]) > 0.1 && theMesh->domain[iElem] == 8)
    //    printf("iElem : %d \t grad %lf\n", iElem, sqrt(dudx[map[i]]*dudx[map[i]] + dudy[map[i]]*dudy[map[i]]));
    /*int countBelow = 0;
    int countAbove = 0;
    double mu[3];
    double m = 0.;
    double mx = 0.;
    double mxd = 0.;
    for (i = 0 ; i < theMotor->size ; i++){
        getMu(theMotor->hystereticCurveB, theMotor->hystereticCurveH, theMotor->nHystereticCurve, dudx[i]*dudx[i], mu);
        double dmudb2 = (mu[2] - mu[0])/(4.* dudx[i] * 1e-4);
        mx = fmax(mx, mu[1]);
        mxd = fmax(mxd, dmudb2);
        m = fmax(dudx[i]*dudx[i], m);
        //printf("1/mu = %10.6le\t %10.6le\n", mu[1], dmudb2);
    }
    printf("b2 = %lf\t 1/mu = %lf\t d1/mu/db2 = %lf\n", m, mx, mxd);*/
    free(dudy);
    return dudx;
}


int interpolate(double *X, double *Y, double *U, double u, int n, GLfloat coord[][2]){
    int i, j, init=n;
    double xsi, eta;
    if ((U[0] - u) * (U[1] - u) < 0.) {
        U[0] -= 1e-10; U[1] -= 1e-10;
        coord[n][0] = X[0] + (X[1] - X[0]) / (U[1] - U[0]) * (u - U[0]);
        coord[n][1] = Y[0] + (Y[1] - Y[0]) / (U[1] - U[0]) * (u - U[0]);
        n++;
    }
    if ((U[0] - u) * (U[2] - u) < 0.) {
        U[0] -= 1e-10; U[2] -= 1e-10;
        coord[n][0] = X[0] + (X[2] - X[0]) / (U[2] - U[0]) * (u - U[0]);
        coord[n][1] = Y[0] + (Y[2] - Y[0]) / (U[2] - U[0]) * (u - U[0]);
        n++;
    }
    if ((U[1] - u) * (U[2] - u) < 0.) {
        U[1] -= 1e-10; U[2] -= 1e-10;
        coord[n][0] = X[1] * (u - U[2]) / (U[1]-U[2]) + X[2] * (U[1] - u) / (U[1] - U[2]);
        coord[n][1] = Y[1] * (u - U[2]) / (U[1]-U[2]) + Y[2] * (U[1] - u) / (U[1] - U[2]);
        n++;
    }
    if ((n - init)%2 != 0) {
        return init;
        printf("n_prev = %d \t n_after = %d\n", init, n);
        printf("U1 = %.3le\t U2 = %.3le\t U3 = %.3le\t u = %.3le\n", U[0], U[1], U[2], u);
    }
    return n;
}

void drawLevels(motor *theMotor, double *f, double level, int subMode){
    
    motorMesh *theMesh = theMotor->mesh;
    
    int map[3], nLoc = theMesh->nLocalNode;
    int shift = 0;
    double Xloc[3], Yloc[3], Uloc[3];
    
    GLfloat (*coord)[2] = malloc(sizeof(coord[0]) * 2 * theMesh->nElem);
	
	int n, iElem;
	if (subMode == 0) {
	    for (iElem = 0, n=0 ; iElem < theMesh->nElemDomain[0] ; iElem++){
            getInfo(theMesh, f, iElem, map, Xloc, Yloc, Uloc);
            n = interpolate(Xloc, Yloc, Uloc, level, n, coord);
	    }
	    for (int i = 0 ; i < 8 ; i++) // 8 for rotor_core
	        shift += theMesh->nElemDomain[i];
	        
	    for (iElem = shift ; iElem < shift+theMesh->nElemDomain[8] ; iElem++){
            getInfo(theMesh, f, iElem, map, Xloc, Yloc, Uloc);
            n = interpolate(Xloc, Yloc, Uloc, level, n, coord);
	    }
	}
	else {
        for (iElem = 0, n=0 ; iElem < theMesh->nElem ; iElem++){
            getInfo(theMesh, f, iElem, map, Xloc, Yloc, Uloc);
            n = interpolate(Xloc, Yloc, Uloc, level, n, coord);
	    }
	}
	
    bov_points_t* pointset = bov_points_new(coord, n, GL_STATIC_DRAW);
    bov_points_set_width(pointset, 0.001/zoom_init);
    bov_points_set_color(pointset, (float[4]) {0., 0., 0., 0.35});
    bov_fast_lines_draw(theCurrentWindow, pointset, 0, BOV_TILL_END);
    bov_points_delete(pointset);    
    
    return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// - TOOL FUNCTIONS - //////////////////////////////////////////////////////

void shadeTheWindow(animation *myAnim, double alpha) {
    GLfloat coordRECT[][2] = {{h/2., w/2.}};
    bov_points_t *pointsRECT = bov_points_new(coordRECT, 1, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsRECT, 2);
    SHADE_TEXT[3] = alpha;
    bov_points_set_color(pointsRECT, SHADE_TEXT);
    bov_points_set_width(pointsRECT, w);
    bov_points_set_marker(pointsRECT, 0.f);
    bov_points_draw(theCurrentWindow, pointsRECT, 0, 1);
    bov_points_delete(pointsRECT);    
}

void makePlay(animation *myAnim) {
    myAnim->startTime = glfwGetTime();
    myAnim->stopTime = myAnim->duration;
    myAnim->discreteTime = 0.;
    myAnim->running = 1;
}

void makePause(animation *myAnim) {
    myAnim->running = 0;
    myAnim->stopTime = 0.;
    myAnim->discreteTime = 0.;
    myAnim->time = 0.;
}

void applyMeshChange(animation *myAnim) {
    int nbc = 8;
    const char meshMessage[] = {"%-*s1 \tFOR\t MESH XS\n" "%-*s2 \tFOR\t MESH  S\n"
                                "%-*s3 \tFOR\t MESH  M\n" "%-*s4 \tFOR\t MESH  L\n"
                                "%-*s5 \tFOR\t MESH XL\n"};
    sprintf(myAnim->messages[3], meshMessage, nbc, "", nbc, "", nbc, "TAP", nbc, "", nbc, "");
    
    myAnim->funStartTime = glfwGetTime();
}

void playButton(animation *myAnim, GLfloat coordPLAY[][2]) {
    bov_points_t *pointsPLAY = bov_points_new(coordPLAY, 3, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsPLAY, 2);
    bov_points_set_color(pointsPLAY, (float[4]) {0.0, 0.0, 0.0, 0.5});
    bov_triangles_draw(theCurrentWindow, pointsPLAY, 0, BOV_TILL_END);
    bov_points_delete(pointsPLAY);
    /*double ratioW = w/w_init;
    double ratioH = h/h_init;
    glfemDrawMsgParams("\x91", (double[2]){44.*ratioW, 529.*ratioH}, BLACK_TEXT, 50.f, PIXEL);*/
}

void pauseButton(animation *myAnim, GLfloat coordPAUSE[][2]) {
    bov_points_t *pointsPAUSE = bov_points_new(coordPAUSE, 2, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsPAUSE, 2);	
    bov_points_set_color(pointsPAUSE, TRSP_TEXT);
    bov_points_set_width(pointsPAUSE, 17.);	
    bov_points_set_outline_color(pointsPAUSE, (float[4]) {0.0, 0.0, 0.0, 0.5});
    bov_points_set_outline_width(pointsPAUSE, 8.);
    bov_points_set_marker(pointsPAUSE, 1.175f);
    bov_points_draw(theCurrentWindow, pointsPAUSE, 0, 2);
    bov_points_delete(pointsPAUSE);
    /*double ratioW = w/w_init;
    double ratioH = h/h_init;
    glfemDrawMsgParams("\x92", (double[2]){38.*ratioW, 529.*ratioH}, BLACK_TEXT, 50.f, PIXEL);*/
}

void circleButton(animation *myAnim, GLfloat coordCIRCLE[][2]) {
	bov_points_t *pointsCIRCLE = bov_points_new(coordCIRCLE, 1, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsCIRCLE, 2);
	bov_points_set_color(pointsCIRCLE, TRSP_TEXT);
	bov_points_set_width(pointsCIRCLE, 40.);
	bov_points_set_outline_color(pointsCIRCLE, (float[4]) {0.0, 0.0, 0.0, 0.5});
	bov_points_set_outline_width(pointsCIRCLE, 8.);
	bov_points_set_marker(pointsCIRCLE, 0.f);
	bov_points_draw(theCurrentWindow, pointsCIRCLE, 0, 1);
    bov_points_delete(pointsCIRCLE); 
}

void easterEgg(animation *myAnim, int n) {
    glfemDrawMsgParams("Easter Eggs found :", (double[2]){654.f, 40.f}, BLACK_TEXT, 20.f, PIXEL);
    float ratioH = h/h_init;
    float ratioW = w/w_init;
        
    int count = 0;
    for (int i = 0 ; i < strlen(myAnim->easterEgg) ; i++)
        if (myAnim->easterEgg[i] == '1') 
            count ++;
            
    float end = 650.f + 200.f*count/n;
    
    GLfloat coordFIXE[][2] = {{750.f*ratioW, 20.*ratioH}};
    bov_points_t *pointsF = bov_points_new(coordFIXE, 1, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsF, 2);
    bov_points_set_color(pointsF, TRSP_TEXT);
    bov_points_set_width(pointsF, 100.*ratioW);
    bov_points_set_outline_color(pointsF, (float[4]) {0.0, 0.0, 0.0, 0.7});
    bov_points_set_outline_width(pointsF, 2.f);
    bov_points_set_marker(pointsF, 2.1f);
    bov_points_draw(theCurrentWindow, pointsF, 0, 2);
    bov_points_delete(pointsF);

    GLfloat coordMOVE[][2] = {{650.f*ratioW, 10*ratioH}, {end*ratioW, 10*ratioH}, {end*ratioW, 30.*ratioH},
                              {end*ratioW, 29.8*ratioH}, {650.f*ratioW, 29.8*ratioH}, {650.f*ratioW, 9.8*ratioH} };
    bov_points_t* pointsM = bov_points_new(coordMOVE, 6, GL_STATIC_DRAW);
    bov_points_set_space_type(pointsM, 2);
    bov_points_set_color(pointsM, (float[4]) {0.0, 0.0, 0.0, 0.7});
    bov_points_set_width(pointsM, 0.*ratioW);
    bov_triangles_draw(theCurrentWindow, pointsM, 0, BOV_TILL_END);
    bov_points_delete(pointsM);
    
    if (count == 6 && myAnim->completed == 0) {
        myAnim->delayTime = glfwGetTime();
        myAnim->completed = 1;
    }
    if (count == 6 && myAnim->delayTime + 3. < glfwGetTime() && myAnim->flag != 4 && myAnim->reallyCompleted == 0) {
        myAnim->flag = 4;
        myAnim->funStartTime = glfwGetTime();
        myAnim->reallyCompleted = 1;
    }
}

void goToMode0(animation *myAnim) {
    myAnim->mode = 0;
    myAnim->subMode = 0;
    myAnim->nbSubModes = 7;
}

void goToMode1(animation *myAnim) {
    myAnim->mode = 1;
    myAnim->subMode = 0;
    myAnim->nbSubModes = 2;
}

void expNotation(double value, double *result) {
    if (fabs(value) < 1e-10) {
        result[1] = 0.;
        result[0] = 0.;
    }else {
        result[1] = 3. * floor(log10(fabs(value)) / 3.);
        result[0] = value / pow(10, result[1]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////// - DISPLAY FUNCTIONS - /////////////////////////////////////////////////////

void displayMode0(animation *myAnim, motor *theMotor, int subMode, femMesh **theDomains) {
    motorMesh *theMotorMesh = theMotor->mesh;
    int indexMsg = 1;
    
    if (subMode == 6) {
        sprintf(myAnim->messages[indexMsg], "Bravo, on voit plus rien maintenant");
        myAnim->easterEgg[4] = '1';
        return;
    }
    
    glfemPlotMesh((femMesh*) theMotorMesh);
    
    femMesh *theStator    = theDomains[0]; femMesh *theGap       = theDomains[11];
    femMesh *theCoil_AP   = theDomains[1]; femMesh *theCoil_AN   = theDomains[2];
    femMesh *theCoil_BP   = theDomains[3]; femMesh *theCoil_BN   = theDomains[4];
    femMesh *theCoil_CP   = theDomains[5]; femMesh *theCoil_CN   = theDomains[6];
    femMesh *theStatorGap = theDomains[7]; femMesh *theRotor     = theDomains[8];
    femMesh *theRotorAir  = theDomains[9]; femMesh *theRotorGap  = theDomains[10];
    
    femMesh *theCoils[6]  = {theCoil_AP , theCoil_AN, theCoil_BP,
                             theCoil_BN, theCoil_CP , theCoil_CN};
                             
    if (subMode == 0) {
        glfemSetColor(GLFEM_BLUE); glfemPlotMesh(theRotor); glfemPlotMesh(theStator);
        glfemSetColor(GLFEM_GREEN); glfemPlotMesh(theGap);
        glfemSetColor(GLFEM_RED); glfemPlotMesh(theRotorGap); glfemPlotMesh(theStatorGap);
        
        for (int i = 0 ; i < 6 ; i++){
            if (theMotor->js[1+i] > 1) {
                glfemSetColor(GLFEM_RED);
                glfemPlotMesh(theCoils[i]);}
            if (theMotor->js[1+i] < -1) {
                glfemSetColor(GLFEM_BLACK);
                glfemPlotMesh(theCoils[i]);}}}
                
    if (subMode == 1) {
        glfemSetColor(GLFEM_BLUE); glfemPlotMesh(theStator); 
        sprintf(myAnim->messages[indexMsg], "Stator ");}
        
    else if (subMode == 2) {
        glfemSetColor(GLFEM_BLUE); glfemPlotMesh(theRotor);
        sprintf(myAnim->messages[indexMsg], "Rotor ");}
        
    else if (subMode == 3) {
        glfemSetColor(GLFEM_RED);
        for (int i = 0 ; i < 6 ; i+=2)
            glfemPlotMesh(theCoils[i]);
        glfemSetColor(GLFEM_BLACK);
        for (int i = 1 ; i < 6 ; i+=2)
            glfemPlotMesh(theCoils[i]);
        sprintf(myAnim->messages[indexMsg], "Bobines ");}
        
    else if (subMode == 4) {
        glfemSetColor(GLFEM_RED); glfemPlotMesh(theRotorGap);
        glfemPlotMesh(theStatorGap); glfemPlotMesh(theRotorAir);
        glfemSetColor(GLFEM_GREEN); glfemPlotMesh(theGap);
        sprintf(myAnim->messages[indexMsg], "Entrefer ");}
        
    else if (subMode == 5) {
        glfemSetColor(GLFEM_GREEN); glfemPlotMesh(theGap);
        sprintf(myAnim->messages[indexMsg], "Zone a maillage variable");
    }
}

void displayMode1(animation *myAnim, motor *theMotor, int subMode, femMesh **theDomains) {
    
    motorMesh *theMotorMesh = theMotor->mesh;
    glfemPlotMesh((femMesh*) theMotorMesh);
    
    femMesh *theStator    = theDomains[0]; femMesh *theGap       = theDomains[11];
    femMesh *theCoil_AP   = theDomains[1]; femMesh *theCoil_AN   = theDomains[2];
    femMesh *theCoil_BP   = theDomains[3]; femMesh *theCoil_BN   = theDomains[4];
    femMesh *theCoil_CP   = theDomains[5]; femMesh *theCoil_CN   = theDomains[6];
    femMesh *theStatorGap = theDomains[7]; femMesh *theRotor     = theDomains[8];
    femMesh *theRotorAir  = theDomains[9]; femMesh *theRotorGap  = theDomains[10];
    
    femMesh *theCoils[6]  = {theCoil_AP , theCoil_AN, theCoil_BP,
                             theCoil_BN, theCoil_CP , theCoil_CN};
    
    double *u = theMotor->a;
    double *grad;
    
    glfemSetScale((femMesh*)theMotorMesh, u);
    glfemPlotSolution(theRotor, u);
    glfemPlotSolution(theStator, u);
    
    if (subMode == 1) {
        for (int i = 0 ; i < 6 ; i++)
            glfemPlotSolution(theCoils[i], theMotor->a);
        glfemPlotSolution(theStatorGap, theMotor->a);
        glfemPlotSolution(theRotorGap, theMotor->a);
        glfemPlotSolution(theGap, theMotor->a);
        glfemPlotSolution(theRotorAir, theMotor->a);
    }
    
    if (myAnim->fieldLines) {
        //grad = getGradientPLOT(theMotor);
        int nLevels = 10;
        for (int i = 1 ; i < nLevels ; i ++) {
            drawLevels(theMotor, theMotor->a, max_colormap * i/((double) nLevels), myAnim->subMode);
            drawLevels(theMotor, theMotor->a, min_colormap * i/((double) nLevels), myAnim->subMode);
        }
    }
}

void displayButton(animation *myAnim) {
    float ratioH = h/h_init;
    float ratioW = w/w_init;
    
    GLfloat coordPLAY[][2] = { {50.*ratioW, 530.*ratioH}, {75.*ratioW, 545.*ratioH}, {50.*ratioW, 560.*ratioH} };
    GLfloat coordPAUSE[][2] = {{59.*ratioW-7, 545.*ratioH}, {59*ratioW+7., 545.*ratioH}};
    GLfloat coordCIRCLE[][2] = { {59.*ratioW, 545.*ratioH}};
    
    circleButton(myAnim, coordCIRCLE);    
    if (myAnim->running) playButton(myAnim, coordPLAY);
    else pauseButton(myAnim, coordPAUSE);
}

void displayNonLinear(animation *myAnim) {
    glfemDrawMsgParams("\x8a", (double[2]){30.f, 350.f}, (GLfloat[4]){1., 0., 0., 0.9}, 50.f, PIXEL);
    glfemDrawMsgParams(" non\nlinear", (double[2]){30.f, 330.f}, (GLfloat[4]){1., 0., 0., 0.9}, 20.f, PIXEL);
}

void displayCurrent(animation *myAnim) {
    glfemDrawMsgParams("Manual", (double[2]){30.f, 265.f}, (GLfloat[4]){0., 0., 0., 0.75}, 18.f, PIXEL);
    glfemDrawMsgParams("Current", (double[2]){27.f, 250.f}, (GLfloat[4]){0., 0., 0., 0.75}, 18.f, PIXEL);
}

void displayCredits(animation *myAnim) {
    float t = bov_window_get_time(theCurrentWindow);
    size_t len1, len2, len3, len=0;
	float a = 0.1;
	float beat = 2*M_PI*70/60;  // heartbeat frequency
    
	float varColor[4] = {sin(beat/8. * t) * 0.1 + 0.25,
                         sin(beat/4. * t + 2.0) * 0.3 + 0.5,
                         sin(beat/2. * t + 5.0) * 0.1 + 0.75,
                         1.0};
	
	float wSerie[3] = {w*.7f, w*.5f, w*.3f};
	float hSerie[3] = {h*.6f, h*.45f, h*.25f};
	float bSerie[3] = {wSerie[0]*2./9., wSerie[1]*2./11., wSerie[2]*2./8.};
	char sSerie[][12] = {"LEPL 1110", "2020 - 2021", "\xb6 pheniXXC"};
	size_t lenSerie[3] = {strlen(sSerie[0]), strlen(sSerie[1]), strlen(sSerie[2])};
	
	for (int i = 1 ; i < 1+lenSerie[0]+lenSerie[1]+lenSerie[2] ; i++)
	    if (i*a < t - myAnim->funStartTime) 
	        len = i;
	
	if (len < 1+lenSerie[0]) {
	    sSerie[0][len]=0;
	    sSerie[1][0]=0;
	    sSerie[2][0]=0;
	}
	else if (len < 1+lenSerie[0]+lenSerie[1]) {
	    len -= lenSerie[0];
	    sSerie[0][lenSerie[0]]=0;
	    sSerie[1][len] = 0;
	    sSerie[2][0]=0;
	}
	else {
	    len -= lenSerie[0] + lenSerie[1];
	    sSerie[0][lenSerie[0]]=0;
	    sSerie[1][lenSerie[1]]=0;
	    sSerie[2][len] = 0;
	}
	
	bov_text_t* text1 = bov_text_new((const GLubyte *)sSerie[0], GL_STATIC_DRAW);
    bov_text_t* text2 = bov_text_new((const GLubyte *)sSerie[1], GL_STATIC_DRAW);
    bov_text_t* text3 = bov_text_new((const GLubyte *)sSerie[2], GL_STATIC_DRAW);
    
    bov_text_t *textSerie[3] = {text1, text2, text3};
    
	bov_text_set_boldness(text1, .4 * sin(beat*t) - .3);
	bov_text_set_boldness(text2, .7 * sin(beat*t-3*M_PI/12.) - .2);
	bov_text_set_boldness(text3, .7 * sin(beat*t-6*M_PI/12.) - .2);
    
    for (int i = 0 ; i < 3 ; i++) {
        if (i < 2)
            bov_text_set_pos(textSerie[i], (GLfloat[2]) {(w - wSerie[i])/2.f, hSerie[i]});
        else 
            bov_text_set_pos(textSerie[i], (GLfloat[2]) {(w - 1.2*wSerie[i])/2.f, hSerie[i]});
        bov_text_set_fontsize(textSerie[i], bSerie[i]);
        bov_text_set_space_type(textSerie[i], 2);
        bov_text_set_outline_width(textSerie[i], 0.f);
        bov_text_set_color(textSerie[i], varColor);
        bov_text_draw(theCurrentWindow, textSerie[i]);
        bov_text_delete(textSerie[i]);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// - MAIN FUNCTIONS - //////////////////////////////////////////////////////

void updateParameters(animation *myAnim) {

    int nb = myAnim->nbSubModes;
    char action = myAnim->action;
    int i, count;
    
    // get the time
    myAnim->time = (glfwGetTime() - myAnim->startTime) * myAnim->fastFactor;
    
    // if stop reached, time is freezed
    if (myAnim->stopTime <= myAnim->time) myAnim->time = myAnim->stopTime;
    if (myAnim->discreteTime < myAnim->stopTime)
        myAnim->running = 1;
    else
        myAnim->running = 0;
    
    // flag = 0 : nothing special
    // flag = 1 : asked for non-existent mesh
    // flag = 2 : 
    // flag = 3 : 
    // mode = 0 : colored domains
    // mode = 1 : magnetic potential
    // mode = 3 : credits
    //////////////////////////
    
    if ( myAnim->meshChange == 0 && 
        ((action == 'u' && myAnim->mode==1) || (action == 'd' && myAnim->mode==0))) {
        action = 'M';
        myAnim->meshChange = 0;
    }
    
    switch (action) {
        case 'r':
            myAnim->subMode = ((myAnim->subMode + 1) % nb + nb) % nb; //always positive
            break;
            
        case 'l':
            myAnim->subMode = ((myAnim->subMode - 1) % nb + nb) % nb; //always positive
            break;
        
        case 'u':
            if (myAnim->meshChange) {
                goToMode0(myAnim);
                myAnim->meshChange = 0;
            }
            else // handled all other possibilities above (+/-10 l.)
                goToMode1(myAnim);
            break;
        
        case 'd':
            if (myAnim->meshChange) {
                goToMode1(myAnim);
                myAnim->meshChange = 0;
            }
            else // handled other possibilities above (+/-10 l.)
                goToMode0(myAnim);
            break;
                
        case 'C':
            if (myAnim->mode == 3) {
                myAnim->mode = 1;
                myAnim->informations = 1;
            }
            else {
                myAnim->mode = 3;
                myAnim->informations = 0;
                myAnim->flag = 0;
                myAnim->funStartTime = glfwGetTime();
                myAnim->easterEgg[3] = '1';
                makePause(myAnim);
            }
            break;
        
        case 'D':
            myAnim->displayMesh = 1 - myAnim->displayMesh;
            break;
        
         case 'F':
            myAnim->fieldLines = 1 - myAnim->fieldLines;
            break;
        
        case 'H':
            if (theCurrentWindow->help_needed)
                myAnim->informations = 0;
            else 
                myAnim->informations = 1;
            break;    
        
        case 'I':
            myAnim->informations = 1 - myAnim->informations;
            break;
        
        case 'J':
            if (myAnim->manCurrent == 0) {
                myAnim->manCurrent = 1;
                break;
            }
            for (i=0, count=0; (i < 6) && (count < 2) ; i++) {
                if (fabs(myAnim->motor->js[1+i]) > 1e-3) {
                    myAnim->motor->js[1+(i+2)%6] = myAnim->motor->js[1+i];
                    myAnim->motor->js[1+i] = 0.;
                    count ++;
                }
            }
            break;
        
        case 'K':
            if ((myAnim->mode == 0) && (myAnim->meshChange == 0)) {
                myAnim->easterEgg[5] = '1';
                myAnim->flag = 3;
                myAnim->funStartTime = glfwGetTime();
            }
            goToMode0(myAnim);
            break;
            
        case 'L':
            myAnim->motor->nonLinearFlag = 1 - myAnim->motor->nonLinearFlag;
            break;
        
        case 'M':
            if (myAnim->meshChange == 0) {
                myAnim->meshChange = 1;
                makePause(myAnim);
                myAnim->informations = 1;
                applyMeshChange(myAnim);
            }
            else {
                myAnim->meshChange = 0;
                myAnim->flag = 0;
            }
            break;
        
        case 'S':
            if ((myAnim->mode == 1) && (myAnim->meshChange == 0)) {
                myAnim->easterEgg[5] = '1';
                myAnim->flag = 3;
                myAnim->funStartTime = glfwGetTime();
            }
            goToMode1(myAnim);
            break;
        
        case 'P':
            if (myAnim->running) makePause(myAnim);
            else if (myAnim->mode == 0 || myAnim->mode == 1) makePlay(myAnim);
            break;
        
        case 'R':
            motorAdaptMesh(myAnim->motor, -1. * myAnim->motor->theta);
            motor *theMotor = myAnim->motor;
            for (int i = 0 ; i < theMotor->size ; i++){
                double th = atan2(theMotor->mesh->Y[i], theMotor->mesh->X[i]);
                double r = hypot(theMotor->mesh->Y[i], theMotor->mesh->X[i]);
                if (theMotor->movingNodes[i]) theMotor->a[i] = sin(3*th);
                else theMotor->a[i] = -sin(3*th);
            }
            myAnim->motor->theta = 0.;
            myAnim->motor->omega = 0.;
            myAnim->running = 0;
            myAnim->stopTime = 0.;
            myAnim->time = 0.;
            myAnim->discreteTime = 0.;
            myAnim->fullTime = 0.;
            myAnim->iteration = 0;
            break;
        
        case 'A':
            //myAnim->fastFactor *= 1.1;
            break;
            
        case 'Z':
            //myAnim->fastFactor /= 1.1;
            break;
    }
}

void updateDisplay(animation *myAnim) {

    if (myAnim->displayMesh)
        current_line_width = 0.0001;
    else
        current_line_width = 0.;
    
    if (myAnim->mode == 0)
        displayMode0(myAnim, myAnim->motor, myAnim->subMode, myAnim->domain);
    
    else if (myAnim->mode == 1)
        displayMode1(myAnim, myAnim->motor, myAnim->subMode, myAnim->domain);
    
    else if (myAnim->mode == 3) {
        displayCredits(myAnim);
        return;
    }
    
    if (myAnim->informations)
        displayButton(myAnim);
        
    if (myAnim->informations && myAnim->motor->nonLinearFlag)
        displayNonLinear(myAnim);
        
    if (myAnim->informations && myAnim->manCurrent)
        displayCurrent(myAnim);
        
    if (theCurrentWindow->help_needed || myAnim->meshChange)
        shadeTheWindow(myAnim, 0.95);
}

void updateMessages(animation *myAnim) {
    
    motor *theMotor = myAnim->motor;
    motorMesh *theMotorMesh = myAnim->motor->mesh;
    
    sprintf(myAnim->messages[0], "# nodes = %5d\n# elems = %5d",
            theMotorMesh->nNode, theMotorMesh->nElem);
       
    if (myAnim->informations == 0) 
        return;
        
    glfemDrawMsgParams(myAnim->messages[0],(double[2]){20.0, 650.0}, BLACK_TEXT, 20.f, PIXEL);
    double res1[2], res2[2];
    
    if ((myAnim->mode == 1) || (myAnim->subMode == 0)) {
        sprintf(myAnim->messages[1], "Time      = %5.2f s\niteration = %5d",
                myAnim->fullTime, myAnim->iteration);
        
        expNotation(myAnim->C, res1);
        expNotation(theMotor->omega * 30 / M_PI, res2);
        sprintf(myAnim->messages[2], "C = %7.2lf E%-2.0f N m\n\xdc = %7.2lf E%-2.0f rpm\n\xcb = %7.2lf     deg\n", 
                res1[0], res1[1], res2[0], res2[1], fmod(180./M_PI * theMotor->theta, 360.));

        glfemSetTextColor(LGRED_TEXT);
        glfemDrawMessage(myAnim->messages[1], (double[2]){20.0, 45.0});
        glfemDrawMessage(myAnim->messages[2], (double[2]){685.0, 672.0});           
    }
    else if ((myAnim->mode == 0) || (myAnim->subMode > 0)) {
        glfemDrawMessage(myAnim->messages[1], (double[2]){20.0, 25.0});
    }
    
    if (myAnim->meshChange) {
        glfemDrawMsgParams(myAnim->messages[3], (double[2]){90.f, 420.f}, RED_TEXT, 40.f, PIXEL);
        glfemDrawMsgParams("\x87", (double[2]){800.f, 350.f}, BLACK_TEXT, 50.f, PIXEL);
        if (myAnim->funStartTime + 2. < glfwGetTime()) {
            myAnim->easterEgg[0] = '1';
            sprintf(myAnim->messages[4], "    Alors, on se decide \n" "j'ai pas tout mon temps ...\n\n" "             \x89");
            glfemDrawMsgParams(myAnim->messages[4], (double[2]){315, 560.f}, BLACK_TEXT, 20.f, PIXEL);            
        }
    }
    else if (myAnim->mode == 0) 
        glfemDrawMsgParams("\x83", (double[2]){800.f, 420.f}, BLACK_TEXT, 50.f, PIXEL);
    else if (myAnim->mode == 1)
        glfemDrawMsgParams("\x80", (double[2]){800.f, 490.f}, BLACK_TEXT, 50.f, PIXEL);
    
    if (myAnim->flag == 1) {
        myAnim->easterEgg[1] = '1';
        sprintf(myAnim->messages[4], "    WOW WOW WOW\n" "     tu es fou\n" "tu veux une segfault \n" "         ?\n");
        glfemDrawMsgParams(myAnim->messages[4], (double[2]){135.f, 210.f}, BLACK_TEXT, 60.f, PIXEL);
        if (myAnim->funStartTime + 3. < glfwGetTime())
            myAnim->flag = 0;
    }
    else if (myAnim->flag == 2) {
        myAnim->easterEgg[2] = '1';
        shadeTheWindow(myAnim, 0.85);
        sprintf(myAnim->messages[4], "  OOOOUFFF,\n" "  ca c'est \n" "BIG BIG BIG !\n");
        glfemDrawMsgParams(myAnim->messages[4], (double[2]){110.f, 450.f}, BLACK_TEXT, 100.f, PIXEL);
        if (myAnim->funStartTime + 3. < glfwGetTime())
            myAnim->flag = 0;
    }
    else if (myAnim->flag == 3) {
        shadeTheWindow(myAnim, 0.85);
        sprintf(myAnim->messages[4], " Oui, 'fin,\n" "on n'y etait\n" " pas deja ?");
        glfemDrawMsgParams(myAnim->messages[4], (double[2]){300.f, 400.f}, BLACK_TEXT, 50.f, PIXEL);
        if (myAnim->funStartTime + 1.5 < glfwGetTime())
            myAnim->flag = 0;
    }
    
    easterEgg(myAnim, 6);
    if (myAnim->flag == 4) {
        shadeTheWindow(myAnim, 1.);
        sprintf(myAnim->messages[4], "\xa1");
        glfemDrawMsgParams(myAnim->messages[4], (double[2]){w_init*0.33, h_init*0.3}, BLACK_TEXT, 300.f, PIXEL);
        if (myAnim->funStartTime + 2. < glfwGetTime())
            myAnim->flag = 0;
    }
    
}

