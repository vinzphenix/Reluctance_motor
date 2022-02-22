/*
 *  glfem.h - BOV version
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include "BOV.h"
#include <stdlib.h>
#define GLFW_INCLUDE_GLU
#include "motor.h"
#include <GLFW/glfw3.h>

static float GLFEM_BLACK[4] = {0.0, 0.0, 0.0, 1.0};
static float GLFEM_BLUE[4] = {0.0, 0.0, 1.0, 1.0};
static float GLFEM_RED[4] = {1.0, 0.0, 0.0, 1.0};
static float GLFEM_GREEN[4] = {0.0, 1.0, 0.0, 1.0};
static float GLFEM_BACKGROUND[4] = {0.9, 0.9, 0.8, 0.0};

void glfemWindowCreate(const char *windowName, int w, int h, int n, double *x, double *y);
void glfemWindowFree();
void glfemWindowUpdate();
void glfemWindowUpdateAndWait();
int glfemWindowShouldClose();
void glfemWindowSetHelpMessage(const char *message);
char glfemGetAction();

void glfemSetScale(femMesh *theMesh, double *u);
void glfemSetColor(float color[4]);
void glfemSetTextColor(float color[4]);
void glfemSetLineWidth(float width);
void glfemSetColorLine(float color[4]);
void glfemDrawMessage(char *message, double pos[2]);
void glfemDrawMsgParams(char *message, double pos[2], float color[4], double size, int type);
void glfemDrawNodes(double *x, double *y, int n);
void glfemDrawElement(double *x, double *y, int n);

static void glfemKeyCallback(GLFWwindow *self, int key, int scancode, int action, int mods);
int glfemGetKey(char theKey);

void glfemPlotMesh(femMesh *theMesh);
double glfemScale(double minimum, double maximum, double value);
void glfemPlotSolution(femMesh *theMesh, double *u);

void makePlay(animation *myAnim);
void updateParameters(animation *myAnim);
void updateDisplay(animation *myAnim); //, bov_points_t* pointsPLAY, bov_points_t* pointsPAUSE);
void updateMessages(animation *myAnim);

#endif
