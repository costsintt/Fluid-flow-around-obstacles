#ifndef FLOWCALC_H
#define FLOWCALC_H

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <stdio.h>

struct BorderValuesCalculator
{
    double bordersVel;
    void(*calc)(struct BorderValuesCalculator *self,
                 size_t xSize, size_t ySize,
                 double** dzeta, double** psi, double** u, double** v,
                 double Re, double dx, double dy, double dt);
};

void calcBorderValues(double velOfMovementBorders, size_t xSize, size_t ySize,
                      double** dzeta, double** psi, double** u, double** v,
                      double Re, double dx, double dy, double dt);

void BorderValuesCalculator_calcBorderValues(
                            struct BorderValuesCalculator *self,
                            size_t xSize, size_t ySize,
                            double** dzeta, double** psi, double** u, double** v,
                            double Re, double dx, double dy, double dt);


void printArr(double** ar, size_t xSize, size_t ySize);

//nabla*nabla(arr(x, y)) = f(x, y)
//this method can calculate only inner cells!
void ellipticExplicitSolve(double **arr, double **f,
                          size_t xSize, size_t ySize, double dx, double dy, double maxErr);

//this method can calculate only inner cells!
//i'm not sure about this statement, but I don't need border cells anyway
void parabolicUpwindSolve(double **arr, double **ddarrxx, double **ddarryy,
                         size_t xSize, size_t ySize, double dx, double dy, double dt, double Re);

void calcNextIter(int8_t** maskForObstacle, size_t xSize, size_t ySize,
                  double** dzeta, double** psi, double** u, double** v,
                  double Re, double dx, double dy, double dt, double maxErr,
                  struct BorderValuesCalculator *borderCalc);


#endif