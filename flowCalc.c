#include "flowCalc.h"

void print(double** ar, size_t xSize, size_t ySize)
{
    for(size_t i = 0; i < ySize; i++)
    {
        for(size_t j = 0; j < xSize; j++)
        {
            printf("%lf\t", ar[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void calcBorderValues(double velOfMovementBorders, size_t xSize, size_t ySize,
                      double** dzeta, double** psi, double** u, double** v,
                      double Re, double dx, double dy, double dt)
{
    //border flow cells
    for(size_t j = 0; j < xSize; j++)
        psi[0][j] = psi[1][j] + velOfMovementBorders * dy;
    for(size_t j = 0; j < xSize; j++)
        psi[ySize-1][j] = psi[ySize-2][j] + velOfMovementBorders * dy;
    for(size_t i = 0; i < ySize; i++)
        psi[i][0] = psi[i][1];
    for(size_t i = 0; i < ySize; i++)
        psi[i][xSize-1] = psi[i][xSize-2];

    //border swirl cell
    for(size_t j = 0; j < xSize; j++)
        dzeta[0][j] = 2 * (psi[1][j] - psi[0][j] + velOfMovementBorders * dy) / dy / dy;
    for(size_t j = 0; j < xSize; j++)
        dzeta[ySize-1][j] = 2 * (psi[ySize-2][j] - psi[ySize-1][j] + velOfMovementBorders * dy) / dy / dy;
    for(size_t i = 0; i < ySize; i++)
        dzeta[i][0] = dzeta[i][1];
    for(size_t i = 0; i < ySize; i++)
        dzeta[i][xSize-1] = dzeta[i][xSize-2];


    //border velocities
    for(size_t j = 0; j < xSize; j++)
    {
        u[0][j] = velOfMovementBorders;
        v[0][j] = 0.0;
    }
    for(size_t j = 0; j < xSize; j++)
    {
        u[ySize-1][j] = velOfMovementBorders;
        v[ySize-1][j] = 0.0;
    }
    for(size_t i = 0; i < ySize; i++)
    {
        u[i][0] = u[i][1];
        v[i][0] = v[i][1];
    }
    for(size_t i = 0; i < ySize; i++)
    {
        u[i][xSize-1] = u[i][xSize-2];
        v[i][xSize-1] = v[i][xSize-2];
    }
}

void BorderValuesCalculator_calcBorderValues(
                            struct BorderValuesCalculator *self,
                            size_t xSize, size_t ySize,
                            double** dzeta, double** psi, double** u, double** v,
                            double Re, double dx, double dy, double dt)
{
    calcBorderValues(self->bordersVel, xSize, ySize, dzeta, psi, u, v, Re, dx, dy, dt);
}

void ellipticExplicitSolve(double **arr, double **f,
                          size_t xSize, size_t ySize, double dx, double dy, double maxErr)
{
    double err;
    double maxErrLocal;
    for(size_t n = 0; n < 1000; n++)
    {
        maxErrLocal = 0.0;
        for(size_t i = 1; i < ySize - 1; i++)
        {
            for(size_t j = 1; j < xSize - 1; j++)
            {
                err = arr[i][j];
                arr[i][j] =
                (
                    dx * dx * (arr[i][j+1] + arr[i][j-1])
                    + dy * dy * (arr[i+1][j] + arr[i-1][j])
                    - dx * dx * dy * dy * f[i][j]
                ) / (2 * (dx * dx + dy * dy));
                err -= arr[i][j];
                err = fabs(err);
                if(maxErrLocal < err) maxErrLocal = err;
            }
        }
        if(maxErrLocal < maxErr) break;
    }
}

void parabolicUpwindSolve(double **arr, double **ddarrxx, double **ddarryy,
                         size_t xSize, size_t ySize, double dx, double dy, double dt, double Re)
{
    double ddtx =  dt / dx;
    double ddty = dt / dy;
    for(size_t i = 1; i < ySize - 1; i++)
    {
        for(size_t j = 1; j < xSize - 1; j++)
        {
            double a = ddarrxx[i][j] > 0 ? 0 : 1;
            double b = ddarryy[i][j] > 0 ? 0 : 1;

            arr[i][j] = arr[i][j]
            - ddtx * ddarrxx[i][j] * (1 - a) * (arr[i][j] - arr[i-1][j])
            - ddtx * a * ddarrxx[i][j] * (arr[i+1][j] - arr[i][j])
            - ddty * ddarryy[i][j] * (1 - b) * (arr[i][j] - arr[i][j-1])
            - ddty * b * ddarryy[i][j] * (arr[i][j+1] - arr[i][j])
            + dt / Re *
            (   (arr[i+1][j] - 2.0 * arr[i][j] + arr[i-1][j]) / dx / dx
              + (arr[i][j+1] - 2.0 * arr[i][j] + arr[i][j-1]) / dy / dy
            );
        }
    }
}

void calcNextIter(int8_t** maskForObstacle, size_t xSize, size_t ySize,
                  double** dzeta, double** psi, double** u, double** v,
                  double Re, double dx, double dy, double dt, double maxErr,
                  struct BorderValuesCalculator *borderCalc)
{
    //inner swirl cells
    parabolicUpwindSolve(dzeta, u, v, xSize, ySize, dx, dy, dt, Re);    

    //inner flow cells
    ellipticExplicitSolve(psi, dzeta, xSize, ySize, dx, dy, maxErr);

    borderCalc->calc(borderCalc, xSize, ySize, dzeta, psi, u, v, Re, dx, dy, dt);

    //inner velocities
    for(size_t i = 1; i < ySize - 1; i++)
    {
        for(size_t j = 1; j < xSize - 1; j++)
        {
            u[i][j] =   (psi[i+1][j] - psi[i-1][j]) / dy / 2;
            v[i][j] = - (psi[i][j+1] - psi[i][j-1]) / dx / 2;
        }
    }
}