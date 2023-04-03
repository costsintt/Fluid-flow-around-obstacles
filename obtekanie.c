#include <stdlib.h>
#include <stdint.h>

#include <stdio.h>
#include <math.h>

// адвекция диффузия решается явной схемой с разностями против потока
//

size_t xSize = 7;
size_t ySize = 9;
double velOfMovementBorders = 0.2;
double dt = 0.01;
double dx = 0.1;
double dy = 0.1;
double maxErr = 0.00001;

//remove
void print(double** ar)
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

//remove
void init(double** dzeta, double** psi, double** u, double** v)
{
    for(size_t j = 0; j < xSize; j++)
        dzeta[0][j] = 2 * (psi[1][j] - psi[0][j] + velOfMovementBorders * dy) / dy / dy;
    for(size_t j = 0; j < xSize; j++)
        dzeta[ySize-1][j] = 2 * (psi[ySize-2][j] - psi[ySize-1][j] + velOfMovementBorders * dy) / dy / dy;

    for(size_t j = 0; j < xSize; j++)
        psi[0][j] = psi[1][j] + velOfMovementBorders * dy;
    for(size_t j = 0; j < xSize; j++)
        psi[ySize-1][j] = psi[ySize-2][j] + velOfMovementBorders * dy;

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
}

int main()
{ 
    double tMax = 0.1;
    double tCurr = 0.0;
    double Re = 5;

    double ddtx =  dt / dx;
    double ddty = dt / dy;

    //remove
    //maskForObstacle: 
    //  1 in the cell means there is obstacle in the cell
    //  0 in the cell means there is no obstacle in the cell
    int8_t** maskForObstacle = calloc(ySize, sizeof(int8_t*));
    for(size_t i = 0; i < ySize; i++) maskForObstacle[i] = calloc(xSize, sizeof(int8_t));

    //remove
    double** dzeta = calloc(ySize, sizeof(double*)); //the swirls
    for(size_t i = 0; i < ySize; i++) dzeta[i] = calloc(xSize, sizeof(double));

    //remove
    double** psi = calloc(ySize, sizeof(double*)); //the flows
    for(size_t i = 0; i < ySize; i++) psi[i] = calloc(xSize, sizeof(double));

    //remove
    double** u = calloc(ySize, sizeof(double*)); //the x-components of the velocities
    for(size_t i = 0; i < ySize; i++) u[i] = calloc(xSize, sizeof(double));

    //remove
    double** v = calloc(ySize, sizeof(double*)); //the y-components of the velocities
    for(size_t i = 0; i < ySize; i++) v[i] = calloc(xSize, sizeof(double));

    //remove //init dzeta, psi, u, v
    init(dzeta, psi, u, v);

    //main loop
    for(; tCurr < tMax; tCurr += dt)
    {
        //inner swirl cells
        for(size_t i = 1; i < ySize - 1; i++)
        {
            for(size_t j = 1; j < xSize - 1; j++)
            {
                double a = u[i][j] > 0 ? 0 : 1;
                double b = v[i][j] > 0 ? 0 : 1;

                dzeta[i][j] = dzeta[i][j]
                - ddtx * u[i][j] * (1 - a) * (dzeta[i][j] - dzeta[i-1][j])
                - ddtx * a * u[i][j] * (dzeta[i+1][j] - dzeta[i][j])
                - ddty * v[i][j] * (1 - b) * (dzeta[i][j] - dzeta[i][j-1])
                - ddty * b * v[i][j] * (dzeta[i][j+1] - dzeta[i][j])
                + dt / Re *
                    (   (dzeta[i+1][j] - 2.0 * dzeta[i][j] + dzeta[i-1][j]) / dx / dx
                      + (dzeta[i][j+1] - 2.0 * dzeta[i][j] + dzeta[i][j-1]) / dy / dy
                    );
            }
        }

        //inner flow cells
        for(size_t n = 0; n < 1000; n++)
        {
            double err;
            double maxErrLocal = 0.0;
            for(size_t i = 1; i < ySize - 1; i++)
            {
                for(size_t j = 1; j < xSize - 1; j++)
                {
                    err = psi[i][j];
                    psi[i][j] =
                    (
                        dx * dx * (psi[i][j+1] + psi[i][j-1])
                        + dy * dy * (psi[i+1][j] + psi[i-1][j])
                        - dx * dx * dy * dy * dzeta[i][j]
                    ) / (2 * (dx * dx + dy * dy));
                    err -= psi[i][j];
                    err = fabs(err);
                    if(maxErrLocal < err) maxErrLocal = err;
                }
            }
            if(maxErrLocal < maxErr) break;
        }
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

    return 0;
}