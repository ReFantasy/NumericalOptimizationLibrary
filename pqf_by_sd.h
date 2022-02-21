/**
 *
 * Solving positive definite quadratic function by steepest descent method
 * f(x) = 1/2 * x * G * x + b * x + c;
 *
 * Author  : ReFantasy
 * Date    : 2022-02-19
 */

#ifndef __PQF_BY_SD_H__
#define __PQF_BY_SD_H__
#include "Eigen/Dense"
#include "global.h"

FLOAT pqf(Vector x, Matrix G, Vector b, FLOAT c);

Vector SteepestDescent(Matrix G, Vector b, FLOAT c, Vector x0, FLOAT gk_norm = 10e-5);

#endif //__QPDF_BY_SD_H__
