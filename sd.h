/**
 *
 * Solving positive definite quadratic function by steepest descent method
 * f(x) = 1/2 * x * G * x + b * x + c;
 *
 * Author  : ReFantasy
 * Date    : 2022-02-19
 */

#ifndef __SD_H__
#define __SD_H__
#include "Eigen/Dense"
#include "global.h"

Vector SteepestDescent(TargetFunctor &fucntor, Vector x0, FLOAT gk_norm = 10e-5);

#endif //__SD_H__