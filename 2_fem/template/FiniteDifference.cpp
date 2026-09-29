#include "../FEMApp.h"

void FEMApp::computeGradientFD(const VectorXF &x, VectorXF &gradient, double eps)
{
    // TODO: Complete the finite difference approximation of the gradient.
    // grad_FD(i) = (f(x + eps * e_i)-f(x - eps * e_i))/(2 * eps);
    // Use the function objectiveFunction to compute f(x).
    gradient = VectorXF::Zero(x.size());
}

void FEMApp::computeHessianFD(const VectorXF &x, MatrixXF &H, double eps)
{
    // TODO: Complete the finite difference computation for the Hessian.
    // hess_FD(i,j) = (grad(x + eps * e_j)(i) - grad(x - eps * e_j)(i)) / (2 * eps);
    // You can use computeGradientFD or objectiveGradient to compute gradients.
    H = MatrixXF::Zero(x.size(), x.size());
}