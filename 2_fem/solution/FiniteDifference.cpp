#include "../FEMApp.h"

void FEMApp::computeGradientFD(const VectorXF &x, VectorXF &gradient, double eps)
{
    int n = x.size();
    gradient = VectorXF::Zero(n);
    double E_plus, E_minus;
    VectorXF x_plus, x_minus;
    for (int i = 0; i < n; ++i)
    {
        x_plus = x;
        x_minus = x;
        x_plus(i) += eps;
        x_minus(i) -= eps;

        E_plus = 0;
        objectiveFunction(x_plus, E_plus);
        E_minus = 0;
        objectiveFunction(x_minus, E_minus);
        gradient(i) = (E_plus - E_minus) / (2 * eps);
    }
}

void FEMApp::computeHessianFD(const VectorXF &x, MatrixXF &H, double eps)
{
    int n = x.size();
    H = MatrixXF::Zero(n, n);
    VectorXF grad_plus, grad_minus;
    VectorXF x_plus, x_minus;
    for (int i = 0; i < n; ++i)
    {
        x_plus = x;
        x_minus = x;
        x_plus(i) += eps;
        x_minus(i) -= eps;

        objectiveGradient(x_plus, grad_plus);
        objectiveGradient(x_minus, grad_minus);
        H.col(i) = (grad_plus - grad_minus) / (2 * eps);
    }
}