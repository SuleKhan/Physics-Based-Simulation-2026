#include "../FEMApp.h"

void FEMApp::computeGradientFD(const VectorXF& x, VectorXF& gradient, double eps)
{
    gradient = VectorXF::Zero(x.size());

    for (int i = 0; i < x.size(); ++i)
    {
        VectorXF x_plus = x;
        VectorXF x_minus = x;

        x_plus(i) += eps;
        x_minus(i) -= eps;

        double E_plus = 0.0;
        double E_minus = 0.0;

        objectiveFunction(x_plus, E_plus);
        objectiveFunction(x_minus, E_minus);

        gradient(i) = (E_plus - E_minus) / (2.0 * eps);
    }
}

void FEMApp::computeHessianFD(const VectorXF& x, MatrixXF& H, double eps)
{
    H = MatrixXF::Zero(x.size(), x.size());

    for (int j = 0; j < x.size(); ++j)
    {
        VectorXF x_plus = x;
        VectorXF x_minus = x;

        x_plus(j) += eps;
        x_minus(j) -= eps;

        VectorXF grad_plus;
        VectorXF grad_minus;

        computeGradientFD(x_plus, grad_plus, eps);
        computeGradientFD(x_minus, grad_minus, eps);

        for (int i = 0; i < x.size(); ++i)
        {
            H(i, j) = (grad_plus(i) - grad_minus(i)) / (2.0 * eps);
        }
    }
}