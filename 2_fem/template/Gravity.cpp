#include "../FEMApp.h"

void FEMApp::addGravityEnergy(const VectorXF& x, double& E)
{
    for (int i = 0; i < x.size() / 2; ++i)
    {
        Vector2F xi = x.segment<2>(2 * i);
        double mass = M(2 * i);

        E += -mass * gravity.dot(xi);
    }
}

void FEMApp::addGravityGradient(const VectorXF& x, VectorXF& gradient)
{
    for (int i = 0; i < x.size() / 2; ++i)
    {
        double mass = M(2 * i);

        gradient.segment<2>(2 * i) += -mass * gravity;
    }
}

void FEMApp::addGravityHessian(const VectorXF& x, TripletListF& hessian_triplets)
{
    (void)x;
    (void)hessian_triplets;
}