#include "../FEMApp.h"

void FEMApp::addSystemEnergy(const VectorXF& x, double& E)
{
    for (int i = 0; i < pins.size(); ++i)
    {
        addPinEnergy(x, i, E);
    }

    for (int i = 0; i < elements.size(); ++i)
    {
        addTriangleEnergy(x, i, E);
    }

    addGravityEnergy(x, E);
}

void FEMApp::addSystemGradient(const VectorXF& x, VectorXF& gradient)
{
    for (int i = 0; i < pins.size(); ++i)
    {
        addPinGradient(x, i, gradient);
    }

    for (int i = 0; i < elements.size(); ++i)
    {
        addTriangleGradient(x, i, gradient);
    }

    addGravityGradient(x, gradient);
}

void FEMApp::addSystemHessian(const VectorXF& x, TripletListF& hessian_triplets)
{
    for (int i = 0; i < pins.size(); ++i)
    {
        addPinHessian(x, i, hessian_triplets);
    }

    for (int i = 0; i < elements.size(); ++i)
    {
        addTriangleHessian(x, i, hessian_triplets);
    }

    addGravityHessian(x, hessian_triplets);
}