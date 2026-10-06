#include "../FEMApp.h"

void FEMApp::addSystemEnergy(const VectorXF &x, double &E)
{
    int num_pins = pins.size();
    int num_elements = elements.size();

    for (int i = 0; i < num_pins; ++i)
    {
        addPinEnergy(x, i, E);
    }

    for (int i = 0; i < num_elements; ++i)
    {
        addTriangleEnergy(x, i, E);
    }

    if (use_gravity)
    {
        addGravityEnergy(x, E);
    }
}

void FEMApp::addSystemGradient(const VectorXF &x, VectorXF &gradient)
{
    int num_pins = pins.size();
    int num_elements = elements.size();

    for (int i = 0; i < num_pins; ++i)
    {
        addPinGradient(x, i, gradient);
    }

    for (int i = 0; i < num_elements; ++i)
    {
        addTriangleGradient(x, i, gradient);
    }

    if (use_gravity)
    {
        addGravityGradient(x, gradient);
    }
}

void FEMApp::addSystemHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    int num_pins = pins.size();
    int num_elements = elements.size();

    for (int i = 0; i < num_pins; ++i)
    {
        addPinHessian(x, i, hessian_triplets);
    }

    for (int i = 0; i < num_elements; ++i)
    {
        addTriangleHessian(x, i, hessian_triplets);
    }

    if (use_gravity)
    {
        addGravityHessian(x, hessian_triplets);
    }
}