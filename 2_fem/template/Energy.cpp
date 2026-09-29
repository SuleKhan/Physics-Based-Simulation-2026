#include "../FEMApp.h"

void FEMApp::addSystemEnergy(const VectorXF &x, double &E)
{
    /// TODO: Compute total potential energy of the system.
    /// HINT: Use addPinEnergy, addTriangleEnergy, and addGravityEnergy.
}

void FEMApp::addSystemGradient(const VectorXF &x, VectorXF &gradient)
{
    /// TODO: Compute gradient of total potential energy.
}

void FEMApp::addSystemHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    /// TODO: Compute Hessian of total potential energy.
}