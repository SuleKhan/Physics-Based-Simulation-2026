#include "../FEMApp.h"

void FEMApp::addGravityEnergy(const VectorXF &x, double &E)
{
    /// TODO: Implement gravitational potential energy.
    /// Make sure not to overwrite E, but add to it!
    E += 0.0;
}

void FEMApp::addGravityGradient(const VectorXF &x, VectorXF &gradient)
{
    /// TODO: Implement gradient of gravitational potential energy.
    /// Make sure not to overwrite gradient, but add to it!
    gradient += VectorXF::Zero(x.size());
}

void FEMApp::addGravityHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    /// TODO: Implement Hessian of gravitational potential energy.
    /// Add hessian entries using hessian_triplets.emplace_back(row, col, value);
    hessian_triplets.emplace_back(0, 0, 0.0);
}