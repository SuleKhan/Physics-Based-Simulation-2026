#include "../FEMApp.h"

void FEMApp::addPinEnergy(const VectorXF &x, int pin_index, double &E)
{
    const Pin &pin = pins[pin_index];

    // suppress unused variable warnings until you implement this
    (void)pin;

    /// TODO: Implement pin spring potential energy.
    /// Make sure not to overwrite E, but add to it!
    E += 0.0;
}

void FEMApp::addPinGradient(const VectorXF &x, int pin_index, VectorXF &gradient)
{
    const Pin &pin = pins[pin_index];

    // suppress unused variable warnings until you implement this
    (void)pin;

    /// TODO: Implement gradient of pin spring potential energy.
    /// Make sure not to overwrite gradient, but add to it!
    gradient += VectorXF::Zero(x.size());
}

void FEMApp::addPinHessian(const VectorXF &x, int pin_index, TripletListF &hessian_triplets)
{
    const Pin &pin = pins[pin_index];

    // suppress unused variable warnings until you implement this
    (void)pin;

    /// TODO: Implement Hessian of pin spring potential energy.
    /// Add hessian entries using hessian_triplets.emplace_back(row, col, value);
    hessian_triplets.emplace_back(0, 0, 0.0);
}