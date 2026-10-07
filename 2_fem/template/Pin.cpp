#include "../FEMApp.h"

void FEMApp::addPinEnergy(const VectorXF& x, int pin_index, double& E)
{
    const Pin& pin = pins[pin_index];

    /// TODO: Implement pin spring potential energy.
    /// Make sure not to overwrite E, but add to it!
    Vector2F p = x.segment<2>(2 * pin.node);
    Vector2F d = p - pin.pos;

    E += 0.5 * pin.k * d.squaredNorm();
}

void FEMApp::addPinGradient(const VectorXF& x, int pin_index, VectorXF& gradient)
{
    const Pin& pin = pins[pin_index];

    /// TODO: Implement gradient of pin spring potential energy.
    /// Make sure not to overwrite gradient, but add to it!

    Vector2F p = x.segment<2>(2 * pin.node);
    Vector2F d = p - pin.pos;

    gradient.segment<2>(2 * pin.node) += pin.k * d;
}

void FEMApp::addPinHessian(const VectorXF& x, int pin_index, TripletListF& hessian_triplets)
{
    const Pin& pin = pins[pin_index];

    // x is not needed for the pin Hessian
    (void)x;

    /// TODO: Implement Hessian of pin spring potential energy.
    /// Add hessian entries using hessian_triplets.emplace_back(row, col, value);

    int i = 2 * pin.node;

    hessian_triplets.emplace_back(i, i, pin.k);
    hessian_triplets.emplace_back(i + 1, i + 1, pin.k);
}