#include "../FEMApp.h"

void FEMApp::addPinEnergy(const VectorXF &x, int pin_index, double &E)
{
    const Pin &pin = pins[pin_index];
    Vector2F pos = x.segment<2>(2 * pin.node);
    Vector2F diff = pos - pin.pos;
    E += 0.5 * pin.k * diff.squaredNorm();
}

void FEMApp::addPinGradient(const VectorXF &x, int pin_index, VectorXF &gradient)
{
    const Pin &pin = pins[pin_index];
    Vector2F pos = x.segment<2>(2 * pin.node);
    Vector2F diff = pos - pin.pos;
    gradient.segment<2>(2 * pin.node) += pin.k * diff;
}

void FEMApp::addPinHessian(const VectorXF &x, int pin_index, TripletListF &hessian_triplets)
{
    const Pin &pin = pins[pin_index];
    hessian_triplets.emplace_back(2 * pin.node, 2 * pin.node, pin.k);
    hessian_triplets.emplace_back(2 * pin.node + 1, 2 * pin.node + 1, pin.k);
}