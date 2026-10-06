#include "../FEMApp.h"

void FEMApp::addGravityEnergy(const VectorXF &x, double &E)
{
    int num_nodes = x.rows() / 2;

    for (int i = 0; i < num_nodes; ++i)
    {
        Vector2F pos = x.segment<2>(2 * i);
        E -= M(2 * i) * gravity(0) * pos(0);
        E -= M(2 * i + 1) * gravity(1) * pos(1);
    }
}

void FEMApp::addGravityGradient(const VectorXF &x, VectorXF &gradient)
{
    int num_nodes = x.rows() / 2;

    for (int i = 0; i < num_nodes; ++i)
    {
        gradient(2 * i) -= M(2 * i) * gravity(0);
        gradient(2 * i + 1) -= M(2 * i + 1) * gravity(1);
    }
}

void FEMApp::addGravityHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    // Gravity does not contribute to the Hessian since it is linear in position.
}