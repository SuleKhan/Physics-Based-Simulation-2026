#include "../FEMApp.h"

void FEMApp::initElement(int triangle_index)
{
    Element &elem = elements[triangle_index];

    /// TODO: Compute rest shape area and inverse of rest shape matrix
    elem.area = 0.0;
    elem.DX_inv = Matrix2F::Zero();
}

Matrix2F FEMApp::deformationGradient(const VectorXF &x, int triangle_index)
{
    const Element &elem = elements[triangle_index];

    // suppress unused variable warnings until you implement this
    (void)elem;

    /// TODO: Compute deformation gradient
    Matrix2F F_mat = Matrix2F::Zero();

    return F_mat;
}

void FEMApp::addTriangleEnergy(const VectorXF &x, int triangle_index, double &E)
{
    Matrix2F F = deformationGradient(x, triangle_index);
    double J = F.determinant();
    double mu = shearModulus;
    double kappa = bulkModulus;

    // suppress unused variable warnings until you implement this
    (void)J;
    (void)mu;
    (void)kappa;

    /// TODO: Compute the energy for this triangle element
    double energy = 0.0;

    E += energy;
}