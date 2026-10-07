#include "../FEMApp.h"

void FEMApp::initElement(int triangle_index)
{
    Element &elem = elements[triangle_index];

    /// TODO: Compute rest shape area and inverse of rest shape matrix
    elem.area = 0.0;
    elem.DX_inv = Matrix2F::Zero();

    int v0 = elem.vertices(0);
    int v1 = elem.vertices(1);
    int v2 = elem.vertices(2);

    Vector2F X0 = X.segment<2>(2 * v0);
    Vector2F X1 = X.segment<2>(2 * v1);
    Vector2F X2 = X.segment<2>(2 * v2);

    Matrix2F DX;
    DX.col(0) = X1 - X0;
    DX.col(1) = X2 - X0;

    elem.area = 0.5 * std::abs(DX.determinant());
    elem.DX_inv = DX.inverse();
}

Matrix2F FEMApp::deformationGradient(const VectorXF &x, int triangle_index)
{
    const Element &elem = elements[triangle_index];

    // suppress unused variable warnings until you implement this
    (void)elem;

    /// TODO: Compute deformation gradient

    int v0 = elem.vertices(0);
    int v1 = elem.vertices(1);
    int v2 = elem.vertices(2);

    Vector2F x0 = x.segment<2>(2 * v0);
    Vector2F x1 = x.segment<2>(2 * v1);
    Vector2F x2 = x.segment<2>(2 * v2);

    Matrix2F dx;
    dx.col(0) = x1 - x0;
    dx.col(1) = x2 - x0;

    Matrix2F F_mat = dx * elem.DX_inv;

    return F_mat;
}

void FEMApp::addTriangleEnergy(const VectorXF &x, int triangle_index, double &E)
{
    const Element& elem = elements[triangle_index];

    Matrix2F F = deformationGradient(x, triangle_index);
    double J = F.determinant();
    double mu = shearModulus;
    double kappa = bulkModulus;

    double logJ = std::log(J);

    double energy_density =
        0.5 * mu * (F.squaredNorm() - 2.0)
        - mu * logJ
        + 0.5 * kappa * logJ * logJ;

    double energy = elem.area * energy_density;

    E += energy;
}