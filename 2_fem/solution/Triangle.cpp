#include "../FEMApp.h"

void FEMApp::initElement(int triangle_index)
{
    Element &elem = elements[triangle_index];
    int v0 = elem.vertices(0);
    int v1 = elem.vertices(1);
    int v2 = elem.vertices(2);

    Vector2F p0 = X.segment<2>(2 * v0);
    Vector2F p1 = X.segment<2>(2 * v1);
    Vector2F p2 = X.segment<2>(2 * v2);

    Matrix2F DX;
    DX.col(0) = p1 - p0;
    DX.col(1) = p2 - p0;

    elem.DX_inv = DX.inverse();
    elem.area = 0.5 * ((p1 - p0).x() * (p2 - p0).y() - (p1 - p0).y() * (p2 - p0).x());
}

Matrix2F FEMApp::deformationGradient(const VectorXF &x, int triangle_index)
{
    const Element &elem = elements[triangle_index];
    int v0 = elem.vertices(0);
    int v1 = elem.vertices(1);
    int v2 = elem.vertices(2);

    Vector2F p0 = x.segment<2>(2 * v0);
    Vector2F p1 = x.segment<2>(2 * v1);
    Vector2F p2 = x.segment<2>(2 * v2);

    Matrix2F Ds;
    Ds.col(0) = p1 - p0;
    Ds.col(1) = p2 - p0;

    Matrix2F F_mat = Ds * elem.DX_inv;
    return F_mat;
}

void FEMApp::addTriangleEnergy(const VectorXF &x, int triangle_index, double &E)
{
    Matrix2F F = deformationGradient(x, triangle_index);
    double J = F.determinant();
    double mu = shearModulus;
    double kappa = bulkModulus;

    double energy_density = 0.5 * mu * (F.squaredNorm() - 2) - mu * std::log(J) + 0.5 * kappa * std::log(J) * std::log(J);
    double rest_area = elements[triangle_index].area;
    E += energy_density * rest_area;
}
