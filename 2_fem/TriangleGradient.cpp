#include "FEMApp.h"

void FEMApp::addTriangleGradient(const VectorXF &x, int triangle_index, VectorXF &gradient)
{
    Matrix2F F = deformationGradient(x, triangle_index);
    double J = F.determinant();
    double mu = shearModulus;
    double kappa = bulkModulus;

    Matrix2F dEdF = mu * F + F.inverse().transpose() * (-mu + kappa * std::log(J));
    Matrix2F m = dEdF * elements[triangle_index].DX_inv.transpose() * elements[triangle_index].area;

    int v0 = elements[triangle_index].vertices(0);
    int v1 = elements[triangle_index].vertices(1);
    int v2 = elements[triangle_index].vertices(2);

    gradient.segment<2>(2 * v0) += -m.col(0) - m.col(1);
    gradient.segment<2>(2 * v1) += m.col(0);
    gradient.segment<2>(2 * v2) += m.col(1);
}

void FEMApp::addTriangleHessian(const VectorXF &x, int triangle_index, TripletListF &hessian_triplets)
{
    Matrix2F F = deformationGradient(x, triangle_index);
    Matrix2F Finv = F.inverse();
    Matrix2F FinvT = Finv.transpose();
    double J = F.determinant();
    double mu = shearModulus;
    double kappa = bulkModulus;

    Matrix2F dXInv = elements[triangle_index].DX_inv;
    double restShapeArea = elements[triangle_index].area;

    Matrix2F dF, dP, tmpM, dH;
    const double dDs[6][4] = {{-1, -1, 0, 0}, {0, 0, -1, -1}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};

    int v0 = elements[triangle_index].vertices(0);
    int v1 = elements[triangle_index].vertices(1);
    int v2 = elements[triangle_index].vertices(2);
    int x_idx[6] = {2 * v0, 2 * v0 + 1, 2 * v1, 2 * v1 + 1, 2 * v2, 2 * v2 + 1};

    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 4; ++j)
            dF(j / 2, j % 2) = dDs[i][j];
        dF = dF * dXInv;
        dP = mu * dF + (mu - kappa * log(J)) * FinvT * dF.transpose() * FinvT;
        tmpM = Finv * dF;
        dP = dP + kappa * (tmpM(0, 0) + tmpM(1, 1)) * FinvT;
        dH = restShapeArea * dP * dXInv.transpose();

        hessian_triplets.emplace_back(x_idx[i], 2 * v0, -dH(0, 0) - dH(0, 1));
        hessian_triplets.emplace_back(x_idx[i], 2 * v0 + 1, -dH(1, 0) - dH(1, 1));
        hessian_triplets.emplace_back(x_idx[i], 2 * v1, dH(0, 0));
        hessian_triplets.emplace_back(x_idx[i], 2 * v1 + 1, dH(1, 0));
        hessian_triplets.emplace_back(x_idx[i], 2 * v2, dH(0, 1));
        hessian_triplets.emplace_back(x_idx[i], 2 * v2 + 1, dH(1, 1));
    }
}