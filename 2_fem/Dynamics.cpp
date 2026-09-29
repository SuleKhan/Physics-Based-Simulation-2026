#include "FEMApp.h"

void FEMApp::addDynamicsPotential(const VectorXF &x, double &E)
{
    for (int i = 0; i < x.size(); ++i)
    {
        double a = 1.0 / pow(m_dt, 2) * (x(i) - x_prev(i) - m_dt * v(i));
        E += 0.5 * pow(m_dt, 2) * a * M(i) * a;
    }
}

void FEMApp::addDynamicsGradient(const VectorXF &x, VectorXF &gradient)
{
    for (int i = 0; i < x.size(); ++i)
    {
        double a = 1.0 / pow(m_dt, 2) * (x(i) - x_prev(i) - m_dt * v(i));
        gradient(i) += a * M(i);
    }
}

void FEMApp::addDynamicsHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    for (int i = 0; i < x.size(); ++i)
    {
        double k = 1.0 / pow(m_dt, 2) * M(i);
        hessian_triplets.emplace_back(i, i, k);
    }
}