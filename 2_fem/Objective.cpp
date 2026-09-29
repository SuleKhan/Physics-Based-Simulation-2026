#include "FEMApp.h"
#include "Optimization.h"

bool FEMApp::step()
{
    Optimization optimization;
    optimization.optimizer = Optimization::NEWTON;
    optimization.objective_function = [&](const VectorXF &y, double &E)
    {
        return objectiveFunction(y, E);
    };
    optimization.gradient_function = [&](const VectorXF &y, double &E, VectorXF &g)
    {
        objectiveFunction(y, E);
        return objectiveGradient(y, g);
    };
    optimization.hessian_function = [&](const VectorXF &y, double &E, VectorXF &g, TripletListF &H)
    {
        objectiveFunction(y, E);
        objectiveGradient(y, g);
        return objectiveHessian(y, H);
    };

    VectorXF x = x_prev;
    Optimization::OptimizationStatus status = optimization.step(x);
    if (status == Optimization::SUCCESS)
    {
        // update position and velocity
        v = (x - x_prev) / m_dt;
        x_prev = x;
    }
    else if (status == Optimization::CONVERGED)
    {
        std::cout << "FEMApp::step: optimization converged." << std::endl;
    }
    else
    {
        // optimization step failed, keep previous state
        std::cerr << "FEMApp::step: optimization step failed." << std::endl;
        m_simulating = false;
    }

    return status == Optimization::SUCCESS;
}

bool FEMApp::objectiveFunction(const VectorXF &x, double &E)
{
    E = 0;
    addSystemEnergy(x, E);
    if (use_dynamics)
    {
        addDynamicsPotential(x, E);
    }
    return true;
}

bool FEMApp::objectiveGradient(const VectorXF &x, VectorXF &gradient)
{
    gradient = VectorXF::Zero(x.size());
    addSystemGradient(x, gradient);
    if (use_dynamics)
    {
        addDynamicsGradient(x, gradient);
    }
    return true;
}

bool FEMApp::objectiveHessian(const VectorXF &x, TripletListF &hessian_triplets)
{
    hessian_triplets.clear();
    addSystemHessian(x, hessian_triplets);
    if (use_dynamics)
    {
        addDynamicsHessian(x, hessian_triplets);
    }
    return true;
}

bool FEMApp::checkGradient()
{
    VectorXF grad, grad_fd;
    objectiveGradient(x_prev, grad);
    computeGradientFD(x_prev, grad_fd, 1e-6);
    double error = (grad - grad_fd).cwiseAbs().maxCoeff();

    bool check = error < 1e-6;
    if (check)
    {
        std::cout << "FEMApp::checkGradient: gradient check passed." << std::endl;
    }
    else
    {
        std::cout << "FEMApp::checkGradient: gradient check failed." << std::endl;
    }
    std::cout << "Max error: " << error << std::endl;

    return check;
}

bool FEMApp::checkHessian()
{
    MatrixXF H, H_fd;
    TripletListF hessian_triplets;
    objectiveHessian(x_prev, hessian_triplets);
    int n = x_prev.size();
    H = MatrixXF::Zero(n, n);
    for (const auto &triplet : hessian_triplets)
    {
        H(triplet.row(), triplet.col()) += triplet.value();
    }
    computeHessianFD(x_prev, H_fd, 1e-6);
    double error = (H - H_fd).cwiseAbs().maxCoeff();

    bool check = error < 1e-3;
    if (check)
    {
        std::cout << "FEMApp::checkHessian: Hessian check passed." << std::endl;
    }
    else
    {
        std::cout << "FEMApp::checkHessian: Hessian check failed." << std::endl;
    }
    std::cout << "Max error: " << error << std::endl;

    return check;
}