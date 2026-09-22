#pragma once

#include "VecMatDef.h"

class Optimization
{
public:
    enum Optimizer
    {
        GRADIENT_DESCENT, //
        NEWTON
    };

    enum OptimizationStatus
    {
        CONVERGED, //
        SUCCESS,
        FAILURE
    };

public:
    Optimizer optimizer = GRADIENT_DESCENT;
    int tolerance_exponent = -16;

public:
    /// Gradient and hessian functions also compute objective value and lower-order derivatives.
    std::function<bool(const VectorXF &, double &)> objective_function;
    std::function<bool(const VectorXF &, double &, VectorXF &)> gradient_function;
    std::function<bool(const VectorXF &, double &, VectorXF &, TripletListF &)> hessian_function;

    Optimization()
    {
        objective_function = [&](const VectorXF &y, double &energy)
        { return false; };
        gradient_function = [&](const VectorXF &y, double &energy, VectorXF &gradient)
        { return false; };
        hessian_function = [&](const VectorXF &y, double &energy, VectorXF &gradient, TripletListF &hessian)
        { return false; };
    }

public:
    /// Solve for x in Ax = b. Return true on success.
    bool linearSolve(const TripletListF &hessian, const VectorXF &b, VectorXF &x);

private:
    OptimizationStatus getDirectionGradientDescent(const VectorXF &y, VectorXF &dy, double &initial_objective_value);

    OptimizationStatus getDirectionNewton(const VectorXF &y, VectorXF &dy, double &initial_objective_value);

public:
    OptimizationStatus step(VectorXF &y);

    /// Find step along search direction dy that decreases the objective value. Update y accordingly.
    bool lineSearch(VectorXF &y, const VectorXF &dy, double initial_objective_value);

    bool lineSearch(VectorXF &y, const VectorXF &dy);

    /// Check gradient for scalar function of vector argument.
    static bool checkGradient(const VectorXF &y, VectorXF &error,
                              const std::function<bool(const VectorXF &, double &)> &func,
                              const std::function<bool(const VectorXF &, VectorXF &)> &grad_func, double epsilon,
                              int print_level = 0);

    /// Check hessian for scalar function of vector argument.
    static bool checkHessian(const VectorXF &y, MatrixXF &error,
                             const std::function<bool(const VectorXF &, VectorXF &)> &grad_func,
                             const std::function<bool(const VectorXF &, MatrixXF &)> &hess_func, double epsilon,
                             int print_level = 0);
};
