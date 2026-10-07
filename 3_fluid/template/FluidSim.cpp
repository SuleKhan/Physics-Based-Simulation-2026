#include "../FluidApp.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

void FluidApp::solvePoisson()
{
    double dx2 = m_dx * m_dx;
    double residual = m_tolerance + 1; // initial residual
    double rho = 1;

    Array2d& p = p_pressure.x();

    for (int it = 0; it < m_max_iter; ++it)
    {
        // Gauss-Seidel sweep (in place, so new values are used as soon as available)
        for (int y = 1; y < m_res_y - 1; ++y)
        {
            for (int x = 1; x < m_res_x - 1; ++x)
            {
                double b = -p_divergence.x()(x, y) / m_dt * rho; // right-hand side

                p(x, y) = (p(x + 1, y) + p(x - 1, y) +
                    p(x, y + 1) + p(x, y - 1) -
                    dx2 * b) / 4.0;
            }
        }

        // Residual: squared L2-norm of (b - A p)
        residual = 0;
        for (int y = 1; y < m_res_y - 1; ++y)
        {
            for (int x = 1; x < m_res_x - 1; ++x)
            {
                double b = -p_divergence.x()(x, y) / m_dt * rho; // right-hand side

                double lhs = (p(x + 1, y) + p(x - 1, y) +
                    p(x, y + 1) + p(x, y - 1) -
                    4.0 * p(x, y)) / dx2;
                double cellResidual = b - lhs;

                residual += cellResidual * cellResidual;
            }
        }

        // L2-norm of the residual
        residual = sqrt(residual);

        // Average per grid cell
        residual /= (m_res_x - 2) * (m_res_y - 2);

        // cout << "Pressure solver: iter=" << it << ", res=" << residual << endl;
    }
}

void FluidApp::correctVelocity()
{
    double rho = 1.0;
    Array2d &p = p_pressure.x();
    Array2d &u = p_velocity.x();
    Array2d &v = p_velocity.y();

    // Note: velocity u_{i+1/2} is practically stored at i+1, hence xV_{i}  -= dt * (p_{i} - p_{i-1}) / dx
    // u stored at x means u_{x-1/2}
    for (int y = 1; y < m_res_y - 1; ++y)
        for (int x = 1; x < m_res_x; ++x)
            u(x, y) -= m_dt / rho * (p(x, y) - p(x - 1, y)) / m_dx;

    for (int y = 1; y < m_res_y; ++y)
        for (int x = 1; x < m_res_x - 1; ++x)
            v(x, y) -= m_dt / rho * (p(x, y) - p(x, y - 1)) / m_dx;
}

void FluidApp::advectValues()
{
    // store original values
    Array2d d(p_density.x());
    Array2d u(p_velocity.x());
    Array2d v(p_velocity.y());
    // move forward
    advectDensitySL(u, v);
    advectVelocitySL(u, v);

    if (m_use_maccormack)
    {
        // store forward advection
        Array2d d_forward(p_density.x());
        Array2d u_forward(p_velocity.x());
        Array2d v_forward(p_velocity.y());
        MacCormackUpdate(d, d_forward, u, u_forward, v, v_forward);
        MacCormackClamp(d, d_forward, u, u_forward, v, v_forward);
    }
}

void FluidApp::advectDensitySL(const Array2d &u, const Array2d &v)
{
    Array2d &d = p_density.x();
    Array2d d_tmp(d.size(0), d.size(1));

    // Densities, grid centers
    for (int y = 1; y < m_res_y - 1; ++y)
    {
        for (int x = 1; x < m_res_x - 1; ++x)
        {
            // TODO: Compute the velocity (at cell center)
            double last_x_velocity = 0.5 * (u(x, y) + u(x + 1, y));
            double last_y_velocity = 0.5 * (v(x, y) + v(x, y + 1));

            // TODO: Find the last position of the particle (in grid coordinates)
            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Densities are known between 1 and res-2
            if (last_x < 1)
                last_x = 1;
            if (last_y < 1)
                last_y = 1;
            if (last_x > m_res_x - 2)
                last_x = m_res_x - 2;
            if (last_y > m_res_y - 2)
                last_y = m_res_y - 2;

            // Determine corners for bilinear interpolation
            int x_low = (int)last_x;
            int y_low = (int)last_y;
            int x_high = x_low + 1;
            int y_high = y_low + 1;

            // Compute the interpolation weights
            double x_weight = last_x - x_low;
            double y_weight = last_y - y_low;

            // TODO: Bilinear interpolation
            // TODO: Bilinear interpolation
            d_tmp(x, y) =
                (1 - x_weight) * (1 - y_weight) * d(x_low, y_low) +
                x_weight * (1 - y_weight) * d(x_high, y_low) +
                (1 - x_weight) * y_weight * d(x_low, y_high) +
                x_weight * y_weight * d(x_high, y_high);
        }
    }

    // Copy the values in temp to the original buffers
    d = d_tmp;
}

void FluidApp::advectVelocitySL(const Array2d &u, const Array2d &v)
{
    Array2d &u_in = p_velocity.x();
    Array2d &v_in = p_velocity.y();

    Array2d u_tmp(u_in.size(0), u_in.size(1));
    Array2d v_tmp(v_in.size(0), v_in.size(1));

    // Velocities (u), MAC grid
    for (int y = 1; y < m_res_y - 1; ++y)
    {
        for (int x = 1; x < m_res_x; ++x)
        {
            // TODO: Compute the velocity
            double last_x_velocity = u(x, y);
            double last_y_velocity = 0.25 * (v(x - 1, y) + v(x, y) +
                v(x - 1, y + 1) + v(x, y + 1));

            // TODO: Find the last position of the particle (in grid coordinates)
            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Being conservative, one can say that the velocities are known between 1.5 and res-2.5
            // (the MAC grid is inside the known densities, which are between 1 and res - 2)
            if (last_x < 1.5)
                last_x = 1.5;
            if (last_y < 1.5)
                last_y = 1.5;
            if (last_x > m_res_x - 1.5)
                last_x = m_res_x - 1.5;
            if (last_y > m_res_y - 2.5)
                last_y = m_res_y - 2.5;

            // Determine corners for bilinear interpolation
            int x_low = (int)last_x;
            int y_low = (int)last_y;
            int x_high = x_low + 1;
            int y_high = y_low + 1;

            // Compute the interpolation weights
            double x_weight = last_x - x_low;
            double y_weight = last_y - y_low;

            // TODO: Bilinear interpolation
            u_tmp(x, y) =
                (1 - x_weight) * (1 - y_weight) * u_in(x_low, y_low) +
                x_weight * (1 - y_weight) * u_in(x_high, y_low) +
                (1 - x_weight) * y_weight * u_in(x_low, y_high) +
                x_weight * y_weight * u_in(x_high, y_high);
        }
    }

    // Velocities (v), MAC grid
    for (int y = 1; y < m_res_y; ++y)
    {
        for (int x = 1; x < m_res_x - 1; ++x)
        {
            // TODO: Compute the velocity
            double last_x_velocity = 0.25 * (u(x, y - 1) + u(x + 1, y - 1) +
                u(x, y) + u(x + 1, y));
            double last_y_velocity = v(x, y);

            // TODO: Find the last position of the particle (in grid coordinates)
            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Being conservative, one can say that the velocities are known between 1.5 and res-2.5
            // (the MAC grid is inside the known densities, which are between 1 and res - 2)
            if (last_x < 1.5)
                last_x = 1.5;
            if (last_y < 1.5)
                last_y = 1.5;
            if (last_x > m_res_x - 2.5)
                last_x = m_res_x - 2.5;
            if (last_y > m_res_y - 1.5)
                last_y = m_res_y - 1.5;

            // Determine corners for bilinear interpolation
            double x_low = (int)last_x;
            double y_low = (int)last_y;
            double x_high = x_low + 1;
            double y_high = y_low + 1;

            // Compute the interpolation weights
            double x_weight = last_x - x_low;
            double y_weight = last_y - y_low;

            // TODO: Bilinear interpolation
            v_tmp(x, y) =
                (1 - x_weight) * (1 - y_weight) * v_in(x_low, y_low) +
                x_weight * (1 - y_weight) * v_in(x_high, y_low) +
                (1 - x_weight) * y_weight * v_in(x_low, y_high) +
                x_weight * y_weight * v_in(x_high, y_high);
        }
    }

    // Copy the values in temp to the original buffers
    u_in = u_tmp;
    v_in = v_tmp;
}

void FluidApp::MacCormackUpdate(const Array2d &d, const Array2d &d_forward, const Array2d &u, const Array2d &u_forward, const Array2d &v, const Array2d &v_forward)
{
    // move backward
    m_dt *= -1;
    advectDensitySL(u, v);
    advectVelocitySL(u, v);

    // store backward advection
    Array2d d_backward(p_density.x());
    Array2d u_backward(p_velocity.x());
    Array2d v_backward(p_velocity.y());
    m_dt *= -1;

    Array2d d_tmp(d_forward);
    Array2d u_tmp(u_forward);
    Array2d v_tmp(v_forward);
    // MacCormack Update
    for (int y = 1; y < m_res_y - 1; ++y)
        for (int x = 1; x < m_res_x - 1; ++x)
            d_tmp(x, y) = d_forward(x, y) + 0.5 * (d(x, y) - d_backward(x, y));

    for (int y = 1; y < m_res_y - 1; ++y)
        for (int x = 1; x < m_res_x; ++x)
            u_tmp(x, y) = u_forward(x, y) + 0.5 * (u(x, y) - u_backward(x, y));

    for (int y = 1; y < m_res_y; ++y)
        for (int x = 1; x < m_res_x - 1; ++x)
            v_tmp(x, y) = v_forward(x, y) + 0.5 * (v(x, y) - v_backward(x, y));

    p_density.x() = d_tmp;
    p_velocity.x() = u_tmp;
    p_velocity.y() = v_tmp;
}

void FluidApp::MacCormackClamp(const Array2d &d, const Array2d &d_forward, const Array2d &u, const Array2d &u_forward, const Array2d &v, const Array2d &v_forward)
{
    Array2d d_tmp(p_density.x());
    Array2d u_tmp(p_velocity.x());
    Array2d v_tmp(p_velocity.y());
    // Clamp density
    for (int y = 1; y < m_res_y - 1; ++y)
    {
        for (int x = 1; x < m_res_x - 1; ++x)
        {
            double last_x_velocity = 0.5 * (u(x, y) + u(x + 1, y));
            double last_y_velocity = 0.5 * (v(x, y) + v(x, y + 1));

            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Densities are known between 1 and res-2
            if (last_x < 1)
                last_x = 1;
            if (last_y < 1)
                last_y = 1;
            if (last_x > m_res_x - 2)
                last_x = m_res_x - 2;
            if (last_y > m_res_y - 2)
                last_y = m_res_y - 2;

            // Determine corners for bilinear interpolation
            int x_low = (int)last_x;
            int y_low = (int)last_y;
            int x_high = x_low + 1;
            int y_high = y_low + 1;

            double d_min = 1e10;
            double d_max = -1e10;

            d_min = std::min(d(x_low, y_low), d_min);
            d_min = std::min(d(x_low, y_high), d_min);
            d_min = std::min(d(x_high, y_low), d_min);
            d_min = std::min(d(x_high, y_high), d_min);

            d_max = std::max(d(x_low, y_low), d_max);
            d_max = std::max(d(x_low, y_high), d_max);
            d_max = std::max(d(x_high, y_low), d_max);
            d_max = std::max(d(x_high, y_high), d_max);

            // TODO: clamp d
            // clamp d
            if (d_tmp(x, y) < d_min || d_tmp(x, y) > d_max)
                d_tmp(x, y) = d_forward(x, y);
        }
    }

    // Clamp velocities (u), MAC grid
    for (int y = 1; y < m_res_y - 1; ++y)
    {
        for (int x = 1; x < m_res_x; ++x)
        {
            // TODO: Compute the velocity
            double last_x_velocity = u(x, y);
            double last_y_velocity = 0.25 * (v(x - 1, y) + v(x, y) +
                v(x - 1, y + 1) + v(x, y + 1));

            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Being conservative, one can say that the velocities are known between 1.5 and res-2.5
            // (the MAC grid is inside the known densities, which are between 1 and res - 2)
            if (last_x < 1.5)
                last_x = 1.5;
            if (last_y < 1.5)
                last_y = 1.5;
            if (last_x > m_res_x - 1.5)
                last_x = m_res_x - 1.5;
            if (last_y > m_res_y - 2.5)
                last_y = m_res_y - 2.5;

            // Determine corners for bilinear interpolation
            int x_low = (int)last_x;
            int y_low = (int)last_y;
            int x_high = x_low + 1;
            int y_high = y_low + 1;

            double u_min = 1e10;
            double u_max = -1e10;

            u_min = std::min(u(x_low, y_low), u_min);
            u_min = std::min(u(x_low, y_high), u_min);
            u_min = std::min(u(x_high, y_low), u_min);
            u_min = std::min(u(x_high, y_high), u_min);

            u_max = std::max(u(x_low, y_low), u_max);
            u_max = std::max(u(x_low, y_high), u_max);
            u_max = std::max(u(x_high, y_low), u_max);
            u_max = std::max(u(x_high, y_high), u_max);

            // TODO: clamp u
            // clamp u
            if (u_tmp(x, y) < u_min || u_tmp(x, y) > u_max)
                u_tmp(x, y) = u_forward(x, y);
        }
    }

    // Clamp velocities (v), MAC grid
    for (int y = 1; y < m_res_y; ++y)
    {
        for (int x = 1; x < m_res_x - 1; ++x)
        {
            // TODO: Compute the velocity
            double last_x_velocity = 0.25 * (u(x, y - 1) + u(x + 1, y - 1) +
                u(x, y) + u(x + 1, y));
            double last_y_velocity = v(x, y);

            double last_x = x - m_dt * last_x_velocity / m_dx;
            double last_y = y - m_dt * last_y_velocity / m_dx;

            // Make sure the coordinates are inside the boundaries
            // Being conservative, one can say that the velocities are known between 1.5 and res-2.5
            // (the MAC grid is inside the known densities, which are between 1 and res - 2)
            if (last_x < 1.5)
                last_x = 1.5;
            if (last_y < 1.5)
                last_y = 1.5;
            if (last_x > m_res_x - 2.5)
                last_x = m_res_x - 2.5;
            if (last_y > m_res_y - 1.5)
                last_y = m_res_y - 1.5;

            // Determine corners for bilinear interpolation
            double x_low = (int)last_x;
            double y_low = (int)last_y;
            double x_high = x_low + 1;
            double y_high = y_low + 1;

            double v_min = 1e10;
            double v_max = -1e10;

            v_min = std::min(v(x_low, y_low), v_min);
            v_min = std::min(v(x_low, y_high), v_min);
            v_min = std::min(v(x_high, y_low), v_min);
            v_min = std::min(v(x_high, y_high), v_min);

            v_max = std::max(v(x_low, y_low), v_max);
            v_max = std::max(v(x_low, y_high), v_max);
            v_max = std::max(v(x_high, y_low), v_max);
            v_max = std::max(v(x_high, y_high), v_max);

            // TODO: clamp v
            // clamp v
            if (v_tmp(x, y) < v_min || v_tmp(x, y) > v_max)
                v_tmp(x, y) = v_forward(x, y);
        }
    }

    p_density.x() = d_tmp;
    p_velocity.x() = u_tmp;
    p_velocity.y() = v_tmp;
}

#pragma GCC diagnostic pop