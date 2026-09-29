#include "../MassSpringApp.h"

#include <iostream>

/////////////////////////////////////
//////// EX 1 ///////////
/////////////////////////////////////
bool MassSpringApp::step()
{
    // perform time integration with different integrators

    // useful functions: normalized(), norm()
    Vector3F spring_dir = (m_position - m_spring_start).normalized();
    double spring_norm = (m_position - m_spring_start).norm();

    /// TODO: Compute the acceleration of the mass.
    /// HINT: use m_spring_length, m_spring_stiffness, m_damping, m_mass, m_gravity to compute forces.
    Vector3F f_int = -m_spring_stiffness * (spring_norm - m_spring_length) * spring_dir;
    Vector3F f_damp = -m_damping * m_velocity;
    Vector3F f_ext = m_mass * m_gravity;
    Vector3F f = f_int + f_damp + f_ext;
    Vector3F a = f / m_mass;

    switch (m_method)
    {
    case 0:
    {
        // analytical solution
        float g = m_gravity.norm();
        float m = m_mass;
        float c1 = m * g / m_spring_stiffness;
        float alpha = -m_damping / (2 * m);
        float beta = sqrt(4 * m_spring_stiffness * m - m_damping * m_damping) /
                     (2 * m);
        float c2 = -c1 * alpha / beta;

        float v_y = exp(alpha * m_time) *
                    (c1 * (alpha * cos(beta * m_time) - beta * sin(beta * m_time)) +
                     c2 * (alpha * sin(beta * m_time) + beta * cos(beta * m_time)));
        float p_y = exp(alpha * m_time) *
                        (c1 * cos(beta * m_time) + c2 * sin(beta * m_time)) -
                    m_spring_length - m * g / m_spring_stiffness;

        m_velocity = Eigen::Vector3d(0, v_y, 0);
        m_position = Eigen::Vector3d(0, p_y, 0);
        break;
    }

    case 1:
    {
        // explicit euler
        m_position += m_dt * m_velocity;
        m_velocity += m_dt * a;
        break;
    }

    case 2:
    {
        // symplectic euler
        m_velocity += m_dt * a;
        m_position += m_dt * m_velocity;
        break;
    }
    case 3:
    {
        // explicit midpoint
        Eigen::Vector3d vm = m_velocity + 0.5 * m_dt * a;
        Eigen::Vector3d pm = m_position + 0.5 * m_dt * m_velocity;
        spring_dir = (pm - m_spring_start).normalized();
        spring_norm = (pm - m_spring_start).norm();
        f_int = -m_spring_stiffness * (spring_norm - m_spring_length) * spring_dir;
        f_damp = -m_damping * vm;
        f_ext = m_mass * m_gravity;
        f = f_int + f_damp + f_ext;
        a = f / m_mass;

        m_position += m_dt * vm;
        m_velocity += m_dt * a;
        break;
    }
    case 4:
    {
        // implicit euler
        float v_y = (m_mass * m_velocity.y() - m_dt *
                                                   (m_spring_stiffness * (m_position.y() + m_spring_length) - m_mass * m_gravity.y())) /
                    (m_spring_stiffness * m_dt * m_dt + m_damping * m_dt + m_mass);
        float p_y = m_position.y() + m_dt * v_y;

        m_velocity = Eigen::Vector3d(0, v_y, 0);
        m_position = Eigen::Vector3d(0, p_y, 0);
        break;
    }
    default:
        std::cerr << m_method << " is not a valid integrator method."
                  << std::endl;
    }

    return false;
}
