#include "../MassSpringApp.h"

#include <iostream>
#include <cassert>

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
    Vector3F f = -m_spring_stiffness * (spring_norm - m_spring_length) * spring_dir - m_damping * m_velocity + m_mass * m_gravity;
	Vector3F a = f / m_mass;

    // Suppress unused variable warnings until you implement the methods.
    (void)a;
    (void)spring_dir;
    (void)spring_norm;

    switch (m_method)
    {
    case 0:
    {
        /// TODO: Implement analytical solution for 1D vertical motion.

		double g = -m_gravity.y();

        if (m_time == 0.0)
		{
		    alpha = - m_damping / (2 * m_mass);
            beta = sqrt((4 * m_spring_stiffness * m_mass) - (m_damping * m_damping)) / (2 * m_mass);

		    c1 = m_pos0.y() + m_spring_length + ((m_mass * g) / m_spring_stiffness);
            c2 = (m_vel0.y() - (alpha * c1)) / beta;
		}


		m_position.y() = exp(alpha * m_time) * (c1 * cos(beta * m_time) + c2 * sin(beta * m_time)) - m_spring_length - ((m_mass * g) / m_spring_stiffness);

        
        break;
    }

    case 1:
    {
        /// TODO: Implement explicit Euler method.
        m_position = m_position + m_dt * m_velocity;
        m_velocity = m_velocity + m_dt * a;
        break;
    }

    case 2:
    {
        /// TODO: Implement symplectic Euler method.
        //spring_dir = (m_position - m_spring_start).normalized();
        //spring_norm = (m_position - m_spring_start).norm();
        //f = -m_spring_stiffness * (spring_norm - m_spring_length) * spring_dir - m_damping * m_velocity + m_mass * m_gravity;
        //a = f / m_mass;
        m_velocity = m_velocity + m_dt * a;
        m_position = m_position + m_dt * m_velocity;
        break;
    }
    case 3:
    {
        /// TODO: Implement explicit midpoint method.
		Vector3F mid_position = m_position + 0.5 * m_dt * m_velocity;
		Vector3F mid_velocity = m_velocity + 0.5 * m_dt * a;
        
        m_position = m_position + m_dt * mid_velocity;

        spring_dir = (mid_position - m_spring_start).normalized();
        spring_norm = (mid_position - m_spring_start).norm();
        f = -m_spring_stiffness * (spring_norm - m_spring_length) * spring_dir - m_damping * mid_velocity + m_mass * m_gravity;
        a = f / m_mass;

        m_velocity = m_velocity + m_dt * a;
        break;
    }
    case 4:
    {

        /// TODO: Implement implicit Euler method.
        /// Assume 1D vertical motion and implement this in closed form.
        
        double h = m_dt;
        double k = m_spring_stiffness;
        double d = m_damping;
        double m = m_mass;
        double L = m_spring_length;

        double y = m_position.y();
        double v = m_velocity.y();
        double ys = m_spring_start.y();
        double gy = m_gravity.y();

        // Solve for v_{n+1}
        double numerator =
            v + (h / m) * (-k * (y - ys + L) + m * gy);

        double denominator =
            1.0 + (h * d / m) + (h * h * k / m);

        double v_new = numerator / denominator;

        // Solve for y_{n+1}
        double y_new = y + h * v_new;

        m_velocity.y() = v_new;
        m_position.y() = y_new;


        /// Implicit Euler using fixed-point iteration

        //double h = m_dt;

        //// Initial guess for v_{n+1}
        //Vector3F v_new = m_velocity;

        //const int max_iterations = 100;
        //const double tolerance = 1e-8;

        //for (int i = 0; i < max_iterations; ++i)
        //{
        //    // x_{n+1} = x_n + h * v_{n+1}
        //    Vector3F x_new = m_position + h * v_new;

        //    // Spring at x_{n+1}
        //    Vector3F spring_vec = x_new - m_spring_start;
        //    double spring_norm = spring_vec.norm();
        //    Vector3F spring_dir = spring_vec.normalized();

        //    // Force at x_{n+1}, v_{n+1}
        //    Vector3F f =
        //        -m_spring_stiffness *
        //        (spring_norm - m_spring_length) *
        //        spring_dir
        //        - m_damping * v_new
        //        + m_mass * m_gravity;

        //    Vector3F a_new = f / m_mass;

        //    // Fixed-point update:
        //    // v_{n+1} = v_n + h * a(x_{n+1}, v_{n+1})
        //    Vector3F v_next = m_velocity + h * a_new;

        //    // Check convergence
        //    double error = (v_next - v_new).norm();

        //    v_new = v_next;

        //    if (error < tolerance)
        //        break;
        //}

        //// Once v_{n+1} has converged, calculate x_{n+1}
        //m_position = m_position + h * v_new;
        //m_velocity = v_new;

        break;
    }
    default:
        assert(false && "Unknown method");
        break;
    }

    return false;
}
