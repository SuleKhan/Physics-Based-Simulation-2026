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
    Vector3F a = Vector3F::Zero();

    // Suppress unused variable warnings until you implement the methods.
    (void)a;
    (void)spring_dir;
    (void)spring_norm;

    switch (m_method)
    {
    case 0:
    {
        /// TODO: Implement analytical solution for 1D vertical motion.
        m_position = Vector3F::Zero();
        m_velocity = Vector3F::Zero();
        break;
    }

    case 1:
    {
        /// TODO: Implement explicit Euler method.
        m_position = Vector3F::Zero();
        m_velocity = Vector3F::Zero();
        break;
    }

    case 2:
    {
        /// TODO: Implement symplectic Euler method.
        m_position = Vector3F::Zero();
        m_velocity = Vector3F::Zero();
        break;
    }
    case 3:
    {
        /// TODO: Implement explicit midpoint method.
        m_position = Vector3F::Zero();
        m_velocity = Vector3F::Zero();
        break;
    }
    case 4:
    {
        /// TODO: Implement implicit Euler method. You can assume 1D vertical motion and implement this in closed form.
        /// For an extra challenge, try implementing the general 3D case using fixed point iteration - then you can drag the mass in 3D!
        m_position = Vector3F::Zero();
        m_velocity = Vector3F::Zero();
        break;
    }
    default:
        assert(false && "Unknown method");
        break;
    }

    return false;
}
