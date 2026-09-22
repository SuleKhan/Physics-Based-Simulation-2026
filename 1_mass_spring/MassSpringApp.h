#pragma once

#include "PBSSubApp.h"
#include "ImGuiHelpers.h"
#include "VecMatDef.h"
#include "CameraHelper.h"
#include <iostream>

#define CAMERA_PAN_SENSITIVITY 0.003f
#define CAMERA_ROTATE_SENSITIVITY 0.003f
#define CAMERA_ZOOM_SENSITIVITY 0.02f

class MassSpringApp : public PBSSubApp
{
public:
    bool m_simulating = false;
    bool m_single_step = false;
    int m_method = 0; // integrator method

    double m_dt = 1e-2;
    double m_mass = 1.0;
    double m_spring_stiffness = 5.0;
    double m_spring_length = 3.0;
    double m_damping = 0.1;

    Vector3F m_gravity = Vector3F(0, -9.81, 0);
    Vector3F m_spring_start = Vector3F::Zero();
    Vector3F m_position = Vector3F(0, -3.0, 0);
    Vector3F m_velocity = Vector3F::Zero();

    double m_time = 0.0;
    std::vector<Vector3F> m_trajectory;

    /// Camera state
    PBSCamera app_camera;
    bool dragging_mass = false;

    /// Overlays
    bool show_axes = true;
    bool show_reference_plane = true;

    /// Implementation correctness checking
    std::string m_check_info;

    void mainLoop() override
    {
        if (m_simulating || m_single_step)
        {
            m_single_step = false;
            step();

            m_time += m_dt;
            m_trajectory.push_back(m_position);
        }
    }

    void resetSimulation()
    {
        m_position = Eigen::Vector3d(0, -3.0, 0);
        m_velocity = Eigen::Vector3d::Zero();
        m_time = 0.0;
        m_trajectory.clear();
    }

    void makeConfigWindow() override
    {
        if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Simulate", &m_simulating);
            if (ImGui::Button("Single Step"))
            {
                m_single_step = true;
            }
            if (ImGui::Button("Reset Simulation"))
            {
                resetSimulation();
            }
            if (ImGui::Button("Clear Trajectory"))
            {
                m_trajectory.clear();
            }
        }
        if (ImGui::CollapsingHeader("Simulation Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Combo("Integrator", &m_method, "Analytical (1D)\0Explicit Euler\0Symplectic Euler\0Explicit Midpoint\0Implicit Euler (1D)\0");
            if (m_method == 0 || m_method == 4)
            {
                ImGui::TextWrapped("Note: Analytical and Implicit Euler methods only work for 1D vertical motion. Dragging the mass in 3D will not work with these methods.");
            }
            ImGui::InputDouble("Timestep", &m_dt);
            ImGui::InputDouble("Mass", &m_mass);
            ImGui::InputDouble("Spring Stiffness", &m_spring_stiffness);
            ImGui::InputDouble("Damping", &m_damping);
        }
        if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Show Axes", &show_axes);
            ImGui::Checkbox("Show Reference Plane", &show_reference_plane);
        }
    }

    void makeAnalysisWindow() override
    {
        if (ImGui::Button("Check Assignment"))
        {
            m_check_info = m_check();
        }
        ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(10, 10, 10, 255));
        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(255, 255, 255, 255));
        ImGui::BeginChild("ScrollingRegion", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);
        ImGui::TextWrapped("%s", m_check_info.c_str());
        ImGui::EndChild();
        ImGui::PopStyleColor(2);
    }

    MassSpringApp()
    {
        /// Initialize camera.
        app_camera.eye = Vector3F(1, 1, 5);
        app_camera.center = Vector3F(0, 0, 0);
        app_camera.set_up_direction(Vector3F(0, 1, 0));
        app_camera.height = 10.0;
    }

    void getViewerData(std::vector<PBSViewerData> &viewer_data, PBSCamera &viewer_camera) override
    {
        /// Copy camera state.
        viewer_camera = app_camera;

        /// Viewer data for spring.
        viewer_data.emplace_back();
        PBSViewerData &viewer_data_spring = viewer_data.back();
        /// Render spring endpoints.
        viewer_data_spring.points = MatrixXF(2, 3);
        viewer_data_spring.points.row(0) = m_spring_start;
        viewer_data_spring.points.row(1) = m_position;
        viewer_data_spring.points_c = MatrixXF::Zero(2, 3);
        viewer_data_spring.point_size = 0.1;
        /// Render spring as line.
        viewer_data_spring.lines_v = MatrixXF(2, 3);
        viewer_data_spring.lines_v.row(0) = m_spring_start;
        viewer_data_spring.lines_v.row(1) = m_position;
        viewer_data_spring.lines_e = MatrixXI(1, 2);
        viewer_data_spring.lines_e.row(0) << 0, 1;
        viewer_data_spring.lines_c = MatrixXF::Zero(1, 3);
        viewer_data_spring.edge_width = 0.035;

        /// Viewer data for mass.
        viewer_data.emplace_back();
        PBSViewerData &viewer_data_mass = viewer_data.back();
        /// Render mass as point.
        viewer_data_mass.points = MatrixXF(1, 3);
        viewer_data_mass.points.row(0) = m_position;
        viewer_data_mass.points_c = MatrixXF::Zero(1, 3);
        viewer_data_mass.points_c.row(0) = Vector3F(0, 0, 1);
        viewer_data_mass.point_size = 0.5;

        /// Viewer data for trajectory.
        viewer_data.emplace_back();
        PBSViewerData &viewer_data_trajectory = viewer_data.back();
        /// Render trajectory as points.
        viewer_data_trajectory.points = MatrixXF(m_trajectory.size(), 3);
        for (size_t i = 0; i < m_trajectory.size(); i++)
        {
            viewer_data_trajectory.points.row(i) = m_trajectory[i];
        }
        viewer_data_trajectory.points_c = MatrixXF::Zero(m_trajectory.size(), 3);
        viewer_data_trajectory.points_c.rowwise() = Vector3F(0, 1, 0).transpose();
        viewer_data_trajectory.point_size = 0.05;

        /// Viewer data for coordinate axes.
        if (show_axes)
        {
            viewer_data.emplace_back();
            PBSViewerData &viewer_data_axes = viewer_data.back();

            viewer_data_axes.points = MatrixXF(3, 3);
            viewer_data_axes.points.row(0) = Vector3F(1, 0, 0);
            viewer_data_axes.points.row(1) = Vector3F(0, 1, 0);
            viewer_data_axes.points.row(2) = Vector3F(0, 0, 1);
            viewer_data_axes.points_c = MatrixXF::Zero(3, 3);
            viewer_data_axes.points_c.row(0) = Vector3F(1, 0, 0);
            viewer_data_axes.points_c.row(1) = Vector3F(0, 1, 0);
            viewer_data_axes.points_c.row(2) = Vector3F(0, 0, 1);
            viewer_data_axes.point_size = 0.2;

            viewer_data_axes.lines_v = MatrixXF(4, 3);
            viewer_data_axes.lines_v.row(0) = Vector3F(0, 0, 0);
            viewer_data_axes.lines_v.row(1) = Vector3F(1, 0, 0);
            viewer_data_axes.lines_v.row(2) = Vector3F(0, 1, 0);
            viewer_data_axes.lines_v.row(3) = Vector3F(0, 0, 1);
            viewer_data_axes.lines_e = MatrixXI(3, 2);
            viewer_data_axes.lines_e.row(0) << 0, 1;
            viewer_data_axes.lines_e.row(1) << 0, 2;
            viewer_data_axes.lines_e.row(2) << 0, 3;
            viewer_data_axes.lines_c = MatrixXF(3, 3);
            viewer_data_axes.lines_c.row(0) = Vector3F(1, 0, 0);
            viewer_data_axes.lines_c.row(1) = Vector3F(0, 1, 0);
            viewer_data_axes.lines_c.row(2) = Vector3F(0, 0, 1);
            viewer_data_axes.edge_width = 0.1;

            viewer_data_axes.alpha = 1.0;
        }

        /// Viewer data for reference plane.
        if (show_reference_plane)
        {
            viewer_data.emplace_back();
            PBSViewerData &viewer_data_plane = viewer_data.back();

            double plane_radius = 5.0;
            int plane_resolution = 10;
            viewer_data_plane.lines_v = MatrixXF((plane_resolution + 1) * 4, 3);
            for (int i = 0; i <= plane_resolution; i++)
            {
                double t = (double)i / (double)plane_resolution;
                double x = (t - 0.5) * 2.0 * plane_radius;
                viewer_data_plane.lines_v.row(i * 4 + 0) = Vector3F(x, 0, -plane_radius);
                viewer_data_plane.lines_v.row(i * 4 + 1) = Vector3F(x, 0, plane_radius);
                viewer_data_plane.lines_v.row(i * 4 + 2) = Vector3F(-plane_radius, 0, x);
                viewer_data_plane.lines_v.row(i * 4 + 3) = Vector3F(plane_radius, 0, x);
            }
            viewer_data_plane.lines_e = MatrixXI((plane_resolution + 1) * 2, 2);
            for (int i = 0; i <= plane_resolution; i++)
            {
                viewer_data_plane.lines_e.row(i * 2 + 0) << i * 4 + 0, i * 4 + 1;
                viewer_data_plane.lines_e.row(i * 2 + 1) << i * 4 + 2, i * 4 + 3;
            }
            viewer_data_plane.lines_c = MatrixXF::Constant((plane_resolution + 1) * 2, 3, 0.5);
            viewer_data_plane.edge_width = 0.02;
            viewer_data_plane.alpha = 1.0;
        }

        /// Copy camera state.
        viewer_camera = app_camera;
    }

    bool callbackKeyPressed(const PBSControlState &control_state, int key) override
    {
        switch (key)
        {
        case ImGuiKey_Space:
            m_simulating = !m_simulating;
            break;
        default:
            break;
        }
        return false;
    }

    bool callbackMouseDown(const PBSControlState &control_state, int button, const Vector2F &mouse_pos) override
    {
        if (button == MOUSE_LEFT)
        {
            Vector2I screen_dims = control_state.screen_dimensions;

            // check if mouse is near mass
            Vector3F camera_dir = app_camera.get_camera_direction();
            Vector3F plane_normal = camera_dir;
            Vector3F plane_point = m_position;
            Vector3F ray_origin = app_camera.eye;
            Vector2F screen_coords(mouse_pos.x() / screen_dims.x(), mouse_pos.y() / screen_dims.y());
            Vector3F ray_dir = ((app_camera.center - app_camera.eye) + (screen_coords.x() - 0.5f) * app_camera.height * 2.0 * (screen_dims.x() * 1.0f / screen_dims.y()) * app_camera.get_right_direction() -
                                (screen_coords.y() - 0.5f) * app_camera.height * 2.0 * app_camera.get_up_direction())
                                   .normalized();

            // intersect ray with plane
            if (plane_normal.dot(ray_dir) > 1e-6)
            {
                double t = (plane_point - ray_origin).dot(plane_normal) / plane_normal.dot(ray_dir);
                Vector3F intersect_point = ray_origin + t * ray_dir;
                double dist_to_mass = (intersect_point - m_position).norm();
                if (dist_to_mass < 0.5f)
                {
                    dragging_mass = true;
                }
            }
        }
        return false;
    }

    bool callbackMouseMove(const PBSControlState &control_state, const Vector2F &mouse_pos,
                           const Vector2F &mouse_delta) override
    {
        if (control_state.modifiers[MOUSE_RIGHT])
        {
            if (control_state.modifiers[KEY_SHIFT])
            {
                CameraHelper::panCameraFromDrag(app_camera, mouse_delta, CAMERA_PAN_SENSITIVITY);
            }
            else
            {
                CameraHelper::rotateCameraFromDrag(app_camera, mouse_delta, CAMERA_ROTATE_SENSITIVITY);
            }
        }
        if (control_state.modifiers[MOUSE_LEFT])
        {
            if (dragging_mass)
            {
                Vector2I screen_dims = control_state.screen_dimensions;

                Vector3F camera_dir = app_camera.get_camera_direction();
                Vector3F plane_normal = camera_dir;
                Vector3F plane_point = m_position;
                Vector3F ray_origin = app_camera.eye;
                Vector2F screen_coords(mouse_pos.x() / screen_dims.x(), mouse_pos.y() / screen_dims.y());
                Vector3F ray_dir = ((app_camera.center - app_camera.eye) + (screen_coords.x() - 0.5f) * app_camera.height * 2.0 * (screen_dims.x() * 1.0f / screen_dims.y()) * app_camera.get_right_direction() -
                                    (screen_coords.y() - 0.5f) * app_camera.height * 2.0 * app_camera.get_up_direction())
                                       .normalized();

                // intersect ray with plane
                if (plane_normal.dot(ray_dir) > 1e-6)
                {
                    double t = (plane_point - ray_origin).dot(plane_normal) / plane_normal.dot(ray_dir);
                    Vector3F intersect_point = ray_origin + t * ray_dir;
                    m_position = intersect_point;
                    m_velocity = Vector3F::Zero(); // reset velocity when dragging
                }
            }
        }

        return false;
    }

    bool callbackMouseUp(const PBSControlState &control_state, int button, const Vector2F &mouse_pos) override
    {
        if (button == MOUSE_LEFT)
        {
            dragging_mass = false;
        }
        return false;
    }

    bool callbackMouseScroll(const PBSControlState &control_state, float t) override
    {
        Vector3F camera_dist = app_camera.eye - app_camera.center;
        app_camera.eye = app_camera.center + camera_dist * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);
        app_camera.height = app_camera.height * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);

        return false;
    }

    bool step();

    void createTests();

    void createRandomTests(int num_tests_per_method);

    std::string m_check();
};