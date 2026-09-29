#pragma once

#include "PBSSubApp.h"
#include "ImGuiHelpers.h"
#include "VecMatDef.h"
#include "CameraHelper.h"
#include <iostream>
#include <numbers>

struct Pin
{
    int node;     // index of the node being pinned
    Vector2F pos; // position to pin to
    double k;     // spring constant
};

struct Element
{
    Vector3I vertices; // vertex indices
    double area;       // undeformed area
    Matrix2F DX_inv;   // inverse of undeformed shape matrix
};

class FEMApp : public PBSSubApp
{
public:
    bool m_simulating = false;

    VectorXF x_prev, X; // deformed and undeformed node positions
    VectorXF v;         // node velocities
    VectorXF M;         // mass matrix diagonal
    std::vector<Element> elements;

    bool use_dynamics = true;
    double m_dt = 1e-2; // timestep

    bool use_gravity = true;
    Vector2F gravity = Vector2F(0, -9.81);

    /// Pins
    int dragging_pin = -1; // index of pin object being dragged
    Vector2F pin_pos_0;
    double pin_angle_0;
    Vector2F pin_pos_1;
    double pin_angle_1;
    std::vector<Pin> pins;

    /// Camera state
    PBSCamera app_camera;

    /// Mesh parameters
    double rx = 0.5;
    double ry = 0.05;
    int nx = 31;
    int ny = 3;

    // Material parameters
    double mass_per_area = 1.0;
    double shearModulus = 1024.;
    double bulkModulus = 1024.;
    double pin_k = 1e5;

    /// Implementation correctness checking
    std::string m_check_info;

    void mainLoop() override
    {
        if (m_simulating)
        {
            step();
        }
    }

    void resetSimulation()
    {
        int n_nodes = nx * ny;
        int n_elements = (nx - 1) * (ny - 1) * 2;

        elements.resize(n_elements);

        X = VectorXF::Zero(2 * n_nodes);
        v = VectorXF::Zero(2 * n_nodes);

        M = VectorXF::Zero(2 * n_nodes);

        // create grid mesh
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int node = j * nx + i;
                double px = rx * (2.0 * i / (nx - 1) - 1.0);
                double py = ry * (2.0 * j / (ny - 1) - 1.0);
                X(2 * node) = px;
                X(2 * node + 1) = py;
            }
        }
        x_prev = X;
        // create elements
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                int elem = (j * (nx - 1) + i) * 2;
                int n0 = j * nx + i;
                int n1 = j * nx + (i + 1);
                int n2 = (j + 1) * nx + i;
                int n3 = (j + 1) * nx + (i + 1);
                elements[elem].vertices = Vector3I(n0, n1, n2);
                elements[elem + 1].vertices = Vector3I(n1, n3, n2);
            }
        }
        for (int i = 0; i < n_elements; ++i)
        {
            initElement(i);
        }

        // nodal masses
        for (int i = 0; i < n_elements; ++i)
        {
            // compute area of triangle element
            int v0 = elements[i].vertices(0);
            int v1 = elements[i].vertices(1);
            int v2 = elements[i].vertices(2);
            Vector2F p0 = X.segment<2>(2 * v0);
            Vector2F p1 = X.segment<2>(2 * v1);
            Vector2F p2 = X.segment<2>(2 * v2);
            double area = 0.5 * ((p1 - p0).x() * (p2 - p0).y() - (p1 - p0).y() * (p2 - p0).x());
            for (int j = 0; j < 3; ++j)
            {
                // each node gets 1/3 of element mass
                int vertex = elements[i].vertices(j);
                M(2 * vertex) += mass_per_area * area / 3.0;
                M(2 * vertex + 1) += mass_per_area * area / 3.0;
            }
        }

        // set up pins
        pin_pos_0 = Vector2F(-0.5, 0.0);
        pin_angle_0 = std::numbers::pi;
        pin_pos_1 = Vector2F(0.5, 0.0);
        pin_angle_1 = 0.0;
        setPins();
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

    void setPins()
    {
        pins.clear();
        // pin left edge
        for (int j = 0; j < ny; ++j)
        {
            int node = j * nx;
            double offset_y = ry * (2.0 * j / (ny - 1) - 1.0);
            Vector2F pos = pin_pos_0 + Vector2F(offset_y * std::sin(pin_angle_0), -offset_y * std::cos(pin_angle_0));
            pins.push_back({node, pos, pin_k});
        }
        // pin right edge
        for (int j = 0; j < ny; ++j)
        {
            int node = j * nx + (nx - 1);
            double offset_y = ry * (2.0 * j / (ny - 1) - 1.0);
            Vector2F pos = pin_pos_1 + Vector2F(-offset_y * std::sin(pin_angle_1), offset_y * std::cos(pin_angle_1));
            pins.push_back({node, pos, pin_k});
        }
    }

    void makeConfigWindow() override
    {
        ImGui::Checkbox("Simulate", &m_simulating);
        if (ImGui::Button("Reset Simulation"))
        {
            resetSimulation();
        }
        ImGui::Checkbox("Use Dynamics", &use_dynamics);
        ImGui::Checkbox("Use Gravity", &use_gravity);
        if (ImGui::Button("Check Gradient"))
        {
            checkGradient();
        }
        if (ImGui::Button("Check Hessian"))
        {
            checkHessian();
        }
    }

    FEMApp()
    {
        /// Initialize camera.
        app_camera.eye = Vector3F(0, 0, 5);
        app_camera.center = Vector3F(0, 0, 0);
        app_camera.set_up_direction(Vector3F(0, 1, 0));
        app_camera.height = 0.8;

        /// Initialize simulation state.
        resetSimulation();
    }

    void getViewerData(std::vector<PBSViewerData> &viewer_data, PBSCamera &viewer_camera) override
    {
        viewer_data.clear();

        /// Draw FEM mesh.
        viewer_data.emplace_back();
        PBSViewerData &data_elements = viewer_data.back();

        int n_vertices = X.rows() / 2;
        int n_triangles = elements.size();
        data_elements.mesh_v = MatrixXF::Zero(n_vertices, 3);
        data_elements.mesh_f = MatrixXI::Zero(n_triangles, 3);
        for (int i = 0; i < n_vertices; ++i)
        {
            data_elements.mesh_v(i, 0) = x_prev(2 * i);
            data_elements.mesh_v(i, 1) = x_prev(2 * i + 1);
            data_elements.mesh_v(i, 2) = 0.0;
        }
        for (int i = 0; i < n_triangles; ++i)
        {
            data_elements.mesh_f.row(i) = elements[i].vertices.transpose();
        }
        data_elements.mesh_c = MatrixXF::Zero(n_triangles, 3);
        data_elements.mesh_flat_material = true;
        VectorXF triangle_energies = VectorXF::Zero(n_triangles);
        for (int i = 0; i < n_triangles; ++i)
        {
            addTriangleEnergy(x_prev, i, triangle_energies(i));
        }
        double max_energy = triangle_energies.maxCoeff();
        if (max_energy < 1e-8)
            max_energy = 1e-8;
        for (int i = 0; i < n_triangles; ++i)
        {
            double t = triangle_energies(i) / max_energy;
            data_elements.mesh_c(i, 0) = t;
            data_elements.mesh_c(i, 1) = 0.0;
            data_elements.mesh_c(i, 2) = 1.0 - t;
        }

        data_elements.lines_v = MatrixXF::Zero(3 * n_triangles, 3);
        data_elements.lines_e = MatrixXI::Zero(3 * n_triangles, 2);
        data_elements.lines_c = MatrixXF::Constant(3 * n_triangles, 3, 0.0);
        data_elements.edge_width = 0.003;
        for (int i = 0; i < n_triangles; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                int v0 = elements[i].vertices[j];
                data_elements.lines_v(3 * i + j, 0) = x_prev(2 * v0);
                data_elements.lines_v(3 * i + j, 1) = x_prev(2 * v0 + 1);
                data_elements.lines_v(3 * i + j, 2) = 0.0;
                data_elements.lines_e(3 * i + j, 0) = 3 * i + j;
                data_elements.lines_e(3 * i + j, 1) = 3 * i + (j + 1) % 3;
            }
        }

        /// Draw pins.
        viewer_data.emplace_back();
        PBSViewerData &data_pins = viewer_data.back();
        int n_pins = static_cast<int>(pins.size());
        data_pins.points = MatrixXF::Zero(n_pins, 3);
        data_pins.points_c = MatrixXF::Zero(n_pins, 3);
        data_pins.points_c.rowwise() = Vector3F(1.0, 0.0, 0.0).transpose();
        data_pins.point_size = 0.01;
        for (int i = 0; i < n_pins; ++i)
        {
            data_pins.points(i, 0) = pins[i].pos.x();
            data_pins.points(i, 1) = pins[i].pos.y();
            data_pins.points(i, 2) = 0.0;
        }
        data_pins.lines_v = MatrixXF::Zero(4, 3);
        data_pins.lines_e = MatrixXI::Zero(2, 2);
        data_pins.lines_c = MatrixXF::Zero(2, 3);
        data_pins.lines_c.rowwise() = Vector3F(1.0, 0.0, 0.0).transpose();
        data_pins.edge_width = 0.005;

        data_pins.lines_v.row(0).head(2) = pin_pos_0.transpose();
        data_pins.lines_v.row(1).head(2) = pin_pos_0.transpose() + 0.1 * Vector2F(cos(pin_angle_0), sin(pin_angle_0)).transpose();
        data_pins.lines_e.row(0) << 0, 1;
        data_pins.lines_v.row(2).head(2) = pin_pos_1.transpose();
        data_pins.lines_v.row(3).head(2) = pin_pos_1.transpose() + 0.1 * Vector2F(cos(pin_angle_1), sin(pin_angle_1)).transpose();
        data_pins.lines_e.row(1) << 2, 3;

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
            Vector2F screen_coords(mouse_pos.x() / screen_dims.x(), mouse_pos.y() / screen_dims.y());
            Vector3F ray_dir = ((app_camera.center - app_camera.eye) + (screen_coords.x() - 0.5f) * app_camera.height * 2.0 * (screen_dims.x() * 1.0f / screen_dims.y()) * app_camera.get_right_direction() -
                                (screen_coords.y() - 0.5f) * app_camera.height * 2.0 * app_camera.get_up_direction())
                                   .normalized();

            // Intersect ray with xy plane
            double t = -app_camera.eye.z() / ray_dir.z();
            Vector2F xy_point = (app_camera.eye + t * ray_dir).head(2);

            std::vector<Vector2F> pin_objects;
            pin_objects.push_back(pin_pos_0);
            pin_objects.push_back(pin_pos_0 + 0.1 * Vector2F(cos(pin_angle_0), sin(pin_angle_0)));
            pin_objects.push_back(pin_pos_1);
            pin_objects.push_back(pin_pos_1 + 0.1 * Vector2F(cos(pin_angle_1), sin(pin_angle_1)));

            dragging_pin = -1;
            double min_dist = 0.02;
            for (int i = 0; i < (int)pin_objects.size(); ++i)
            {
                double dist = (xy_point - pin_objects[i]).norm();
                if (dist < min_dist)
                {
                    min_dist = dist;
                    dragging_pin = i;
                }
            }
            setPins();
        }
        return false;
    }

    bool callbackMouseMove(const PBSControlState &control_state, const Vector2F &mouse_pos,
                           const Vector2F &mouse_delta) override
    {
        if (dragging_pin != -1)
        {
            Vector2I screen_dims = control_state.screen_dimensions;
            Vector2F screen_coords(mouse_pos.x() / screen_dims.x(), mouse_pos.y() / screen_dims.y());
            Vector3F ray_dir = ((app_camera.center - app_camera.eye) + (screen_coords.x() - 0.5f) * app_camera.height * 2.0 * (screen_dims.x() * 1.0f / screen_dims.y()) * app_camera.get_right_direction() -
                                (screen_coords.y() - 0.5f) * app_camera.height * 2.0 * app_camera.get_up_direction())
                                   .normalized();

            // Intersect ray with xy plane
            double t = -app_camera.eye.z() / ray_dir.z();
            Vector2F xy_point = (app_camera.eye + t * ray_dir).head(2);

            // TODO: assign pin pos/angle based on dragging_pin
            if (dragging_pin == 0)
            {
                pin_pos_0 = xy_point;
            }
            else if (dragging_pin == 1)
            {
                Vector2F dir = xy_point - pin_pos_0;
                pin_angle_0 = std::atan2(dir.y(), dir.x());
            }
            else if (dragging_pin == 2)
            {
                pin_pos_1 = xy_point;
            }
            else if (dragging_pin == 3)
            {
                Vector2F dir = xy_point - pin_pos_1;
                pin_angle_1 = std::atan2(dir.y(), dir.x());
            }
            setPins();
        }
        return false;
    }

    bool callbackMouseUp(const PBSControlState &control_state, int button, const Vector2F &mouse_pos) override
    {
        if (button == MOUSE_LEFT)
        {
            dragging_pin = -1;
        }
        return false;
    }

    bool callbackMouseScroll(const PBSControlState &control_state, float t) override
    {
        return false;
    }

    /// Objective.cpp
public:
    bool step();

    bool objectiveFunction(const VectorXF &x, double &E);
    bool objectiveGradient(const VectorXF &x, VectorXF &gradient);
    bool objectiveHessian(const VectorXF &x, TripletListF &hessian_triplets);

    bool checkGradient();
    bool checkHessian();

    /// FiniteDifference.cpp
public:
    void computeGradientFD(const VectorXF &x, VectorXF &gradient, double eps = 1e-6);
    void computeHessianFD(const VectorXF &x, MatrixXF &H, double eps = 1e-6);

    /// Dynamics.cpp
public:
    void addDynamicsPotential(const VectorXF &x, double &E);
    void addDynamicsGradient(const VectorXF &x, VectorXF &gradient);
    void addDynamicsHessian(const VectorXF &x, TripletListF &hessian_triplets);

    /// Energy.cpp
public:
    void addSystemEnergy(const VectorXF &x, double &E);
    void addSystemGradient(const VectorXF &x, VectorXF &gradient);
    void addSystemHessian(const VectorXF &x, TripletListF &hessian_triplets);

    /// Pin.cpp
public:
    void addPinEnergy(const VectorXF &x, int pin_index, double &E);
    void addPinGradient(const VectorXF &x, int pin_index, VectorXF &gradient);
    void addPinHessian(const VectorXF &x, int pin_index, TripletListF &hessian_triplets);

public:
    /// Triangle.cpp
    void initElement(int triangle_index);
    Matrix2F deformationGradient(const VectorXF &x, int triangle_index);
    void addTriangleEnergy(const VectorXF &x, int triangle_index, double &E);
    /// TriangleGradient.cpp
    void addTriangleGradient(const VectorXF &x, int triangle_index, VectorXF &gradient);
    void addTriangleHessian(const VectorXF &x, int triangle_index, TripletListF &hessian_triplets);

    /// Gravity.cpp
public:
    void addGravityEnergy(const VectorXF &x, double &E);
    void addGravityGradient(const VectorXF &x, VectorXF &gradient);
    void addGravityHessian(const VectorXF &x, TripletListF &hessian_triplets);

public:
    std::string m_check();

};
