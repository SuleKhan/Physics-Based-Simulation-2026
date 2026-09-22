#pragma once

#include <Eigen/Geometry>
#include <GLFW/glfw3.h> // For key
#include "VecMatDef.h"

enum PBSModifier
{
    MOUSE_LEFT,
    MOUSE_MIDDLE,
    MOUSE_RIGHT,
    KEY_SHIFT,
    KEY_CTRL,
    KEY_ALT,
    NUM_MODIFIERS /// The integer value of this will be the number of modifiers.
};

struct PBSControlState
{
    bool modifiers[PBSModifier::NUM_MODIFIERS];
    Vector2I screen_dimensions;

    PBSControlState() : modifiers{} { std::fill(std::begin(modifiers), std::end(modifiers), false); }
};

struct PBSViewerData
{
    MatrixXF points;
    MatrixXF points_c;
    MatrixXF lines_v;
    MatrixXI lines_e;
    MatrixXF lines_c;
    MatrixXF mesh_v;
    MatrixXI mesh_f;
    MatrixXF mesh_c;
    MatrixXF mesh_tex_uv;
    MatrixXF mesh_tex_r;
    MatrixXF mesh_tex_g;
    MatrixXF mesh_tex_b;

    bool mesh_flat_material = false;
    double alpha = 1.0;
    bool show_lines = false;
    bool points_quad = false;
    double point_size = 0.01;
    double edge_width = 0.003;

    bool write_png = false;
    std::string png_file_name;
};

struct PBSCamera
{
    Vector3F eye;
    Vector3F center;

    /// Half the window height, i.e. the distance (normal to camera direction) from center to top of window.
    double height;

private:
    /// Private to require calling get_up_direction, which ensures a vector perpendicular to camera direction.
    Vector3F up;

public:
    [[nodiscard]] Vector3F get_camera_direction() const { return (center - eye).normalized(); }

    [[nodiscard]] Vector3F get_right_direction() const { return get_camera_direction().cross(up).normalized(); }

    [[nodiscard]] Vector3F get_up_direction() const
    {
        return get_right_direction().cross(get_camera_direction()).normalized();
    }

    void set_up_direction(const Vector3F &up_direction) { up = up_direction; }
};

class PBSSubApp
{
public:
    /// Get data to display in 3D graphics window.
    virtual void getViewerData(std::vector<PBSViewerData> &viewer_data, PBSCamera &camera) = 0;

    /// Construct ImGui config menu on left of screen.
    virtual void makeConfigWindow() {}

    /// Construct ImGui analysis window on right of screen.
    virtual void makeAnalysisWindow() {}

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackKeyPressed(const PBSControlState &control_state, int key) { return true; }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseScroll(const PBSControlState &control_state, float t) { return true; }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseDown(const PBSControlState &control_state, int button, const Eigen::Vector2d &mouse_pos)
    {
        return true;
    }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseUp(const PBSControlState &control_state, int button, const Eigen::Vector2d &mouse_pos)
    {
        return true;
    }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseMove(const PBSControlState &control_state, const Eigen::Vector2d &mouse_pos,
                                   const Eigen::Vector2d &mouse_delta)
    {
        return true;
    }

    /// Main loop function.
    virtual void mainLoop() {}
};
