#pragma once

#include "PBSSubApp.h"

namespace CameraHelper
{
    void rotateCameraFromDrag(PBSCamera &camera, const Vector2F &drag, double sensitivity);

    void panCameraFromDrag(PBSCamera &camera, const Vector2F &drag, double sensitivity);
} // namespace CameraHelper
