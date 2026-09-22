#pragma once

#include <memory>
#include <vector>

class PBSSubApp;

class PBSApp
{
public:
    std::shared_ptr<PBSSubApp> subapp;

    PBSApp(std::shared_ptr<PBSSubApp> subapp_ref) : subapp(subapp_ref) {}

public:
    void launch();
};
