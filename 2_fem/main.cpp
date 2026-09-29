#include "PBSApp.h"
#include "FEMApp.h"

int main(int argc, char *argv[])
{
    FEMApp app;
    PBSApp pbsapp(std::make_shared<FEMApp>(app));
    pbsapp.launch();

    return 0;
}