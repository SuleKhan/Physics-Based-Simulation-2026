#include "PBSApp.h"
#include "MassSpringApp.h"

int main(int argc, char *argv[])
{
    MassSpringApp app;
    PBSApp pbsapp(std::make_shared<MassSpringApp>(app));
    pbsapp.launch();

    return 0;
}