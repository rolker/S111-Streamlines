#include "s111_streamlines.h"

int main(int argc, char *argv[])
{
    s111::FlowField flowField;
    for(int i = 1; i < argc; i++)
        flowField.addFile(argv[i]);
    return 0;
}
