#include <iostream>
#include "parameters.h"
#include "worker.h"

int main(int argc, char** argv)
{
    if (argc<2)
        {
        std::cout<<" Usage: " << argv[0] << " [input file name] " << std::endl;
        std::exit(0);
        }
    parameters par(argv[1]);
    par.analyse_subspaces(); // maybe move to parmeter constructor 
    par.print();    // maybe move to parmeter constructor

    worker sys;
    if (par.TCSA)
        {
        sys.run(par);
        }

    return 0;

}
