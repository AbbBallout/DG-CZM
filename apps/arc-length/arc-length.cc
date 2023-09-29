
#include "../../include/arc-length/arc_truss.h"
#include "../../include/arc-length/ARC.h"

#include <iostream>

  int main(int argc, char** argv)
{
  try
    {
      // TrussProblem<1,2> prob;
      // prob.run() ; 

        ARC::ProblemParameters<2> par;
        ParameterAcceptor::initialize(argv[1]);
  
      
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        ARC::ARC<2> problem(par.degree,par);
        problem.run(); 

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
 
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
