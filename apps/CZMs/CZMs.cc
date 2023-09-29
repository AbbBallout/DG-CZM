#include "../../include/CZMs/Exponential_TSL.h"  // Namespace Exponential_TSL
#include "../../include/CZMs/fenics_delamitation.h"  // Namespace fenics_delamitation
#include "../../include/CZMs/TSL.h" // Namespace TSL
#include "../../include/CZMs/TSL_parallel.h" // Namespace TSL_prllel
#include "../../include/CZMs/TSL_prallel_old.h" // Namespace TSL_prllel_old
#include "../../include/CZMs/TSL_everywhere.h"  // TSL_everywhere



 int main(int argc, char** argv)
{   
   using namespace dealii; 
  
    try
     {
      if (argc == 1)
      { 
        std::cout << "argc = " << argc << "  No input paramter passed \n " ;  
        exit(1); 
      }          

        TSL_everywhere::ProblemParameters<2> par;
        ParameterAcceptor::initialize(argv[1]);
  
      
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        TSL_everywhere::TSL_everywhere<2> problem(par.degree,par);
        problem.run(); 
     }

    
    catch (std::exception& exc)
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
    return 0;
}


