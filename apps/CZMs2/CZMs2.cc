
#include "../../include/CZMs2/TSL_everywhere.h"
#include "../../include/CZMs2/TSL-Freddi.h" // TSL_Freddi
#include "../../include/CZMs2/TSL-Friction-Freddi.h" //TSL_Friction_Freddii
#include "../../include/CZMs2/Friction-extrinsic.h"  //Extrinsice_friction

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

        Extrinsice_friction::ProblemParameters<2> par;
        ParameterAcceptor::initialize(argv[1]);
  
      
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        Extrinsice_friction::Extrinsice_friction<2> problem(par.degree,par);
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


