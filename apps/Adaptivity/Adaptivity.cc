
#include "../../include/Adaptivity/Friction-adaptivity.h" //Friction_adaptivity
#include "../../include/Adaptivity/Friction_adaptivity_everywhere.h" //Friction_adaptivity_everywhere

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

        Friction_adaptivity_everywhere::ProblemParameters<2> par;
        ParameterAcceptor::initialize(argv[1]);
  
      
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        Friction_adaptivity_everywhere::Friction_adaptivity_everywhere<2> problem(par);
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


