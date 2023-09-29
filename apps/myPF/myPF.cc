#include "../../include/myPF/my_phase_field.h"


  int main(int argc, char** argv)
{   
    using namespace dealii; 
   // deallog.depth_console(1);
   
    try
     {
     Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      
      if (argc == 1)
      { 
        std::cout << "argc = " << argc << "  No input paramter passed \n " ;  
        exit(1); 
      }          
          ProblemParameters<2> par;
          ParameterAcceptor::initialize("input.prm");
          par.setMoreParameters() ; 
          PhaseField<2> elastic_problem(par);
          elastic_problem.run();
         
        //  using namespace phasefieldSNES;

         // Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        // TopLevel<2> elastic_problem;
        // elastic_problem.run();
         

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

