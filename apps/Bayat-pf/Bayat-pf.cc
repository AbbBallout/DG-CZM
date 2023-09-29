#include "../../include/Bayat-pf/Newton-DG-15.h"
#include "../../include/Bayat-pf/Newton-DG-elasticity.h"
#include "../../include/Bayat-pf/Newton-DG-elasticity-CZM.h"
#include "../../include/Bayat-pf/Newton-CG-elasticity-CZM.h"
#include "../../include/Bayat-pf/Laplace-DG.h"
#include "../../include/Bayat-pf/Newton-DG-elasticity-CZM-parallel.h"
#include "../../include/Bayat-pf/Newton-DG-elasticity-CZM-Kinsol.h"
#include "../../include/Bayat-pf/Newton-DG-elasticity-CZM-simplified.h"



// for $Q_1$, $Q_2$

 int main(int argc, char** argv)
{   
    using namespace dealii; 
   // deallog.depth_console(1);


    try
     {
      if (argc == 1)
      { 
        std::cout << "argc = " << argc << "  No input paramter passed \n " ;  
        exit(1); 
      }          
          DG_Elasticity_CZM::ProblemParameters<2> par;
          ParameterAcceptor::initialize(argv[1]);
  
      
       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        DG_Elasticity_CZM::DG_Elasticity_CZM<2> problem(par.degree,par);
        problem.run(); 


        // DG_Elasticity_CZM_simplified::DG_Elasticity_CZM<2> problem(1);
        //  problem.run(); 

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


