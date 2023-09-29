#include "../../include/Nguyen/CG_Elasticity_Linear.h"   // namespace  CG_Elasticity_Linear
#include "../../include/Nguyen/CG_Elasticity_Nonlinear.h" //namespace  CG_Elasticity_Nonlinear
#include "../../include/Nguyen/DG_Elasticity_Linear.h"  // namespace  DG_Elasticity_Linear
#include "../../include/Nguyen/DG_Elasticity_Nonlinear.h" //namespace  DG_Elasticity_Nonlinear
#include "../../include/Nguyen/Nguyen.h" //namespace  Nguyen
#include  "../../include/Nguyen/ModeIdelamination_Bayat.h"  //napespace \\BayatIDela

// for $Q_1$, $Q_2$

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
        DG_Elasticity_Linear::ProblemParameters<2> par;
        ParameterAcceptor::initialize(argv[1]);
  
      
         Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        DG_Elasticity_Linear::DG_Elasticity_Linear<2> problem(par.degree,par);
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


