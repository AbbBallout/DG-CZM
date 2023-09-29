#ifndef myPF
#define myPF

// Parallel output is bad fix it
// Step 18 holds the entire mesh. use step 40 
// Adaptive mesh
// Write documentation
// Assembly with threads
// iterative solver 
// Steps 15, 18, 20 and 55 were used
// Step 44 for extension ( parallel::shared::Triangulation to  parallel::distributed::Triangulation )

// Start from Steps 40, 44, 55 


#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/numerics/data_out.h>


#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>


namespace myPF
{
    using namespace dealii;

    template <int dim>
    class ProblemParameters : public ParameterAcceptor
    {
    public:
        ProblemParameters() : ParameterAcceptor("main")
        {
            add_parameter("Poisson", Poisson, " ");
            add_parameter("Young", Young, "");
            add_parameter("k_res", k_res, " ");
            add_parameter("alpha_0", alpha_0, " ");
            add_parameter("ell", ell, " ", this->prm, Patterns::Integer());
            add_parameter("Gc", Gc, " ");
            add_parameter("error_tol", error_tol, " ");
            add_parameter("damage_type", damage_type, " ", this->prm, Patterns::Integer(0, 2));
        }

        void setMoreParameters()
        {
            if (damage_type == 1)
                c_w = 2.;
            else if (damage_type == 2)
                c_w = 8. / 3.;
        }

        double Poisson = 0.0;
        double Young = 0.0;
        double k_res = 0.0;
        double alpha_0 = 0.0;
        int ell = 0;
        double Gc = 0.0;
        double error_tol = 0.0;
        int damage_type = 0;

        double c_w = 0.69;
    };

    // step 28
    class Data
    {
    };

    template <int dim>
    class PhaseField
    {
    public:
        PhaseField(const ProblemParameters<dim> &par);
        ~PhaseField();
        void run();

    private:
        void setup_system();
        void assemble_system();
        void solve();
        void refine_grid();
        void output_results(const unsigned int cycle) const;

        parallel::distributed::Triangulation<dim> triangulation;
        FESystem<dim> fe;
        DoFHandler<dim> dof_handler;

        AffineConstraints<double> constraints;

        const QGauss<dim> quadrature_formula;

        PETScWrappers::MPI::BlockSparseMatrix system_matrix;
        PETScWrappers::MPI::BlockVector system_rhs;
        PETScWrappers::MPI::BlockVector solution_update;
        PETScWrappers::MPI::BlockVector current_solution;


        // Timestep parameters
        double real_time = 0.;
        double real_time_step = 1;
        double real_time_final = 30;
        unsigned int real_timestep_no = 0;

        std::vector<IndexSet> owned_partitioning;
        std::vector<IndexSet> relevant_partitioning;

        MPI_Comm mpi_communicator;
        const unsigned int n_mpi_processes;
        const unsigned int this_mpi_process;

        ConditionalOStream pcout;
        TimerOutput computing_timer;

        const ProblemParameters<dim> &par;
    };

    template <int dim>
    PhaseField<dim>::PhaseField(const ProblemParameters<dim> &par)
        : par(par),
          triangulation(MPI_COMM_WORLD),  // With adaptive mesh refinement check stepp 55 
          fe(FE_Q<dim>(1), dim, FE_Q<dim>(1), 1),
          dof_handler(triangulation),
          quadrature_formula(fe.degree + 1),
          mpi_communicator(MPI_COMM_WORLD),
          n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)),
          this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)),
          pcout(std::cout, this_mpi_process == 0),
          computing_timer(pcout, TimerOutput::never, TimerOutput::wall_times)
    {
    }

    template <int dim>
    PhaseField<dim>::~PhaseField()
    {
        dof_handler.clear();
    }

    template <int dim>
    void PhaseField<dim>::setup_system()
    {
        TimerOutput::Scope t(computing_timer, "Setup system");

        dof_handler.distribute_dofs(fe);

    
    // see how much of this code should be in initial step 
    std::vector<unsigned int> phasefield_sub_blocks(dim + 1, 0);
    phasefield_sub_blocks[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, phasefield_sub_blocks);
    
    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, phasefield_sub_blocks);
    
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_s = dofs_per_block[1];
    
    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
      << n_u << '+' << n_s << ')' << std::endl;

    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u);
    owned_partitioning[1] = dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_s);  
    const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_s); 


  // For applying BC 
  // Why are you using brackets? 
   {
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler,constraints);
    constraints.close();
   }

 {
     // Coupling table taken form step 55 as is     

    //   Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    //   for (unsigned int c = 0; c < dim + 1; ++c)
    //     for (unsigned int d = 0; d < dim + 1; ++d)
    //       if (c == dim && d == dim)
    //         coupling[c][d] = DoFTools::none;
    //       else if (c == dim || d == dim || c == d)
    //         coupling[c][d] = DoFTools::always;
    //       else
    //         coupling[c][d] = DoFTools::none;
    
     BlockDynamicSparsityPattern   dsp(relevant_partitioning);

    DoFTools::make_sparsity_pattern(dof_handler,
                                   dsp,
                                  constraints,
                                /*keep constrained dofs*/ false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                           dof_handler.locally_owned_dofs(),
                                           mpi_communicator,
                                           locally_relevant_dofs);

    system_matrix.reinit(owned_partitioning, dsp,mpi_communicator);
  }
    current_solution.reinit(owned_partitioning,relevant_partitioning,mpi_communicator);
    solution_update.reinit(owned_partitioning,relevant_partitioning,mpi_communicator);
    system_rhs.reinit(owned_partitioning, mpi_communicator);
 

    }

    template <int dim>
    void PhaseField<dim>::assemble_system()
    { 
        TimerOutput::Scope t(computing_timer, "assembly");
        system_matrix = 0;
        system_rhs    = 0;

        FEValues<dim> fe_values(fe,
                       quadrature_formula,
                       update_values | update_gradients |
                        update_quadrature_points | update_JxW_values);

        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
        const unsigned int n_q_points    = quadrature_formula.size();

        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     cell_rhs(dofs_per_cell);
        
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        const FEValuesExtractors::Vector displacements(0);
        const FEValuesExtractors::Scalar phase_field(dim);
    }

      template <int dim>
    void PhaseField<dim>::refine_grid()
    {
        // Adaptive mesh refinement goes here 
        pcout<< "No adaptive was executed refinement  \n" ;   

    }

    template <int dim>
    void PhaseField<dim>::run()
    { 
        for (unsigned int cycle = 0; cycle < 2; ++cycle)
        {
            pcout << "Cycle " << cycle << ':' << std::endl;
            if (cycle == 0)
            {
                pcout << "Starting with n of processes = " << n_mpi_processes << " \n";
                GridGenerator::hyper_cube(triangulation);
                triangulation.refine_global(3); 
            }
            else
            refine_grid();

            pcout << "Number of active cells:       "
                  << triangulation.n_active_cells() << " (by partition:";
            for (unsigned int p = 0; p < n_mpi_processes; ++p)
                pcout << (p == 0 ? ' ' : '+') << (GridTools::count_cells_with_subdomain_association(triangulation, p));
            pcout << ")" << std::endl;

            setup_system();
            assemble_system();
            
        }
        computing_timer.print_summary();
    }

} // close namespace

#endif
