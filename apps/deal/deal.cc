
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;


struct CopyData
{
};

int main(int argc, char** argv)
{   

    MPI_Comm mpi_communicator(MPI_COMM_WORLD);
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
    ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0));

     parallel::distributed::Triangulation<2>  tri(mpi_communicator ) ;
    Point<2> P1, P2(1, 2);
    std::vector<unsigned int> repetitions{1, 2};
    GridGenerator::subdivided_hyper_rectangle(tri, repetitions, P1, P2, true);
    tri.refine_global(3);
    for (const auto &cell : tri.active_cell_iterators())
    {
        Point<2> p = cell->center();
        if (p[1] > 1.0)
            cell->set_material_id(1);
        else
            cell->set_material_id(2);
    }
    
    FE_DGQ<2> fe(1);
    const QGauss<2> quadrature(2);
    const QGauss<1> face_quadrature(2);
    const MappingQ1<2, 2> mapping;
    DoFHandler<2> dof_handler(tri);
    dof_handler.distribute_dofs(fe);
    const UpdateFlags cell_flags =  update_quadrature_points;
    const UpdateFlags face_flags = update_quadrature_points;

    MeshWorker::ScratchData<2> scratch_data(mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData copy_data;


    const auto face_worker = [&](const auto &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 MeshWorker::ScratchData<2> &scratch_data,
                                 CopyData &copy_data)
    {  
       
       const FEInterfaceValues<2> &fe_iv = scratch_data.reinit(cell, f, sf, ncell, nf, nsf);


      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();


      if(cell->material_id() != ncell->material_id())
      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          cell->set_refine_flag();
          ncell->set_refine_flag();
        }

    };

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          nullptr,
                          nullptr,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                          nullptr,
                          face_worker);

      tri.execute_coarsening_and_refinement();
      DataOut<2> data_out;
      data_out.attach_dof_handler(dof_handler);

      Vector<double> subdomain(tri.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = tri.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      data_out.build_patches();

      data_out.write_vtu_with_pvtu_record("", "solution", 1, mpi_communicator, 3, 0);
     
}


#if 0

struct CopyData
{
};

int main(int argc, char** argv)
{   

    MPI_Comm mpi_communicator(MPI_COMM_WORLD);
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
    ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0));

    parallel::distributed::Triangulation<2> tri(mpi_communicator);
    Point<2> P1, P2(1, 2);
    std::vector<unsigned int> repetitions{1, 2};
    GridGenerator::subdivided_hyper_rectangle(tri, repetitions, P1, P2, true);
    tri.refine_global(3);
    for (const auto &cell : tri.active_cell_iterators())
    {
        Point<2> p = cell->center();
        if (p[1] > 1.0)
            cell->set_material_id(1);
        else
            cell->set_material_id(2);
    }
    FE_DGQ<2> fe(1);
    const QGauss<2> quadrature(2);
    const QGauss<1> face_quadrature(2);
    const MappingQ1<2, 2> mapping;
    DoFHandler<2> dof_handler(tri);
    dof_handler.distribute_dofs(fe);
    double face_counter=0;
    const UpdateFlags cell_flags = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;
    const UpdateFlags face_flags = update_values | update_gradients |
                                   update_quadrature_points |
                                   update_normal_vectors | update_JxW_values;

    MeshWorker::ScratchData<2> scratch_data(mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData copy_data;

    const auto cell_worker =
        [&](const auto &cell, auto &scratch_data, auto &copy_data) {};

    const auto boundary_worker = [&](const auto &cell,
                                     const unsigned int &face_no,
                                     MeshWorker::ScratchData<2> &scratch_data,
                                     CopyData &copy_data) {};
    const auto copier = [&](const auto &c) {};

    const auto face_worker = [&](const auto &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 MeshWorker::ScratchData<2> &scratch_data,
                                 CopyData &copy_data)
    {  
       
       face_counter++; 
    };

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                          boundary_worker,
                          face_worker);
   face_counter = Utilities::MPI::sum(face_counter, mpi_communicator);                       
   pcout<< "face_counter " << face_counter << "\n" ;
}
#endif


#if 0
template <int dim>
class LaplaceProblem
{
public:
    LaplaceProblem();

    void run();

private:
    void setup_system();

    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
    TimerOutput computing_timer;

};

template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
    : mpi_communicator(MPI_COMM_WORLD), pcout(std::cout,
                                              (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
{
}

template <int dim>
void LaplaceProblem<dim>::setup_system()
{
    TimerOutput::Scope t(computing_timer, "setup");
}
template <int dim>
void LaplaceProblem<dim>::run()
{
    setup_system();
    computing_timer.print_summary();
    computing_timer.reset();
}

int main(int argc, char *argv[])
{

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    LaplaceProblem<2> laplace_problem_2d;
    laplace_problem_2d.run();
}
#endif

#if 0
  template <int dim>
    class ProblemParameters : public ParameterAcceptor
    {
    public:
        ProblemParameters() : ParameterAcceptor("main")
        {
  
            add_parameter("para_i", para_i, " ",this->prm, Patterns::Integer());
            add_parameter("para_f", para_f, " ",this->prm, Patterns::Double());
            add_parameter("cycles", cycles, " ",this->prm, Patterns::Integer());
            add_parameter("fracture_at", fracture_at, " ",this->prm, Patterns::Integer());
            add_parameter("t_0", t_0, " ",this->prm, Patterns::Double());
            add_parameter("step_length", step_length, " ",this->prm, Patterns::Double());
    
        }

      
        int para_i ;
        double para_f;
        int cycles;
        int fracture_at;
        double t_0;
        double step_length;



    };


int main()
{          ProblemParameters<2> par;
          ParameterAcceptor::initialize("input.prm");
         std::cout<<par.step_length<< " \n";
}
#endif

#if 0
// Can't interpolate DG to tthe boundary
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
using namespace dealii;



int main()
{
    // create an feq obkect and interpolate a value of 1 at boundary 0
    {
        Triangulation<2> tria;
        GridGenerator::hyper_cube(tria, 0, 1, true);
        tria.refine_global(4);
        FE_Q<2> fe(1);
        DoFHandler<2> dof(tria);

        Vector<double> solution;

        dof.distribute_dofs(fe);
        solution.reinit(dof.n_dofs());

        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof,
                                                 0,
                                                 Functions::ConstantFunction<2>(1),
                                                 boundary_values);
        for (auto &boundary_value : boundary_values)
            solution(boundary_value.first) = boundary_value.second;

        DataOut<2> data_out;

        data_out.attach_dof_handler(dof);
        data_out.add_data_vector(solution, "solution");
        data_out.build_patches();

        const std::string filename = "solution-feq.vtu";
        std::ofstream output(filename);
        data_out.write_vtu(output);
    }

    // Same as above but
    // we create an feGQq obkect and interpolate
    {
        Triangulation<2> tria;
        GridGenerator::hyper_cube(tria, 0, 1, true);
        tria.refine_global(4);
        FE_DGQ<2> fe(1);                 // only difference between this and the above code 
        DoFHandler<2> dof(tria);

        Vector<double> solution;

        dof.distribute_dofs(fe);
        solution.reinit(dof.n_dofs());

        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof,
                                                 0,
                                                 Functions::ConstantFunction<2>(1),
                                                 boundary_values);
        for (auto &boundary_value : boundary_values)
            solution(boundary_value.first) = boundary_value.second;

        DataOut<2> data_out;

        data_out.attach_dof_handler(dof);
        data_out.add_data_vector(solution, "solution");
        data_out.build_patches();

        const std::string filename = "solution-fedgq.vtu";
        std::ofstream output(filename);
        data_out.write_vtu(output);
    }
}

#endif