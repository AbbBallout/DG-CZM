#ifndef ARCTRUSS_H
#define ARCTRUSS_H


#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
 
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
 
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_values.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
 
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>

using namespace dealii; 

  template <int spacedim>
  void right_hand_side(const std::vector<Point<spacedim>> &points,
                       std::vector<Tensor<1, spacedim>> &  values)
  {
    // AssertDimension(values.size(), points.size());
    // Assert(dim >= 2, ExcNotImplemented());
 
    Point<spacedim> point;
    point(0) = 0.5;
    const double load = 100.0 , tol=1e-3;
     
          
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
     {
         if ( std::abs(points[point_n][0]-0.5)<tol)
          {values[point_n][1] = -load/2 ;
           std::cout<<"\called \n" ; }
        else
          values[point_n][1] = 0.0;
      }

  }
  
template <int dim , int spacedim>
  class TrussProblem
  {
  public:
    TrussProblem();
    void run();
 
  private:
    void create_mesh();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results(const unsigned int cycle) const;
 
    Triangulation<dim,spacedim> triangulation;
    DoFHandler<dim,spacedim>    dof_handler;
 
    FESystem<dim,spacedim> fe;
 
    AffineConstraints<double> constraints;
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> solution;
    Vector<double> system_rhs;
  };

 template <int dim , int spacedim>
  TrussProblem<dim,spacedim>::TrussProblem()
    : dof_handler(triangulation)
    , fe(FE_Q<dim,spacedim>(1), spacedim)
  {}


template <int dim , int spacedim>
void TrussProblem<dim,spacedim>::create_mesh()
{

  GridIn<dim , spacedim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream reading("truss-arc.msh");
  gridin.read_msh(reading);
 

  triangulation.refine_global(3);
   

     std::cout << "Mesh info:" << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;

    // std::map<types::boundary_id, unsigned int> boundary_count;
    // for (const auto &cell : triangulation.active_cell_iterators())
     

   std::ofstream out("grid-1.vtu");
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << "grid-1.vtu" << std::endl << std::endl;
}

  template <int dim , int spacedim>
  void TrussProblem<dim, spacedim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
 
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<spacedim>(spacedim),
                                             constraints);

        VectorTools::interpolate_boundary_values(dof_handler,
                                             1,
                                             Functions::ZeroFunction<spacedim>(spacedim),
                                             constraints);
    constraints.close();
 
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);
 
    system_matrix.reinit(sparsity_pattern);
  }

  template <int dim , int spacedim>
  void TrussProblem<dim,spacedim>::assemble_system()
  {
    QGauss<dim> quadrature_formula(fe.degree + 1);
 
    FEValues<dim, spacedim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);
 
    Functions::ConstantFunction<spacedim> lambda(1.), mu(1.);
 
    std::vector<Tensor<1, spacedim>> rhs_values(n_q_points);
 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;
 
        fe_values.reinit(cell);
 
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);
 
        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
 
            for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j =
                  fe.system_to_component_index(j).first;
 
                for (const unsigned int q_point :
                     fe_values.quadrature_point_indices())
                  {
                    cell_matrix(i, j) +=
                      (                                                  
                        (fe_values.shape_grad(i, q_point)[component_i] * 
                         fe_values.shape_grad(j, q_point)[component_j] * 
                         lambda_values[q_point])                         
                        +                                                
                        (fe_values.shape_grad(i, q_point)[component_j] * 
                         fe_values.shape_grad(j, q_point)[component_i] * 
                         mu_values[q_point])                             
                        +                                                
                        ((component_i == component_j) ?        
                           (fe_values.shape_grad(i, q_point) * 
                            fe_values.shape_grad(j, q_point) * 
                            mu_values[q_point]) :              
                           0)                                  
                        ) *                                    
                      fe_values.JxW(q_point);                  
                  }
              }
          }
 
        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
 
            for (const unsigned int q_point :
                 fe_values.quadrature_point_indices())
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             rhs_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }
 
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }

  template <int dim , int spacedim>
  void TrussProblem<dim,spacedim>::solve()
  {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);
 
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
 
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
 
    constraints.distribute(solution);
  }

  template <int dim , int spacedim>
  void TrussProblem<dim,spacedim>::output_results(const unsigned int cycle) const
  {  



         std::vector<std::string> solution_names(spacedim, "u");
 
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(spacedim,
                   DataComponentInterpretation::component_is_part_of_vector);
  DataOut<dim,spacedim> data_out;
  data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           interpretation);

  data_out.build_patches();
  std::ofstream output("solution" + std::to_string(cycle) + ".vtu");
  data_out.write_vtu(output);

  }


template <int dim , int spacedim>
void TrussProblem<dim, spacedim>::run()
{
  create_mesh();
  setup_system();
  assemble_system();
  solve();
  output_results(0); 
}






#endif
