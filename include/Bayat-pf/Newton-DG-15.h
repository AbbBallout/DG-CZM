#ifndef STEP15DG
#define STEP15DG


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
 
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
#include <deal.II/fe/fe_q.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
 
 #include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
 
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
 
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>
#include <fstream>
#include <iostream>
 #include <deal.II/lac/sparse_direct.h>
 
#include <deal.II/numerics/solution_transfer.h>
 
namespace Step15
{
  using namespace dealii;
 
 
 
 double get_penalty_factor(const unsigned int fe_degree,
                          const double       cell_extent_left,
                          const double       cell_extent_right)
{
  const unsigned int degree = std::max(1U, fe_degree);
  return degree * (degree + 1.) * 0.5 *
         (1. / cell_extent_left + 1. / cell_extent_right);
}

struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> joint_dof_indices;
  std::array<double, 2>                values;
  std::array<unsigned int, 2>          cell_indices;
};
 
 
 
struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;
  double                               value;
  unsigned int                         cell_index;
 
 
  template <class Iterator>
  void reinit(const Iterator &cell, const unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};




  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();
    void run();
 
  private:
    void   setup_system(const bool initial_step);
    void   assemble_system();
    void   solve();
    void   interpolate_solution();
    void   set_boundary_values();
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
    void   output_results(const unsigned int refinement_cycle) const;
 
    Triangulation<dim> triangulation;
 
    DoFHandler<dim> dof_handler;
    FE_DGQ<dim>       fe;
 
    AffineConstraints<double> hanging_node_constraints;
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;

  const QGauss<dim>     quadrature;
  const QGauss<dim - 1> face_quadrature;
  const QGauss<dim>     quadrature_overintegration;
  const QGauss<dim - 1> face_quadrature_overintegration;
  const unsigned int    degree;
 const double diffusion_coefficient = 1.;
  const MappingQ1<dim>  mapping;
  using ScratchData = MeshWorker::ScratchData<dim>;
  };
 
 
 
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };
 
 
  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p,
                                    const unsigned int /*component*/) const
  {
    if(p[1]<0)
    return 0 ;
    
    return 1;
  }
 
 
 
 
  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
  : degree(1)
  , quadrature(degree + 1)
  , face_quadrature(degree + 1)
  , quadrature_overintegration(degree + 2)
  , face_quadrature_overintegration(degree + 2)
  , fe(degree)
  , dof_handler(triangulation)
  , mapping()
  {}
 
 
 
 
  template <int dim>
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
  { 



  



    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(dof_handler.n_dofs());
 
        hanging_node_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        hanging_node_constraints.close();
      }
 
 
 
    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
 
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
 
    hanging_node_constraints.condense(dsp);
 
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }
 
 
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system()
{
    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v          = scratch_data.reinit(cell);
        const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
        copy_data.reinit(cell, dofs_per_cell);
 
        const auto &       q_points    = scratch_data.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const std::vector<double> &JxW = scratch_data.get_JxW_values();

     std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);
        fe_v.get_function_gradients(current_solution,
                                         old_solution_gradients);
 
 
        for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                copy_data.cell_matrix(i, j) +=
                  diffusion_coefficient *     // nu
                  fe_v.shape_grad(i, point) * // grad v_h
                  fe_v.shape_grad(j, point) * // grad u_h
                  JxW[point];                 // dx
 
              copy_data.cell_rhs(i) += -                      
                  diffusion_coefficient *     // nu
                  fe_v.shape_grad(i, point) * // grad v_h
                  old_solution_gradients[point] * // grad u_h
                  JxW[point];                 // dx
            }
      };
 
    const auto boundary_worker = [&](const auto &        cell,
                                     const unsigned int &face_no,
                                     auto &              scratch_data,
                                     auto &              copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);
 
      const auto &       q_points      = scratch_data.get_quadrature_points();
      const unsigned int n_q_points    = q_points.size();
      const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;
 
      const std::vector<double> &        JxW = scratch_data.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        scratch_data.get_normal_vectors();
 
      std::vector<double> g(n_q_points);
       BoundaryValues<dim> boundary;
       boundary.value_list(q_points, g);
      
      const double extent1 = cell->measure() / cell->face(face_no)->measure();
      const double penalty = get_penalty_factor(degree, extent1, extent1);

      std::vector<double> old_solutions(n_q_points);
        fe_fv.get_function_values(current_solution,
                                         old_solutions);

      std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);
        fe_fv.get_function_gradients(current_solution,
                                         old_solution_gradients);

      for (const auto &face : cell->face_iterators())
      if ( (face->boundary_id() == 2) || (face->boundary_id() == 3) )
      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=
                (-diffusion_coefficient *        // - nu
                   fe_fv.shape_value(i, point) * // v_h
                   (fe_fv.shape_grad(j, point) * // (grad u_h .
                    normals[point])              //  n)
 
                 - diffusion_coefficient *         // - nu
                     (fe_fv.shape_grad(i, point) * // (grad v_h .
                      normals[point]) *            //  n)
                     fe_fv.shape_value(j, point)   // u_h
 
                 + diffusion_coefficient * penalty * // + nu sigma
                     fe_fv.shape_value(i, point) *   // v_h
                     fe_fv.shape_value(j, point)     // u_h
 
                 ) *
                JxW[point]; // dx
 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            copy_data.cell_rhs(i) += 
            -(
              -diffusion_coefficient *        // - nu
                   fe_fv.shape_value(i, point) * // v_h
                   (old_solution_gradients[point] * // (grad u_h .
                    normals[point])              //  n)
 
                 - diffusion_coefficient *         // - nu
                     (fe_fv.shape_grad(i, point) * // (grad v_h .
                      normals[point]) *            //  n)
                     old_solutions[point]   // u_h
 
                 + diffusion_coefficient * penalty * // + nu sigma
                     fe_fv.shape_value(i, point) *   // v_h
                     old_solutions[point]     // u_h

             )*JxW[point] 


              +(-diffusion_coefficient *        // - nu
                 (fe_fv.shape_grad(i, point) * // (grad v_h .
                  normals[point]) *            //  n)
                 g[point]                      // g
 
 
               + diffusion_coefficient * penalty *        // + nu sigma
                   fe_fv.shape_value(i, point) * g[point] // v_h g
 
               ) *
              JxW[point]; // dx
        }
    };
 
    const auto face_worker = [&](const auto &        cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &        ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 auto &              scratch_data,
                                 auto &              copy_data) {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);
 
      copy_data.face_data.emplace_back();
      CopyDataFace &     copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      copy_data_face.cell_rhs.reinit(n_dofs_face);

 
      const std::vector<double> &        JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
       const auto &       q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = get_penalty_factor(degree, extent1, extent2);

       std::vector<Tensor<1,dim>> average(n_q_points);
        fe_iv.get_average_of_function_gradients(current_solution,
                                         average);

       std::vector<double> jumps(n_q_points);
        fe_iv.get_jump_in_function_values(current_solution,
                                         jumps);

      for (const unsigned int point : fe_iv.quadrature_point_indices())
        {
          for (const unsigned int i : fe_iv.dof_indices())
           { for (const unsigned int j : fe_iv.dof_indices())
              copy_data_face.cell_matrix(i, j) += 
                (-diffusion_coefficient *                     // - nu
                   fe_iv.jump_in_shape_values(i, point) *     // [v_h]
                   (fe_iv.average_of_shape_gradients(j,       
                                                     point) * // ({grad u_h} .
                    normals[point])                           //  n)
 
                 -
                 diffusion_coefficient *                         // - nu
                   (fe_iv.average_of_shape_gradients(i, point) * // (grad v_h .
                    normals[point]) *                            //  n)
                   fe_iv.jump_in_shape_values(j, point)          // [u_h]
 
                 + diffusion_coefficient * penalty *        // + nu sigma
                     fe_iv.jump_in_shape_values(i, point) * // [v_h]
                     fe_iv.jump_in_shape_values(j, point)   // [u_h]
 
                 ) *
                JxW[point]; // dx

            copy_data_face.cell_rhs(i)= 
                  -(-diffusion_coefficient *                     // - nu
                   fe_iv.jump_in_shape_values(i, point) *     // [v_h]
                   (average[point] * // ({grad u_h} .
                    normals[point])                           //  n)
 
                 -
                 diffusion_coefficient *                         // - nu
                   (fe_iv.average_of_shape_gradients(i, point) * // (grad v_h .
                    normals[point]) *                            //  n)
                   jumps[point]         // [u_h]
 
                 + diffusion_coefficient * penalty *        // + nu sigma
                     fe_iv.jump_in_shape_values(i, point) * // [v_h]
                     jumps[point]   // [u_h]
 
                 ) *
                JxW[point]; // dx                                // - nu
      
           }
        }
    };
 
    AffineConstraints<double> constraints;
    constraints.close();
    const auto copier = [&](const auto &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);
 
      for (auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };
 
 
    const UpdateFlags cell_flags = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;
    const UpdateFlags face_flags = update_values | update_gradients |
                                   update_quadrature_points |
                                   update_normal_vectors | update_JxW_values;
 
    ScratchData scratch_data(
     mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData copy_data;
 
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);


  }
 
  
 
 
 
 
  template <int dim>
  void MinimalSurfaceProblem<dim>::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 system_rhs.l2_norm() * 1e-6);
   SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update,system_rhs);
 
    hanging_node_constraints.distribute(newton_update);
 
    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);
  }


 
  template <int dim>
  double MinimalSurfaceProblem<dim>::determine_step_length() const
  {
    return 1.0;
  }
 
 
 
 
  template <int dim>
  void MinimalSurfaceProblem<dim>::output_results(
    const unsigned int refinement_cycle) const
  {




    const std::string filename = "sol_Q" + Utilities::int_to_string(degree, 1) +
                               "-" + Utilities::int_to_string(refinement_cycle, 2) +
                               ".vtu";
  std::ofstream output(filename);
 
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(current_solution, "u", DataOut<dim>::type_dof_data);
     data_out.add_data_vector(newton_update, "update", DataOut<dim>::type_dof_data);
  data_out.build_patches(mapping);
  data_out.write_vtu(output);
  }
 
 
 
  template <int dim>
  void MinimalSurfaceProblem<dim>::run()
  {
     GridGenerator::hyper_cube(triangulation,-1,1,true);
     triangulation.refine_global(4);
 
     setup_system(/*first time=*/true);
     
 
    double       last_residual_norm = std::numeric_limits<double>::max();
    unsigned int refinement_cycle   = 0;
    do
      {
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl;
 
        if (refinement_cycle != 0)
             setup_system(false);

 
        
 
        for (unsigned int inner_iteration = 0; inner_iteration < 5;
             ++inner_iteration)
          {
            assemble_system();
            last_residual_norm = system_rhs.l2_norm();
 
            solve();
 
            std::cout << "  Residual: " << last_residual_norm << std::endl;

          }
 
        output_results(refinement_cycle);
 
        ++refinement_cycle;
        std::cout << std::endl;
      }
    while (last_residual_norm > 1e-3 && refinement_cycle<10);
   }
} // namespace Step15
#endif