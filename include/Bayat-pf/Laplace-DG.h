
#ifndef  STEP74_H
#define  STEP74_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
 
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
 
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>
 
 #include <deal.II/lac/constraint_matrix.h>
namespace Step74
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
  class SIPGLaplace
  {
  public:
    SIPGLaplace();
    void run();
 
  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;
 
    void   compute_errors();
    void   compute_error_estimate();
    double compute_energy_norm_error();
 
    Triangulation<dim>    triangulation;
    const unsigned int    degree;
    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> face_quadrature;
    const MappingQ1<dim>  mapping;
 
    using ScratchData = MeshWorker::ScratchData<dim>;
 
    const FE_DGQ<dim> fe;
    DoFHandler<dim>   dof_handler;
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;

 
    ConvergenceTable convergence_table;
 
    const double diffusion_coefficient = 1.;
 
    
  };
 
  template <int dim>
  SIPGLaplace<dim>::SIPGLaplace()
    : degree(1)
    , quadrature(degree + 1)
    , face_quadrature(degree + 1)
    , mapping()
    , fe(degree)
    , dof_handler(triangulation)
  {
  }
 
 
 
  template <int dim>
  void SIPGLaplace<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
 
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }
 
 
 
  template <int dim>
  void SIPGLaplace<dim>::assemble_system()
  {
    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v          = scratch_data.reinit(cell);
        const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
        copy_data.reinit(cell, dofs_per_cell);
 
        const auto &       q_points    = scratch_data.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const std::vector<double> &JxW = scratch_data.get_JxW_values();
 

 
        for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                copy_data.cell_matrix(i, j) +=
                  diffusion_coefficient *     // nu
                  fe_v.shape_grad(i, point) * // grad v_h
                  fe_v.shape_grad(j, point) * // grad u_h
                  JxW[point];                 // dx
 
              copy_data.cell_rhs(i) += 0.0;
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
 


 
      const double extent1 = cell->measure() / cell->face(face_no)->measure();
      const double penalty = get_penalty_factor(degree, extent1, extent1);
     
      
  if((cell->face(face_no)->boundary_id() == 0) || (cell->face(face_no)->boundary_id() == 1)) 
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
              (-diffusion_coefficient *        // - nu
                 (fe_fv.shape_grad(i, point) * // (grad v_h .
                  normals[point]) *            //  n)
                 (cell->face(face_no)->boundary_id())                     // g
 
 
               + diffusion_coefficient * penalty *        // + nu sigma
                   fe_fv.shape_value(i, point) *(cell->face(face_no)->boundary_id()) // v_h g
 
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
                                 auto &              copy_data)
      {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);
 
      copy_data.face_data.emplace_back();
      CopyDataFace &     copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
 
      const std::vector<double> &        JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
 
      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = get_penalty_factor(degree, extent1, extent2);
 
      for (const unsigned int point : fe_iv.quadrature_point_indices())
        {
          for (const unsigned int i : fe_iv.dof_indices())
            for (const unsigned int j : fe_iv.dof_indices())
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
        }
    };
 
    AffineConstraints<double> constraints;
    constraints.close();
    const auto copier = [&](const auto &c)
    {
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
  void SIPGLaplace<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
     
  }
 
 
 
  template <int dim>
  void SIPGLaplace<dim>::output_results(const unsigned int cycle) const
  {
    const std::string filename = "sol_Q" + Utilities::int_to_string(degree, 1) +
                                 "-" + Utilities::int_to_string(cycle, 2) +
                                 ".vtu";
    std::ofstream output(filename);
 
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);
    data_out.build_patches(mapping);
    data_out.write_vtu(output);
  }
 
 
 
 
  template <int dim>
  void SIPGLaplace<dim>::run()
  {



        GridGenerator::hyper_cube(triangulation,0,1,true);


        triangulation.refine_global(4);

    {   
         typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						                               endc = dof_handler.end();
                                           
        for (; cell != endc; ++cell)
        {
          if ((cell->center()[1]) > 0.9 )
            {
              if ((cell->center()[0] > 0.9)  || (cell->center()[0] < 0.1))
                cell->set_refine_flag();
            }
        }
     }
  {
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						                               endc = dof_handler.end();
      triangulation.execute_coarsening_and_refinement();

        for (; cell != endc; ++cell)
        {
          if ((cell->center()[1]) > 0.9 )
            {
              if ((cell->center()[0] > 0.9)  || (cell->center()[0] < 0.1))
                cell->set_refine_flag();
            }
        }
      triangulation.execute_coarsening_and_refinement();
  }

 
        std::cout << "  Number of active cells       : "
                  << triangulation.n_active_cells() << std::endl;
        setup_system();
 
        std::cout << "  Number of degrees of freedom : " << dof_handler.n_dofs()
                  << std::endl;
 
        assemble_system();
        solve();
        output_results(1);
       
        std::cout << std::endl;
      
 
  }
} // namespace Step74
 
 




#endif 