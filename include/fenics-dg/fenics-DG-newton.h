#ifndef PAHSEFILEDDGNEWTON_H
#define PAHSEFILEDDGNEWTON_H

// Inspired from step 15

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/timer.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>


#include <iostream>
#include <fstream>
#include <deal.II/numerics/matrix_tools.h>
#include <map>

namespace pfdgnewton
{
 using namespace dealii;
  
  template <int dim>
class Lambda : public Function<dim>
{
public:
  Lambda(int material_Id)
  {
    material_Id= material_Id ;
  }
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;

  private:
  int material_Id ; 
};

template <int dim>
double Lambda<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
 
  if(material_Id != 1)
  return 1236.0 ; 
   else
   return 14285.714285;
}

template <int dim>
class Mu : public Function<dim>
{
public:
  Mu(int material_Id)
  {
    material_Id= material_Id ;
  }

  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;

    private:
  int material_Id ; 
};

template <int dim>
double Mu<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{ 
  if(material_Id != 1)
  return 1236.0 ;  
   else
   return 3571.428571428572;

}

template <int dim>
class GC : public Function<dim>
{
public:
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double GC<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
  if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <1e-3 ) 
      || 
      (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <1e-3 )   )
  return 0.5 ;
  else
  return 1.0 ; 
}

template <int dim>
class SIG : public Function<dim>
{
public:
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double SIG<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
  if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <1e-3  ) 
      || 
      (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <1e-3  )   )
  return 50.0 ;
  else
  return 400.0 ; 
}


 template <int dim>
    struct ScratchData
    {
        ScratchData(const Mapping<dim> &mapping,
                    const FiniteElement<dim> &fe,
                    const Quadrature<dim> &quadrature,
                    const Quadrature<dim - 1> &quadrature_face,
                    const UpdateFlags update_flags = update_values |
                                                     update_gradients |
                                                     update_quadrature_points |
                                                     update_JxW_values,
                    const UpdateFlags interface_update_flags =
                        update_values | update_gradients | update_quadrature_points |
                        update_JxW_values | update_normal_vectors)
            : fe_values(mapping, fe, quadrature, update_flags), fe_interface_values(mapping,
                                                                                    fe,
                                                                                    quadrature_face,
                                                                                    interface_update_flags)
        {
        }
        ScratchData(const ScratchData<dim> &scratch_data)
            : fe_values(scratch_data.fe_values.get_mapping(),
                        scratch_data.fe_values.get_fe(),
                        scratch_data.fe_values.get_quadrature(),
                        scratch_data.fe_values.get_update_flags()),
              fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                                  scratch_data.fe_interface_values.get_fe(),
                                  scratch_data.fe_interface_values.get_quadrature(),
                                  scratch_data.fe_interface_values.get_update_flags())
        {
        }
        FEValues<dim> fe_values;
        FEInterfaceValues<dim> fe_interface_values;
    };

    struct CopyDataFace
    {
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;
        std::vector<types::global_dof_index> joint_dof_indices;
    };
    struct CopyData
    {
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;
        std::vector<types::global_dof_index> local_dof_indices;
        std::vector<CopyDataFace> face_data;

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
    class BoundaryValues : public Function<dim>
    {
    public:
     BoundaryValues(double disp)
        : Function<dim>(dim),
         disp(disp)
      {}
      virtual void vector_value(const Point<dim> &p,
                                Vector<double> &  value) const override;
      private:
      double disp ; 
    };
    
    
    template <int dim>
    void BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                           Vector<double> &  values) const
    {
           Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
            
            if(p[0]>0.99)
            {values(0) = disp;
            values(1) = 0.0;
            }
            else
            {
              values(0) = 0.0;
            values(1) = 0.0;
            }
    }


  template <int dim>
  class PF_DG
  {
  public:
    PF_DG(const unsigned int degree);
    void run();
 
  private:
    void   import_mesh();
    void   setup_system(const bool initial_step);
    void   assemble_system();
    void   solve();
   void    set_boundary_values(); 
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
    void   output_results(const unsigned int cycle) const;
 
    Triangulation<dim> triangulation;
 
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe;
    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;
    const MappingQ1<dim> mapping;
    AffineConstraints<double> constraints;
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> current_solution;
    //Vector<double> previous_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;

    const unsigned int degree;

    ConvergenceTable convergence_table;
    TimerOutput computing_timer;
    unsigned int timestep_number;

    Tensor<1,dim> disp;
  };
 
 

 
 
 
 
  template <int dim>
  PF_DG<dim>::PF_DG(const unsigned int degree)
    : degree(degree), 
    quadrature(degree + 1), 
    face_quadrature(degree + 1),
    mapping(), 
    fe(FE_DGQ<dim>(degree), dim ), 
    dof_handler(triangulation), 
    computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)

{
      disp[0]=0.0;
      disp[1]=0.0;
}
 
 
   template <int dim>
  void PF_DG<dim>::import_mesh()
{
  GridIn<2> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("fenicsr2.5.msh");
  gridin.read_msh(f);

  //    std::ofstream out("grid-1.vtu");
  // GridOut       grid_out;
  // grid_out.write_vtu(triangulation, out);
  // std::cout << " written to " << "grid-1.vtu" << std::endl << std::endl;

}
 
  template <int dim>
  void PF_DG<dim>::setup_system(const bool initial_step)
  { 
    TimerOutput::Scope t(computing_timer, "set up");

    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(dof_handler.n_dofs());
        //previous_solution.reinit(dof_handler.n_dofs());
      }
 
 
    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
   
   

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
 
    sparsity_pattern.copy_from(dsp);
   
    system_matrix.reinit(sparsity_pattern);

    std::cout << "dofs = " << dof_handler.n_dofs() << std::endl ; 
  }
 
 
  template <int dim>
  void PF_DG<dim>::assemble_system()
  { 
    TimerOutput::Scope t(computing_timer, "Assemble");

     const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) 
      {  
        scratch_data.fe_values.reinit(cell);
        const FEValues<dim>   &fe_v = scratch_data.fe_values;

        const unsigned int   dofs_per_cell = scratch_data.fe_values.get_fe().n_dofs_per_cell();
        copy_data.reinit(cell, dofs_per_cell);
        const auto &       q_points    = scratch_data.fe_values.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const std::vector<double> &JxW = fe_v.get_JxW_values();




        std::vector<double> lambda_values(n_q_points);
        std::vector<double> mu_values(n_q_points);       

        Lambda<dim> lambda(cell->material_id()); 
        Mu<dim> mu(cell->material_id()); 

        lambda.value_list(scratch_data.fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(scratch_data.fe_values.get_quadrature_points(), mu_values);


         const FEValuesExtractors::Vector displacements(0);


        std::vector<Tensor<2,dim>> old_solution_gradients(n_q_points);
        fe_v[displacements].get_function_gradients(current_solution,old_solution_gradients);

        std::vector<double> old_solution_divergences(n_q_points);
        fe_v[displacements].get_function_divergences(current_solution,old_solution_divergences);

        for (unsigned int point = 0; point < n_q_points; ++point)            
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
           {

                    copy_data.cell_matrix(i, j) += 
                 (  
                (
                fe_v[displacements].divergence(i,point) * 
                   fe_v[displacements].divergence(j,point) * 
                   lambda_values[point])                         
                  +                                                
                  (scalar_product(fe_v[displacements].gradient(i,point),fe_v[displacements].gradient (j,point)) * 
                   mu_values[point])                             
                  +  
                  ( scalar_product(fe_v[displacements].gradient (i,point), transpose(fe_v[displacements].gradient (j,point)))*
                  mu_values[point])                                 
                  )*JxW[point];                 // dx                                                // dx
           }
             copy_data.cell_rhs(i) +=  
                 -(  
                (
                fe_v[displacements].divergence(i,point) * 
                   old_solution_divergences[point] * 
                   lambda_values[point])                         
                  +                                                
                  (scalar_product(fe_v[displacements].gradient(i,point),old_solution_gradients[point]) * 
                   mu_values[point])                             
                  +  
                 ( scalar_product(fe_v[displacements].gradient (i,point), transpose(old_solution_gradients[point]))*
                 mu_values[point])                                 
                  )*JxW[point];                   
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
      
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
       fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

      const auto &       q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      copy_data.face_data.emplace_back();



      CopyDataFace &     copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      copy_data_face.cell_rhs.reinit(n_dofs_face) ; 
      const std::vector<double> &        JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
         const FEValuesExtractors::Vector displacements(0);

         std::vector<Tensor<1,dim>> old_solution(n_q_points);
         fe_iv[displacements].get_jump_in_function_values(current_solution,old_solution);


        std::vector<double> GC_values(n_q_points);
        std::vector<double> sig_values(n_q_points);       

        GC<dim> gc;
        SIG<dim> sig;  

          Tensor<1,dim> tangential;
       

         gc.value_list(fe_iv.get_quadrature_points(), GC_values);
         sig.value_list(fe_iv.get_quadrature_points(), sig_values);
         double law, delta_0 , gmax;
         
         for (unsigned int point = 0; point < n_q_points; ++point)
         {
           delta_0=  GC_values[point]/sig_values[point]/std::exp(1);
           tangential = cross_product_2d(normals[point]);
           gmax= std::sqrt(old_solution[point]*normals[point]*old_solution[point]*normals[point]+ 
           2.0*2.0*old_solution[point]*tangential*old_solution[point]*tangential);
           law= GC_values[point]/(delta_0*delta_0) * std::exp(-gmax/delta_0); 

          for (unsigned int i = 0; i < n_dofs_face; ++i)  
          {   

             for (unsigned int j = 0; j < n_dofs_face; ++j) 
             {

                copy_data_face.cell_matrix(i, j) +=
                GC_values[point]/(delta_0*delta_0)*std::exp(-gmax*delta_0)*
                fe_iv[displacements].jump_in_values(i, point)* 
                   ( 
                      (-1/(delta_0*gmax+1e-10))*
                     (fe_iv[displacements].jump_in_values(j, point)*normals[point]*old_solution[point]*normals[point]+
                     2.0*2.0*fe_iv[displacements].jump_in_values(j, point)*tangential*old_solution[point]*tangential)*
                     old_solution[point]
                     +fe_iv[displacements].jump_in_values(j, point)
                   )          
                     *JxW[point] 
                 
                //   +law*fe_iv[displacements].jump_in_values(i, point)  *            
                //    fe_iv[displacements].jump_in_values(j, point)  *
                //  JxW[point]
                ; 

            }
            copy_data_face.cell_rhs(i)+=-law*                              
                   fe_iv[displacements].jump_in_values(i, point)              
                    *old_solution[point]*
                JxW[point]; 
          }
         } 

    };

        const auto boundary_worker = [&](const auto &        cell,
                                     const unsigned int &face_no,
                                     auto &              scratch_data,
                                     auto &              copy_data)
                                     
      {  scratch_data.fe_interface_values.reinit(cell, face_no);
      const FEFaceValuesBase<dim> &fe_fv =
        scratch_data.fe_interface_values.get_fe_face_values(0);

      const auto &       q_points      = fe_fv.get_quadrature_points();
      const unsigned int n_q_points    = q_points.size();
      const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;
      
      const std::vector<double> &        JxW = fe_fv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =fe_fv.get_normal_vectors();
             const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1,dim>> old_solution(n_q_points);
      fe_fv[displacements].get_function_values(current_solution,old_solution);


      for (const auto &face : cell->face_iterators())
      if ( (face->boundary_id() == 1) || (face->boundary_id() == 2) )
      for (unsigned int point = 0; point < n_q_points; ++point)
        { 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=                                                                                      // - nu            
                    +1e+10 *                       
                      fe_fv[displacements].value(i, point) *   // v_h
                      fe_fv[displacements].value(j, point)     // u_h
                   *
                 JxW[point]; // dx

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {   
                copy_data.cell_rhs(i) += 

               //old_solution[point]                                              
                -  1e+10 *                 
                   fe_fv[displacements].value(i, point) 
                   * (-(face->boundary_id()-1)*disp)  *
                   JxW[point]; // dx

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


    ScratchData<dim> scratch_data(mapping, fe, quadrature,  face_quadrature);
      

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
  void PF_DG<dim>::solve()
  {

  TimerOutput::Scope t(computing_timer, "solve");

   std::cout<< "solving" << std::endl ;

   SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update,system_rhs);
 
   // constraints.distribute(newton_update);
 
    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);
  }
 
 // PROBLEM with boundary values 
 // PROBELM why does printing in the loop works when printing outside the loop exists? 
   template <int dim>
  void PF_DG<dim>::set_boundary_values()
  { 
   
    // std::map<types::global_dof_index, double> boundary_values;
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                          1,                               //boundary ID                                              
    //                                          Functions::ZeroFunction<dim>(dim),
    //                                          boundary_values);
    
    // for (auto &boundary_value : boundary_values)
    // { 
    //   current_solution(boundary_value.first) = boundary_value.second; 
    // }



    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                          2,                               //boundary ID 
    //                                          BoundaryValues<dim>(disp[0]),
    //                                          boundary_values);


    // for (auto &boundary_value : boundary_values)
    // { 
    //   current_solution(boundary_value.first) = boundary_value.second; 
    // }
   
    //constraints.distribute(current_solution);


  }
 
  template <int dim>
  double PF_DG<dim>::compute_residual(const double alpha) const
  {
    Vector<double> residual(dof_handler.n_dofs());
 
    Vector<double> evaluation_point(dof_handler.n_dofs());
    evaluation_point = current_solution;
    evaluation_point.add(alpha, newton_update);
 
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> gradients(n_q_points);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_residual = 0;
        fe_values.reinit(cell);
 
        fe_values.get_function_gradients(evaluation_point, gradients);
 
 
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1. / std::sqrt(1 + gradients[q] * gradients[q]);
 
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
                                   * coeff                    // * a_n
                                   * gradients[q]             // * \nabla u_n
                                   * fe_values.JxW(q));       // * dx
          }
 
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          residual(local_dof_indices[i]) += cell_residual(i);
      }
 
    constraints.condense(residual);
 
    for (types::global_dof_index i :
         DoFTools::extract_boundary_dofs(dof_handler))
      residual(i) = 0;
 
    return residual.l2_norm();
  }
 
 
 
  template <int dim>
  double PF_DG<dim>::determine_step_length() const
  {
    return 0.05;
  }
 
 
 
 
  template <int dim>
  void PF_DG<dim>::output_results(
    const unsigned int cycle) const
  {
     std::vector<std::string> solution_names(dim, "u");
 
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler,
                           current_solution,
                           solution_names,
                           interpretation);

  data_out.build_patches();
  std::ofstream output("solution" + std::to_string(cycle) + ".vtu");
  // std::ofstream output("solution" + std::to_string(timestep_number) + ".vtu");
  data_out.write_vtu(output);
  }
 
 
 
  template <int dim>
  void PF_DG<dim>::run()
  {
    import_mesh();
    triangulation.refine_global();

 
    setup_system(true);
    disp[0]+=0.0002;

    //  BoundaryValues<dim> BC(disp[0]);

    // VectorTools::project(dof_handler,
    //                      constraints,
    //                      QGauss<dim>(fe.degree + 1),
    //                      BC,
    //                      current_solution);


 
    double       last_residual_norm = std::numeric_limits<double>::max();
    
    unsigned int cycle   = 0;
    do
      {
        std::cout << "cycle " << cycle << std::endl;
        
        if(cycle==0)
        disp[0]=0.0002;
        else
        disp[0]=0.0;
 
            assemble_system();
            last_residual_norm = system_rhs.l2_norm();
 
            solve();
 
          //  std::cout << "  Residual: " << compute_residual(0) << std::endl;
            std::cout << "  system_rhs l2_norm " << last_residual_norm << std::endl;
            std::cout << "  update l2_norm " << newton_update.l2_norm() << std::endl;
 
   // if(cycle%1 == 0)
    {
        TimerOutput::Scope t(computing_timer, "output");
       output_results(cycle);
    //    output_stress(cycle);
    }
 
        cycle++;
        std::cout << std::endl;
      }
    while (last_residual_norm > 1e-1 && cycle<3);

         computing_timer.print_summary();
       computing_timer.reset();
  }



}//closing napespace 


#endif