#ifndef PAHSEFILEDDG4_H
#define PAHSEFILEDDG4_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/base/symmetric_tensor.h>
 
#include <deal.II/physics/transformations.h>
#include <fstream>
#include <iostream>
 
namespace pfdgm
{
  using namespace dealii;

template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;

  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;

    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << '(' << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
    // std::ofstream out(filename);
  // GridOut       grid_out;
  // grid_out.write_vtu(triangulation, out);
  // std::cout << " written to " << filename << std::endl << std::endl;
}

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
  if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <0.0003 ) 
      || 
      (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <0.0003 )   )
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
  if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <0.0003 ) 
      || 
      (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <0.0003 )   )
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
  struct PointHistory
  {
    SymmetricTensor<2, dim> old_stress;
  };


  template <int dim>
  class PF_DGM
  {
  public:
    PF_DGM(unsigned int degree);
 
    void run();
 
  private:
    void import_mesh();
    void setup_system(const bool initial_step);
    void assemble_system();
    void solve_time_step();
    void output_results(const unsigned int cycle) const;
    void output_stress(const unsigned int cycle) const;

   void move_mesh(); 
   void reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id);
    unsigned int degree;
    MPI_Comm     mpi_communicator;
 
    FESystem<dim>                             fe;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;
    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;
    const MappingQ1<dim> mapping;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
 
    AffineConstraints<double> constraints;
 
    PETScWrappers::MPI::SparseMatrix system_matrix;
    PETScWrappers::MPI::Vector      solution_np1;
    PETScWrappers::MPI::Vector      solution_np0;
    PETScWrappers::MPI::Vector        system_rhs;
 
     unsigned int timestep_number;
    Tensor<1,dim> disp;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };
 
 
 
  template <int dim>
  PF_DGM<dim>::PF_DGM(unsigned int degree)
    : degree(degree)
    ,quadrature(degree + 1) 
    ,face_quadrature(degree + 1)
    ,mapping()
    , mpi_communicator(MPI_COMM_WORLD)
    , fe(FE_DGQ<dim>(degree), dim)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {
          disp[0]=0.0;
      disp[1]=0.0;
  }
 
 
  template <int dim>
  void PF_DGM<dim>::import_mesh()
  {
  // const Point<2> bottom_hole_origin(0.4, 0);
  // const Point<2> top_hole_origin(0.6, 0.5);


  // const SphericalManifold<2> bottom_manifold(bottom_hole_origin);
  // const SphericalManifold<2> top_manifold(bottom_hole_origin);


  GridIn<2> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("fenics.msh");
  gridin.read_msh(f);

  //print_mesh_info(triangulation, "grid-1.vtu");
  }
 
  template <int dim>
  void PF_DGM<dim>::setup_system(const bool initial_step)
  {
    TimerOutput::Scope t(computing_timer, "setup");
     
  dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
 
    solution_np1.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
   
     if (initial_step)
    solution_np0.reinit(locally_owned_dofs, mpi_communicator);


   constraints.clear();
   //constraints.reinit(locally_relevant_dofs);
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                          0,
    //                                          Functions::ZeroFunction<dim>(dim),
    //                                          constraints);
    // std::vector<double> disp(2);
    // disp[0]=0.01;
    // disp[1]=0; 
    //     VectorTools::interpolate_boundary_values(dof_handler,
    //                                          1,
    //                                          Functions::ConstantFunction<dim>(disp),
    //                                          constraints);

    constraints.close();
   
    DynamicSparsityPattern dsp(locally_relevant_dofs);
     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);
 
    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    pcout << "dofs = " << dof_handler.n_dofs() << std::endl ; 
  }
 
 
 
  template <int dim>
  void PF_DGM<dim>::assemble_system()
  {
    TimerOutput::Scope t(computing_timer, "assembly");
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

      if (cell->is_locally_owned())
        for (unsigned int point = 0; point < n_q_points; ++point)            
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
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
         fe_iv[displacements].get_jump_in_function_values(solution_np0,old_solution);


        std::vector<double> GC_values(n_q_points);
        std::vector<double> sig_values(n_q_points);       

        GC<dim> gc;
        SIG<dim> sig;  

          Tensor<1,dim> tangential;
       

         gc.value_list(fe_iv.get_quadrature_points(), GC_values);
         sig.value_list(fe_iv.get_quadrature_points(), sig_values);
         double law, delta ;
         
       if (cell->is_locally_owned())
         for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < n_dofs_face; ++i)  
             for (unsigned int j = 0; j < n_dofs_face; ++j) 
             { 
              delta=  GC_values[point]/sig_values[point]/std::exp(1);
              tangential = cross_product_2d(normals[point]);
              law= GC_values[point]/(delta*delta) * std::exp(- std::sqrt
              (old_solution[point]*normals[point]*old_solution[point]*normals[point]+ 2*2*old_solution[point]*tangential*old_solution[point]*tangential)
              /delta);
              copy_data_face.cell_matrix(i, j) += law*                              
                   fe_iv[displacements].jump_in_values(i, point)  *            
                   fe_iv[displacements].jump_in_values(j, point)  *
                JxW[point]; 

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
       
      if (cell->is_locally_owned())
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
                                                                             
                +  1e+10 *                 
                   fe_fv[displacements].value(i, point) 
                   * (face->boundary_id()-1)*disp *
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

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }
 
 
 
  template <int dim>
  void PF_DGM<dim>::solve_time_step()
  {
    TimerOutput::Scope              t(computing_timer, "solve");
   PETScWrappers::MPI::Vector completely_distributed_solution(
      locally_owned_dofs, mpi_communicator);
    SolverControl                    solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.solve(system_matrix,  completely_distributed_solution , system_rhs);
    solution_np1 = completely_distributed_solution; 

   
  }
 
 

 
 
 
  template <int dim>
  void PF_DGM<dim>::output_results(const unsigned int cycle) const
  {
    const Vector<double> localized_solution(solution_np1);
      DataOut<dim>  data_out;
 
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);
 
      data_out.attach_dof_handler(dof_handler);
      

     std::vector<std::string> solution_names(dim, "u");
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
 
           data_out.add_data_vector(localized_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

 
      Vector<double> mpi_owner(triangulation.n_active_cells());
      mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      data_out.add_data_vector(mpi_owner, "owner");
 
      data_out.build_patches(mapping,
                             fe.degree);
 
      const std::string filename =
        "solution_" + Utilities::int_to_string(cycle, 3) + ".vtu";
      data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
    
  }

    template <int dim>
  void PF_DGM<dim>::output_stress(const unsigned int cycle) const
  {   

     PETScWrappers::MPI::Vector strains;
     strains.reinit(locally_owned_dofs, mpi_communicator);
       FEValues<dim> fe_values(fe,
                                quadrature,
                                update_values | update_gradients |
                                update_quadrature_points | update_JxW_values); 
   

   const unsigned int q_points = quadrature.size();

   std::vector<Tensor<2,dim>> gradu(q_points);
   std::vector< std::vector< Tensor< 1, dim > > > previous_gradient (q_points,std::vector<Tensor<1,dim> >(dim));

  std::vector<double> lambda_values(q_points);
  std::vector<double> mu_values(q_points);     

  typename DoFHandler<dim>::active_cell_iterator 
  cell = dof_handler.begin_active(), endc = dof_handler.end();

  Lambda<dim> lambda(cell->material_id()); 
  Mu<dim> mu(cell->material_id()); 

  const FEValuesExtractors::Vector displacements(0);
    Tensor<2,dim> Identity;
    Identity[0][0]=1; Identity[1][0]=0; Identity[0][1]=0; Identity[1][1]=1; 

   
  }
  
 
  template <int dim>
 void PF_DGM<dim>::move_mesh()
{
  
 
  std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
  for (auto &cell : dof_handler.active_cell_iterators())
    for (const auto v : cell->vertex_indices())
      if (vertex_touched[cell->vertex_index(v)] == false)
        {
          vertex_touched[cell->vertex_index(v)] = true;
 
          Point<dim> vertex_displacement;
          for (unsigned int d = 0; d < dim; ++d)
            vertex_displacement[d] =
              solution_np1(cell->vertex_dof_index(v, d));
 
          cell->vertex(v) += vertex_displacement;
        }
}

 template <int dim>
  void PF_DGM<dim>::reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id)
 {  
          reaction_stress = 0; 
        QGauss<dim-1> face_quadrature_formula(5);

        FEFaceValues<dim> fe_face_values (fe, 
                                          face_quadrature_formula,
                                          UpdateFlags(update_values             | 				
                                                      update_gradients          | 
                                                      update_quadrature_points  |
                                                      update_normal_vectors     |
                                                      update_JxW_values));

       const unsigned int n_face_q_points = face_quadrature_formula.size();

        std::vector<Tensor<2,dim>> gradu(n_face_q_points);
        std::vector< std::vector< Tensor< 1, dim > > > previous_gradient (n_face_q_points,std::vector<Tensor<1,dim> >(dim));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        
        std::vector<double> lambda_values(n_face_q_points);
        std::vector<double> mu_values(n_face_q_points);       

        Lambda<dim> lambda(cell->material_id()); 
        Mu<dim> mu(cell->material_id()); 


     
       const FEValuesExtractors::Vector displacements(0);
       Tensor<2,dim> Identity;
       Identity[0][0]=1; Identity[1][0]=0; Identity[0][1]=0; Identity[1][1]=1; 


            for (; cell!=endc; ++cell)      
                for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    	        if (cell->face(face)->boundary_id() == boundary_id)
	            {
                    fe_face_values.reinit (cell, face);      
                    fe_face_values[displacements].get_function_gradients(solution_np1,gradu) ;
                    lambda.value_list(fe_face_values.get_quadrature_points(), lambda_values);
                    mu.value_list(fe_face_values.get_quadrature_points(), mu_values); 
                    const std::vector<double>         &JxW = fe_face_values.get_JxW_values();
                    const std::vector<Tensor<1, dim>> &normals =fe_face_values.get_normal_vectors();                         
       
	                 for (unsigned int point=0; point<n_face_q_points; ++point)
                    {    auto strain = 0.5*(gradu[point]+transpose(gradu[point]));                                               
                         
                         reaction_stress += (lambda_values[point]*trace(strain)*Identity  + 2*mu_values[point]*strain )

                                           *normals[point]*JxW[point] ;
                     }                        
	            }
            
      //    ;      
 }
 
  template <int dim>
  void PF_DGM<dim>::run()
  {


    timestep_number = 0;
    AffineConstraints<double> constraints;
    constraints.close();

    import_mesh();
    triangulation.refine_global(2);
      double error ;
     unsigned nonlin = 0 ;
     std::ofstream forces("forces.txt");


    for (unsigned int cycle = 0; cycle <450; ++cycle)
    { error=1;
      nonlin = 0 ;
     disp[0]+=0.0002;
       pcout<< "cycle = " << cycle << " and dislacement = " << disp[0] << "\n" ; 
      while(error>1e-4 &&  nonlin<200)
    { 
       pcout<< "   nonlin  = " << nonlin << "\n" ; 
        if (cycle==0 & nonlin==0)
    // if (cycle==0)
     setup_system(true);
     else
     setup_system(false);
     
     nonlin++ ; 

     assemble_system();
     solve_time_step();
     
     //stress

     //error
      PETScWrappers::MPI::Vector  diff = solution_np0 ; 
      diff.add(-1,  solution_np1);
      error =  diff.l2_norm();
       pcout << "    error = " << error << std::endl ; 
 
      //solution update
      solution_np0 =  solution_np1;


    }
    //output
    
    //std::cout << solution_np1[0] << " " <<  solution_np1[50000] << " " << solution_np1[25000]  <<"\n" ; 
  //if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
    if(cycle%2 == 0)
    {
        TimerOutput::Scope t(computing_timer, "output");
        output_results(cycle);
    }
    


     forces.open("forces.txt", std::ios_base::app) ;
      Tensor<1,dim> ten;
     reaction(ten,1);
      pcout<< "   reaction force = " << ten[0] << "\t" << ten[1] << "\n" ; 

    forces << disp[0] <<  "\t"  << ten[0]  << "\t" << ten[1] << "\n" ;
     forces.close();
  }
    

   computing_timer.print_summary();
       computing_timer.reset();
    
  }


 

  
} // namespace 
 
 















#endif