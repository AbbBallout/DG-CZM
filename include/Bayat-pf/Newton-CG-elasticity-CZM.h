#ifndef NEWTON_CG_ELASTICITY_CZM_CZM_H
#define NEWTON_CG_ELASTICITY_CZM_CZM_H

/**
 * @brief Based on the paper: Locking-free interface failure modeling by a cohesivediscontinuous Galerkin method for matching andnonmatching meshes 2019
 *  Bayat
 */

// @sect3{Include files}

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
#include <deal.II/numerics/data_postprocessor.h>
#include <iostream>
#include <fstream>

namespace CG_Elasticity_CZM
{ 
using namespace dealii;


template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues(double displacement)
    : Function<dim>(dim),
     displacement(displacement)
  {}
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &  value) const override;

  private:
  double displacement ; 
};


template <int dim>
void BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                       Vector<double> &  values) const
{
       Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));


        values(0) = 0 ;
        values(1) = displacement; 
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
  // return 1236.0 ; 
  // return 0;
  //  else
  //  return 14285.714285;
  return 0 ;
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




/**
 * @brief DG scratch data
 * 
 * @tparam dim 
 */
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
class StressPostprocessor : public DataPostprocessorTensor<dim>
{
public:
  StressPostprocessor ()
    :
    DataPostprocessorTensor<dim> ("stress",
                                  update_gradients)
  {}
 
  virtual
  void evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data, std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (input_data.solution_gradients.size(),
                     computed_quantities.size());
    Tensor<2,dim> strain;
    Tensor<2,dim> stress;
    Tensor<2,dim> Identity;
    Identity[0][0]=1.0; Identity[1][0]=0.0; Identity[0][1]=0.0; Identity[1][1]=1.0; 
    double lambda, mu ; 
    const typename DoFHandler<dim>::cell_iterator current_cell = input_data.template get_cell<dim>();
     
    if(current_cell->material_id() != 1)  
    {lambda =  1236.0 ; 
    mu   =    1236.0 ;  }
   else 
    {lambda =  14285.714285;
    mu    = 3571.428571428572;}


    for (unsigned int p=0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension (computed_quantities[p].size(),
                         (Tensor<2,dim>::n_independent_components));

        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
          { 
            strain[d][e] = (input_data.solution_gradients[p][d][e] + input_data.solution_gradients[p][e][d]) / 2; 
          }
        
         stress = lambda*trace(strain)*Identity +  2* mu*strain ; 
        
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
              =stress[d][e];    
      }
  }
};

template <int dim>
class StrainPostprocessor : public DataPostprocessorTensor<dim>
{
public:
  StrainPostprocessor ()
    :
    DataPostprocessorTensor<dim> ("strain",
                                  update_gradients)
  {}
 
  virtual
  void evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data, std::vector<Vector<double> > &computed_quantities) const override
  {
    AssertDimension (input_data.solution_gradients.size(),
                     computed_quantities.size());

    for (unsigned int p=0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension (computed_quantities[p].size(),
                         (Tensor<2,dim>::n_independent_components));

        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
              =(input_data.solution_gradients[p][d][e]
                 +
                 input_data.solution_gradients[p][e][d]) / 2;
      }
  }
};

/**
 * @brief  Main class
 * 
 * @tparam dim 
*/

 template <int dim>
  class CG_Elasticity_CZM
  {
  public:
    CG_Elasticity_CZM(const unsigned int degree);
    void run();
 
  private:
    void   setup_system(const bool initial_step);
    void   assemble_system();
    void   solve();
   void    set_boundary_values(); 
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
   void output_results(const unsigned int cycle,const unsigned int refinement_level) const;


    void reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id);
 
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
    Vector<double> newton_update;
    Vector<double> system_rhs;

    const unsigned int degree;

    ConvergenceTable convergence_table;
    TimerOutput computing_timer;
    unsigned int timestep_number;

    Tensor<1,dim> disp;
  };



  template <int dim>
  CG_Elasticity_CZM<dim>::CG_Elasticity_CZM(const unsigned int degree)
    : degree(degree), 
    quadrature(degree + 1), 
    face_quadrature(degree + 1),
    mapping(), 
    fe(FE_Q<dim>(degree), dim ), 
    dof_handler(triangulation), 
    computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)

{
      disp[0]=0.0;
      disp[1]=0.0;
}
 
  template <int dim>
  void CG_Elasticity_CZM<dim>::setup_system(const bool initial_step)
  { 
    TimerOutput::Scope t(computing_timer, "set up");

    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(dof_handler.n_dofs());
      }
 
 
    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
   
   

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
 
    sparsity_pattern.copy_from(dsp);
   
    system_matrix.reinit(sparsity_pattern);

    std::cout << "dofs = " << dof_handler.n_dofs() << std::endl ; 
  }


  template <int dim>
  void CG_Elasticity_CZM<dim>::assemble_system()
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

          Tensor<dim,dim> Identity;
       Identity[0][0]=1; Identity[1][1]=1; 

        for (unsigned int point = 0; point < n_q_points; ++point)            
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            { 
                Tensor<dim,dim> straini = 0.5*(fe_v[displacements].gradient(i, point)+transpose(fe_v[displacements].gradient(i, point))) ;

              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
             { 
                Tensor<dim,dim> strainj = 0.5*(fe_v[displacements].gradient(j, point)+transpose(fe_v[displacements].gradient(j, point))) ;
                Tensor<dim,dim> stress = lambda_values[point]*trace(strainj)*Identity + 2*mu_values[point]*strainj ;

                copy_data.cell_matrix(i, j) += scalar_product(straini,stress)*JxW[point];
                                                       
              } 
            
            Tensor<dim,dim> old_strain = 0.5*(old_solution_gradients[point]+transpose(old_solution_gradients[point])) ;
            Tensor<dim,dim> old_stress = lambda_values[point]*trace(old_strain)*Identity + 2*mu_values[point]*old_strain ;

             
             copy_data.cell_rhs(i) +=  - scalar_product(straini,old_stress)*JxW[point];                  
            
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
      
      
    };

        const auto boundary_worker = [&](const auto &        cell,
                                     const unsigned int &face_no,
                                     auto &              scratch_data,
                                     auto &              copy_data)
                                     
      {  


                                       
        
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

  constraints.condense(system_matrix);
  constraints.condense(system_rhs);
 
 {
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           2,
                                           Functions::ZeroFunction<dim>(dim),
                                           boundary_values);

  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     newton_update,
                                     system_rhs);  
 }         

  {
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           3,
                                           Functions::ZeroFunction<dim>(dim),
                                           boundary_values);

  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     newton_update,
                                     system_rhs);  
 }     


  }

  template <int dim>
  void CG_Elasticity_CZM<dim>::solve()
  {

  TimerOutput::Scope t(computing_timer, "solve");

   std::cout<< "solving" << std::endl ;

   SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update,system_rhs);
 
    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);

  }


  template <int dim>
  double CG_Elasticity_CZM<dim>::determine_step_length() const
  {
    return 1;
  }
 
 
 
 
 template <int dim>
  void CG_Elasticity_CZM<dim>::output_results(const unsigned int cycle,const unsigned int refinement_level) const
  {

     const StrainPostprocessor<dim> strain;
     const StressPostprocessor<dim> stress;
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
     std::vector<std::string> solution_names(dim, "u");
  
  DataOut<dim> data_out;
  
  data_out.add_data_vector(dof_handler,current_solution,solution_names,interpretation);
   data_out.add_data_vector (dof_handler,current_solution, strain);
   data_out.add_data_vector (dof_handler,current_solution, stress);
  data_out.build_patches();


  std::ofstream output("solution"+ std::to_string(refinement_level) +"_" + std::to_string(cycle)  + ".vtu");
  // std::ofstream output("solution" + std::to_string(timestep_number) + ".vtu");
  data_out.write_vtu(output);
  
  }


    template <int dim>
  void CG_Elasticity_CZM<dim>::reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id)
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
                    fe_face_values[displacements].get_function_gradients(current_solution,gradu) ;
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
                                               
}

  template <int dim>
  void CG_Elasticity_CZM<dim>::run()
  {
    
    for(unsigned int refinement_level = 0 ; refinement_level<3 ; refinement_level++  )
   { 
    triangulation.clear();
    //Mesh geenration 
    const Point<dim> P1, P2(1,2);
    const std::vector<unsigned int> repetitions{2, 4};  
    GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,P1,P2,true);
    triangulation.refine_global(3+refinement_level);
 
    Point<dim> p ;
    for (const auto &cell : triangulation.active_cell_iterators())
    {p=cell->center() ; 
    // if(std::abs(p[1]-1)<0.05 )
    cell->set_material_id(1) ;
    // else
    // cell->set_material_id(2) ;
    }
    
 


    double error ;
    unsigned nonlin = 0 ;
    disp[0]=0;
    disp[1]=0;
     std::ofstream forces("forces" + std::to_string(refinement_level) +   ".txt");

   for (unsigned int cycle = 0; cycle <20; ++cycle)
      {
     error=1;
      nonlin = 0 ;
     disp[1]+=1e-4;
       std::cout<< "cycle = " << cycle << " and dislacement = " << disp[1] << "\n" ; 
        
       std::cout<< "   nonlin  = " << nonlin << "\n" ; 
       while(error>1e-4 &&  nonlin<200)
       {
        if (cycle==0 & nonlin==0)
    // if (cycle==0)
     setup_system(true);
     else
     setup_system(false);

  {std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           2,
                                           BoundaryValues<dim>(0.0),
                                           boundary_values);

  for (auto &boundary_value : boundary_values)
    current_solution(boundary_value.first) = boundary_value.second;
 
 // constraints.distribute(current_solution);
  }

  {std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           3,
                                           BoundaryValues<dim>(disp[1]),
                                           boundary_values);
  for (auto &boundary_value : boundary_values)
    current_solution(boundary_value.first) = boundary_value.second;
 
  //constraints.distribute(current_solution);
  }

     
     nonlin++ ; 

     assemble_system();
     solve();
     
     error = newton_update.l2_norm();
          //  std::cout << "  Residual: " << compute_residual(0) << std::endl;
            std::cout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;
            std::cout << "  update l2_norm " << newton_update.l2_norm() << std::endl;
      }
    if(cycle%1 == 0)
    {
        TimerOutput::Scope t(computing_timer, "output");
        output_results(cycle,refinement_level);
    }
     forces.open("forces" + std::to_string(refinement_level) +   ".txt", std::ios_base::app) ;
      Tensor<1,dim> ten; 
      
     reaction(ten,3);
      std::cout<<  "   reaction force = " << ten[0] << "\t" << ten[1] << "\n" ; 
      
        std::cout << std::endl;
         forces << disp[1] <<  "\t"  << ten[0]  << "\t" << ten[1] << "\n" ; 
     forces.close();
      }

         computing_timer.print_summary();
       computing_timer.reset();
  }

  }





} //closing namespace 











#endif 
