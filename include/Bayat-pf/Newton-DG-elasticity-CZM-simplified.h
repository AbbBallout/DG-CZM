#ifndef NEWTON_DG_ELASTICITY_CZM_CZM_SIMPLIFIED_H
#define NEWTON_DG_ELASTICITY_CZM_CZM_SIMPLIFIED_H


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
#include <deal.II/base/parameter_acceptor.h>

#include <iostream>
#include <fstream>
#include <deal.II/base/hdf5.h>
namespace DG_Elasticity_CZM_simplified
{ 
using namespace dealii;


  template <int dim>
class Lambda : public Function<dim>
{
public:
  Lambda()
  {}
   
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double Lambda<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
 // return 16644.64353758137;
 // return  57.6923155396654; 
  return 0;
}

template <int dim>
class Mu : public Function<dim>
{
public:
  Mu()
  {}

  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;

};

template <int dim>
double Mu<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{ 
  
  return 38.46153987229933;
 // return 33.35556876745394;

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
     
//lambda = 57.6923155396654;
    //lambda=16644.64353758137; 
    lambda=0;  
    mu   =  38.46153987229933;  
   // mu   = 33.35556876745394;
     

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

template <int dim>
struct PointHistory
{
  bool is_damaged; 
};


/**
 * @brief  Main class
 * 
 * @tparam dim 
*/

 template <int dim>
  class DG_Elasticity_CZM
  {
  public:
    DG_Elasticity_CZM(const unsigned int degree);
    void run();
 
  private:
    void   setup_system(const bool initial_step);
    void   assemble_system(const int non_lin,const unsigned int cycle);
    void   solve(const unsigned int cycle);
    //double compute_residual(const double alpha) const;
    double determine_step_length() ;
   void output_results(const unsigned int cycle,const unsigned int refinement_level) const;


    void reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id);
    void interface_traction (Tensor<1,dim> &interface_stress);

 
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


    std::vector<PointHistory<dim>> quadrature_point_history;
    void setup_quadrature_point_history();

  };

// TODO: maybe switch rows and columns quadrature_point_history[face][history] 
  template <int dim>
  void DG_Elasticity_CZM<dim>::setup_quadrature_point_history()
  {
 
    triangulation.clear_user_data();
 
    {
      std::vector<PointHistory<dim>> tmp;
      quadrature_point_history.swap(tmp);
    }

  quadrature_point_history.resize(triangulation.n_active_faces()*face_quadrature.size()) ;


     unsigned int history_index = 0;
     for (auto &face : triangulation.active_face_iterators())
     {  
        face->set_user_pointer(&quadrature_point_history[history_index]);
         history_index += face_quadrature.size();
     }
 
    Assert(history_index == quadrature_point_history.size(),      ExcInternalError());
  }

  template <int dim>
  DG_Elasticity_CZM<dim>::DG_Elasticity_CZM(const unsigned int degree)
    :degree(degree), 
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
  void DG_Elasticity_CZM<dim>::setup_system(const bool initial_step)
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
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
 
    sparsity_pattern.copy_from(dsp);
   
    system_matrix.reinit(sparsity_pattern);

    std::cout << "dofs = " << dof_handler.n_dofs() << std::endl ; 
  }


  template <int dim>
  void DG_Elasticity_CZM<dim>::assemble_system(const int non_lin,const unsigned int cycle)
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

        Lambda<dim> lambda; 
        Mu<dim> mu; 

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
      
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
       fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
        copy_data.face_data.emplace_back();
          CopyDataFace &     copy_data_face = copy_data.face_data.back();

      const auto &       q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      

      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      copy_data_face.cell_rhs.reinit(n_dofs_face) ; 
      const std::vector<double> &        JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

        std::vector<double> lambda_values(n_q_points);
        std::vector<double> mu_values(n_q_points);       

        Lambda<dim> lambda; 
        Mu<dim> mu; 

        lambda.value_list(q_points, lambda_values);
        mu.value_list(q_points, mu_values);


         const FEValuesExtractors::Vector displacements(0);

         std::vector<Tensor<1,dim>> old_solution_jumps(n_q_points);
         fe_iv[displacements].get_jump_in_function_values(current_solution,old_solution_jumps);
          

        std::vector<Tensor<2,dim>> gradu(n_q_points);
        fe_iv[displacements].get_average_of_function_gradients(current_solution,gradu);

          Tensor<dim,dim> Identity;
       Identity[0][0]=1; Identity[1][1]=1; 

          
          const double penalty = (1e+3)* ((2*mu_values[0])*(3*lambda_values[0]+2*mu_values[0])/(lambda_values[0]+mu_values[0]))
          /(cell->measure()); 
          PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());

          const double beta=1.0; // do you multiply this with the other thingies? 
          const double t_0 =3.0;          // maximim insterface strngth 
          const double lambda_f = 0.001 ;  //elongation at full faliure 
          const double m=1, n=2.0 ;   
          double seperation; 
          Tensor<1,dim> g,g_eff ;  
          Tensor<1,dim> tr, tr_eff;  //traction   

          Tensor<2,dim> TCZ; 
          Tensor<1,dim> TCZ_res; 
           Tensor<1,dim> tangential; 


         for (unsigned int point = 0; point < n_q_points; ++point)
         { 
            Tensor<2,dim> old_strain = 0.5*(gradu[point]+transpose(gradu[point]));
            Tensor<2,dim> old_stress = (lambda_values[point]*trace(old_strain)*Identity +2*mu_values[point]*old_strain);

             tangential= cross_product_2d(normals[point]);
             tr[0]= beta*scalar_product(old_stress,outer_product(normals[point],tangential));
             tr[1]= scalar_product(old_stress,outer_product(normals[point],normals[point]));
             tr_eff[0]= tr[0];
             tr_eff[1]= std::abs(tr[1])/2 + tr[1]/2 ;
            

             if(tr_eff.norm()>t_0)
               quadrature_points_history[point].is_damaged=true ;

            
            g=old_solution_jumps[point];
            g=-g;
            g_eff[0]=beta*g[0]; 
            g_eff[1]=std::abs(g[1])/2 + g[1]/2;
            seperation = g_eff.norm();


          for (unsigned int i = 0; i < n_dofs_face; ++i)  
          {   
          
            Tensor<2,dim> straini = 0.5*(fe_iv[displacements].average_of_gradients(i, point)+transpose(fe_iv[displacements].average_of_gradients(i, point))) ;


             for (unsigned int j = 0; j < n_dofs_face; ++j) 
             {  
               Tensor<2,dim> strainj = 0.5*(fe_iv[displacements].average_of_gradients(j, point)+transpose(fe_iv[displacements].average_of_gradients(j, point))) ;

                if(quadrature_points_history[point].is_damaged==false  || cell->face(f)->user_flag_set () ==false) 
                {
                copy_data_face.cell_matrix(i, j) += 
                (
                -fe_iv[displacements].jump_in_values(i, point)*
                (lambda_values[point]*trace(strainj)*Identity +2*mu_values[point]*strainj )*normals[point]

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini)*normals[point]*
                fe_iv[displacements].jump_in_values(j, point)  //Symetry term

                + penalty*fe_iv[displacements].jump_in_values(i, point)*fe_iv[displacements].jump_in_values(j, point)
                )
                *JxW[point] ;
                }
                else
                {

                //////traction law //////////////
                TCZ[0][0]= std::pow(beta,3)*std::pow(g[0],2); 
                TCZ[0][1]= beta*g[0]*(std::abs(g[1])/2 - g[1]/2); 
                TCZ[1][0]= beta*beta*g[0]*(std::abs(g[1])/2 - g[1]/2); 
                TCZ[1][1]= std::pow((std::abs(g[1])/2 - g[1]/2),2); 

                TCZ= - TCZ* ((lambda_f+seperation*(m-1))/(seperation*seperation*(lambda_f-seperation)));

               TCZ[0][0]+= beta;
               TCZ[1][1]+= 0.5*(1+((g[1] > 0) ? 1 : ((g[1] < 0) ? -1 : 0)));
  
               TCZ=TCZ*(t_0/seperation) *std::pow((lambda_f-seperation) /lambda_f ,m) ; 
 
              TCZ[1][1]+= penalty*std::pow(2,-n)*n*std::pow(std::abs(g[1])-g[1],n)/std::abs(g[1]);
              
               
                   copy_data_face.cell_matrix(i, j) +=
                   ( fe_iv[displacements].jump_in_values(i, point))*(TCZ*fe_iv[displacements].jump_in_values(j, point))*JxW[point] ; 

              }
            
            }

            if(quadrature_points_history[point].is_damaged==false ||  cell->face(f)->user_flag_set () ==false) 
            {
            copy_data_face.cell_rhs(i)+=               
              -(
                -fe_iv[displacements].jump_in_values(i, point)*
                old_stress *normals[point]

                -(lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini)*normals[point]*
                 old_solution_jumps[point]  //Symetry term

                + penalty*fe_iv[displacements].jump_in_values(i, point)*old_solution_jumps[point]
                )
                *JxW[point] ; 

            }
            else{

                
                TCZ_res = g_eff*(t_0/seperation)*std::pow( (lambda_f-seperation) / lambda_f ,m);
                TCZ_res[1]+= -penalty* std::pow(std::abs(g[1])/2-g[1]/2,n);
                  TCZ_res=-TCZ_res;

              copy_data_face.cell_rhs(i)+= 
              -fe_iv[displacements].jump_in_values(i, point)*TCZ_res*JxW[point];
           

          }


        }
      } 

    };

          const auto boundary_worker = [&](const auto &        cell,
                                     const unsigned int &face_no,
                                     auto &              scratch_data,
                                     auto &              copy_data)
                                     
       { 
        scratch_data.fe_interface_values.reinit(cell, face_no);
      const FEFaceValuesBase<dim> &fe_fv =
        scratch_data.fe_interface_values.get_fe_face_values(0);

      const auto &       q_points      = fe_fv.get_quadrature_points();
      const unsigned int n_q_points    = q_points.size();
      const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;
      
      const std::vector<double> &        JxW = fe_fv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =fe_fv.get_normal_vectors();
             const FEValuesExtractors::Vector displacements(0);
             const FEValuesExtractors::Scalar displacements_x(0);
             const FEValuesExtractors::Scalar displacements_y(1);


      std::vector<Tensor<1,dim>> old_solution_values(n_q_points);
      fe_fv[displacements].get_function_values(current_solution,old_solution_values);

        std::vector<Tensor<2,dim>> gradu(n_q_points);
        std::vector< std::vector< Tensor< 1, dim > > > previous_gradient (n_q_points,std::vector<Tensor<1,dim> >(dim));

      std::vector<double> lambda_values(n_q_points);
        std::vector<double> mu_values(n_q_points);       

        Lambda<dim> lambda; 
        Mu<dim> mu; 

        lambda.value_list(q_points, lambda_values);
        mu.value_list(q_points, mu_values);

          const double penalty = (1e+3)* ((2*mu_values[0])*(3*lambda_values[0]+2*mu_values[0])/(lambda_values[0]+mu_values[0]))
          /(cell->measure()); 

          Tensor<dim,dim> Identity;
       Identity[0][0]=1; Identity[1][1]=1; 

        Tensor<1,dim> tangential;
      
    if((cell->face(face_no)->boundary_id() == 2) ) 
      for (unsigned int point = 0; point < n_q_points; ++point)
        { 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {     
          
          Tensor<2,dim> straini = 0.5*(fe_fv[displacements].gradient(i, point)+transpose(fe_fv[displacements].gradient(i, point))) ;

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {                                                                                    
               
               Tensor<2,dim> strainj = 0.5*(fe_fv[displacements].gradient(j, point)+transpose(fe_fv[displacements].gradient(j, point))) ;

                copy_data.cell_matrix(i, j) += 
                (
                -fe_fv[displacements].value(i, point)*
                (lambda_values[point]*trace(strainj)*Identity +2*mu_values[point]*strainj )*normals[point]

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*
                fe_fv[displacements].value(j, point)  //Symetry term

                + penalty*fe_fv[displacements].value(i, point)*fe_fv[displacements].value(j, point)
                )
                *JxW[point] ;
           }
          


            Tensor<2,dim> strain = 0.5*(gradu[point]+transpose(gradu[point])); 

              
                copy_data.cell_rhs(i) += 

              -(
                -fe_fv[displacements].value(i, point)*
                (lambda_values[point]*trace(strain)*Identity +2*mu_values[point]*strain )*normals[point]

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*
                old_solution_values[point]                                                                   //Symetry term

                + penalty*fe_fv[displacements].value(i, point)*old_solution_values[point]
                )
                *JxW[point] 

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*
                ((cell->face(face_no)->boundary_id()-2)*disp)*JxW[point]               //Symetry term

                  +penalty *                 
                   fe_fv[displacements].value(i, point) 
                   * ((cell->face(face_no)->boundary_id()-2)*disp)  *
                   JxW[point]; // dx

          }
                                      
        }

      //    if(cell->face(face_no)->boundary_id() == 3 ) 
      // for (unsigned int point = 0; point < n_q_points; ++point)
      //   { 
      //     for (unsigned int i = 0; i < dofs_per_cell; ++i)
      //   {     
      //     tangential= cross_product_2d(normals[point]);
      //     Tensor<2,dim> straini = 0.5*(fe_fv[displacements].gradient(i, point)+transpose(fe_fv[displacements].gradient(i, point))) ;

      //       for (unsigned int j = 0; j < dofs_per_cell; ++j)
      //       {                                                                                    
               
      //          Tensor<2,dim> strainj = 0.5*(fe_fv[displacements].gradient(j, point)+transpose(fe_fv[displacements].gradient(j, point))) ;

      //           copy_data.cell_matrix(i, j) += 
      //           (
      //           -fe_fv[displacements_x].value(i, point)*
      //           (lambda_values[point]*trace(strainj)*Identity +2*mu_values[point]*strainj )tangential*tangential

      //           -
      //           (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*tangential*tangential*
      //           fe_fv[displacements_y].value(j, point)  //Symetry term

      //           + penalty*fe_fv[displacements_y].value(i, point)*fe_fv[displacements_y].value(j, point)
      //           )
      //           *JxW[point] ;
      //      }
          


      //       Tensor<2,dim> strain = 0.5*(gradu[point]+transpose(gradu[point])); 

              
      //           copy_data.cell_rhs(i) += 

      //         -(
      //           -fe_fv[displacements_y].value(i, point)*
      //           (lambda_values[point]*trace(strain)*Identity +2*mu_values[point]*strain )*normals[point]*normals[point]

      //           -
      //           (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*normals[point]*
      //           (old_solution_values[point] *normals[point] )                                                                 //Symetry term

      //           + penalty*fe_fv[displacements_y].value(i, point)*old_solution_values[point]*normals[point]
      //           )
      //           *JxW[point] 

      //           -
      //           (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*normals[point]*
      //           ((cell->face(face_no)->boundary_id()-2)*disp[1])*JxW[point]               //Symetry term

      //             +penalty *                 
      //              fe_fv[displacements_y].value(i, point) 
      //              * ((cell->face(face_no)->boundary_id()-2)*disp[1])  *
      //              JxW[point]; // dx

      //     }
                           
      //   } 



         if(cell->face(face_no)->boundary_id() == 3 ) 
      for (unsigned int point = 0; point < n_q_points; ++point)
        { 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {     
          
          Tensor<2,dim> straini = 0.5*(fe_fv[displacements].gradient(i, point)+transpose(fe_fv[displacements].gradient(i, point))) ;

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {                                                                                    
               
               Tensor<2,dim> strainj = 0.5*(fe_fv[displacements].gradient(j, point)+transpose(fe_fv[displacements].gradient(j, point))) ;

                copy_data.cell_matrix(i, j) += 
                (
                -fe_fv[displacements_y].value(i, point)*
                (lambda_values[point]*trace(strainj)*Identity +2*mu_values[point]*strainj )*normals[point]*normals[point]

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*normals[point]*
                fe_fv[displacements_y].value(j, point)  //Symetry term

                + penalty*fe_fv[displacements_y].value(i, point)*fe_fv[displacements_y].value(j, point)
                )
                *JxW[point] ;
           }
          


            Tensor<2,dim> strain = 0.5*(gradu[point]+transpose(gradu[point])); 

              
                copy_data.cell_rhs(i) += 

              -(
                -fe_fv[displacements_y].value(i, point)*
                (lambda_values[point]*trace(strain)*Identity +2*mu_values[point]*strain )*normals[point]*normals[point]

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*normals[point]*
                (old_solution_values[point] *normals[point] )                                                                 //Symetry term

                + penalty*fe_fv[displacements_y].value(i, point)*old_solution_values[point]*normals[point]
                )
                *JxW[point] 

                -
                (lambda_values[point]*trace(straini)*Identity +2*mu_values[point]*straini )*normals[point]*normals[point]*
                ((cell->face(face_no)->boundary_id()-2)*disp[1])*JxW[point]               //Symetry term

                  +penalty *                 
                   fe_fv[displacements_y].value(i, point) 
                   * ((cell->face(face_no)->boundary_id()-2)*disp[1])  *
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
                                                 cdf.cell_rhs,
                                                 cdf.joint_dof_indices,
                                                system_matrix,
                                                system_rhs);
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
  void DG_Elasticity_CZM<dim>::solve(const unsigned int cycle)
  {

  TimerOutput::Scope t(computing_timer, "solve");

   std::cout<< "solving" << std::endl ;

   SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update,system_rhs);

      double alpha;
    if(cycle<55|| cycle > 70)
    alpha = 1.0;
    else
    alpha = determine_step_length();

    current_solution.add(alpha, newton_update);

  }


  template <int dim>
  double DG_Elasticity_CZM<dim>::determine_step_length() 
  {
    return 0.4;
  }
 
 
 template <int dim>
  void DG_Elasticity_CZM<dim>::output_results(const unsigned int cycle,const unsigned int refinement_level) const
  {
     

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
     std::vector<std::string> solution_names(dim, "u");
  
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler,current_solution,solution_names,interpretation);
      const StrainPostprocessor<dim> strain;
     const StressPostprocessor<dim> stress;
   data_out.add_data_vector (dof_handler,current_solution, strain);
   data_out.add_data_vector (dof_handler,current_solution, stress);
  data_out.build_patches();

  std::ofstream output("solution"+ std::to_string(refinement_level) +"_" + std::to_string(cycle)  + ".vtu");
 data_out.write_vtu(output);
  
  }


    template <int dim>
  void DG_Elasticity_CZM<dim>::reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id)
{       
        reaction_stress = 0; 
        QGauss<dim-1> face_quadrature_formula(degree + 2);

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

        Lambda<dim> lambda; 
        Mu<dim> mu; 


     
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
void DG_Elasticity_CZM<dim>::interface_traction (Tensor<1,dim> &interface_stress)
{
          interface_stress = 0; 
        QGauss<dim-1> face_quadrature_formula(degree + 2);

        FEInterfaceValues<dim> fe_iv (fe, 
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

        Lambda<dim> lambda; 
        Mu<dim> mu; 


     
       const FEValuesExtractors::Vector displacements(0);
       Tensor<2,dim> Identity;
       Identity[0][0]=1; Identity[1][0]=0; Identity[0][1]=0; Identity[1][1]=1; 

            
            for (; cell!=endc; ++cell)      
                for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    	        if (cell->face(face)->user_flag_set())
	            {  
                //I want to visit the face only once 
                cell->face(face)->clear_user_flag();

                    fe_iv.reinit (cell, face);      
                    fe_iv[displacements].get_average_of_function_gradients(current_solution,gradu) ;
                    lambda.value_list(fe_iv.get_quadrature_points(), lambda_values);
                    mu.value_list(fe_iv.get_quadrature_points(), mu_values); 
                    const std::vector<double>         &JxW = fe_iv.get_JxW_values();
                    const std::vector<Tensor<1, dim>> &normals =fe_iv.get_normal_vectors();                         
       
	                 for (unsigned int point=0; point<n_face_q_points; ++point)
                    {    auto strain = 0.5*(gradu[point]+transpose(gradu[point]));                                               
                         
                                      
                          interface_stress += (lambda_values[point]*trace(strain)*Identity  + 2*mu_values[point]*strain )

                                           *normals[point]*JxW[point] ;

                     }                        
	            }

        for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face : cell->face_indices())
      if (std::abs(cell->face(face)->center()[1]-1)<1e-3)
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();

}


  template <int dim>
  void DG_Elasticity_CZM<dim>::run()
  {
    
    for(unsigned int refinement_level = 0 ; refinement_level<1 ; refinement_level++  )
   { 
    triangulation.clear(); 
    //Mesh geenration 
    const Point<dim> P1, P2(1,2);
    const std::vector<unsigned int> repetitions{1, 2};  
    GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,P1,P2,true);
    triangulation.refine_global(0+refinement_level);
 

    for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face : cell->face_indices())
      if (std::abs(cell->face(face)->center()[1]-1)<1e-3)
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();
 

   setup_quadrature_point_history(); 
    double error ;
    unsigned nonlin = 0 ;
    

     std::ofstream forces("forces" + std::to_string(refinement_level) +   ".txt");
    unsigned int max_nonlin=0;
   for (unsigned int cycle = 0; cycle <1; ++cycle)
      {
     error=1;
      nonlin = 0 ; 
      
      disp[1]+=1e-3;
       std::cout<< "cycle = " << cycle << " and dislacement = " << disp[1] << "\n" ; 
      

       while((error>1e-4 && nonlin<10000)  )
       {
        if (cycle==0 & nonlin==0)
     
     setup_system(true);
     else
     setup_system(false);

   output_results(nonlin,refinement_level);
     nonlin++ ; 
     
     max_nonlin=std::max(max_nonlin,nonlin);

     std::cout<< "  non lin = " << nonlin << std::endl;
     assemble_system(nonlin,cycle);

     solve(cycle);
     
     error = system_rhs.l2_norm();
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
      Tensor<1,dim> inter;
     reaction(ten,3);
     interface_traction(inter);
      std::cout<<  "   reaction force = " << ten[0] << "\t" << ten[1] << "\n" ; 
      
        std::cout << std::endl;
         forces << disp[0] << "\t"<<disp[1] <<  "\t"  
         << ten[0]  << "\t" << ten[1] << "\t"
          << inter[0]  << "\t" << inter[1] << "\n" ;

     forces.close();
            std::cout<<"max nonlin iterations " << max_nonlin << " \n" ;
      }

         computing_timer.print_summary();
       computing_timer.reset();

  }

  }


} //closing namespace 



#endif 
