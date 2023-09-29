#ifndef PAHSEFILEDDG2_H
#define PAHSEFILEDDG2_H



//using namespace dealii;



// template <int dim>
// void print_mesh_info(const Triangulation<dim> &triangulation,
//                      const std::string &       filename)
// {
//   std::cout << "Mesh info:" << std::endl
//             << " dimension: " << dim << std::endl
//             << " no. of cells: " << triangulation.n_active_cells() << std::endl;

//   {
//     std::map<types::boundary_id, unsigned int> boundary_count;
//     for (const auto &face : triangulation.active_face_iterators())
//       if (face->at_boundary())
//         boundary_count[face->boundary_id()]++;

//     std::cout << " boundary indicators: ";
//     for (const std::pair<const types::boundary_id, unsigned int> &pair :
//          boundary_count)
//       {
//         std::cout << pair.first << '(' << pair.second << " times) ";
//       }
//     std::cout << std::endl;
//   }

//   // std::ofstream out(filename);
//   // GridOut       grid_out;
//   // grid_out.write_vtu(triangulation, out);
//   // std::cout << " written to " << filename << std::endl << std::endl;
// }





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

#include <map>

namespace pfdg
{
using namespace dealii;

template <int dim>
class Lambda : public Function<dim>
{
public:
  Lambda(types::material_id material_Id)
  {
    material_Id= material_Id ;
  }
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;

  private:
  types::material_id material_Id ; 
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
  Mu(types::material_id material_Id)
  {
    material_Id= material_Id ;
  }

  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;

    private:
  types::material_id material_Id ; 
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
 GC(bool is_at_interface)
  {
    this->is_at_interface= is_at_interface ;
  }
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
  private:
  bool is_at_interface; 
};

template <int dim>
double GC<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
  // if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <1e-3 ) 
  //     || 
  //     (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <1e-3 ))
  if(is_at_interface==true)
  return 0.5 ;
  else
  return 1.0 ; 
}

template <int dim>
class SIG : public Function<dim>
{
public:
 SIG(bool is_at_interface)
  {
    this->is_at_interface = is_at_interface ;
   
  }
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
    private:
  bool is_at_interface;
};

template <int dim>
double SIG<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
  // if( (std::abs(std::sqrt((p[0]-0.4)*(p[0]-0.4) + (p[1]-0)*(p[1]-0)) - 0.25)  <1e-3  ) 
  //     || 
  //     (std::abs(std::sqrt((p[0]-0.6)*(p[0]-0.6) + (p[1]-0.5)*(p[1]-0.5)) - 0.25)  <1e-3  )   )
  if(is_at_interface==true)
  return 50.0 ;
  else
  return 400.0 ; 
}

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
  double g_max;
};

template <int dim>
class PF_DG
{
public:
    PF_DG(const unsigned int degree);
    void run();

private:
    void import_mesh();
    void setup_system(const bool initial_step);
    void assemble_system();
    void solve_time_step();
    void output_results(const unsigned int cycle,const unsigned int refinement_level) const;
    void output_strain(const unsigned int cycle,const unsigned int refinement_level) const;
    void output_stress(const unsigned int cycle,const unsigned int refinement_level) const;


    void move_mesh(); 
    void reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id);
    void commulative_opening (double &opening, const bool &interface_or_bulk);


    Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;
    const MappingQ1<dim> mapping;
    AffineConstraints<double> constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    
    Vector<double> solution_np1;
    Vector<double> solution_np0;
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
  void PF_DG<dim>::setup_quadrature_point_history()
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
 
double get_penalty_factor(const unsigned int fe_degree,
                          const double       cell_extent_left,
                          const double       cell_extent_right)
{
  const unsigned int degree = std::max(1U, fe_degree);
  return degree * (degree + 1.) * 0.5 *
         (1. / cell_extent_left + 1. / cell_extent_right);
}


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
  // const Point<2> bottom_hole_origin(0.4, 0);
  // const Point<2> top_hole_origin(0.6, 0.5);


  // const SphericalManifold<2> bottom_manifold(bottom_hole_origin);
  // const SphericalManifold<2> top_manifold(bottom_hole_origin);


  GridIn<2> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("fenicsr2.5.2.msh");
  gridin.read_msh(f);

      std::ofstream out("grid-1.vtu");
    GridOut       grid_out;
    grid_out.write_vtu(triangulation, out);
    std::cout << " written to " << "grid-1.vtu" << std::endl << std::endl;

  //print_mesh_info(triangulation, "grid-1.vtu");


}


  template <int dim>
  void PF_DG<dim>::setup_system(const bool initial_step)
  {
     TimerOutput::Scope t(computing_timer, "set up");

     if (initial_step)
     {dof_handler.distribute_dofs(fe);
      solution_np0.reinit(dof_handler.n_dofs());
      solution_np1.reinit(dof_handler.n_dofs());
     }
     
     system_rhs.reinit(dof_handler.n_dofs());
     
   



   constraints.clear();
    // DoFTools::make_hanging_node_constraints(dof_handler, constraints);
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
   
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
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

        GC<dim> gc(cell->face(f)->user_flag_set());
        SIG<dim> sig(cell->face(f)->user_flag_set());

         gc.value_list(fe_iv.get_quadrature_points(), GC_values);
         sig.value_list(fe_iv.get_quadrature_points(), sig_values);
         




          PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());
          
              
       //   Assert(quadrature_points_history >=&quadrature_point_history.front(),ExcInternalError());
       //   Assert(quadrature_points_history <=&quadrature_point_history.back(),ExcInternalError());
        

          Tensor<1,dim> tangential;
         
         double law, delta,g ; // g is effective displacement 
         
         for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < n_dofs_face; ++i)  
             for (unsigned int j = 0; j < n_dofs_face; ++j) 
             {
             

              delta=  GC_values[point]/sig_values[point]/std::exp(1);
              tangential = cross_product_2d(normals[point]);
              g = std::sqrt(old_solution[point]*normals[point]*old_solution[point]*normals[point]+ 2.0*2.0*old_solution[point]*tangential*old_solution[point]*tangential);
               
             quadrature_points_history[point].g_max = std::max(quadrature_points_history[point].g_max,g); 
              g = quadrature_points_history[point].g_max ;
               

              law= GC_values[point]/(delta*delta) * std::exp(-g / delta ) ;
              
             


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
       
        Lambda<dim> lambda(cell->material_id()); 
        Mu<dim> mu(cell->material_id()); 

        std::vector<double> lambda_values(n_q_points);
        std::vector<double> mu_values(n_q_points);  

    const double extent1 = cell->measure() / cell->face(face_no)->measure();
    const double penalty = get_penalty_factor(degree, extent1, extent1);

      for (const auto &face : cell->face_iterators())
    if((cell->face(face_no)->boundary_id() == 1) || (cell->face(face_no)->boundary_id() == 2)) 
      { 
        lambda.value_list(q_points, lambda_values);
        mu.value_list(q_points, mu_values);

       for (unsigned int point = 0; point < n_q_points; ++point)
        { 
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=                                                                                      // - nu            
                    +(2*mu_values[point]+3*lambda_values[point])*penalty *                       
                      fe_fv[displacements].value(i, point) *   // v_h
                      fe_fv[displacements].value(j, point)     // u_h
                   *
                 JxW[point]; // dx
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {   
      
                copy_data.cell_rhs(i) += 
                                                                             
               +(2*mu_values[point]+3*lambda_values[point])*penalty *              
                   fe_fv[displacements].value(i, point) 
                   * (cell->face(face_no)->boundary_id()-1)*disp *
               JxW[point]; // dx

              }


                                       
        }
      } 
    };
   // 
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
  void PF_DG<dim>::solve_time_step()
  { 
    TimerOutput::Scope t(computing_timer, "solve");

   std::cout<< "solving" << std::endl ;

   SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution_np1,system_rhs);


  // SolverControl            solver_control(2500, 1e-6);
  // SolverCG<Vector<double>> cg(solver_control);
 
  // PreconditionSSOR<SparseMatrix<double>> preconditioner;
  // preconditioner.initialize(system_matrix, 1.2);
 
  // cg.solve(system_matrix, solution_np1, system_rhs, preconditioner);
  // std::cout<<"  CG converged after "<< solver_control.last_step() << std::endl; 
  
  }

  template <int dim>
  void PF_DG<dim>::output_results(const unsigned int cycle,const unsigned int refinement_level) const
  {
     
     
     std::vector<std::string> solution_names(dim, "u");
 
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler,
                           solution_np0,
                           solution_names,
                           interpretation);

  data_out.build_patches();
  std::ofstream output("solution"+ std::to_string(refinement_level) +"_" + std::to_string(cycle)  + ".vtu");
  // std::ofstream output("solution" + std::to_string(timestep_number) + ".vtu");
  data_out.write_vtu(output);
  
  }

    template <int dim>
    void PF_DG<dim>::output_strain(const unsigned int cycle,const unsigned int refinement_level) const
    { 
       StrainPostprocessor<dim> strain;
        
       DataOut<dim> data_out;
       data_out.attach_dof_handler (dof_handler);
        
                                 
       data_out.add_data_vector (solution_np0, strain);
       data_out.build_patches ();
       std::ofstream output("strain"+ std::to_string(refinement_level) +"_" + std::to_string(cycle)  + ".vtu");
       data_out.write_vtu (output);

    }


    template <int dim>
    void PF_DG<dim>::output_stress(const unsigned int cycle,const unsigned int refinement_level) const
    {
       
         StressPostprocessor<dim> stress;
        
       DataOut<dim> data_out;
       data_out.attach_dof_handler (dof_handler);
        
                                 
       data_out.add_data_vector (solution_np0, stress);
       data_out.build_patches ();
       std::ofstream output("stress"+ std::to_string(refinement_level) +"_" + std::to_string(cycle)  + ".vtu");
       data_out.write_vtu (output);


    }

 template <int dim>
 void PF_DG<dim>::move_mesh()
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
  void PF_DG<dim>::reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id)
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
                    fe_face_values[displacements].get_function_gradients(solution_np0,gradu) ;
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
 void PF_DG<dim>::commulative_opening (double &opening, const bool &interface_or_bulk)
{       
        opening = 0; 
        QGauss<dim-1> face_quadrature_formula(4);

        FEFaceValues<dim> fe_face_values (fe, 
                                          face_quadrature_formula,
                                          UpdateFlags(update_JxW_values| 
                                                      update_quadrature_points ));
        const unsigned int n_face_q_points = face_quadrature_formula.size();

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        


        std::vector<double> GC_values(n_face_q_points);
        std::vector<double> sig_values(n_face_q_points);       


           for (; cell!=endc; ++cell)      
                for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    	        if ((cell->face(face)->at_boundary() == false) && (cell->face(face)->user_flag_set() ==interface_or_bulk))
	            {    
                    fe_face_values.reinit (cell, face);      
                    const std::vector<double>         &JxW = fe_face_values.get_JxW_values(); 
                     PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(face)->user_pointer()); 
                      
                     GC<dim> gc(cell->face(face)->user_flag_set());
                      SIG<dim> sig(cell->face(face)->user_flag_set());
                               gc.value_list(fe_face_values.get_quadrature_points(), GC_values);
                         sig.value_list(fe_face_values.get_quadrature_points(), sig_values);
       
	                 for (unsigned int point=0; point<n_face_q_points; ++point)
                    {        
                         if(quadrature_points_history[point].g_max > GC_values[point]/sig_values[point]/std::exp(1)/10)
                         opening+=quadrature_points_history[point].g_max*JxW[point]; 
                     }                        
	            }                                         
}


// it takes 5 hours amount of time to run with tol 1e-4, dsip=0.0002 and 300 cylces , and 200,000 dofs per iteration 
template <int dim>
void  PF_DG<dim>::run()
{   
for(unsigned int refinement_level = 0 ; refinement_level<2 ; refinement_level++  )
   { 
    disp[0] = 0.0 ; 
    timestep_number = 0;
    AffineConstraints<double> constraints;
    constraints.close();

    import_mesh();
    triangulation.refine_global(refinement_level);

   // flag the material interface edges. Weird that this should be done after refinement 
    for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face : cell->face_indices())
      if ((cell->face(face)->at_boundary() == false)&& (cell->material_id() != cell->neighbor(face)->material_id()))
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();
   
    setup_quadrature_point_history(); 
      double error ;
     unsigned nonlin = 0 ;
     std::ofstream forces("forces" + std::to_string(refinement_level) +   ".txt");


    for (unsigned int cycle = 0; cycle <600; ++cycle)
    { error=1;
      nonlin = 0 ;
     disp[0]+=1e-4;
       std::cout<< "cycle = " << cycle << " and dislacement = " << disp[0] << "\n" ; 
      while(error>1e-4 &&  nonlin<200)
    { 
       std::cout<< "   nonlin  = " << nonlin << "\n" ; 
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
      Vector<double>  diff = solution_np0 ; 
      diff.add(-1,  solution_np1);
      error =  diff.l2_norm();
       std::cout << "    error = " << error << std::endl ; 
 
      //solution update
      solution_np0 =  solution_np1;


    }
    //output
    
    //std::cout << solution_np1[0] << " " <<  solution_np1[50000] << " " << solution_np1[25000]  <<"\n" ; 
  //if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
    if(cycle%2 == 0)
      {
        TimerOutput::Scope t(computing_timer, "output");
        output_results(cycle,refinement_level);
        output_strain(cycle,refinement_level);
        output_stress(cycle,refinement_level);
      }
    
   
     
     forces.open("forces" + std::to_string(refinement_level) +   ".txt", std::ios_base::app) ;

      Tensor<1,dim> ten; 
      double interface_opening, bulk_opening;
     reaction(ten,1);
     commulative_opening(interface_opening,0);
     commulative_opening(bulk_opening,1);

      std::cout<<  "   reaction force = " << ten[0] << "\t" << ten[1] << "\n" ; 

     forces << disp[0] <<  "\t"  << ten[0]  << "\t" << ten[1] << "\t" << interface_opening<< "\t" << bulk_opening << "\n" ; 
     forces.close();
  }
    

   computing_timer.print_summary();
       computing_timer.reset();

 }
}


}//closing napespace 


#endif
