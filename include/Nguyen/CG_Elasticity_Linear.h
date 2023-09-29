#ifndef CG_LINEAR_H
#define CG_LINEAR_H

// Contnious Galerkin series code for linear FEM with prm  

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
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <iostream>
#include <fstream>

namespace CG_Elasticity_Linear
{
  using namespace dealii;

  template <int dim>
  class ProblemParameters : public ParameterAcceptor
  {
  public:
    ProblemParameters() : ParameterAcceptor("main")
    {

      add_parameter("E", E, " ", this->prm, Patterns::Double());
      add_parameter("nu", nu, " ", this->prm, Patterns::Double());
      add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
      add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
      add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
      add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
      add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
    }

    int cycles = 0, initial_refinement = 0, degree = 0;
    double E = 0, nu = 0, displacementx = 0, displacementy = 0;
  };


  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues(Tensor<1, dim> displacement)
        : Function<dim>(dim),
          displacement(displacement)
    {
    }
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &value) const override;

  private:
    Tensor<1, dim> displacement;
  };

  template <int dim>
  void BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const
  {
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));

    values(0) = displacement[0];
    values(1) = displacement[1];
  }


  template <int dim>
  class Lambda : public Function<dim>
  {
  public:
    Lambda(int material_Id, const ProblemParameters<dim> &par)
    {
      this->material_Id = material_Id;
      this->E = par.E;
      this->nu = par.nu;
    }
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

  private:
    int material_Id;
    double E, nu;
  };

  template <int dim>
  double Lambda<dim>::value(const Point<dim> &p,
                            const unsigned int /*component*/) const
  { 
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    if (material_Id == 1)
      return lambda;
    else
    {
      std::cout << "exit at lambda value \n";
      exit(1);
    }
  }

  template <int dim>
  class Mu : public Function<dim>
  {
  public:
    Mu(int material_Id, const ProblemParameters<dim> &par)
    {
      this->material_Id = material_Id;
      this->E = par.E;
      this->nu = par.nu;
    }

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

  private:
    int material_Id;
    double E, nu;
  };

  template <int dim>
  double Mu<dim>::value(const Point<dim> &p,
                        const unsigned int /*component*/) const
  {
    double mu= E / (2 * (1 + nu));
    if (material_Id == 1)    
      return mu;
    else
    {
      std::cout << "exit at mu value \n";
      exit(1);
    }
  }

  template <int dim>
  class StressPostprocessor : public DataPostprocessorTensor<dim>
  {
  public:
    StressPostprocessor(const ProblemParameters<dim> &par)
        : DataPostprocessorTensor<dim>("stress",
                                       update_gradients)
    { 
      this->E = par.E;
      this->nu = par.nu;
    }

    virtual void evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data, 
                                        std::vector<Vector<double>> &computed_quantities) const override
    {
      AssertDimension(input_data.solution_gradients.size(),
                      computed_quantities.size());

      Tensor<2, dim> strain;
      Tensor<2, dim> stress;
      Tensor<2, dim> Identity;
      Identity[0][0] = 1.0; Identity[1][0] = 0.0;
      Identity[0][1] = 0.0; Identity[1][1] = 1.0;
      double lambda, mu;
      const typename DoFHandler<dim>::cell_iterator current_cell = input_data.template get_cell<dim>();
      
      if (current_cell->material_id() == 1)
      { 
        lambda =  E * nu / ((1 + nu) * (1 - 2 * nu));
        mu = E / (2 * (1 + nu));
      }
      else
      {
        std::cout << "exit at stress output \n";
        exit(1);
      }

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(),
                        (Tensor<2, dim>::n_independent_components));

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
          {
            strain[d][e] = (input_data.solution_gradients[p][d][e] + input_data.solution_gradients[p][e][d]) / 2;
          }

        stress = lambda * trace(strain) * Identity + 2 * mu * strain;


       for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = stress[d][e];

      }

      Tensor<2,dim> str_avg; 
       for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

       for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            str_avg[d][e]+=computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))];
      }


      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

       for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))]
            =str_avg[d][e]/(input_data.solution_gradients.size());
      }

      
    }
      private:
    double E, nu;
  };

  template <int dim>
  class StrainPostprocessor : public DataPostprocessorTensor<dim>
  {
  public:
    StrainPostprocessor()
        : DataPostprocessorTensor<dim>("strain",
                                       update_gradients)
    {
    }

    virtual void evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data, std::vector<Vector<double>> &computed_quantities) const override
    {
      AssertDimension(input_data.solution_gradients.size(),
                      computed_quantities.size());

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(),
                        (Tensor<2, dim>::n_independent_components));

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = 
            (input_data.solution_gradients[p][d][e] +input_data.solution_gradients[p][e][d]) /2;
      }
       Tensor<2,dim> str_avg; 
       for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

       for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            str_avg[d][e]+=computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))];
      }


      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

       for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))]
            =str_avg[d][e]/(input_data.solution_gradients.size());
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

  /**
   * @brief  Main class
   *
   * @tparam dim
   */

  template <int dim>
  class CG_Elasticity_Linear
  {
  public:
    CG_Elasticity_Linear(const unsigned int degree, const ProblemParameters<dim> &par);
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void output_results(unsigned int refinement_level, unsigned int cycle) const;
    void reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id);
    void interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress);
    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FESystem<dim> fe;

    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;

    const MappingQ1<dim> mapping;
    AffineConstraints<double> constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    const unsigned int degree;

    ConvergenceTable convergence_table;
    TimerOutput computing_timer;

    Tensor<1, dim> disp;

    const ProblemParameters<dim> &par;
  };

  template <int dim>
  CG_Elasticity_Linear<dim>::CG_Elasticity_Linear(const unsigned int degree, const ProblemParameters<dim> &par)
      : par(par),
        degree(degree),
        quadrature(degree +1),
        face_quadrature(degree + 1),
        mapping(),
        fe(FE_Q<dim>(degree), dim),
        dof_handler(triangulation),
        computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)

  {
    disp[0] = 0.0;
    disp[1] = 0.0;
  }

  template <int dim>
  void CG_Elasticity_Linear<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "set up");

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::component_wise(dof_handler);
    
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    std::cout << "dofs = " << dof_handler.n_dofs() << std::endl ; 
  }

  template <int dim>
  void CG_Elasticity_Linear<dim>::assemble_system()
  {
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][1] = 1;

    TimerOutput::Scope t(computing_timer, "Assemble");

    const auto cell_worker =
        [&](const auto &cell, auto &scratch_data, auto &copy_data)
    {
      scratch_data.fe_values.reinit(cell);
      const FEValues<dim> &fe_v = scratch_data.fe_values;

      const unsigned int dofs_per_cell = scratch_data.fe_values.get_fe().n_dofs_per_cell();
      copy_data.reinit(cell, dofs_per_cell);
      const auto &q_points = scratch_data.fe_values.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      const std::vector<double> &JxW = fe_v.get_JxW_values();

      std::vector<double> lambda_values(n_q_points);
      std::vector<double> mu_values(n_q_points);

      Lambda<dim> lambda(cell->material_id(), par);
      Mu<dim> mu(cell->material_id(), par);

      lambda.value_list(scratch_data.fe_values.get_quadrature_points(), lambda_values);
      mu.value_list(scratch_data.fe_values.get_quadrature_points(), mu_values);

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<2, dim>> old_solution_gradients(n_q_points);
      fe_v[displacements].get_function_gradients(solution, old_solution_gradients);

      for (unsigned int point = 0; point < n_q_points; ++point)
        for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        {
          Tensor<2, dim> straini = 0.5 * (fe_v[displacements].gradient(i, point) + transpose(fe_v[displacements].gradient(i, point)));

          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
          {
            Tensor<2, dim> strainj = 0.5 * (fe_v[displacements].gradient(j, point) + transpose(fe_v[displacements].gradient(j, point)));
            Tensor<2, dim> stressj = lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj;

            copy_data.cell_matrix(i, j) += scalar_product(straini, stressj) * JxW[point];
          }

          copy_data.cell_rhs(i) += 0.0;
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
                                               cdf.cell_rhs,
                                               cdf.joint_dof_indices,
                                               system_matrix,
                                               system_rhs);
      }
    };

    ScratchData<dim> scratch_data(mapping, fe, quadrature, face_quadrature);

    CopyData copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells);

     std::map<types::global_dof_index, double> boundary_values;

    const FEValuesExtractors::Scalar x_componenet(0);
    const FEValuesExtractors::Scalar y_componenet(1);

  //    boundary_values.clear();
  //   VectorTools::interpolate_boundary_values(dof_handler,
  //                                            0,
  //                                            Functions::ZeroFunction<dim>(2),
  //                                            boundary_values,
  //                                            fe.component_mask(x_componenet));
  //   MatrixTools::apply_boundary_values(boundary_values,
  //                                      system_matrix,
  //                                      solution,
  //                                      system_rhs);


    boundary_values.clear();
    VectorTools::interpolate_boundary_values(dof_handler,
                                             2,
                                             Functions::ZeroFunction<dim>(2),
                                             boundary_values
                                            // ,fe.component_mask(y_componenet)
                                             );
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

  boundary_values.clear();
    VectorTools::interpolate_boundary_values(dof_handler,
                                             3,
                                             BoundaryValues<dim>(disp),
                                             boundary_values
                                 //            , fe.component_mask(y_componenet)
                                             );

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                     system_rhs);
   }

  template <int dim>
  void CG_Elasticity_Linear<dim>::solve()
  {

    TimerOutput::Scope t(computing_timer, "solve");

    std::cout << "solving" << std::endl;

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
  }


  template <int dim>
  void CG_Elasticity_Linear<dim>::output_results(unsigned int refinement_level, unsigned int cycle) const
  {

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(dim,
                       DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> solution_names(dim, "u");

    DataOut<dim> data_out;

    data_out.add_data_vector(dof_handler, solution, solution_names, interpretation);
    const StrainPostprocessor<dim> strain;
    const StressPostprocessor<dim> stress(par);
    data_out.add_data_vector(dof_handler, solution, strain);
    data_out.add_data_vector(dof_handler, solution, stress);
    data_out.build_patches();

    std::ofstream output("solution" + std::to_string(refinement_level) + "_" + std::to_string(cycle) + ".vtu");

    data_out.write_vtu(output);



  }




  template <int dim>
  void CG_Elasticity_Linear<dim>::reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id)
  {
    reaction_stress = 0;
    QGauss<dim - 1> face_quadrature_formula(degree + 2 );

    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     UpdateFlags(update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_normal_vectors |
                                                 update_JxW_values));

    const unsigned int n_face_q_points = face_quadrature_formula.size();

    std::vector<Tensor<2, dim>> gradu(n_face_q_points);
    std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_face_q_points, std::vector<Tensor<1, dim>>(dim));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    std::vector<double> lambda_values(n_face_q_points);
    std::vector<double> mu_values(n_face_q_points);

    Lambda<dim> lambda(cell->material_id(), par);
    Mu<dim> mu(cell->material_id(), par);

    const FEValuesExtractors::Vector displacements(0);
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][0] = 0;
    Identity[0][1] = 0;
    Identity[1][1] = 1;
    Tensor<1,dim> tangential;
    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->boundary_id() == boundary_id)
        {
          fe_face_values.reinit(cell, face);
          fe_face_values[displacements].get_function_gradients(solution, gradu);
          lambda.value_list(fe_face_values.get_quadrature_points(), lambda_values);
          mu.value_list(fe_face_values.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_face_values.get_JxW_values();
          const std::vector<Tensor<1, dim>> &normals = fe_face_values.get_normal_vectors();

          for (unsigned int point = 0; point < n_face_q_points; ++point)
          {
            Tensor<2,dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2,dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) ;

            tangential = cross_product_2d(normals[point]);
            reaction_stress[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
            reaction_stress[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];
          }
        }
  }




  template <int dim>
  void CG_Elasticity_Linear<dim>::interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
  {
    interface_stress = 0;
    interface_jump = 0;
    Tensor<1, dim> interface_jump_temp;
    interface_jump_temp = 0;

    QGauss<dim - 1> face_quadrature_formula(degree + 2);

    FEInterfaceValues<dim> fe_iv(fe,
                                 face_quadrature_formula,
                                 UpdateFlags(update_values |
                                             update_gradients |
                                             update_quadrature_points |
                                             update_normal_vectors |
                                             update_JxW_values));

    const unsigned int n_face_q_points = face_quadrature_formula.size();

    std::vector<Tensor<2, dim>> gradu(n_face_q_points);
    std::vector<Tensor<1, dim>> jumpu(n_face_q_points);

    std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_face_q_points, std::vector<Tensor<1, dim>>(dim));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    std::vector<double> lambda_values(n_face_q_points);
    std::vector<double> mu_values(n_face_q_points);

    Lambda<dim> lambda(cell->material_id(), par);
    Mu<dim> mu(cell->material_id(), par);

    const FEValuesExtractors::Vector displacements(0);
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][0] = 0;
    Identity[0][1] = 0;
    Identity[1][1] = 1;
    Tensor<1,dim> tangential;
    int counter = 0;

    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->user_flag_set())
        {
          counter++;
          // I want to visit the face only once
          cell->face(face)->clear_user_flag();

          auto cellnb = cell->neighbor(face);
          int nface_number = cell->neighbor_face_no(face);

          fe_iv.reinit(cell, face, numbers::invalid_unsigned_int, cellnb, nface_number, numbers::invalid_unsigned_int);

          fe_iv[displacements].get_jump_in_function_values(solution, jumpu);
          fe_iv[displacements].get_average_of_function_gradients(solution, gradu);
          lambda.value_list(fe_iv.get_quadrature_points(), lambda_values);
          mu.value_list(fe_iv.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_iv.get_JxW_values();
          const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

          for (unsigned int point = 0; point < n_face_q_points; ++point)
          {
             Tensor<2,dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2,dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) ;
            tangential = cross_product_2d(normals[point]);

            interface_stress[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
            interface_stress[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];

            interface_jump_temp += jumpu[point];
          }

          interface_jump_temp[0] = interface_jump_temp[0] / n_face_q_points;
          interface_jump_temp[1] = interface_jump_temp[1] / n_face_q_points;

          interface_jump += interface_jump_temp / (std::pow(2, par.initial_refinement));
          interface_jump_temp = 0;
        }

    //If the geomtery is flipped you have a problem 
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const unsigned int face : cell->face_indices())
        if (std::abs(cell->face(face)->center()[1] - 1) < 1e-10)
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();
  }

  template <int dim>
  void CG_Elasticity_Linear<dim>::run()
  {

    for (unsigned int refinement_level = 0; refinement_level < 1; refinement_level++)
    {
      triangulation.clear();
      // Mesh geenration
      const Point<dim> P1, P2(1, 2);
      const std::vector<unsigned int> repetitions{1, 2};
      GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
      triangulation.refine_global(par.initial_refinement + refinement_level);
        GridTools::distort_random(0.3, triangulation, true, 556);

      Point<dim> p;
      for (const auto &cell : triangulation.active_cell_iterators())
      {
        p = cell->center();
        //  if(std::abs(p[1]-1)<0.02 )
        cell->set_material_id(1);
        //  else
        //  cell->set_material_id(2) ;
      }

      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int face : cell->face_indices())
          if (std::abs(cell->face(face)->center()[1] -1) < 1e-10)
            cell->face(face)->set_user_flag();
          else
            cell->face(face)->clear_user_flag();

      GridTools::rotate(0, triangulation);

      std::ofstream forces("forces" + std::to_string(refinement_level) + ".txt");

      for (unsigned int cycle = 0; cycle < par.cycles; ++cycle)
      {

        disp[0] += par.displacementx;
        disp[1] += par.displacementy;

        std::cout << " ####### cycle = " << cycle << " and displacement = " << disp[1] << " ###### \n";

        setup_system();

        assemble_system();

        solve();

        forces.open("forces" + std::to_string(refinement_level) + ".txt", std::ios_base::app);
        Tensor<1, dim> ten;
        Tensor<1, dim> jump;
        Tensor<1, dim> inter;

        reaction(ten, 3);
        interface_jump_and_traction(jump, inter);
        std::cout << "   reaction force = " << ten[0] << "\t" << ten[1] << "\n";

        std::cout << std::endl;
        forces << disp[0] << "\t" << disp[1] << "\t"
               << ten[0] << "\t" << ten[1] << "\t"
               << jump[0] << "\t" << jump[1] << "\t"
               << inter[0] << "\t" << inter[1] << "\n";

        forces.close();

        if (cycle % 1 == 0)
        {
          TimerOutput::Scope t(computing_timer, "output");
          output_results(refinement_level, cycle);
        }
      }
    }

    computing_timer.print_summary();
    computing_timer.reset();
  }

} // closing namespace

#endif