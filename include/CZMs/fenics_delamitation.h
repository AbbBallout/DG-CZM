#ifndef fenics_delamitation_H
#define fenics_delamitation_H

// Mode I delametation Bayat example

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

namespace fenics_delamitation
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
      add_parameter("sig_c", sig_c, " ", this->prm, Patterns::Double());
      add_parameter("delta_c", delta_c, " ", this->prm, Patterns::Double());
      add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
      add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
      add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
      add_parameter("symmetry", symmetry, " ", this->prm, Patterns::Integer());
      add_parameter("step_length", step_length, " ", this->prm, Patterns::Double());
      add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
      add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
      add_parameter("error", error, " ", this->prm, Patterns::Double());
      add_parameter("is_user_penalty", is_user_penalty, " ", this->prm, Patterns::Bool());
          add_parameter("is_implicit", is_implicit, " ", this->prm, Patterns::Bool());
          add_parameter("n", n, " ", this->prm, Patterns::Double());
          add_parameter("penalty_factor", penalty_factor, " ", this->prm, Patterns::Double());
      add_parameter("penalty", penalty, " ", this->prm, Patterns::Double());
      add_parameter("unloading", unloading, " ", this->prm, Patterns::Integer());
      add_parameter("reloading", reloading, " ", this->prm, Patterns::Integer());
         add_parameter("output_frequency", output_frequency, " ", this->prm, Patterns::Integer());
    }
    

    int cycles = 0, initial_refinement = 0, degree = 0, symmetry = 0, unloading=0, reloading=0,output_frequency=0;
    double E = 0, nu = 0, sig_c=0, delta_c=0, step_length = 0, displacementx = 0, n=0, displacementy = 0, error = 0, penalty = 0 , penalty_factor=0;
    bool is_user_penalty = 0 , is_implicit=0;
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
    double mu = E / (2 * (1 + nu));
    if (material_Id == 1)
      return mu;
    else
    {
      std::cout << "exit at mu value \n";
      exit(1);
    }
  }

  /**
   * @brief DG scratch data
   *
   * @tparam dim
   */

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
      Identity[0][0] = 1.0;
      Identity[1][0] = 0.0;
      Identity[0][1] = 0.0;
      Identity[1][1] = 1.0;
      double lambda, mu;
      const typename DoFHandler<dim>::cell_iterator current_cell = input_data.template get_cell<dim>();

      if (current_cell->material_id() == 1)
      {
        lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
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

      // Tensor<2, dim> str_avg;
      // for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      // {

      //   for (unsigned int d = 0; d < dim; ++d)
      //     for (unsigned int e = 0; e < dim; ++e)
      //       str_avg[d][e] += computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))];
      // }

      // for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      // {

      //   for (unsigned int d = 0; d < dim; ++d)
      //     for (unsigned int e = 0; e < dim; ++e)
      //       computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = str_avg[d][e] / (input_data.solution_gradients.size());
      // }
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
                (input_data.solution_gradients[p][d][e] + input_data.solution_gradients[p][e][d]) / 2;
      }
      Tensor<2, dim> str_avg;
      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            str_avg[d][e] += computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))];
      }

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = str_avg[d][e] / (input_data.solution_gradients.size());
      }
    }
  };

  template <int dim>
  struct PointHistory
  {
    bool is_damaged;
    double max_seperation;
    bool is_reunload;

    Tensor<1, dim> slope;
    Tensor<1, dim> traction_at_reunload;
    Tensor<1, dim> seperation_at_reunload;
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
  class fenics_delamitation
  {
  public:
    fenics_delamitation(const unsigned int degree, const ProblemParameters<dim> &par);
    void run();

  private:
    void setup_system(const bool initial_step);
     void assemble_system(const unsigned int cycle, const unsigned int non_lin, double &error);
    void solve(const unsigned int cycle);
    void output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const;

    void reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id);
    void interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress);
    void interface_jump_and_traction_at_qpoints(const unsigned int cycle);
    void plot_over_line();

    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FESystem<dim> fe;

    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;
        const QTrapezoid<dim - 1> damaged_face_quadrature;
    



    const MappingQ1<dim> mapping;
    AffineConstraints<double> constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;

    const unsigned int degree;

    ConvergenceTable convergence_table;
    TimerOutput computing_timer;

    Tensor<1, dim> disp;
    std::vector<PointHistory<dim>> quadrature_point_history;
    void setup_quadrature_point_history();
    const ProblemParameters<dim> &par;


    bool is_relaxed;
  };

  template <int dim>
  fenics_delamitation<dim>::fenics_delamitation(const unsigned int degree, const ProblemParameters<dim> &par)
      : par(par),
        degree(degree),
        quadrature(degree + 1),
        face_quadrature(degree+1),
        damaged_face_quadrature(),
        mapping(),
        fe(FE_DGQ<dim>(degree), dim),
        dof_handler(triangulation),
        computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)

  {
    disp[0] = 0.0;
    disp[1] = 0.0;
  }

  template <int dim>
  void fenics_delamitation<dim>::setup_quadrature_point_history()
  {

    triangulation.clear_user_data();

    {
      std::vector<PointHistory<dim>> tmp;
      quadrature_point_history.swap(tmp);
    }

    quadrature_point_history.resize(triangulation.n_active_faces() * damaged_face_quadrature.size());

    unsigned int history_index = 0;
    for (auto &face : triangulation.active_face_iterators())
    {
      face->set_user_pointer(&quadrature_point_history[history_index]);
      history_index += damaged_face_quadrature.size();
    }

    Assert(history_index == quadrature_point_history.size(), ExcInternalError());
  }

  template <int dim>
  void fenics_delamitation<dim>::setup_system(const bool initial_step)
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

    std::cout << "dofs = " << dof_handler.n_dofs() << std::endl;
  }

  template <int dim>
  void fenics_delamitation<dim>::assemble_system(const unsigned int cycle, const unsigned int non_lin,double &error)
  {

    TimerOutput::Scope t(computing_timer, "Assemble");

    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][1] = 1;

    using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    const auto cell_worker =
        [&](const Iterator &cell, auto &scratch_data, auto &copy_data)
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
      fe_v[displacements].get_function_gradients(current_solution, old_solution_gradients);

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

          Tensor<2, dim> old_strain = 0.5 * (old_solution_gradients[point] + transpose(old_solution_gradients[point]));
          Tensor<2, dim> old_stress = lambda_values[point] * trace(old_strain) * Identity + 2 * mu_values[point] * old_strain;
          copy_data.cell_rhs(i) += -scalar_product(straini, old_stress) * JxW[point];
        }
    };

    const auto boundary_worker = [&](const Iterator &cell,
                                     const unsigned int &face_no,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data)
    {
      scratch_data.fe_interface_values.reinit(cell, face_no);
      const FEFaceValuesBase<dim> &fe_fv =
          scratch_data.fe_interface_values.get_fe_face_values(0);

      const auto &q_points = fe_fv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      const unsigned int n_facet_dofs = fe_fv.get_fe().n_dofs_per_cell();
      const std::vector<double> &JxW = fe_fv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();
      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1, dim>> old_solution_values(n_q_points);
      fe_fv[displacements].get_function_values(current_solution, old_solution_values);

      std::vector<Tensor<2, dim>> gradu(n_q_points);
      std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_q_points, std::vector<Tensor<1, dim>>(dim));

      std::vector<double> lambda_values(n_q_points);
      std::vector<double> mu_values(n_q_points);

      Lambda<dim> lambda(cell->material_id(), par);
      Mu<dim> mu(cell->material_id(), par);

      lambda.value_list(q_points, lambda_values);
      mu.value_list(q_points, mu_values);

      Tensor<1, dim> zero;
      zero[0] = 0;
      zero[1] = 0;
      // ------- restricing x and y on boundaries 2 adn 3 -------

      if ((cell->face(face_no)->boundary_id() == 0))
        for (unsigned int point = 0; point < n_q_points; ++point)
        {
          double penalty = (par.penalty_factor) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
          if (par.is_user_penalty)
            penalty = par.penalty;

          Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

          for (unsigned int i = 0; i < n_facet_dofs; ++i)
          {
            Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

            for (unsigned int j = 0; j < n_facet_dofs; ++j)
            {
              Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

              copy_data.cell_matrix(i, j) +=
                  (-fe_fv[displacements].value(i, point) *
                       (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]

                   + par.symmetry *
                         (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                         fe_fv[displacements].value(j, point) // Symetry term

                   + penalty * fe_fv[displacements].value(i, point) * fe_fv[displacements].value(j, point)) *
                  JxW[point];
            }

            copy_data.cell_rhs(i) +=

                -(
                    -fe_fv[displacements].value(i, point) *
                        (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point]

                    + par.symmetry *
                          (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                          old_solution_values[point] // Symetry term

                    + penalty * fe_fv[displacements].value(i, point) * old_solution_values[point]) *
                    JxW[point]

                + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                      zero * JxW[point] // Symetry term

                + penalty *
                      fe_fv[displacements].value(i, point) * zero *
                      JxW[point]; // dx
          }
        }
        

      // Applying y BC only on the top and bottom walls

  const FEValuesExtractors::Scalar displacements_x(0);
      const FEValuesExtractors::Scalar displacements_y(1);

 bool is_at_corner=false;
      if(cell->face(face_no)->boundary_id() == 3)
      for(int i=0;i<4 ;i<i++) 
       if(cell->face(i)->boundary_id() == 1)
       is_at_corner=true;

      if ( cell->face(face_no)->boundary_id() == 3)
      if (is_at_corner)
     for (unsigned int point = 0; point < n_q_points; ++point)
            {
              double penalty = (1e+3) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
              if (par.is_user_penalty)
                penalty = par.penalty;

              Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

              for (unsigned int i = 0; i < n_facet_dofs; ++i)
              {
                Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

                for (unsigned int j = 0; j < n_facet_dofs; ++j)
                {
                  Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                  copy_data.cell_matrix(i, j) +=
                      (-fe_fv[displacements_y].value(i, point) *
                           (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point] * normals[point]

                       + par.symmetry *
                             (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point] * normals[point]

                        + par.symmetry *
                              (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                              (old_solution_values[point] * normals[point]) // Symetry term

                         +penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point] * normals[point]) *
                        JxW[point]

                    + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                           disp[1] * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * disp[1] *
                          JxW[point]; // dx
              }
            }


            if(cell->face(face_no)->boundary_id() == 2)
      for(int i=0;i<4 ;i<i++) 
       if(cell->face(i)->boundary_id() == 1)
       is_at_corner=true;
        
        if ( cell->face(face_no)->boundary_id() == 2)
      if (is_at_corner)
     for (unsigned int point = 0; point < n_q_points; ++point)
            { 
              

              double penalty = (1e+3) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
              if (par.is_user_penalty)
                penalty = par.penalty;

              Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

              for (unsigned int i = 0; i < n_facet_dofs; ++i)
              {
                Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

                for (unsigned int j = 0; j < n_facet_dofs; ++j)
                {
                  Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                  copy_data.cell_matrix(i, j) +=
                      (-fe_fv[displacements_y].value(i, point) *
                           (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point] * normals[point]

                       + par.symmetry *
                             (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point] * normals[point]

                        - par.symmetry *
                              (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                              (old_solution_values[point] * normals[point]) // Symetry term
                              //flipped normal

                         -penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point] * normals[point]) *
                        JxW[point] //flipped normal

                    + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point] *
                           (-disp[1]) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * (-disp[1]) *
                          JxW[point]; // dx
              }
            }



     };

   const auto face_worker = [&](const Iterator &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const Iterator &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim> &scratch_data,
                                 CopyData &copy_data)
    {
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      const unsigned int n_dofs_face = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      copy_data_face.cell_rhs.reinit(n_dofs_face);
      const std::vector<double> &JxW = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      std::vector<double> lambda_values(n_q_points);
      std::vector<double> mu_values(n_q_points);

      Lambda<dim> lambda(cell->material_id(), par);
      Mu<dim> mu(cell->material_id(), par);

      lambda.value_list(q_points, lambda_values);
      mu.value_list(q_points, mu_values);

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1, dim>> old_solution_jumps(n_q_points);
      fe_iv[displacements].get_jump_in_function_values(current_solution, old_solution_jumps);

      std::vector<Tensor<2, dim>> gradu(n_q_points);
      fe_iv[displacements].get_average_of_function_gradients(current_solution, gradu);



      if (cell->face(f)->user_flag_set()==false)
      for (unsigned int point = 0; point < n_q_points; ++point)
      {
        Tensor<2, dim> old_strain = 0.5 * (gradu[point] + transpose(gradu[point]));
        Tensor<2, dim> old_stress = (lambda_values[point] * trace(old_strain) * Identity + 2 * mu_values[point] * old_strain);


      double penalty = (par.penalty_factor) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
        if (par.is_user_penalty)
          penalty = par.penalty;

           
       

        for (unsigned int i = 0; i < n_dofs_face; ++i)
        {
          Tensor<2, dim> straini = 0.5 * (fe_iv[displacements].average_of_gradients(i, point) + transpose(fe_iv[displacements].average_of_gradients(i, point)));

          for (unsigned int j = 0; j < n_dofs_face; ++j)
          {
            Tensor<2, dim> strainj = 0.5 * (fe_iv[displacements].average_of_gradients(j, point) + transpose(fe_iv[displacements].average_of_gradients(j, point)));
            
              copy_data_face.cell_matrix(i, j) +=
                  (-fe_iv[displacements].jump_in_values(i, point) *
                       (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]

                   + par.symmetry *
                         (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                         fe_iv[displacements].jump_in_values(j, point) // Symetry term

                   + penalty * fe_iv[displacements].jump_in_values(i, point) * fe_iv[displacements].jump_in_values(j, point)) *
                  JxW[point]; // dx
          }
         
            copy_data_face.cell_rhs(i) +=
                -(
                    -fe_iv[displacements].jump_in_values(i, point) *
                        old_stress * normals[point]

                    + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                          old_solution_jumps[point]

                    + penalty * fe_iv[displacements].jump_in_values(i, point) * old_solution_jumps[point]) *
                JxW[point];
        }
      }
    };

    {
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

    ScratchData<dim> scratch_data(mapping, fe, quadrature,face_quadrature);

    CopyData copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);

   }



  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    const auto cell_worker2 =
        [&](const Iterator &cell, auto &scratch_data, auto &copy_data)
    {};

     const auto boundary_worker2 = [&](const Iterator &cell,
                                     const unsigned int &face_no,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data)
    {};

    const auto face_worker2 = [&](const Iterator &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const Iterator &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim> &scratch_data,
                                 CopyData &copy_data)
    { 
      
 FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      const unsigned int n_dofs_face = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      copy_data_face.cell_rhs.reinit(n_dofs_face);
      const std::vector<double> &JxW = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      std::vector<double> lambda_values(n_q_points);
      std::vector<double> mu_values(n_q_points);

      Lambda<dim> lambda(cell->material_id(), par);
      Mu<dim> mu(cell->material_id(), par);

      lambda.value_list(q_points, lambda_values);
      mu.value_list(q_points, mu_values);

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1, dim>> old_solution_jumps(n_q_points);
      fe_iv[displacements].get_jump_in_function_values(current_solution, old_solution_jumps);

      std::vector<Tensor<2, dim>> gradu(n_q_points);
      fe_iv[displacements].get_average_of_function_gradients(current_solution, gradu);

      PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());
      Tensor<1, dim> g, g_eff;
      Tensor<1, dim> tr, tr_eff; // traction
      Tensor<1, dim> tangential;
      Tensor<2, dim> TCZ;
      Tensor<1, dim> TCZ_res;
      double seperation;
      double n = par.n;
      bool check;

      if (cell->face(f)->user_flag_set())
        for (unsigned int point = 0; point < n_q_points; ++point)
        { 

              double penalty = (par.penalty_factor) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
          if (par.is_user_penalty)
            penalty = par.penalty;

  
          g = old_solution_jumps[point];
          g = -g;
          g_eff[0] = g[0];
          g_eff[1] = std::abs(g[1]) / 2 + g[1] / 2;
          seperation = g_eff.norm();

          if (non_lin == 1)
            if (quadrature_points_history[point].is_damaged)
              quadrature_points_history[point].max_seperation = std::max(seperation, quadrature_points_history[point].max_seperation);

          check = quadrature_points_history[point].is_reunload;

          // if(error<par.error)
          {
            if ((quadrature_points_history[point].max_seperation - seperation) > 1e-5)
            {
              quadrature_points_history[point].is_reunload = true;
            }
            else
              quadrature_points_history[point].is_reunload = false;
          }

          if (non_lin == 1 && quadrature_points_history[point].is_reunload == false)
          {

            quadrature_points_history[point].traction_at_reunload = tr_eff;
            quadrature_points_history[point].seperation_at_reunload = g_eff;
          }
                    double sig_c, delta_c;

          sig_c = par.sig_c;
          delta_c = par.delta_c;

          if (quadrature_points_history[point].is_reunload)
          {
           
            if (par.is_implicit)
              // fenics
              TCZ_res = sig_c * std::exp(1 - quadrature_points_history[point].max_seperation / delta_c) * g_eff / delta_c;
            else
            // shifted
            {
              TCZ_res[0] = g_eff[0] * quadrature_points_history[point].traction_at_reunload[0] / quadrature_points_history[point].seperation_at_reunload[0];
              TCZ_res[1] = g_eff[1] * quadrature_points_history[point].traction_at_reunload[1] / quadrature_points_history[point].seperation_at_reunload[1];
            }
            TCZ = 0.0;

            TCZ[0][0] = 1.0;

            if (g[1] >= 0)
              TCZ[1][1] = 1.0;

            if (par.is_implicit)
              // fenics
              TCZ = TCZ * sig_c * std::exp(1 - quadrature_points_history[point].max_seperation / delta_c) / delta_c;
            // shiftyed
            else
            {
              TCZ[0][0] = TCZ[0][0] * quadrature_points_history[point].traction_at_reunload[0] / quadrature_points_history[point].seperation_at_reunload[0];
              TCZ[1][1] = TCZ[1][1] * quadrature_points_history[point].traction_at_reunload[1] / quadrature_points_history[point].seperation_at_reunload[1];
            }
            if (std::abs(quadrature_points_history[point].seperation_at_reunload[0]) < 1e-12)
            {
              TCZ[0][0] = 0;
              TCZ_res[0] = 0;
            }

            if (std::abs(quadrature_points_history[point].seperation_at_reunload[1]) < 1e-12)
            {
              TCZ[1][1] = 0;
              TCZ_res[1] = 0;
            }

            if (g[1] < 0)
            {
              TCZ_res[1] += -penalty * std::pow(-g[1], n);
              TCZ[1][1] += penalty * n * std::pow(-g[1], n - 1);
            }
          }

          else
          {
            if (par.is_implicit)
              // fenics
              TCZ_res = sig_c * std::exp(1 - seperation / delta_c) * g_eff / delta_c;
            else
              // shifted
              TCZ_res = sig_c * std::exp(-seperation / delta_c) * (g_eff + par.delta_c * quadrature_points_history[point].slope) / delta_c;

            if (g[1] >= 0)
            {
              if (std::abs(seperation) < 1e-12)
                seperation = 1e-12;

              if (par.is_implicit)
              {
                // fenics
                TCZ[0][0] = 1 - std::pow(g_eff[0], 2) / (delta_c * seperation);
                TCZ[0][1] = -g_eff[0] * g_eff[1] / (delta_c * seperation);
                TCZ[1][0] = -g_eff[0] * g_eff[1] / (delta_c * seperation);
                TCZ[1][1] = 1 - std::pow(g_eff[1], 2) / (delta_c * seperation);
              }
              else
              {
                // shifted
                TCZ[0][0] = 1 - (g_eff[0]) * (g_eff[0] + par.delta_c * quadrature_points_history[point].slope[0]) / (delta_c * seperation);
                TCZ[0][1] = -(g_eff[1]) * (g_eff[0] + par.delta_c * quadrature_points_history[point].slope[0]) / (delta_c * seperation);
                TCZ[1][0] = -(g_eff[0]) * (g_eff[1] + par.delta_c * quadrature_points_history[point].slope[1]) / (delta_c * seperation);
                TCZ[1][1] = 1 - (g_eff[1]) * (g_eff[1] + par.delta_c * quadrature_points_history[point].slope[1]) / (delta_c * seperation);
              }
            }
            else
            {
              if (std::abs(seperation) < 1e-12)
                seperation = -1e-12;

              if (par.is_implicit)
              // fenics
              {
                TCZ[0][0] = 1 - std::pow(g_eff[0], 2) / (delta_c * seperation);
                TCZ[0][1] = 0.0;
                TCZ[1][0] = 0.0;
                TCZ[1][1] = 0.0;
              }
              else
              {
                // shited
                TCZ[0][0] = 1 - (g_eff[0]) * (g_eff[0] + par.delta_c * quadrature_points_history[point].slope[0]) / (delta_c * seperation);
                TCZ[0][1] = 0.0;
                TCZ[1][0] = 0.0;
                TCZ[1][1] = 0.0;
              }
            }

            if (par.is_implicit)
              // fenics
              TCZ = TCZ * sig_c * std::exp(1 - seperation / delta_c) / delta_c;
            else
              // shifted
              TCZ = TCZ * sig_c * std::exp(-seperation / delta_c) / delta_c;

            if (g[1] < 0)
            {
              TCZ_res[1] += -penalty * std::pow(-g[1], n);
              TCZ[1][1] += penalty * n * std::pow(-g[1], n - 1);
              
            }
          }

          Tensor<2, dim> old_strain = 0.5 * (gradu[point] + transpose(gradu[point]));
          Tensor<2, dim> old_stress = (lambda_values[point] * trace(old_strain) * Identity + 2 * mu_values[point] * old_strain);

          tangential = cross_product_2d(normals[point]);
          tr[0] = scalar_product(old_stress, outer_product(normals[point], tangential));
          tr[1] = scalar_product(old_stress, outer_product(normals[point], normals[point]));
          tr_eff[0] = tr[0];
          tr_eff[1] = std::abs(tr[1]) / 2 + tr[1] / 2;

          check = quadrature_points_history[point].is_damaged;
          if (cell->face(f)->user_flag_set())
            if (TCZ_res.norm() > par.sig_c)
              quadrature_points_history[point].is_damaged = true;

          if (check == false && quadrature_points_history[point].is_damaged == true)
            quadrature_points_history[point].slope = TCZ_res / TCZ_res.norm();

        



          if (cell->face(f)->user_flag_set() == true)
          {
            std::cout << "Is damaged? " << quadrature_points_history[point].is_damaged << "\n";
            std::cout << "at q point  " << q_points[point][0] << "\n";
            std::cout << "jump\t" << g[0] << "\t" << g[1] << "\n"
                      << "separation " << seperation << "\n";
            std::cout << "slope " << quadrature_points_history[point].slope << "\n";
            std::cout << "T_eff\t" << tr_eff[0] << "\t" << tr_eff[1] << "\n"
                      << "T_eff_norm " << tr_eff.norm() << "\n";
            std::cout << "TCZ\t" << TCZ[0][0] << "\t" << TCZ[0][1] << "\t" << TCZ[1][0] << "\t" << TCZ[1][1] << "\n";
            std::cout << "TCZ_res\t" << TCZ_res[0] << "\t" << TCZ_res[1] << "\n \n";
          }

          for (unsigned int i = 0; i < n_dofs_face; ++i)
          {
            Tensor<2, dim> straini = 0.5 * (fe_iv[displacements].average_of_gradients(i, point) + transpose(fe_iv[displacements].average_of_gradients(i, point)));

            for (unsigned int j = 0; j < n_dofs_face; ++j)
            {
              Tensor<2, dim> strainj = 0.5 * (fe_iv[displacements].average_of_gradients(j, point) + transpose(fe_iv[displacements].average_of_gradients(j, point)));
              if (par.is_implicit || quadrature_points_history[point].is_damaged)
                copy_data_face.cell_matrix(i, j) += (fe_iv[displacements].jump_in_values(i, point)) * (TCZ * fe_iv[displacements].jump_in_values(j, point)) * JxW[point];
              else
                copy_data_face.cell_matrix(i, j) +=
                    (-fe_iv[displacements].jump_in_values(i, point) *
                         (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]

                     + par.symmetry *
                           (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                           fe_iv[displacements].jump_in_values(j, point) // Symetry term

                     + penalty * fe_iv[displacements].jump_in_values(i, point) * fe_iv[displacements].jump_in_values(j, point)) *
                    JxW[point]; // dx
            }
            if (par.is_implicit || quadrature_points_history[point].is_damaged)
              copy_data_face.cell_rhs(i) +=
                  -(-1) * fe_iv[displacements].jump_in_values(i, point) * TCZ_res * JxW[point];
            else
              copy_data_face.cell_rhs(i) +=
                  -(
                      -fe_iv[displacements].jump_in_values(i, point) *
                          old_stress * normals[point]

                      + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
                            old_solution_jumps[point]

                      + penalty * fe_iv[displacements].jump_in_values(i, point) * old_solution_jumps[point]) *
                  JxW[point];
          }
        }
    };


        {
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

    ScratchData<dim> scratch_data(mapping, fe, quadrature, damaged_face_quadrature);

    CopyData copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker2,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker2,
                          face_worker2);

   }

  }


  template <int dim>
  void fenics_delamitation<dim>::solve(const unsigned int cycle)
  {

    TimerOutput::Scope t(computing_timer, "solve");

    std::cout << "solving" << std::endl;

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update, system_rhs);

    double alpha;


    if(is_relaxed)
    alpha = par.step_length;
    else
    alpha = 1;


    current_solution.add(alpha, newton_update);
  }


  template <int dim>
  void fenics_delamitation<dim>::output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const
  {

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(dim,
                       DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> solution_names(dim, "u");

    DataOut<dim> data_out;

    data_out.add_data_vector(dof_handler, current_solution, solution_names, interpretation);
    const StrainPostprocessor<dim> strain;
    const StressPostprocessor<dim> stress(par);
    data_out.add_data_vector(dof_handler, current_solution, strain);
    data_out.add_data_vector(dof_handler, current_solution, stress);
    data_out.build_patches();

    std::ofstream output("solution" + std::to_string(refinement_level) + "_" + std::to_string(cycle) + ".vtu");

  //  std::ofstream output("solution" + std::to_string(refinement_level) + "_" + std::to_string(cycle) + "_" + std::to_string(non_lin) + ".vtu");

    data_out.write_vtu(output);
  }

  template <int dim>
  void fenics_delamitation<dim>::reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id)
  {
    reaction_stress = 0;
    QGauss<dim - 1> face_quadrature_formula(degree + 2);

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
    Tensor<1, dim> tangential;
    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->boundary_id() == boundary_id)
        {
          fe_face_values.reinit(cell, face);
          fe_face_values[displacements].get_function_gradients(current_solution, gradu);
          lambda.value_list(fe_face_values.get_quadrature_points(), lambda_values);
          mu.value_list(fe_face_values.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_face_values.get_JxW_values();
          const std::vector<Tensor<1, dim>> &normals = fe_face_values.get_normal_vectors();

          for (unsigned int point = 0; point < n_face_q_points; ++point)
          {
            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);

            tangential = cross_product_2d(normals[point]);
            reaction_stress[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
            reaction_stress[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];
          }
        }
  }

  template <int dim>
  void fenics_delamitation<dim>::interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
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
    Tensor<1, dim> tangential;
    int counter = 0;

    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->user_flag_set() &&cell->face(face)->center()[0]>50-5*std::pow(2,-(float)par.initial_refinement)-1e-10 )
        {
          counter++;
          // I want to visit the face only once
          cell->face(face)->clear_user_flag();

          auto cellnb = cell->neighbor(face);
          int nface_number = cell->neighbor_face_no(face);

          fe_iv.reinit(cell, face, numbers::invalid_unsigned_int, cellnb, nface_number, numbers::invalid_unsigned_int);

          fe_iv[displacements].get_jump_in_function_values(current_solution, jumpu);
          fe_iv[displacements].get_average_of_function_gradients(current_solution, gradu);
          lambda.value_list(fe_iv.get_quadrature_points(), lambda_values);
          mu.value_list(fe_iv.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_iv.get_JxW_values();
          const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

          for (unsigned int point = 0; point < n_face_q_points; ++point)
          {
            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);
            tangential = cross_product_2d(normals[point]);

            interface_stress[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
            interface_stress[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];

            interface_jump_temp += jumpu[point];
          }

          interface_jump_temp = interface_jump_temp / n_face_q_points;


          interface_jump += interface_jump_temp;
          interface_jump_temp = 0;
        }
    interface_jump=interface_jump/counter;
   
    // If the geomtery is flipped you have a problem
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const unsigned int face : cell->face_indices())
            if (std::abs(cell->face(face)->center()[1] - 1.0) < 1e-5)
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();
  }

template <int dim>
void fenics_delamitation<dim>::interface_jump_and_traction_at_qpoints(const unsigned int cycle)
{
     Tensor<1, dim>  interface_stress;
    Tensor<1, dim> interface_jump;
   

    

    FEInterfaceValues<dim> fe_iv(fe,
                                 damaged_face_quadrature,
                                 UpdateFlags(update_values |
                                             update_gradients |
                                             update_quadrature_points |
                                             update_normal_vectors |
                                             update_JxW_values));

    const unsigned int n_face_q_points = damaged_face_quadrature.size();

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
    Tensor<1, dim> tangential;
   int counter = 0;
 

    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->user_flag_set() )
        {
          
          // I want to visit the face only once
          cell->face(face)->clear_user_flag();

          auto cellnb = cell->neighbor(face);
          int nface_number = cell->neighbor_face_no(face);

          fe_iv.reinit(cell, face, numbers::invalid_unsigned_int, cellnb, nface_number, numbers::invalid_unsigned_int);
          auto& q_points = fe_iv.get_quadrature_points();
          fe_iv[displacements].get_jump_in_function_values(current_solution, jumpu);
          fe_iv[displacements].get_average_of_function_gradients(current_solution, gradu);
          lambda.value_list(fe_iv.get_quadrature_points(), lambda_values);
          mu.value_list(fe_iv.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_iv.get_JxW_values();
          const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
          
          for (unsigned int point = 0; point < n_face_q_points; ++point)
          {  
            counter++;
            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);
            tangential = cross_product_2d(normals[point]);

            interface_stress[0] = scalar_product(stress, outer_product(normals[point], tangential));
            interface_stress[1] = scalar_product(stress, outer_product(normals[point], normals[point]));

            interface_jump = -jumpu[point];

            
            if(cycle==0)
            {
              std::ofstream q_point("q_points/q_point_" + std::to_string(counter) + ".txt");
            }
            else
            {
            std::ofstream q_point;
            q_point.open("q_points/q_point_" + std::to_string(counter) + ".txt", std::ios_base::app);
            q_point<<q_points[point][0]<< "\t"<< interface_jump[0] << "\t" << interface_jump[1]<<"\t" << interface_stress[0] << "\t" <<interface_stress[1]<<"\n" ;
            q_point.close();
            } 

          }

        }
   
   
    // If the geomtery is flipped you have a problem
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const unsigned int face : cell->face_indices())
            if (std::abs(cell->face(face)->center()[1] - 1.0) < 1e-5)
          cell->face(face)->set_user_flag();
        else
          cell->face(face)->clear_user_flag();
}


 template <int dim>
  void fenics_delamitation<dim>::plot_over_line()
  {
    QGauss<dim> quadrature_formula(degree + 2);

    FEValues<dim> fe_v(fe,quadrature_formula,
                                 UpdateFlags(update_values |
                                             update_gradients |
                                             update_quadrature_points |
                                             update_JxW_values));

     const unsigned int n_q_points = quadrature_formula.size();
      std::vector<Tensor<2, dim>> gradu(n_q_points);   
      std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_q_points, std::vector<Tensor<1, dim>>(dim));
     
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
       std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);

    Lambda<dim> lambda(cell->material_id(), par);
    Mu<dim> mu(cell->material_id(), par);

     const FEValuesExtractors::Vector displacements(0);
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][0] = 0;
    Identity[0][1] = 0;
    Identity[1][1] = 1;

    for (; cell != endc; ++cell)
    if (std::abs(cell->center()[1]-1.17)<1e-2)
  {
          fe_v.reinit(cell);
          fe_v[displacements].get_function_gradients(current_solution, gradu);
          lambda.value_list(fe_v.get_quadrature_points(), lambda_values);
          mu.value_list(fe_v.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_v.get_JxW_values();
         Tensor<2, dim> avg_stress;
         for (unsigned int point = 0; point < n_q_points; ++point)
          {
            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);
            avg_stress+=stress;
          }
      
      std::cout<<cell->center()[0]<<" \n" 
      << avg_stress[0][0]<<"\t"<< avg_stress[0][1] <<"\n"
      << avg_stress[1][0]<<"\t"<< avg_stress[1][1] <<" \n \n";


 
  }


  }


  template <int dim>
  void fenics_delamitation<dim>::run()
  {

    for (unsigned int refinement_level = 0; refinement_level < 1; refinement_level++)
    {
      triangulation.clear();
      // Mesh geenration
      const Point<dim> P1(0,0), P2(50, 2);
      const std::vector<unsigned int> repetitions{10, 4};
      GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
      triangulation.refine_global(par.initial_refinement + refinement_level);

       //GridTools::distort_random(0.3, triangulation, true);

  //  {   
    
    //   for(unsigned int i=0; i<1 ; i++)                                    
    //    { 
    //     typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
		// 				            endc = dof_handler.end();

    //     for (; cell != endc; ++cell)
    //       {
    //         // if ((cell->center()[1]) > 1.9 || (cell->center()[1]) < 0.1 )
    //         //   if ((cell->center()[0] > 0.9)  || (cell->center()[0] < 0.1))
    //         //     cell->set_refine_flag();
           
    //       // std::cout << std::pow(2.0,-(float)(i))<<"\n" ; 

    //        if (std::abs(cell->center()[1] -2.5 )<0.5)
    //        cell->set_refine_flag();
    //       }

    //        triangulation.execute_coarsening_and_refinement();
    //    }
    //  }
  
  

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
          if (std::abs(cell->face(face)->center()[1] - 1.0) < 1e-5)
            cell->face(face)->set_user_flag();
          else
            cell->face(face)->clear_user_flag();




      GridTools::rotate(0, triangulation);
      setup_quadrature_point_history();
      double error;
      double error2;

      unsigned int non_lin = 0;

      std::ofstream forces("forces" + std::to_string(refinement_level) + ".txt");

      unsigned int max_nonlin = 0;
      unsigned int at_cycle=0;
      for (unsigned int cycle = 0; cycle < par.cycles; ++cycle)
      {
        error = 1;
        error2=1;
        non_lin = 0;
 

          if(cycle<par.unloading || cycle > par.reloading)
      {disp[0]+= par.displacementx;
     disp[1]+=par.displacementy;
      }
     else
     {disp[0]= disp[0]- par.displacementx;
     disp[1]= disp[1]-par.displacementy;
      }

        std::cout << " ####### cycle = " << cycle << " and displacement = " << disp[1] << " ###### \n";

         while (error > par.error && error2 > par.error && non_lin < 30)
        {

          if (cycle == 0 && non_lin == 0)
            setup_system(true);
          else
            setup_system(false);

          // if (cycle == 0 && non_lin == 0)
          // {
          // assemble_system_fixed_point();
          // solve_fixed_point();
          // }

          max_nonlin = std::max(max_nonlin, non_lin);

          if (max_nonlin == non_lin)
            at_cycle = cycle;

          non_lin++;

          assemble_system(cycle, non_lin, error);

          error = system_rhs.l2_norm();

          //  std::cout << "  Residual: " << compute_residual(0) << std::endl;
          std::cout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;

          if (error > par.error)
          {
            std::cout << "  non lin = " << non_lin << std::endl;
            solve(cycle);
            std::cout << "  update l2_norm " << newton_update.l2_norm() << std::endl;
          }
          else
          {
            assemble_system(cycle, non_lin, error);
            error2 = system_rhs.l2_norm();
          }

          std::cout << "max nonlin iterations " << max_nonlin << "  at_cycle " << at_cycle << "\n------- \n ";


        }



          if (cycle % par.output_frequency == 0)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results(cycle, non_lin, refinement_level);
          }

        forces.open("forces" + std::to_string(refinement_level) + ".txt", std::ios_base::app);
        Tensor<1, dim> ten;
        Tensor<1, dim> jump;
        Tensor<1, dim> inter;

        reaction(ten, 2);
        interface_jump_and_traction(jump, inter);
        std::cout << "   reaction force = " << ten[0] << "\t" << ten[1] << "\n";

        std::cout << std::endl;
        forces << disp[0] << "\t" << disp[1] << "\t"
               << ten[0] << "\t" << ten[1] << "\t"
               << jump[0] << "\t" << jump[1] << "\t"
               << inter[0] << "\t" << inter[1] << "\n";

        forces.close();

        
        interface_jump_and_traction_at_qpoints(cycle) ; 

        std::cout << "max nonlin iterations " << max_nonlin <<  "  at_cycle " << at_cycle << "\n------- \n ";
      }

     // plot_over_line();

      computing_timer.print_summary();
      computing_timer.reset();
    }
  }

} // closing namespace

#endif
