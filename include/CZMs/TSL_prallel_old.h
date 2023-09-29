#ifndef TSL_prllel_old_OLD_H
#define TSL_prllel_old_OLD_H

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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

// #define FORCE_USE_OF_TRILINOS
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
    !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

namespace TSL_prllel_old
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
      add_parameter("Ep", Ep, " ", this->prm, Patterns::Double());
      add_parameter("nup", nup, " ", this->prm, Patterns::Double());
      add_parameter("geometry", geometry, "", this->prm, Patterns::Selection("reactangle|sheer|peeling|delamitation|signle-edge|matrix"));
      add_parameter("is_implicit", is_implicit, " ", this->prm, Patterns::Bool());
      add_parameter("tsl_law", tsl_law, " ", this->prm, Patterns::Selection("exponential|linear"));
      add_parameter("sig_c", sig_c, " ", this->prm, Patterns::Double());
      add_parameter("delta_c", delta_c, " ", this->prm, Patterns::Double());
      add_parameter("lambda_f", lambda_f, " ", this->prm, Patterns::Double());
      add_parameter("regularization", regularization, " ", this->prm, Patterns::Double());
      add_parameter("m", m, " ", this->prm, Patterns::Double());
      add_parameter("beta", beta, " ", this->prm, Patterns::Double());
      add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
      add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
      add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
      add_parameter("quadrature", quadrature, " ", this->prm, Patterns::Selection("trap|gauss|lobatto|gauss-reduced|lobatto-over"));
      add_parameter("symmetry", symmetry, " ", this->prm, Patterns::Integer());
      add_parameter("step_length", step_length, " ", this->prm, Patterns::Double());
      add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
      add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
      add_parameter("error", error, " ", this->prm, Patterns::Double());
      add_parameter("n", n, " ", this->prm, Patterns::Double());
      add_parameter("is_user_penalty", is_user_penalty, " ", this->prm, Patterns::Bool());
      add_parameter("is_special_BC", is_special_BC, " ", this->prm, Patterns::Bool());
      add_parameter("penalty", penalty, " ", this->prm, Patterns::Double());
      add_parameter("penalty_factor", penalty_factor, " ", this->prm, Patterns::Double());
      add_parameter("unloading", unloading, " ", this->prm, Patterns::Integer());
      add_parameter("reloading", reloading, " ", this->prm, Patterns::Integer());
      add_parameter("output_frequency", output_frequency, " ", this->prm, Patterns::Integer());
    }

    int cycles = 0, initial_refinement = 0, degree = 0, symmetry = 0, unloading = 0, reloading = 0, output_frequency = 0;
    double E = 0, nu = 0, Ep = 0, nup = 0, regularization, n = 0, m = 0, sig_c = 0, delta_c = 0, lambda_f = 0, beta = 0, step_length = 0, displacementx = 0, displacementy = 0, error = 0, penalty = 0, penalty_factor = 0;
    bool is_user_penalty = 0, is_implicit = 0, is_special_BC = 0;
    std::string tsl_law = "exponential", quadrature = "trap", geometry = "reactangle";
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
      this->Ep = par.Ep;
      this->nup = par.nup;
    }
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

  private:
    int material_Id;
    double E, nu, Ep, nup;
  };

  template <int dim>
  double Lambda<dim>::value(const Point<dim> &p,
                            const unsigned int /*component*/) const
  {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    double lambdap = (Ep * nup) / ((1 + nup) * (1 - 2 * nup));

    if (material_Id == 1)
      return lambda;
    else
    {
      return lambdap;
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
      this->Ep = par.Ep;
      this->nup = par.nup;
    }

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

  private:
    int material_Id;
    double E, nu, Ep, nup;
  };

  template <int dim>
  double Mu<dim>::value(const Point<dim> &p,
                        const unsigned int /*component*/) const
  {
    double mu = E / (2 * (1 + nu));
    double mup = Ep / (2 * (1 + nup));

    if (material_Id == 1)
      return mu;
    else
    {
      return mup;
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
      this->Ep = par.Ep;
      this->nup = par.nup;
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
        lambda = Ep * nup / ((1 + nup) * (1 - 2 * nup));
        mu = Ep / (2 * (1 + nup));
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
    double E, nu, Ep, nup;
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

      // for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      // {

      //   for (unsigned int d = 0; d < dim; ++d)
      //     for (unsigned int e = 0; e < dim; ++e)
      //       computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = str_avg[d][e] / (input_data.solution_gradients.size());
      // }
    }
  };

  template <int dim>
  struct PointHistory
  {
    bool is_damaged;
    bool is_fully_damaged;
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
  class TSL_prllel_old
  {
  public:
    TSL_prllel_old(const unsigned int degree, const ProblemParameters<dim> &par);
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

    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FESystem<dim> fe;

    const QGauss<dim> quadrature;
    const QGauss<dim - 1> face_quadrature;
    std::unique_ptr<const Quadrature<dim - 1>> damaged_face_quadrature;
    const MappingQ1<dim> mapping;
    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix system_matrix;

    LA::MPI::Vector ghosted_solution;
    LA::MPI::Vector distributed_solution;

    LA::MPI::Vector ghosted_newton_update;
    LA::MPI::Vector distributed_newton_update;

    LA::MPI::Vector system_rhs;

    const unsigned int degree;

    ConvergenceTable convergence_table;

    Tensor<1, dim> disp;
    std::vector<PointHistory<dim>> quadrature_point_history;
    void setup_quadrature_point_history();

    const ProblemParameters<dim> &par;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    ConditionalOStream pcout;
    TimerOutput computing_timer;
    bool is_relaxed;
  };

  template <int dim>
  TSL_prllel_old<dim>::TSL_prllel_old(const unsigned int degree, const ProblemParameters<dim> &par)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        triangulation(mpi_communicator),
        degree(degree),
        quadrature(degree + 1),
        face_quadrature(degree + 1),
        mapping(),
        fe(FE_DGQ<dim>(degree), dim),
        dof_handler(triangulation),
        pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
        computing_timer(mpi_communicator,
                        pcout,
                        TimerOutput::never,
                        TimerOutput::wall_times)

  {
    disp = 0.0;

    if (par.quadrature == "trap")
      damaged_face_quadrature = std::make_unique<const QTrapezoid<dim - 1>>();
    if (par.quadrature == "gauss")
      damaged_face_quadrature = std::make_unique<const QGauss<dim - 1>>(degree + 1);
    if (par.quadrature == "lobatto")
      damaged_face_quadrature = std::make_unique<const QGaussLobatto<dim - 1>>(degree + 1);
    if (par.quadrature == "gauss-reduced")
      damaged_face_quadrature = std::make_unique<const QGauss<dim - 1>>(degree);
    if (par.quadrature == "lobatto-over")
      damaged_face_quadrature = std::make_unique<const QGaussLobatto<dim - 1>>(degree + 2);
  }

  template <int dim>
  void TSL_prllel_old<dim>::setup_quadrature_point_history()
  {

    triangulation.clear_user_data();

    {
      std::vector<PointHistory<dim>> tmp;
      quadrature_point_history.swap(tmp);
    }

    quadrature_point_history.resize(triangulation.n_active_faces() * damaged_face_quadrature->size());

    unsigned int history_index = 0;
    for (auto &face : triangulation.active_face_iterators())
    {
      face->set_user_pointer(&quadrature_point_history[history_index]);
      history_index += damaged_face_quadrature->size();
    }

    Assert(history_index == quadrature_point_history.size(), ExcInternalError());
  }

  template <int dim>
  void TSL_prllel_old<dim>::setup_system(const bool initial_step)
  {
    TimerOutput::Scope t(computing_timer, "set up");
    if (initial_step)
    {
      dof_handler.distribute_dofs(fe);
      locally_owned_dofs = dof_handler.locally_owned_dofs();
      locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
      ghosted_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
      distributed_solution.reinit(locally_owned_dofs, mpi_communicator);
    }
    distributed_solution = ghosted_solution;

    ghosted_newton_update.reinit(locally_owned_dofs,
                                 locally_relevant_dofs,
                                 mpi_communicator);
    distributed_newton_update.reinit(locally_owned_dofs, mpi_communicator);

    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         constraints,
                                         /*keep_constrained_dofs = */ false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);
    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    pcout << "dofs = " << dof_handler.n_dofs() << std::endl;
  }

  template <int dim>
  void TSL_prllel_old<dim>::assemble_system(const unsigned int cycle, const unsigned int non_lin, double &error)
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
      fe_v[displacements].get_function_gradients(ghosted_solution, old_solution_gradients);

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
      fe_fv[displacements].get_function_values(ghosted_solution, old_solution_values);

      std::vector<Tensor<2, dim>> gradu(n_q_points);
      std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_q_points, std::vector<Tensor<1, dim>>(dim));

      std::vector<double> lambda_values(n_q_points);
      std::vector<double> mu_values(n_q_points);

      Lambda<dim> lambda(cell->material_id(), par);
      Mu<dim> mu(cell->material_id(), par);

      lambda.value_list(q_points, lambda_values);
      mu.value_list(q_points, mu_values);

      // ------- restricing x and y on boundaries 2 adn 3 -------

      bool is_at_corner = false;
      if (cell->face(face_no)->boundary_id() == 3 || cell->face(face_no)->boundary_id() == 2)
        for (int i = 0; i < 4; i < i++)
          if (cell->face(i)->boundary_id() == 1)
            is_at_corner = true;

      if (cell->face(face_no)->boundary_id() == 1)
        for (int i = 0; i < 4; i < i++)
          if (cell->face(i)->boundary_id() == 3)
            is_at_corner = true;

      if (par.geometry == "reactangle")
      {
        if (par.is_special_BC == false)
        {
          // if ((cell->face(face_no)->boundary_id() == 2))
          if ((cell->face(face_no)->boundary_id() == 2) || cell->face(face_no)->boundary_id() == 3 )
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
                          ((cell->face(face_no)->boundary_id() - 2) * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * ((cell->face(face_no)->boundary_id() - 2) * disp) *
                          JxW[point]; // dx
              }
            }

          // if (cell->face(face_no)->boundary_id() == 1 && is_at_corner)
          //   for (unsigned int point = 0; point < n_q_points; ++point)
          //   {

          //     double penalty = (par.penalty_factor) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
          //     if (par.is_user_penalty)
          //       penalty = par.penalty;

          //     Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

          //     for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //     {
          //       Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

          //       for (unsigned int j = 0; j < n_facet_dofs; ++j)
          //       {
          //         Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

          //         copy_data.cell_matrix(i, j) +=
          //             (-fe_fv[displacements].value(i, point) *
          //                  (lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]

          //              + par.symmetry *
          //                    (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
          //                    fe_fv[displacements].value(j, point) // Symetry term

          //              + penalty * fe_fv[displacements].value(i, point) * fe_fv[displacements].value(j, point)) *
          //             JxW[point];
          //       }

          //       copy_data.cell_rhs(i) +=

          //           -(
          //               -fe_fv[displacements].value(i, point) *
          //                   (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point]

          //               + par.symmetry *
          //                     (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] *
          //                     old_solution_values[point] // Symetry term

          //               + penalty * fe_fv[displacements].value(i, point) * old_solution_values[point]) *
          //               JxW[point]

          //           + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point]
          //                  * disp * JxW[point] // Symetry term

          //           + penalty *
          //                 fe_fv[displacements].value(i, point)  * disp *
          //                 JxW[point]; // dx
          //     }
          //   }
        }

        else
        {
          const FEValuesExtractors::Scalar displacements_x(0);
          const FEValuesExtractors::Scalar displacements_y(1);

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
                      (-fe_fv[displacements_x].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[0]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                             fe_fv[displacements_x].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_x].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[0]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                              old_solution_values[point][0] // Symetry term

                        + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                          (0.0) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_x].value(i, point) * (0.0) *
                          JxW[point]; // dx
              }
            }

          if ((cell->face(face_no)->boundary_id() == 2))
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[1]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[1]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                              old_solution_values[point][1] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                          (0.0) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * (0.0) *
                          JxW[point]; // dx
              }
            }

          if (cell->face(face_no)->boundary_id() == 3 && is_at_corner)
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]) * normals[point]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point]) * normals[point] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point]) * normals[point]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point]) *
                              old_solution_values[point] * normals[point] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point] * normals[point]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point]) * normals[point] *
                          ((cell->face(face_no)->boundary_id() - 2) * disp * normals[point]) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * ((cell->face(face_no)->boundary_id() - 2) * disp * normals[point]) *
                          JxW[point]; // dx
              }
            }

          // if(cell->face(face_no)->boundary_id() == 3 && is_at_corner)
          //  for (unsigned int point = 0; point < n_q_points; ++point)
          //   for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //    copy_data.cell_rhs(i)+=fe_fv[displacements].value(i, point)*disp*JxW[point];
        }
      }

      if (par.geometry == "sheer")
      {
        if (par.is_special_BC == false)
        {
          if ((cell->face(face_no)->boundary_id() == 0) && cell->material_id() != 1)
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
                          -disp * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * -disp *
                          JxW[point]; // dx
              }
            }

          if ((cell->face(face_no)->boundary_id() == 1) && cell->material_id() == 1)
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
                          disp * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * disp *
                          JxW[point]; // dx
              }
            }
        }

        else
        {
          std::cout << "No special BC were contructed for sheer exapmle\n";
          exit(1);
        }
      }
      if (par.geometry == "peeling")
      {
        if (par.is_special_BC == false)
        {
          if (cell->face(face_no)->boundary_id() == 2)
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
                          ((cell->face(face_no)->boundary_id() - 2) * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * ((cell->face(face_no)->boundary_id() - 2) * disp) *
                          JxW[point]; // dx
              }
            }
          const FEValuesExtractors::Scalar displacements_x(0);
          const FEValuesExtractors::Scalar displacements_y(1);

          if (cell->face(face_no)->boundary_id() == 1 && cell->material_id() == 1)
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

                    + par.symmetry * (lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * disp * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * disp *
                          JxW[point]; // dx
              }
            }
        }

        else
        {
          const FEValuesExtractors::Scalar displacements_x(0);
          const FEValuesExtractors::Scalar displacements_y(1);

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
                      (-fe_fv[displacements_x].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[0]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                             fe_fv[displacements_x].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_x].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[0]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                              old_solution_values[point][0] // Symetry term

                        + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[0] *
                          (0.0) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_x].value(i, point) * (0.0) *
                          JxW[point]; // dx
              }
            }

          if ((cell->face(face_no)->boundary_id() == 2))
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[1]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[1]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                              old_solution_values[point][1] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                          (0.0) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * (0.0) *
                          JxW[point]; // dx
              }
            }

          // if (cell->face(face_no)->boundary_id() == 3 && is_at_corner)
          //   for (unsigned int point = 0; point < n_q_points; ++point)
          //   {

          //     double penalty = (par.penalty_factor) * ((2 * mu_values[point]) * (3 * lambda_values[point] + 2 * mu_values[point]) / (lambda_values[point] + mu_values[point])) / (cell->measure());
          //     if (par.is_user_penalty)
          //       penalty = par.penalty;

          //     Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

          //     for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //     {
          //       Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

          //       for (unsigned int j = 0; j < n_facet_dofs; ++j)
          //       {
          //         Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

          //         copy_data.cell_matrix(i, j) +=
          //             (-fe_fv[displacements_y].value(i, point) *
          //                  ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point]) * normals[point]

          //              + par.symmetry *
          //                    ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point]) * normals[point] *
          //                    fe_fv[displacements_y].value(j, point) // Symetry term

          //              + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
          //             JxW[point];
          //       }

          //       copy_data.cell_rhs(i) +=

          //           -(
          //               -fe_fv[displacements_y].value(i, point) *
          //                   ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point]) * normals[point]

          //               + par.symmetry *
          //                     ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point] * normals[point]) *
          //                     old_solution_values[point] * normals[point] // Symetry term

          //               + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point] * normals[point]) *
          //               JxW[point]

          //           + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point]) * normals[point] *
          //                 disp * normals[point] * JxW[point] // Symetry term

          //           + penalty *
          //                 fe_fv[displacements_y].value(i, point) * disp * normals[point] *
          //                 JxW[point]; // dx
          //     }
          //   }

          if (cell->face(face_no)->boundary_id() == 1 && cell->material_id() == 1)
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[1]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[1]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                              old_solution_values[point][1] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                          disp[1] * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * disp[1] *
                          JxW[point]; // dx
              }
            }

          // if(cell->face(face_no)->boundary_id() == 3 && is_at_corner)
          //  for (unsigned int point = 0; point < n_q_points; ++point)
          //   for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //    copy_data.cell_rhs(i)+=fe_fv[displacements].value(i, point)*disp*JxW[point];
        }
      }

      if (par.geometry == "delamitation")
      {
        if (par.is_special_BC == false)
        {
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
                          ((cell->face(face_no)->boundary_id()) * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * ((cell->face(face_no)->boundary_id()) * disp) *
                          JxW[point]; // dx
              }
            }

          if (cell->face(face_no)->boundary_id() == 3 && is_at_corner)
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
                          (disp)*JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * disp *
                          JxW[point]; // dx
              }
            }

          if (cell->face(face_no)->boundary_id() == 2 && is_at_corner)
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
                          (-disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * -disp *
                          JxW[point]; // dx
              }
            }

          //  if(cell->face(face_no)->boundary_id() == 3 && is_at_corner)
          //     for (unsigned int point = 0; point < n_q_points; ++point)
          //      for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //       copy_data.cell_rhs(i)+=fe_fv[displacements].value(i, point)*disp*JxW[point];

          //    if(cell->face(face_no)->boundary_id() == 2 && is_at_corner)
          //     for (unsigned int point = 0; point < n_q_points; ++point)
          //      for (unsigned int i = 0; i < n_facet_dofs; ++i)
          //       copy_data.cell_rhs(i)+=-fe_fv[displacements].value(i, point)*disp*JxW[point];
        }
        else
        {
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
                          (0.0 * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * (0.0 * disp) *
                          JxW[point]; // dx
              }
            }

          const FEValuesExtractors::Scalar displacements_x(0);
          const FEValuesExtractors::Scalar displacements_y(1);

          if (cell->face(face_no)->boundary_id() == 1 && cell->material_id() == 1)
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[1]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[1]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                              old_solution_values[point][1] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                          disp[1] * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * disp[1] *
                          JxW[point]; // dx
              }
            }

          if (cell->face(face_no)->boundary_id() == 1 && cell->material_id() != 1)
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
                      (-fe_fv[displacements_y].value(i, point) *
                           ((lambda_values[point] * trace(strainj) * Identity + 2 * mu_values[point] * strainj) * normals[point])[1]

                       + par.symmetry *
                             ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                             fe_fv[displacements_y].value(j, point) // Symetry term

                       + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                      JxW[point];
                }

                copy_data.cell_rhs(i) +=

                    -(
                        -fe_fv[displacements_y].value(i, point) *
                            ((lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain) * normals[point])[1]

                        + par.symmetry *
                              ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                              old_solution_values[point][1] // Symetry term

                        + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                        JxW[point]

                    + par.symmetry * ((lambda_values[point] * trace(straini) * Identity + 2 * mu_values[point] * straini) * normals[point])[1] *
                          -disp[1] * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements_y].value(i, point) * -disp[1] *
                          JxW[point]; // dx
              }
            }
        }
      }

      if (par.geometry == "signle-edge")
      {
        if (par.is_special_BC == false)
        {
          if ((cell->face(face_no)->boundary_id() == 2))
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
                          (0.0 * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * (0.0 * disp) *
                          JxW[point]; // dx
              }
            }

          if (cell->face(face_no)->boundary_id() == 3)
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
                          (disp)*JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * disp *
                          JxW[point]; // dx
              }
            }
        }
        else
        {
          pcout << "no special BC set for single edge geometry \n";
          exit(1);
        }
      }

      if (par.geometry == "matrix")
      {

           if (par.is_special_BC == false)
        {
            if ((cell->face(face_no)->boundary_id() == 2) || cell->face(face_no)->boundary_id() == 3)
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
                          ((cell->face(face_no)->boundary_id() - 2) * disp) * JxW[point] // Symetry term

                    + penalty *
                          fe_fv[displacements].value(i, point) * ((cell->face(face_no)->boundary_id() - 2) * disp) *
                          JxW[point]; // dx
              }
            }
        }
        else
        {
          pcout << "no specail BC set for matrix geometry \n";
          exit(1);
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
      fe_iv[displacements].get_jump_in_function_values(ghosted_solution, old_solution_jumps);

      std::vector<Tensor<2, dim>> gradu(n_q_points);
      fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu);

      if (cell->face(f)->user_flag_set() == false)
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

    using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    const auto cell_worker_empty =
        [&](const Iterator &cell, auto &scratch_data, auto &copy_data) {};

    const auto boundary_worker_empty = [&](const Iterator &cell,
                                           const unsigned int &face_no,
                                           ScratchData<dim> &scratch_data,
                                           CopyData &copy_data) {};

    const auto face_worker_tsl = [&](const Iterator &cell,
                                     const unsigned int &f,
                                     const unsigned int &sf,
                                     const Iterator &ncell,
                                     const unsigned int &nf,
                                     const unsigned int &nsf,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data)
    {
      if (cell->face(f)->user_flag_set())
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

        // double lambda_values;
        // double mu_values;
        // double E, nu;
        // E = par.E / 2 + par.Ep / 2;
        // nu = par.nu / 2 + par.nup / 2;
        // lambda_values = (E * nu) / ((1 + nu) * (1 - 2 * nu));
        // mu_values = E / (2 * (1 + nu));

        const FEValuesExtractors::Vector displacements(0);

        std::vector<Tensor<1, dim>> old_solution_jumps(n_q_points);
        fe_iv[displacements].get_jump_in_function_values(ghosted_solution, old_solution_jumps);

        std::vector<Tensor<2, dim>> gradu(n_q_points);
        fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu);

        PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());
        Tensor<1, dim> g, g_eff;
        Tensor<1, dim> tr, tr_eff; // traction
        Tensor<1, dim> tangential;
        Tensor<2, dim> TCZ;
        Tensor<1, dim> TCZ_res;
        double seperation;
        double n = par.n;
        double beta = par.beta;
        bool check;
        Tensor<2, dim> Rot; // rotation

        for (unsigned int point = 0; point < n_q_points; ++point)
        {

          Tensor<2, dim> old_strain = 0.5 * (gradu[point] + transpose(gradu[point]));
          Tensor<2, dim> old_stress = (lambda_values[point] * trace(old_strain) * Identity + 2 * mu_values[point] * old_strain);
          double penalty = (par.penalty_factor) * (par.E/2+par.Ep/2) / (cell->measure());
          if (par.is_user_penalty)
            penalty = par.penalty;

          tangential = cross_product_2d(normals[point]);
          tr[0] = scalar_product(old_stress, outer_product(normals[point], tangential));
          tr[1] = scalar_product(old_stress, outer_product(normals[point], normals[point]));
          tr_eff[0] = tr[0] / beta;
          tr_eff[1] = std::abs(tr[1]) / 2 + tr[1] / 2;

          const double normaly = normals[point][1];
          const int sign_normalx = ((normals[point][0] > 0) ? 1 : ((normals[point][0] < 0) ? -1 : 0));

          Rot[0][0] = normaly;
          Rot[0][1] = -sign_normalx * std::sqrt(1 - normaly * normaly);
          Rot[1][0] = sign_normalx * std::sqrt(1 - normaly * normaly);
          Rot[1][1] = normaly;
          Rot = -Rot;
          g = Rot * old_solution_jumps[point];
          g_eff[0] = g[0] * beta;
          g_eff[1] = std::abs(g[1]) / 2 + g[1] / 2;
          seperation = g_eff.norm();

          check = quadrature_points_history[point].is_damaged;

          if (tr_eff.norm() > par.sig_c || g_eff.norm() > par.delta_c)
            quadrature_points_history[point].is_damaged = true;

          if (check == false && quadrature_points_history[point].is_damaged == true)
            quadrature_points_history[point].slope = tr_eff / tr_eff.norm();

          // quadrature_points_history[point].slope[0] = std::abs(quadrature_points_history[point].slope[0]);

          if (non_lin == 1)
            if (quadrature_points_history[point].is_damaged)
              quadrature_points_history[point].max_seperation = std::max(seperation, quadrature_points_history[point].max_seperation);

          check = quadrature_points_history[point].is_reunload;

          // if(error<par.error)
          // if(non_lin==1)
          {
            if ((quadrature_points_history[point].max_seperation - seperation) > 1e-8)
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

          if (par.is_implicit)
            if (seperation < 1e-12)
              seperation = 1e-12;

          double regularization = par.regularization;

          if (par.tsl_law == "exponential")
          {
            if (quadrature_points_history[point].is_reunload)
            {

              if (par.is_implicit)
                // fenics
                TCZ_res = sig_c * std::exp(1 - quadrature_points_history[point].max_seperation / delta_c) * g_eff / delta_c;
              else
                // explicit
                TCZ_res = sig_c * std::exp(-std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization) / delta_c) * g_eff / std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization);

              TCZ = 0.0;

              TCZ[0][0] = 1.0;

              if (g[1] >= 0)
                TCZ[1][1] = 1.0;

              if (par.is_implicit)
                // implicit
                TCZ = TCZ * sig_c * std::exp(1 - quadrature_points_history[point].max_seperation / delta_c) / delta_c;
              // explicit
              else
                TCZ = TCZ * sig_c * std::exp(-std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization) / delta_c) / std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization);

              // if (std::abs(quadrature_points_history[point].seperation_at_reunload[0]) < 1e-12)
              // {
              //   TCZ[0][0] = 0;
              //   TCZ_res[0] = 0;
              // }

              // if (std::abs(quadrature_points_history[point].seperation_at_reunload[1]) < 1e-12)
              // {
              //   TCZ[1][1] = 0;
              //   TCZ_res[1] = 0;
              // }

              if (g[1] < 0)
              {
                TCZ_res[1] += -penalty * std::pow(-g[1], n);
                TCZ[1][1] += penalty * n * std::pow(-g[1], n - 1);
              }
            }

            else
            { // we are laoding
              if (par.is_implicit)
                // fenics
                TCZ_res = sig_c * std::exp(1 - seperation / delta_c) * g_eff / delta_c;
              else

                TCZ_res = sig_c * std::exp(-std::sqrt(seperation * seperation + regularization) / delta_c) * g_eff / std::sqrt(seperation * seperation + regularization);

              if (g[1] >= 0)
              {

                if (par.is_implicit)
                {
                  // fenics
                  TCZ[0][0] = 1 - std::pow(g_eff[0], 2) / (delta_c * seperation);
                  TCZ[0][1] = -g_eff[0] * g_eff[1] / (delta_c * seperation);
                  TCZ[1][0] = -g_eff[0] * g_eff[1] / (delta_c * seperation);
                  TCZ[1][1] = 1 - std::pow(g_eff[1], 2) / (delta_c * seperation);
                }
                else
                { // explicit
                  TCZ[0][0] = 1 - (g_eff[0]) * (g_eff[0]) * (1 / (par.delta_c * std::sqrt((seperation * seperation + regularization))) + 1 / (seperation * seperation + regularization));
                  TCZ[0][1] = -(g_eff[1]) * (g_eff[0]) * (1 / (par.delta_c * std::sqrt((seperation * seperation + regularization))) + 1 / (seperation * seperation + regularization));
                  TCZ[1][0] = -(g_eff[0]) * (g_eff[1]) * (1 / (par.delta_c * std::sqrt((seperation * seperation + regularization))) + 1 / (seperation * seperation + regularization));
                  TCZ[1][1] = 1 - (g_eff[1]) * (g_eff[1]) * (1 / (par.delta_c * std::sqrt((seperation * seperation + regularization))) + 1 / (seperation * seperation + regularization));
                }
              }
              else
              {

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

                  TCZ[0][0] = 1 - (g_eff[0]) * (g_eff[0]) * (1 / (par.delta_c * std::sqrt((seperation * seperation + regularization))) + 1 / (seperation * seperation + regularization));
                  TCZ[0][1] = 0.0;
                  TCZ[1][0] = 0.0;
                  TCZ[1][1] = 0.0;
                }
              }

              if (par.is_implicit)
                // fenics
                TCZ = TCZ * sig_c * std::exp(1 - seperation / delta_c) / delta_c;
              else
                TCZ = TCZ * sig_c * std::exp(-std::sqrt(seperation * seperation + regularization) / delta_c) / std::sqrt(seperation * seperation + regularization);

              if (g[1] < 0)
              {
                TCZ_res[1] += -penalty * std::pow(-g[1], n);
                TCZ[1][1] += penalty * n * std::pow(-g[1], n - 1);
              }
            }

          } // end if exponential TSL

          if (par.tsl_law == "linear")
          {

            if (quadrature_points_history[point].is_reunload)
            {

              if (par.is_implicit)
                // implcit
                if (quadrature_points_history[point].is_damaged)
                  TCZ_res = sig_c * std::pow((par.lambda_f + par.delta_c - quadrature_points_history[point].max_seperation) / par.lambda_f, par.m) * g_eff / quadrature_points_history[point].max_seperation;
                else
                  TCZ_res = sig_c * g_eff / delta_c;
              else
              // explicit
              {
                TCZ_res = g_eff * (sig_c / std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization)) * std::pow((par.lambda_f - std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization)) / par.lambda_f, par.m);
              }
              TCZ = 0.0;

              TCZ[0][0] = 1.0;

              if (g[1] >= 0)
                TCZ[1][1] = 1.0;

              if (par.is_implicit)
                // implicit
                TCZ = TCZ * (sig_c / quadrature_points_history[point].max_seperation) * std::pow((par.lambda_f + par.delta_c - quadrature_points_history[point].max_seperation) / par.lambda_f, par.m);
              // explicit
              else
              {
                TCZ = TCZ * (sig_c / std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization)) * std::pow((par.lambda_f - std::sqrt(quadrature_points_history[point].max_seperation * quadrature_points_history[point].max_seperation + regularization)) / par.lambda_f, par.m);
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
            /////////we are laoding/////////
            else
            {
              if (par.is_implicit)
              { // implcit
                if (quadrature_points_history[point].is_damaged)
                  TCZ_res = sig_c * std::pow((par.lambda_f + par.delta_c - seperation) / par.lambda_f, par.m) * g_eff / seperation;
                else
                  TCZ_res = sig_c * g_eff / delta_c;
              }
              else
              {
                // explicit
                TCZ_res = g_eff * (sig_c / std::sqrt(seperation * seperation + regularization)) * std::pow((par.lambda_f - std::sqrt(seperation * seperation + regularization)) / par.lambda_f, par.m);
              }
              if (g[1] >= 0)
              {

                if (par.is_implicit)
                {
                  // implicit
                  if (quadrature_points_history[point].is_damaged)
                  {
                    TCZ[0][0] = std::pow(beta, 3) * std::pow(g_eff[0], 2);
                    TCZ[0][1] = beta * beta * g_eff[0] * g_eff[1];
                    TCZ[1][0] = beta * g_eff[0] * g_eff[1];
                    TCZ[1][1] = std::pow(g_eff[1], 2);
                    TCZ = -TCZ * ((par.lambda_f + par.delta_c + seperation * (par.m - 1)) / (seperation * seperation * (par.lambda_f + par.delta_c - seperation)));
                    TCZ[0][0] += 1.0;
                    TCZ[1][1] += 1.0;
                    TCZ = TCZ * (par.sig_c / seperation) * std::pow((par.lambda_f + par.delta_c - seperation) / par.lambda_f, par.m);
                  }

                  else
                  {
                    TCZ[0][0] = sig_c / delta_c;
                    TCZ[1][1] = sig_c / delta_c;
                  }
                }
                else
                {

                  TCZ[0][0] = std::pow(beta, 3) * std::pow(g_eff[0], 2);
                  TCZ[0][1] = beta * beta * g_eff[0] * g_eff[1];
                  TCZ[1][0] = beta * g_eff[0] * g_eff[1];
                  TCZ[1][1] = std::pow(g_eff[1], 2);
                  TCZ = -TCZ * ((par.lambda_f + std::sqrt(seperation * seperation + regularization) * (par.m - 1)) / ((seperation * seperation + regularization) * (par.lambda_f - std::sqrt(seperation * seperation + regularization))));
                  TCZ[0][0] += 1.0;
                  TCZ[1][1] += 1.0;
                  TCZ = TCZ * (sig_c / std::sqrt(seperation * seperation + regularization)) * std::pow((par.lambda_f - std::sqrt(seperation * seperation + regularization)) / par.lambda_f, par.m);
                }
              }
              else
              {

                if (par.is_implicit)
                { // implicit
                  if (quadrature_points_history[point].is_damaged)
                  {
                    TCZ[0][0] = std::pow(beta, 3) * std::pow(g_eff[0], 2);
                    TCZ[0][1] = 0.0;
                    TCZ[1][0] = 0.0;
                    TCZ[1][1] = 0.0;
                    TCZ = -TCZ * ((par.lambda_f + par.delta_c + seperation * (par.m - 1)) / (seperation * seperation * (par.lambda_f + par.delta_c - seperation)));
                    TCZ[0][0] += 1.0;
                    TCZ[1][1] += 0.0;
                    TCZ = TCZ * (par.sig_c / (seperation + par.delta_c)) * std::pow((par.lambda_f - seperation) / par.lambda_f, par.m);
                  }
                  else
                  {
                    TCZ[0][0] = sig_c / delta_c;
                    TCZ[1][1] = sig_c / delta_c;
                  }
                }
                else
                {
                  // explicit
                  TCZ[0][0] = std::pow(beta, 3) * std::pow(g_eff[0], 2);
                  TCZ[0][1] = 0.0;
                  TCZ[1][0] = 0.0;
                  TCZ[1][1] = 0.0;
                  TCZ = -TCZ * ((par.lambda_f + std::sqrt(seperation * seperation + regularization) * (par.m - 1)) / ((seperation * seperation + regularization) * (par.lambda_f - std::sqrt(seperation * seperation + regularization))));

                  TCZ[0][0] += 1.0;
                  TCZ[1][1] += 0.0;
                  TCZ = TCZ * (sig_c / std::sqrt(seperation * seperation + regularization)) * std::pow((par.lambda_f - std::sqrt(seperation * seperation + regularization)) / par.lambda_f, par.m);
                }
              }

              if (g[1] < 0)
              {
                TCZ_res[1] += -penalty * std::pow(-g[1], n);
                TCZ[1][1] += penalty * n * std::pow(-g[1], n - 1);
              }
            }

            if (par.is_implicit == true)
              if (seperation > par.lambda_f + delta_c)
                quadrature_points_history[point].is_fully_damaged = true;

            if (par.is_implicit == false)
              if (seperation > par.lambda_f)
                quadrature_points_history[point].is_fully_damaged = true;

          } // end if linear TSL

          if (par.geometry == "delamitation")
            if (q_points[point][0] > 45)
            {
              quadrature_points_history[point].is_damaged = true;
              quadrature_points_history[point].is_fully_damaged = true;
            }

          if (par.geometry == "signle-edge")
            if (q_points[point][0] < 50)
            {
              quadrature_points_history[point].is_damaged = true;
              quadrature_points_history[point].is_fully_damaged = true;
            }

          //   if(par.geometry=="reactangle")
          //   if ( q_points[point][0]>0.5)
          //   {
          //      quadrature_points_history[point].is_damaged=true;
          // quadrature_points_history[point].is_fully_damaged=true;
          //   }

          if (quadrature_points_history[point].is_fully_damaged)
          {
            TCZ = 0.0;
            TCZ_res = 0.0;
          }

          {
            pcout << "Is damaged? " << quadrature_points_history[point].is_damaged << "\n";
            pcout << "at q point  " << q_points[point][0] << "\n";
            pcout << "Rotation " << Rot << "\n";
            pcout << "jump\t" << g[0] << "\t" << g[1] << "\n"
                  << "separation " << seperation << " max seperation  " << quadrature_points_history[point].max_seperation << "\n";
            pcout << "slope " << quadrature_points_history[point].slope << "\n";
            pcout << "T_eff\t" << tr_eff[0] << "\t" << tr_eff[1] << "\n"
                  << "T_eff_norm " << tr_eff.norm() << "\n";
            pcout << "TCZ\t" << TCZ[0][0] << "\t" << TCZ[0][1] << "\t" << TCZ[1][0] << "\t" << TCZ[1][1] << "\n";
            pcout << "TCZ_res\t" << TCZ_res[0] << "\t" << TCZ_res[1] << "\n \n";
          }

          for (unsigned int i = 0; i < n_dofs_face; ++i)
          {
            Tensor<2, dim> straini = 0.5 * (fe_iv[displacements].average_of_gradients(i, point) + transpose(fe_iv[displacements].average_of_gradients(i, point)));

            for (unsigned int j = 0; j < n_dofs_face; ++j)
            {
              Tensor<2, dim> strainj = 0.5 * (fe_iv[displacements].average_of_gradients(j, point) + transpose(fe_iv[displacements].average_of_gradients(j, point)));
              if (par.is_implicit || quadrature_points_history[point].is_damaged)
                copy_data_face.cell_matrix(i, j) += fe_iv[displacements].jump_in_values(i, point) * (transpose(Rot) * TCZ * Rot) * fe_iv[displacements].jump_in_values(j, point) * JxW[point];
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
                  -fe_iv[displacements].jump_in_values(i, point) * (transpose(Rot) * TCZ_res) * JxW[point];
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

      ScratchData<dim> scratch_data(mapping, fe, quadrature, face_quadrature);

      CopyData copy_data;

      MeshWorker::mesh_loop(dof_handler.begin_active(),
                            dof_handler.end(),
                            cell_worker,
                            copier,
                            scratch_data,
                            copy_data,
                            MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                            boundary_worker,
                            face_worker);
    }

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

      ScratchData<dim> scratch_data(mapping, fe, quadrature, *damaged_face_quadrature);

      CopyData copy_data;

      MeshWorker::mesh_loop(dof_handler.begin_active(),
                            dof_handler.end(),
                            cell_worker_empty,
                            copier,
                            scratch_data,
                            copy_data,
                            MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                            boundary_worker_empty,
                            face_worker_tsl);
    }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  template <int dim>
  void TSL_prllel_old<dim>::solve(const unsigned int cycle)
  {

    TimerOutput::Scope t(computing_timer, "solve");

    pcout << "solving" << std::endl;

    SolverControl solver_control;
#ifdef USE_PETSC_LA
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
#else
    TrilqinosWrappers::SolverDirect::AdditionalData add_data(false, "Amesos_Klu");
    TrilinosWrappers::SolverDirect solver(solver_control, add_data);
#endif

    solver.solve(system_matrix, distributed_newton_update, system_rhs);
    ghosted_newton_update = distributed_newton_update;

    double alpha;

    if (is_relaxed)
      alpha = par.step_length;
    else
      alpha = 1.0;

    distributed_solution.add(alpha, distributed_newton_update);
    ghosted_solution = distributed_solution;
  }

  template <int dim>
  void TSL_prllel_old<dim>::output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const
  {

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(dim,
                       DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> solution_names(dim, "u");

    DataOut<dim> data_out;

    data_out.add_data_vector(dof_handler, ghosted_solution, solution_names, interpretation);
    const StrainPostprocessor<dim> strain;
    const StressPostprocessor<dim> stress(par);
    data_out.add_data_vector(dof_handler, ghosted_solution, strain);
    data_out.add_data_vector(dof_handler, ghosted_solution, stress);
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();

    // std::ofstream output("solution" + std::to_string(refinement_level) + "_" + std::to_string(cycle) + ".vtu");

    // std::ofstream output("solution" + std::to_string(refinement_level) + "_" + std::to_string(cycle) + "_" + std::to_string(non_lin) + ".vtu");

    // data_out.write_vtu(output);

    data_out.write_vtu_with_pvtu_record("./output/", "solution", cycle, mpi_communicator, 2, 0);
  }

  // Note: degree+1 works but we miss outputs for >degree+1
  template <int dim>
  void TSL_prllel_old<dim>::reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id)
  {
    reaction_stress = 0.0;
    QGauss<dim - 1> face_quadrature_formula(degree + 1);

    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature,
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

    Tensor<1, dim> local_reaction;

    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->boundary_id() == boundary_id)
          {
            fe_face_values.reinit(cell, face);
            fe_face_values[displacements].get_function_gradients(ghosted_solution, gradu);
            lambda.value_list(fe_face_values.get_quadrature_points(), lambda_values);
            mu.value_list(fe_face_values.get_quadrature_points(), mu_values);
            const std::vector<double> &JxW = fe_face_values.get_JxW_values();
            const std::vector<Tensor<1, dim>> &normals = fe_face_values.get_normal_vectors();

            for (unsigned int point = 0; point < n_face_q_points; ++point)
            {
              Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
              Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);

              tangential = cross_product_2d(normals[point]);
              local_reaction[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
              local_reaction[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];
            }
          }

    reaction_stress = Utilities::MPI::sum(local_reaction, mpi_communicator);
  }

  template <int dim>
  void TSL_prllel_old<dim>::interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
  {
    interface_stress = 0;
    interface_jump = 0;

    QGauss<dim - 1> face_quadrature_formula(degree + 1);

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

    Tensor<1, dim> local_interface_traction;
    Tensor<1, dim> local_interface_jump;
    unsigned int local_counter = 0;
    unsigned int counter = 0;

    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->user_flag_set())
          {

            auto cellnb = cell->neighbor(face);
            int nface_number = cell->neighbor_face_no(face);

            fe_iv.reinit(cell, face, numbers::invalid_unsigned_int, cellnb, nface_number, numbers::invalid_unsigned_int);

            fe_iv[displacements].get_jump_in_function_values(ghosted_solution, jumpu);
            fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu);
            lambda.value_list(fe_iv.get_quadrature_points(), lambda_values);
            mu.value_list(fe_iv.get_quadrature_points(), mu_values);
            const std::vector<double> &JxW = fe_iv.get_JxW_values();
            const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

            if (cell->material_id() == 1)
            {
              local_counter++;
              for (unsigned int point = 0; point < n_face_q_points; ++point)
              {
                Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
                Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);
                tangential = cross_product_2d(normals[point]);

                local_interface_traction[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
                local_interface_traction[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];

                local_interface_jump += jumpu[point] / n_face_q_points;
              }
            }
          }

    interface_jump = Utilities::MPI::sum(local_interface_jump, mpi_communicator);
    counter = Utilities::MPI::sum(local_counter, mpi_communicator);
    interface_jump = interface_jump / (double)counter;
    interface_stress = Utilities::MPI::sum(local_interface_traction, mpi_communicator);
  }

  // Only works in series I guess
  template <int dim>
  void TSL_prllel_old<dim>::interface_jump_and_traction_at_qpoints(const unsigned int cycle)
  {
    Tensor<1, dim> interface_stress;
    Tensor<1, dim> interface_jump;

    FEInterfaceValues<dim> fe_iv(fe,
                                 *damaged_face_quadrature,
                                 UpdateFlags(update_values |
                                             update_gradients |
                                             update_quadrature_points |
                                             update_normal_vectors |
                                             update_JxW_values));

    const unsigned int n_face_q_points = damaged_face_quadrature->size();

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
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->user_flag_set())
          {

            // I want to visit the face only once
            cell->face(face)->clear_user_flag();

            auto cellnb = cell->neighbor(face);
            int nface_number = cell->neighbor_face_no(face);

            fe_iv.reinit(cell, face, numbers::invalid_unsigned_int, cellnb, nface_number, numbers::invalid_unsigned_int);
            auto &q_points = fe_iv.get_quadrature_points();
            fe_iv[displacements].get_jump_in_function_values(ghosted_solution, jumpu);
            fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu);
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

              if (cycle == 0)
              {
                std::ofstream q_point("q_points/q_point_" + std::to_string(counter) + ".txt");
                if (q_point.is_open() == 0)
                {
                  pcout << "q_points folder doesn't exist \n";
                  exit(1);
                }
              }
              else
              {
                std::ofstream q_point;
                q_point.open("q_points/q_point_" + std::to_string(counter) + ".txt", std::ios_base::app);
                q_point << q_points[point][0] << "\t" << interface_jump[0] << "\t" << interface_jump[1] << "\t" << interface_stress[0] << "\t" << interface_stress[1] << "\n";
                q_point.close();
              }
            }
          }

    // If the geomtery is flipped you have a problem
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const unsigned int face : cell->face_indices())
          if ((cell->face(face)->at_boundary() == false) && (cell->material_id() != cell->neighbor(face)->material_id()))
            cell->face(face)->set_user_flag();
          else
            cell->face(face)->clear_user_flag();
  }

  // Probs Only works in series too
  template <int dim>
  void TSL_prllel_old<dim>::plot_over_line()
  {
    QGauss<dim> quadrature_formula(degree + 1);

    FEValues<dim> fe_v(fe, quadrature_formula,
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
      if (cell->is_locally_owned())
        if (std::abs(cell->center()[1] - 1.17) < 1e-2)
        {
          fe_v.reinit(cell);
          fe_v[displacements].get_function_gradients(ghosted_solution, gradu);
          lambda.value_list(fe_v.get_quadrature_points(), lambda_values);
          mu.value_list(fe_v.get_quadrature_points(), mu_values);
          const std::vector<double> &JxW = fe_v.get_JxW_values();
          Tensor<2, dim> avg_stress;
          for (unsigned int point = 0; point < n_q_points; ++point)
          {
            /////// * JxW
            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
            Tensor<2, dim> stress = (lambda_values[point] * trace(strain) * Identity + 2 * mu_values[point] * strain);
            avg_stress += stress;
          }

          pcout << cell->center()[0] << " \n"
                << avg_stress[0][0] << "\t" << avg_stress[0][1] << "\n"
                << avg_stress[1][0] << "\t" << avg_stress[1][1] << " \n \n";
        }
  }

  template <int dim>
  void TSL_prllel_old<dim>::run()
  {
    pcout << "Running with "
#ifdef USE_PETSC_LA
          << "PETSc"
#else
          << "Trilinos"
#endif
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    for (unsigned int refinement_level = 0; refinement_level < 1; refinement_level++)
    {
      triangulation.clear();
      // Mesh geenration
      Point<dim> P1, P2;
      std::vector<unsigned int> repetitions;
      repetitions.resize(2);

      if (par.geometry == "reactangle" || par.geometry == "sheer")
      {
        P2[0] = 1.0;
        P2[1] = 2.0;
        repetitions[0] = 1;
        repetitions[1] = 2;

        GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
        triangulation.refine_global(par.initial_refinement + refinement_level);

        Point<dim> p;
        for (const auto &cell : triangulation.active_cell_iterators())
        {
          p = cell->center();
          if (p[1] > 1.0)
            cell->set_material_id(1);
          else
            cell->set_material_id(2);
        }
        // GridTools::rotate(-M_PI/4, triangulation);
        //   GridTools::distort_random(0.3, triangulation, true,556);
      }

      if (par.geometry == "peeling")
      {
        P1[0] = 0.0;
        P1[1] = -0.5;
        P2[0] = 10.0;
        P2[1] = 2.0;
        repetitions[0] = 5;
        repetitions[1] = 6;

        GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
        triangulation.refine_global(par.initial_refinement + refinement_level);

        Point<dim> p;
        for (const auto &cell : triangulation.active_cell_iterators())
        {
          p = cell->center();
          if (p[1] > 1.0)
            cell->set_material_id(1);
          else
            cell->set_material_id(2);
        }
      }

      if (par.geometry == "delamitation")
      {
        P1[0] = 0.0;
        P1[1] = 0.0;
        P2[0] = 50.0;
        P2[1] = 2.0;
        repetitions[0] = 10;
        repetitions[1] = 4;

        GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
        triangulation.refine_global(par.initial_refinement + refinement_level);

        Point<dim> p;
        for (const auto &cell : triangulation.active_cell_iterators())
        {
          p = cell->center();
          if (p[1] > 1.0)
            cell->set_material_id(1);
          else
            cell->set_material_id(2);
        }
      }

      if (par.geometry == "signle-edge")
      {
        P2[0] = 100.0;
        P2[1] = 100.0;
        repetitions[0] = 1;
        repetitions[1] = 1;

        GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
        triangulation.refine_global(par.initial_refinement + refinement_level);

        Point<dim> p;
        for (const auto &cell : triangulation.active_cell_iterators())
        {
          p = cell->center();
          if (p[1] > 50.0)
            cell->set_material_id(1);
          else
            cell->set_material_id(2);
        }
      }

      if (par.geometry == "matrix")
      {
        GridIn<2> gridin;
        gridin.attach_triangulation(triangulation);
        std::ifstream f("matrix.msh");

        gridin.read_msh(f);
        triangulation.refine_global(par.initial_refinement);
        std::ofstream out("grid-1.vtu");
        GridOut grid_out;
        grid_out.write_vtu(triangulation, out);
        pcout << " written to "
              << "grid-1.vtu" << std::endl
              << std::endl;
      }

      // {

      //   for (unsigned int i = 0; i < 3; i++)
      //   {
      //     typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
      //                                                    endc = dof_handler.end();

      //     for (; cell != endc; ++cell)
      //     {
      //       // if ((cell->center()[1]) > 1.9 || (cell->center()[1]) < 0.1 )
      //       //   if ((cell->center()[0] > 0.9)  || (cell->center()[0] < 0.1))
      //       //     cell->set_refine_flag();

      //       // pcout << std::pow(2.0,-(float)(i))<<"\n" ;

      //       if (std::abs(cell->center()[1] - 1) < 0.1)
      //         cell->set_refine_flag();
      //     }

      //     triangulation.execute_coarsening_and_refinement();
      //   }
      // }

      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int face : cell->face_indices())
          if ((cell->face(face)->at_boundary() == false) && (cell->material_id() != cell->neighbor(face)->material_id()))
            cell->face(face)->set_user_flag();
          else
            cell->face(face)->clear_user_flag();

      setup_quadrature_point_history();
      double error;
      double error2;
      unsigned int non_lin = 0;

      std::ofstream forces("forces" + std::to_string(refinement_level) + ".txt");
      unsigned int max_nonlin = 0;
      unsigned int at_cycle = 0;
      for (unsigned int cycle = 0; cycle < par.cycles; ++cycle)
      {
        error = 1;
        error2 = 1;
        non_lin = 0;

        if (cycle < par.unloading || cycle > par.reloading)
        {
          disp[0] += par.displacementx;
          disp[1] += par.displacementy;
        }
        else
        {
          disp[0] = disp[0] - par.displacementx;
          disp[1] = disp[1] - par.displacementy;
        }

        pcout << " ####### cycle = " << cycle << " and displacement = " << disp[1] << " ###### \n";
        unsigned int max_nonlin_iter = 100;
        while (error > par.error && error2 > par.error && non_lin < max_nonlin_iter)
        {

          if (cycle == 0 && non_lin == 0)
            setup_system(true);
          else
            setup_system(false);

          max_nonlin = std::max(max_nonlin, non_lin);

          if (max_nonlin == non_lin)
            at_cycle = cycle;

          if (max_nonlin == max_nonlin_iter - 1)
          {
            pcout << "max non_lin iterations reached b4 convergence.  \n";
            exit(1);
          }

          non_lin++;

          assemble_system(cycle, non_lin, error);

          error = system_rhs.l2_norm();

          //  pcout << "  Residual: " << compute_residual(0) << std::endl;
          pcout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;

          if (error > par.error)
          {
            pcout << "  non lin = " << non_lin << std::endl;
            solve(cycle);
            pcout << "  update l2_norm " << distributed_newton_update.l2_norm() << std::endl;
          }
          else
          {
            assemble_system(cycle, non_lin, error);
            error2 = system_rhs.l2_norm();
          }
          pcout << "cycle " << cycle << "\n------- \n ";
          pcout << "max nonlin iterations " << max_nonlin << "  at_cycle " << at_cycle << "\n------- \n ";
        }

        if (cycle % par.output_frequency == 0)
        {
          TimerOutput::Scope t(computing_timer, "output");
          output_results(cycle, non_lin, refinement_level);
        }

        interface_jump_and_traction_at_qpoints(cycle);

        forces.open("forces" + std::to_string(refinement_level) + ".txt", std::ios_base::app);
        if (forces.is_open() == false)
        {
          pcout << "File didn't open\n";
          exit(1);
        }
        Tensor<1, dim> ten;
        Tensor<1, dim> jump;
        Tensor<1, dim> inter;

        reaction(ten, 3);
        interface_jump_and_traction(jump, inter);
        pcout << "   reaction force = " << ten[0] << "\t" << ten[1] << "\n";

        pcout << std::endl;

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          forces << disp[0] << "\t" << disp[1] << "\t"
                 << ten[0] << "\t" << ten[1] << "\t"
                 << jump[0] << "\t" << jump[1] << "\t"
                 << inter[0] << "\t" << inter[1] << "\n";

        forces.close();
      }

      // plot_over_line();

      computing_timer.print_summary();
      computing_timer.reset();
    }
  }

} // closing namespace

#endif
