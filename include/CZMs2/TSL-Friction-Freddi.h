#ifndef TSL_FRICTION_FREDDI_H
#define TSL_FRICTION_FREDDI_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
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
#include <deal.II/grid/grid_in.h>
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
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/base/parameter_acceptor.h>

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

template <int dim>
dealii::Point<dim> grid_y_transform(const dealii::Point<dim> &pt_in)
{
  const double &x = pt_in[0];
  const double &y = pt_in[1];

  const double slope = 0.2;

  dealii::Point<dim> pt_out = pt_in;
  if (std::abs(pt_in[1] - 1) < 0.99)
    pt_out[1] = pt_out[1] + slope * (pt_out[0] - 0.5);

  return pt_out;
}

double heaviside(double x)
{
  return (x < 0) ? 0.0 : 1.0;
}

// #define FORCE_USE_OF_TRILINOS
namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
} // namespace LA

namespace TSL_Friction_Freddi
{
  using namespace dealii;

  template <int dim>
  class ProblemParameters : public ParameterAcceptor
  {
  public:
    ProblemParameters() : ParameterAcceptor("main")
    {
      add_parameter("testing", testing, " ", this->prm, Patterns::Double());
      add_parameter("E", E, " ", this->prm, Patterns::Double());
      add_parameter("nu", nu, " ", this->prm, Patterns::Double());
      add_parameter("Ep", Ep, " ", this->prm, Patterns::Double());
      add_parameter("nup", nup, " ", this->prm, Patterns::Double());
      add_parameter("sig_c", sig_c, " ", this->prm, Patterns::Double());
      add_parameter("sig_ci", sig_ci, " ", this->prm, Patterns::Double());
      add_parameter("delta_c", delta_c, " ", this->prm, Patterns::Double());
      add_parameter("delta_ci", delta_ci, " ", this->prm, Patterns::Double());
      add_parameter("ku", ku, " ", this->prm, Patterns::Double());
      add_parameter("kui", kui, " ", this->prm, Patterns::Double());
      add_parameter("regularization", regularization, " ", this->prm, Patterns::Double());
      add_parameter("m", m, " ", this->prm, Patterns::Double());
      add_parameter("coff", coff, " ", this->prm, Patterns::Double());
      add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
      add_parameter("geometry", geometry, "", this->prm, Patterns::Selection("reactangle|shear|peeling||hole|half_matrix|full_matrix"));
      add_parameter("TSL_law", TSL_law, "", this->prm, Patterns::Selection("damage_freddi"));
      add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
      add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
      add_parameter("quadrature", quadrature, " ", this->prm, Patterns::Selection("trap|gauss|lobatto|gauss-reduced|lobatto-over"));
      add_parameter("symmetry", symmetry, " ", this->prm, Patterns::Integer());
      add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
      add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
      add_parameter("error", error, " ", this->prm, Patterns::Double());
      add_parameter("is_user_penalty", is_user_penalty, " ", this->prm, Patterns::Bool());
      add_parameter("is_TSL_everywhere", is_TSL_everywhere, " ", this->prm, Patterns::Bool());
      add_parameter("is_distorted", is_distorted, " ", this->prm, Patterns::Bool());
      add_parameter("is_implicit", is_implicit, " ", this->prm, Patterns::Bool());
      add_parameter("penalty", penalty, " ", this->prm, Patterns::Double());
      add_parameter("penalty_factor", penalty_factor, " ", this->prm, Patterns::Double());
      add_parameter("unloading", unloading, " ", this->prm, Patterns::Integer());
      add_parameter("reloading", reloading, " ", this->prm, Patterns::Integer());
      add_parameter("output_frequency", output_frequency, " ", this->prm, Patterns::Integer());
    }

    int cycles = 0, initial_refinement = 0, degree = 0, symmetry = 0, unloading = 0, reloading = 0, output_frequency = 0;
    double testing = 0, E = 0, nu = 0, Ep = 0, nup = 0, regularization = 0, m = 0, sig_c = 0, sig_ci = 0,
           delta_c = 0, delta_ci = 0, ku = 0, kui = 0, coff = 0,
           displacementx = 0, displacementy = 0, error = 0, penalty = 0, penalty_factor = 0;
    bool is_implicit = 0, is_user_penalty = 0, is_TSL_everywhere = 0, is_distorted = 0;
    std::string quadrature = "trap", geometry = "reactangle", TSL_law = "damage_freddi";
  };

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

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(),
                        (Tensor<2, dim>::n_independent_components));

        const typename DoFHandler<dim>::cell_iterator current_cell =
            input_data.template get_cell<dim>();

        double elasticity = 0, poisson = 0;

        if (current_cell->material_id() == 1)
        {
          elasticity = E;
          poisson = nu;
        }
        else
        {
          elasticity = Ep;
          poisson = nup;
        }

        double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
        double mu = elasticity / (2 * (1 + poisson));

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
    }
  };

  // # clean here.
  // # maybe have history at nonlin==1 and at every _ non_lin ietration
  template <int dim>
  struct PointHistory
  {
    bool is_damaged;
    double max_seperation;
    Tensor<1, dim> jump;
    double seperation;
    double damage;
    bool is_reunload;

    Tensor<1, dim> output_cords;
    Tensor<1, dim> output_normals;
    Tensor<1, dim> output_jump;
    Tensor<1, dim> output_numerical_traction;
    Tensor<1, dim> output_analytical_traction;
    Tensor<2, dim> output_stress;

    double Freddi_g;

    Tensor<1, dim> traction_init;
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
  class TSL_Friction_Freddi
  {
  public:
    TSL_Friction_Freddi(const unsigned int degree, const ProblemParameters<dim> &par);
    void run();

  private:
    void setup_system(const bool initial_step);
    void assemble_system(const unsigned int cycle, const unsigned int non_lin, const double error);
    void solve(const unsigned int cycle);
    void output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const;
    void q_point_PPV(const unsigned int cycle, const unsigned int non_lin);

    void interface_jump_and_traction_at_qpoints(const unsigned int cycle);
    void reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id);
    void interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress);
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FESystem<dim> fe;

    const QGauss<dim> quadrature;
    std::unique_ptr<const Quadrature<dim - 1>> face_quadrature;
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
  TSL_Friction_Freddi<dim>::TSL_Friction_Freddi(const unsigned int degree, const ProblemParameters<dim> &par)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        triangulation(mpi_communicator),
        degree(degree),
        quadrature(degree + 1),
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

    if (par.quadrature == "gauss")
      face_quadrature = std::make_unique<const QGauss<dim - 1>>(degree + 1);
    if (par.quadrature == "lobatto")
      face_quadrature = std::make_unique<const QGaussLobatto<dim - 1>>(degree + 1);
    if (par.quadrature == "gauss-reduced")
      face_quadrature = std::make_unique<const QGauss<dim - 1>>(degree);
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::setup_quadrature_point_history()
  {

    triangulation.clear_user_data();

    {
      std::vector<PointHistory<dim>> tmp;
      quadrature_point_history.swap(tmp);
    }

    quadrature_point_history.resize(triangulation.n_active_faces() * face_quadrature->size());

    unsigned int history_index = 0;
    for (auto &face : triangulation.active_face_iterators())
    {
      face->set_user_pointer(&quadrature_point_history[history_index]);
      history_index += face_quadrature->size();
    }

    Assert(history_index == quadrature_point_history.size(), ExcInternalError());
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::setup_system(const bool initial_step)
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
  void TSL_Friction_Freddi<dim>::assemble_system(const unsigned int cycle, const unsigned int non_lin, const double error)
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
      double elasticity = 0, poisson = 0;

      if (cell->material_id() == 1)
      {
        elasticity = par.E;
        poisson = par.nu;
      }
      else
      {
        elasticity = par.Ep;
        poisson = par.nup;
      }

      double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
      double mu = elasticity / (2 * (1 + poisson));

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<2, dim>> old_solution_gradients(n_q_points);
      fe_v[displacements].get_function_gradients(ghosted_solution, old_solution_gradients);

      for (unsigned int point = 0; point < n_q_points; ++point)
        for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        {
          Tensor<2, dim> strain_i = 0.5 * (fe_v[displacements].gradient(i, point) + transpose(fe_v[displacements].gradient(i, point)));

          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
          {
            Tensor<2, dim> strain_j = 0.5 * (fe_v[displacements].gradient(j, point) + transpose(fe_v[displacements].gradient(j, point)));
            Tensor<2, dim> stressj = lambda * trace(strain_j) * Identity + 2 * mu * strain_j;

            copy_data.cell_matrix(i, j) += scalar_product(strain_i, stressj) * JxW[point];
          }

          Tensor<2, dim> old_strain = 0.5 * (old_solution_gradients[point] + transpose(old_solution_gradients[point]));
          Tensor<2, dim> old_stress = lambda * trace(old_strain) * Identity + 2 * mu * old_strain;
          copy_data.cell_rhs(i) += -scalar_product(strain_i, old_stress) * JxW[point];
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

      double elasticity = 0, poisson = 0;

      if (cell->material_id() == 1)
      {
        elasticity = par.E;
        poisson = par.nu;
      }
      else
      {
        elasticity = par.Ep;
        poisson = par.nup;
      }

      double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
      double mu = elasticity / (2 * (1 + poisson));

      // ------- restricing x and y on boundaries 2 adn 3 -------

      bool is_at_corner = false;
      if (cell->face(face_no)->boundary_id() == 3 || cell->face(face_no)->boundary_id() == 2)
        for (int i = 0; i < 4; i < i++)
          if (cell->face(i)->boundary_id() == 1)
            is_at_corner = true;

      double penalty = (par.penalty_factor) * elasticity / (cell->measure());
      if (par.is_user_penalty)
        penalty = par.penalty;

      const FEValuesExtractors::Scalar displacements_x(0);
      const FEValuesExtractors::Scalar displacements_y(1);

      if (par.geometry == "reactangle" || par.geometry == "hole" || par.geometry == "full_matrix")
      {
        if ((cell->face(face_no)->boundary_id() == 2) || cell->face(face_no)->boundary_id() == 3)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> strain_i = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strain_j = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements].value(i, point) *
                         (lambda * trace(strain_j) * Identity + 2 * mu * strain_j) * normals[point]

                     + par.symmetry *
                           (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                           fe_fv[displacements].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements].value(i, point) * fe_fv[displacements].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements].value(i, point) *
                          (lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]

                      + par.symmetry *
                            (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                            old_solution_values[point] // Symetry term

                      + penalty * fe_fv[displacements].value(i, point) * old_solution_values[point]) *
                      JxW[point]

                  + par.symmetry * (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                        ((cell->face(face_no)->boundary_id() - 2) * disp) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements].value(i, point) * ((cell->face(face_no)->boundary_id() - 2) * disp) *
                        JxW[point]; // dx
            }
          }
      }

      if (par.geometry == "shear")
      {

        if ((cell->face(face_no)->boundary_id() == 0) && cell->material_id() != 1)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_x].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[0]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                           fe_fv[displacements_x].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_x].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[0]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                            old_solution_values[point][0] // Symetry term

                      + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                        (-disp[0]) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_x].value(i, point) * (-disp[0]) *
                        JxW[point]; // dx
            }
          }

        if ((cell->face(face_no)->boundary_id() == 1) && cell->material_id() == 1)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_x].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[0]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                           fe_fv[displacements_x].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_x].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[0]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                            old_solution_values[point][0] // Symetry term

                      + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                        (disp[0]) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_x].value(i, point) * (disp[0]) *
                        JxW[point]; // dx
            }
          }

        if (cell->face(face_no)->boundary_id() == 2)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_y].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[1]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                           fe_fv[displacements_y].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_y].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[1]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                            old_solution_values[point][1] // Symetry term

                      + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                        -disp[1] * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_y].value(i, point) * -disp[1] *
                        JxW[point]; // dx
            }
          }

        if (cell->face(face_no)->boundary_id() == 3)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_y].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[1]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                           fe_fv[displacements_y].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_y].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[1]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                            old_solution_values[point][1] // Symetry term

                      + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                        disp[1] * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_y].value(i, point) * disp[1] *
                        JxW[point]; // dx
            }
          }
      }

      if (par.geometry == "peeling")
      {
        if ((cell->face(face_no)->boundary_id() == 2))
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> strain_i = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strain_j = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements].value(i, point) *
                         (lambda * trace(strain_j) * Identity + 2 * mu * strain_j) * normals[point]

                     + par.symmetry *
                           (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                           fe_fv[displacements].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements].value(i, point) * fe_fv[displacements].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements].value(i, point) *
                          (lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]

                      + par.symmetry *
                            (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                            old_solution_values[point] // Symetry term

                      + penalty * fe_fv[displacements].value(i, point) * old_solution_values[point]) *
                      JxW[point]

                  + par.symmetry * (lambda * trace(strain_i) * Identity + 2 * mu * strain_i) * normals[point] *
                        (0 * disp) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements].value(i, point) * (0) * disp *
                        JxW[point]; // dx
            }
          }

        if (cell->face(face_no)->boundary_id() == 1 && cell->material_id() == 1)
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_y].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[1]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                           fe_fv[displacements_y].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_y].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[1]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                            old_solution_values[point][1] // Symetry term

                      + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                        disp[1] * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_y].value(i, point) * disp[1] *
                        JxW[point]; // dx
            }
          }
      }

      if (par.geometry == "half_matrix")
      {

        if ((cell->face(face_no)->boundary_id() == 1))
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_x].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[0]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                           fe_fv[displacements_x].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_x].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[0]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                            old_solution_values[point][0] // Symetry term

                      + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                        (disp[0]) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_x].value(i, point) * (disp[0]) *
                        JxW[point]; // dx
            }
          }

        if ((cell->face(face_no)->boundary_id() == 0))
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_x].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[0]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                           fe_fv[displacements_x].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_x].value(i, point) * fe_fv[displacements_x].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_x].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[0]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                            old_solution_values[point][0] // Symetry term

                      + penalty * fe_fv[displacements_x].value(i, point) * old_solution_values[point][0]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[0] *
                        (0.0) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_x].value(i, point) * (0.0) *
                        JxW[point]; // dx
            }
          }

        if ((cell->face(face_no)->boundary_id() == 2))
          for (unsigned int point = 0; point < n_q_points; ++point)
          {

            Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));

            for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              Tensor<2, dim> straini = 0.5 * (fe_fv[displacements].gradient(i, point) + transpose(fe_fv[displacements].gradient(i, point)));

              for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                Tensor<2, dim> strainj = 0.5 * (fe_fv[displacements].gradient(j, point) + transpose(fe_fv[displacements].gradient(j, point)));

                copy_data.cell_matrix(i, j) +=
                    (-fe_fv[displacements_y].value(i, point) *
                         ((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point])[1]

                     + par.symmetry *
                           ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                           fe_fv[displacements_y].value(j, point) // Symetry term

                     + penalty * fe_fv[displacements_y].value(i, point) * fe_fv[displacements_y].value(j, point)) *
                    JxW[point];
              }

              copy_data.cell_rhs(i) +=

                  -(
                      -fe_fv[displacements_y].value(i, point) *
                          ((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point])[1]

                      + par.symmetry *
                            ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                            old_solution_values[point][1] // Symetry term

                      + penalty * fe_fv[displacements_y].value(i, point) * old_solution_values[point][1]) *
                      JxW[point]

                  + par.symmetry * ((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point])[1] *
                        (0.0) * JxW[point] // Symetry term

                  + penalty *
                        fe_fv[displacements_y].value(i, point) * (0.0) *
                        JxW[point]; // dx
            }
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
      scratch_data.fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
      const FEFaceValuesBase<dim> &fe_fv1 = scratch_data.fe_interface_values.get_fe_face_values(0);
      const FEFaceValuesBase<dim> &fe_fv2 = scratch_data.fe_interface_values.get_fe_face_values(1);

      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      // fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
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

      double elasticity = 0, poisson = 0;
      double nelasticity = 0, npoisson = 0;

      if (cell->material_id() == 1)
      {
        elasticity = par.E;
        poisson = par.nu;
      }
      else
      {
        elasticity = par.Ep;
        poisson = par.nup;
      }

      if (ncell->material_id() == 1)
      {
        nelasticity = par.E;
        npoisson = par.nu;
      }
      else
      {
        nelasticity = par.Ep;
        npoisson = par.nup;
      }

      double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
      double mu = elasticity / (2 * (1 + poisson));
      double nlambda = nelasticity * npoisson / ((1 + npoisson) * (1 - 2 * npoisson));
      double nmu = nelasticity / (2 * (1 + npoisson));

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1, dim>> old_solution_jumps(n_q_points);
      fe_iv[displacements].get_jump_in_function_values(ghosted_solution, old_solution_jumps);

      std::vector<Tensor<2, dim>> gradu_avg(n_q_points);
      fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu_avg);

      PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());

      Tensor<1, dim> g, g_eff;
      Tensor<1, dim> q, p;
      Tensor<1, dim> tr, tr_eff;
      Tensor<1, dim> sig_u, sig_trial;
      Tensor<1, dim> tangential;
      Tensor<2, dim> TCZ;
      Tensor<1, dim> TCZ_res;
      double seperation;

      Tensor<2, dim> Rot; // rotation

      for (unsigned int point = 0; point < n_q_points; ++point)
      {

        double sig_c = 0, delta_c = 0, ku = 0;
        if (cell->material_id() == ncell->material_id())
        {
          sig_c = par.sig_c;
          ku = par.ku;
          delta_c = par.delta_c;
        }
        else
        {
          sig_c = par.sig_ci;
          ku = par.kui;
          delta_c = par.delta_ci;
        }

        double penalty = (par.penalty_factor) * ((elasticity < nelasticity) ? nelasticity : elasticity) / (cell->measure());
        if (par.is_user_penalty)
          penalty = par.penalty;

        Tensor<2, dim> strain = 0.5 * (gradu_avg[point] + transpose(gradu_avg[point]));
        Tensor<2, dim> old_stress_avg = (lambda * trace(strain) * Identity + 2 * mu * strain);

        const double normaly = normals[point][1];
        const int sign_normalx = ((normals[point][0] > 0) ? 1 : ((normals[point][0] < 0) ? -1 : 0));
        Rot[0][0] = normaly;
        Rot[0][1] = -sign_normalx * std::sqrt(1 - normaly * normaly);
        Rot[1][0] = sign_normalx * std::sqrt(1 - normaly * normaly);
        Rot[1][1] = normaly;
        Rot = -Rot;

        tangential = cross_product_2d(normals[point]);
        tr[0] = tangential * (old_stress_avg * normals[point]);
        tr[1] = normals[point] * (old_stress_avg * normals[point]);
        tr_eff[0] = tr[0];
        tr_eff[1] = std::abs(tr[1]) / 2 + tr[1] / 2;

        g = Rot * old_solution_jumps[point];
        // p[0] = g[0];
        q[1] = heaviside(g[1]) * g[1];
        g_eff[0] = g[0];
        g_eff[1] = std::abs(g[1]) / 2 + g[1] / 2;
        seperation = g.norm();
        if (par.is_implicit == true)
        {
          quadrature_points_history[point].max_seperation = std::max(quadrature_points_history[point].seperation, quadrature_points_history[point].max_seperation);
          quadrature_points_history[point].seperation = seperation;
        }

        Tensor<2, dim> K;
        K[0][0] = sig_c / delta_c;
        K[1][1] = sig_c / delta_c;

        sig_u = K * g;
        double damage;
        if (par.is_implicit == true)
          damage = ku * (sig_c - sig_u.norm()) / (sig_u.norm() * (sig_c - ku));
        else
          damage = ku * (sig_c - sig_u.norm()) / ((sig_u.norm()) * (sig_c - ku));

        if (par.is_implicit == true)
        {
          damage = std::max(0.0, std::min(1.0, damage));
          quadrature_points_history[point].damage = std::max(quadrature_points_history[point].damage, damage);

          if (g[1] > 0)
            p[0] = g[0];
          else
            p[0] = quadrature_points_history[point].Freddi_g;
        }
        sig_trial = K * (g - (p + q));

        double coff = par.coff;
        double Coulomb = std::abs(sig_trial[0]) + coff * sig_trial[1];

        if (par.is_implicit == true)
        {
          if (g[1] < 0)
          {
            if (Coulomb <= 0)
            {
              // do nothing to p
            }
            else
            {
              quadrature_points_history[point].Freddi_g += (1 / K[0][0]) * (coff * sig_trial[1] + std::abs(sig_trial[0])) * (std::abs(sig_trial[0]) / sig_trial[0]);
              p[0] = quadrature_points_history[point].Freddi_g;
            }
          }
        }

        sig_trial = K * (g - (p + q));
        // double Coulomb2 = std::abs(sig_trial[0]) + coff * sig_trial[1];
        // if (cell->face(f)->user_flag_set() == true)
        //   if (q_points[point][0] < 0.001)
        //     pcout << "columb aft " << Coulomb2 << "\n \n";


    
        //      if(non_lin==1 )
        if (par.is_implicit == true)
        {
          if (sig_u.norm() > sig_c)
            quadrature_points_history[point].is_damaged = true;
        }
        else
        {
          if (tr_eff.norm() > sig_c)
          {
            quadrature_points_history[point].is_damaged = true;
          }
        }

        if (quadrature_points_history[point].is_damaged == false)
          quadrature_points_history[point].traction_init = tr_eff;

     
        //   if(non_lin==1 )
        if (par.is_implicit == true)
        {
          if ((quadrature_points_history[point].damage - damage) > 1e-8)
          {
            quadrature_points_history[point].is_reunload = true;
          }
          else
            quadrature_points_history[point].is_reunload = false;
        }
        if (par.is_implicit == false)
        {
          // if ((quadrature_points_history[point].damage - damage) > 1e-8)
          // {
          //   quadrature_points_history[point].is_reunload = true;
          // }
          // else
          //   quadrature_points_history[point].is_reunload = false;
        }

        double regularization = par.regularization;

        if (par.TSL_law == "damage_freddi")
        {
          if (par.is_implicit == true)
          {
            if (cell->face(f)->user_flag_set() == true)
            // if (q_points[point][0] < 0.001)
            {
              pcout << "damage " << quadrature_points_history[point].damage << " \n";
              pcout << "seperation " << seperation << "\n";
              pcout << "Freddi_g " << quadrature_points_history[point].Freddi_g << "\n";
              pcout << "columb b4 " << Coulomb << "\n";
            }

             
            if (quadrature_points_history[point].is_reunload == true)
              damage = quadrature_points_history[point].damage;

            if (damage > 0 + 1e-8)
            {

              TCZ_res = K * (g - damage * (p + q));

              // dTCZ / du
              TCZ = K;

              // dTCZ /ddamage x ddamage/du
              if (damage < 1 - 1e-8 && quadrature_points_history[point].is_reunload == false)
              {
                Tensor<1, dim> ddamage_du;
                ddamage_du[0] = -K[0][0] * K[0][0] * g[0];
                ddamage_du[1] = -K[1][1] * K[1][1] * g[1];
                ddamage_du *= (1 / sig_u.norm() + 1 / (sig_c - sig_u.norm())) * damage / sig_u.norm();
                TCZ += outer_product(-K * (p + q), ddamage_du);
              }

              // dTCZ/dp x dp/du
              if (g[1] > 0)
              {
                Tensor<2, dim> dp_du;
                dp_du[0][0] = 1.0;
                TCZ += -K * damage * dp_du;
              }
              else
              {
                if (Coulomb < 0)
                {
                  // do nothing
                }
                else
                {
                  Tensor<2, dim> dp_du;
                  dp_du[0][0] = 1.0;
                  dp_du[0][1] = coff * (K[1][1] / K[0][0]) * std::abs(g[0] - quadrature_points_history[point].Freddi_g) / (g[0] - quadrature_points_history[point].Freddi_g);
                  TCZ += -K * damage * dp_du;
                }
              }

              //  dTCZ/dq * dq/du
              Tensor<2, dim> dq_du;
              dq_du[1][1] = heaviside(g[1]);
              TCZ += -K * damage * dq_du;
            }
            else
            {
              TCZ_res = K * g;
              TCZ = K;
            }
          }
          else // expliict
          {

            TCZ_res = 0;
            Tensor<1, dim> shift = invert(K) * quadrature_points_history[point].traction_init;

            g += shift;
            Tensor<1,dim> g2=g; 
            //g[1] = heaviside(g[1]) * g[1];
            sig_u = K * (g);
            damage = ku * (sig_c - sig_u.norm()) / ((sig_u.norm()) * (sig_c - ku));
            damage = std::max(0.0, std::min(1.0, damage));

            if (non_lin == 1)
              quadrature_points_history[point].damage = std::max(damage, quadrature_points_history[point].damage);

            if ((quadrature_points_history[point].damage - damage) > 1e-8)
              quadrature_points_history[point].is_reunload = true;
            else
              quadrature_points_history[point].is_reunload = false;

            if (quadrature_points_history[point].is_reunload == false)
            {
              quadrature_points_history[point].max_seperation = std::max(g.norm(), quadrature_points_history[point].max_seperation);
              if (non_lin == 1)
                quadrature_points_history[point].jump = g;
            }
            
            p[0] = g[0];
            q[1] = heaviside(g[1]) * g[1];
            if (g[1] < 0)
            {
              p[0] = quadrature_points_history[point].Freddi_g;
            }

            g2[0]=0;
            g2[1]= heaviside(-g[1]) * g[1];
            sig_trial = K * (g - (p + q)) + penalty*g2;

            double coff = par.coff;
            double Coulomb = std::abs(sig_trial[0]) + coff * sig_trial[1];
           
            if (g[1] < 0)
            {
              if (Coulomb <= 0)
              {
                // do nothing to p
                
              }
              else
              {

                quadrature_points_history[point].Freddi_g += (1 / K[0][0]) * (coff * sig_trial[1] + std::abs(sig_trial[0])) * (std::abs(sig_trial[0]) / sig_trial[0]);
                p[0] = quadrature_points_history[point].Freddi_g;
              }
            }

            TCZ_res = K * (g - damage * (p + q));

            TCZ = 0;
            TCZ = K;

            // dTCZ /ddamage x ddamage/du
            if (damage < 1 - 1e-10 && quadrature_points_history[point].is_reunload == false)
            {
              Tensor<1, dim> ddamage_du;
              ddamage_du[0] = -K[0][0] * K[0][0] * g[0];
              ddamage_du[1] = -K[1][1] * K[1][1] * g[1];
              ddamage_du *= (1 / sig_u.norm() + 1 / (sig_c - sig_u.norm())) * damage / sig_u.norm();
              TCZ += outer_product(-K * (p + q), ddamage_du);
            }

            // dTCZ/dp x dp/du
            if (g[1] > 0)
            {
              Tensor<2, dim> dp_du;
              dp_du[0][0] = 1.0;
              TCZ += -K * damage * dp_du;
            }
            else
            {
              if (Coulomb < 0)
              { 
                

                // do nothing
              }
              else
              {
                Tensor<2, dim> dp_du;
                dp_du[0][0] = 1.0;
                dp_du[0][1] = coff * (K[1][1] / K[0][0]) * std::abs(g[0] - quadrature_points_history[point].Freddi_g) / (g[0] - quadrature_points_history[point].Freddi_g);
                TCZ += -K * damage * dp_du;
              }
            }

            //  dTCZ/dq * dq/du
            Tensor<2, dim> dq_du;
            dq_du[1][1] = heaviside(g[1]);
            TCZ += -K * damage * dq_du;

            if (quadrature_points_history[point].is_reunload)
            {

              damage = quadrature_points_history[point].damage;
              g -= shift;

              TCZ_res = K * (1 - damage) * g;
              TCZ_res[0] *= 1 + shift[0] / (quadrature_points_history[point].jump[0] - shift[0]);
              TCZ_res[1] *= 1 + shift[1] / (quadrature_points_history[point].jump[1] - shift[1]);

              TCZ = K * (1 - damage);
              TCZ[0][0] *= 1 + shift[0] / (quadrature_points_history[point].jump[0] - shift[0]);
              TCZ[1][1] *= 1 + shift[1] / (quadrature_points_history[point].jump[1] - shift[1]);
            }

            if (g[1] < 0)
            {
              TCZ_res[1] = 0;
              TCZ[1][1] = 0;
              TCZ[0][1] = 0;
              TCZ[1][0] = 0;

              TCZ_res[1] += penalty * g[1];
              TCZ[1][1] += penalty;
            }

            if (cell->face(f)->user_flag_set() == true)
            {
              pcout << "is_damaged " << quadrature_points_history[point].is_damaged << "\n";
              pcout << "K " << K << "\n";
              pcout << "shift " << shift << "\n";
              pcout << "g " << g << "\n";
              pcout << "p " << p << "\n";
              pcout << "q " << q << "\n";
              pcout << "damage " << damage << "\n";
              pcout << "sig_trial " << sig_trial << "\n";
              pcout << "Coulomb " << Coulomb << "\n";
              pcout << "TCZ_res " << TCZ_res << "\n \n";
            }
          }

        } // closing damage freddi

        quadrature_points_history[point].output_cords = q_points[point];
        quadrature_points_history[point].output_normals = normals[point];
        // quadrature_points_history[point].output_jump = invert(Rot) * g_eff;
        // quadrature_points_history[point].output_numerical_traction = tr_eff;
        quadrature_points_history[point].output_analytical_traction = invert(Rot) * TCZ_res;
        quadrature_points_history[point].output_stress = old_stress_avg;

        // if (cell->face(f)->user_flag_set() == true)
        // // if (q_points[point][0] > 0.3 && q_points[point][0] < 0.4)
        // {
        //   pcout << "Is damaged? " << quadrature_points_history[point].is_damaged << "\n";
        //   pcout << "at q point  " << q_points[point][0] << "\n";
        //   pcout << "Criterion  " << disp[1] / lambda_f << "\n";
        //   pcout << "normals\t" << normals[point][0] << "\t" << normals[point][1] << "\n";
        //   pcout << "Rotation " << Rot << "\n";
        //   pcout << "jump\t" << g[0] << "\t" << g[1] << "\n"
        //         << "separation " << seperation << " max seperation  " << quadrature_points_history[point].max_seperation << "\n";
        //   pcout << "old_solution_jumps\t" << old_solution_jumps[point][0] << "\t" << old_solution_jumps[point][1] << "\n";
        //   pcout << "T_eff\t" << tr_eff[0] << "\t" << tr_eff[1] << "\n"
        //         << "T_eff_norm " << tr_eff.norm() << "\n";
        //   pcout << "TCZ\t" << TCZ[0][0] << "\t" << TCZ[0][1] << "\t" << TCZ[1][0] << "\t" << TCZ[1][1] << "\n";
        //   pcout << "TCZ_glob\t" << (invert(Rot) * TCZ * Rot) << "\n";
        //   pcout << "TCZ_res\t" << TCZ_res[0] << "\t" << TCZ_res[1] << " \n";
        //   pcout << "TCZ_res glob \t" << invert(Rot) * TCZ_res << "\n";
        //   pcout << "old_stress_avg\t" << old_stress_avg << "\n \n";
        // }

        if (par.is_implicit == true)
          quadrature_points_history[point].is_damaged = true;

        if (par.is_TSL_everywhere == false)
          if (cell->face(f)->user_flag_set() == false)
            quadrature_points_history[point].is_damaged = false;

        if (par.geometry == "full_matrix")
          if (cell->material_id() != 1 && ncell->material_id() != 1)
            quadrature_points_history[point].is_damaged = false;

        // Protect against over damage for reactangle
        // if (par.geometry == "reactangle" && par.is_TSL_everywhere == true)
        //   if (std::abs((cell->face(f)->center()[1] - 1.0)) > 0.5)
        //     quadrature_points_history[point].is_damaged = false;

        for (unsigned int i = 0; i < n_dofs_face; ++i)
        {

          Tensor<2, dim> straini = 0.5 * (fe_iv[displacements].average_of_gradients(i, point) + transpose(fe_iv[displacements].average_of_gradients(i, point)));
          Tensor<2, dim> stress_i = lambda * trace(straini) * Identity + 2 * mu * straini;

          for (unsigned int j = 0; j < n_dofs_face; ++j)
          {

            Tensor<2, dim> strainj = 0.5 * (fe_iv[displacements].average_of_gradients(j, point) + transpose(fe_iv[displacements].average_of_gradients(j, point)));
            Tensor<2, dim> stress_j = lambda * trace(strainj) * Identity + 2 * mu * strainj;

            if ((quadrature_points_history[point].is_damaged))
              copy_data_face.cell_matrix(i, j) += fe_iv[displacements].jump_in_values(i, point) * (invert(Rot) * TCZ * Rot) * fe_iv[displacements].jump_in_values(j, point) * JxW[point];
            else
              copy_data_face.cell_matrix(i, j) +=
                  (-fe_iv[displacements].jump_in_values(i, point) *
                       (stress_j)*normals[point]

                   + par.symmetry *
                         ((stress_i)) * normals[point] *
                         fe_iv[displacements].jump_in_values(j, point) // Symetry term

                   + penalty * fe_iv[displacements].jump_in_values(i, point) * fe_iv[displacements].jump_in_values(j, point)) *
                  JxW[point]; // dx
          }
          if (quadrature_points_history[point].is_damaged)
            copy_data_face.cell_rhs(i) +=
                -fe_iv[displacements].jump_in_values(i, point) * (invert(Rot) * TCZ_res) * JxW[point];
          else
            copy_data_face.cell_rhs(i) +=
                -(
                    -fe_iv[displacements].jump_in_values(i, point) *
                        old_stress_avg * normals[point]

                    + par.symmetry * stress_i * normals[point] *
                          old_solution_jumps[point]

                    + penalty * fe_iv[displacements].jump_in_values(i, point) * old_solution_jumps[point]) *
                JxW[point];
        }
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

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);

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

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::solve(const unsigned int cycle)
  {

    TimerOutput::Scope t(computing_timer, "solve");
    pcout << "solving" << std::endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.solve(system_matrix, distributed_newton_update, system_rhs);
    ghosted_newton_update = distributed_newton_update;

    double alpha = 1.0;

    distributed_solution.add(alpha, distributed_newton_update);
    ghosted_solution = distributed_solution;
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const
  {
    if (cycle != 0)
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
      Vector<double> subdomain(triangulation.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      data_out.build_patches();

      data_out.write_vtu_with_pvtu_record("./output/", "solution", cycle, mpi_communicator, 3, 0);
      data_out.write_vtu_with_pvtu_record("./output/", "solution_non_lin" + std::to_string(cycle), non_lin, mpi_communicator, 3, 0);
    }
  }

  // Note: degree+1 works but we miss outputs for >degree+1
  template <int dim>
  void TSL_Friction_Freddi<dim>::reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id)
  {
    TimerOutput::Scope t(computing_timer, "other");
    reaction_stress = 0.0;
    QGauss<dim - 1> face_quadrature_formula(degree + 1);

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
    double elasticity = 0, poisson = 0;
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

            if (cell->material_id() == 1)
            {
              elasticity = par.E;
              poisson = par.nu;
            }
            else
            {
              elasticity = par.Ep;
              poisson = par.nup;
            }

            double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
            double mu = elasticity / (2 * (1 + poisson));
            const FEValuesExtractors::Vector displacements(0);

            fe_face_values.reinit(cell, face);
            fe_face_values[displacements].get_function_gradients(ghosted_solution, gradu);

            const std::vector<double> &JxW = fe_face_values.get_JxW_values();
            const std::vector<Tensor<1, dim>> &normals = fe_face_values.get_normal_vectors();

            for (unsigned int point = 0; point < n_face_q_points; ++point)
            {
              Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
              Tensor<2, dim> stress = (lambda * trace(strain) * Identity + 2 * mu * strain);

              tangential = cross_product_2d(normals[point]);
              local_reaction[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
              local_reaction[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];
            }
          }

    reaction_stress = Utilities::MPI::sum(local_reaction, mpi_communicator);
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
  {
    TimerOutput::Scope t(computing_timer, "other");

    interface_stress = 0;
    interface_jump = 0;
    double face_counter = 0;
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][0] = 0;
    Identity[0][1] = 0;
    Identity[1][1] = 1;
    unsigned int counter = 0;

    const auto face_worker = [&](const auto &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim> &scratch_data,
                                 CopyData &copy_data)
    {
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

      const FEFaceValuesBase<dim> &fe_fv1 = fe_iv.get_fe_face_values(0);
      const FEFaceValuesBase<dim> &fe_fv2 = fe_iv.get_fe_face_values(1);

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
      const std::vector<double> &JxW = fe_iv.get_JxW_values();
      double elasticity = 0, poisson = 0;
      double nelasticity = 0, npoisson = 0;

      if (cell->material_id() == 1)
      {
        elasticity = par.E;
        poisson = par.nu;
      }
      else
      {
        elasticity = par.Ep;
        poisson = par.nup;
      }

      if (ncell->material_id() == 1)
      {
        nelasticity = par.E;
        npoisson = par.nu;
      }
      else
      {
        nelasticity = par.Ep;
        npoisson = par.nup;
      }

      double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
      double mu = elasticity / (2 * (1 + poisson));
      double nlambda = nelasticity * npoisson / ((1 + npoisson) * (1 - 2 * npoisson));
      double nmu = nelasticity / (2 * (1 + npoisson));

      const FEValuesExtractors::Vector displacements(0);

      std::vector<Tensor<1, dim>> jumpu(n_q_points);
      fe_iv[displacements].get_jump_in_function_values(ghosted_solution, jumpu);

      std::vector<Tensor<2, dim>> grads_1(n_q_points);
      fe_fv1[displacements].get_function_gradients(ghosted_solution, grads_1);
      std::vector<Tensor<2, dim>> grads_2(n_q_points);
      fe_fv2[displacements].get_function_gradients(ghosted_solution, grads_2);

      if (cell->face(f)->user_flag_set())
        for (unsigned int point = 0; point < n_q_points; ++point)
        {
          counter++;

          Tensor<2, dim> strain1 = 0.5 * (grads_1[point] + transpose(grads_1[point]));
          Tensor<2, dim> strain2 = 0.5 * (grads_2[point] + transpose(grads_2[point]));
          Tensor<2, dim> stress1 = (lambda * trace(strain1) * Identity + 2 * mu * strain1);
          Tensor<2, dim> stress2 = (nlambda * trace(strain2) * Identity + 2 * nmu * strain2);
          Tensor<2, dim> stress = stress1 / 2 + stress2 / 2;
          Tensor<1, dim> tangential = cross_product_2d(normals[point]);

          interface_stress[0] += scalar_product(stress, outer_product(normals[point], tangential)) * JxW[point];
          interface_stress[1] += scalar_product(stress, outer_product(normals[point], normals[point])) * JxW[point];
          interface_jump[0] += jumpu[point] * normals[point];
          interface_jump[1] += jumpu[point] * tangential;
        }
    };

    const auto cell_worker =
        [&](const auto &cell, auto &scratch_data, auto &copy_data) {};

    const auto boundary_worker = [&](const auto &cell,
                                     const unsigned int &face_no,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data) {};
    const auto copier = [&](const auto &c) {};

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
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

    counter = Utilities::MPI::sum(counter, mpi_communicator);
    interface_jump = Utilities::MPI::sum(interface_jump, mpi_communicator);
    interface_stress = Utilities::MPI::sum(interface_stress, mpi_communicator);
    interface_jump = interface_jump / counter;
  }

  template <int dim>
  void TSL_Friction_Freddi<dim>::interface_jump_and_traction_at_qpoints(const unsigned int cycle)
  {
    TimerOutput::Scope t(computing_timer, "other");

    const auto face_worker = [&](const auto &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim> &scratch_data,
                                 CopyData &copy_data)
    {
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());
      if (cell->face(f)->user_flag_set())
        for (unsigned int point = 0; point < n_q_points; ++point)
        {
          std::ofstream q_point;
          if (cycle == 0)
          {
            q_point.open("q_points/q_point_" + std::to_string(quadrature_points_history[point].output_cords[0]) + "_" + std::to_string(point) + ".txt");
            q_point << "disp\t2=cords_x\tcords_y\t4=g_eff_x\tg_eff_y\tg_eff_norm\t7=tr_x\ttr_y\ttr_norm()\t10TSL_x\tTSL_y\tTSL_norm"
                    << "\n";
            q_point.close();
          }
          else
          {
            q_point.open("q_points/q_point_" + std::to_string(quadrature_points_history[point].output_cords[0]) + "_" + std::to_string(point) + ".txt", std::ios_base::app);
            q_point << disp.norm() << "\t"
                    << quadrature_points_history[point].output_cords[0] << "\t" << quadrature_points_history[point].output_cords[1] << "\t"
                    << quadrature_points_history[point].output_jump[0] << "\t" << quadrature_points_history[point].output_jump[1] << "\t" << quadrature_points_history[point].output_jump.norm() << "\t"
                    << quadrature_points_history[point].output_numerical_traction[0] << "\t" << quadrature_points_history[point].output_numerical_traction[1] << "\t" << quadrature_points_history[point].output_numerical_traction.norm() << "\t"
                    << quadrature_points_history[point].output_analytical_traction[0] << "\t" << quadrature_points_history[point].output_analytical_traction[1] << "\t" << quadrature_points_history[point].output_analytical_traction.norm()
                    << "\n";
            q_point.close();
          }
        }
    };
    const auto cell_worker =
        [&](const auto &cell, auto &scratch_data, auto &copy_data) {};

    const auto boundary_worker = [&](const auto &cell,
                                     const unsigned int &face_no,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data) {};
    const auto copier = [&](const auto &c) {};

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
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

  template <int dim>
  void TSL_Friction_Freddi<dim>::q_point_PPV(const unsigned int cycle, const unsigned int non_lin)
  {

    TimerOutput::Scope t(computing_timer, "other");

    const auto face_worker = [&](const auto &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim> &scratch_data,
                                 CopyData &copy_data)
    {
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());
      // unsigned int counter = 0;
      for (unsigned int point = 0; point < n_q_points; ++point)
      {
        // counter++;

        std::ofstream q_point;

        q_point.open("output/cycle_" + std::to_string(cycle) + ".txt", std::ios_base::app);
        q_point << quadrature_points_history[point].output_cords[0] << "," << quadrature_points_history[point].output_cords[1] << ","
                << quadrature_points_history[point].output_jump[0] << "," << quadrature_points_history[point].output_jump[1] << "," << quadrature_points_history[point].output_jump.norm() << ","
                << quadrature_points_history[point].output_numerical_traction[0] << "," << quadrature_points_history[point].output_numerical_traction[1] << "," << quadrature_points_history[point].output_numerical_traction.norm() << ","
                << quadrature_points_history[point].output_analytical_traction[0] << "," << quadrature_points_history[point].output_analytical_traction[1] << "," << quadrature_points_history[point].output_analytical_traction.norm() << ","
                << quadrature_points_history[point].output_stress[0][0] << "," << quadrature_points_history[point].output_stress[1][0] << "," << quadrature_points_history[point].output_stress[1][1]
                << "\n";
        q_point.close();
      }
    };

    const auto cell_worker =
        [&](const auto &cell, auto &scratch_data, auto &copy_data) {};

    const auto boundary_worker = [&](const auto &cell,
                                     const unsigned int &face_no,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data) {};
    const auto copier = [&](const auto &c) {};

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
    CopyData copy_data;

    std::ofstream q_point;
    if (cycle == 0)
    {
      // create files
      for (int i = 1; i < par.cycles; i++)
      {
        q_point.open("output/cycle_" + std::to_string(i) + ".txt");
        q_point << "x"
                << ","
                << "y"
                << ","
                << "g_eff_x"
                << ","
                << "g_eff_y"
                << ","
                << "g_eff_norm"
                << ","
                << "traction_x"
                << ","
                << "traction_y"
                << ","
                << "traction_norm"
                << ","
                << "TSL_x"
                << ","
                << "TSL_y"
                << ","
                << "TSL_norm"
                << ","
                << "sigma_xx"
                << ","
                << "sigma_xy"
                << ","
                << "sigma_yy"
                << "\n";
        q_point.close();
      }
    }
    else
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

  template <int dim>
  void TSL_Friction_Freddi<dim>::run()
  {

    pcout << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    for (unsigned int refinement_level = 0; refinement_level < 1; refinement_level++)
    {
      triangulation.clear();

      if (par.geometry == "reactangle" || par.geometry == "shear")
      {
        Point<dim> P1, P2(1, 2);
        std::vector<unsigned int> repetitions{1, 2};
        GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
        triangulation.refine_global(par.initial_refinement + refinement_level);
        Point<dim> p;
        for (const auto &cell : triangulation.active_cell_iterators())
        {
          p = cell->center();
          // if (p[0] < 0.5)
          // {
          if (p[1] > 1.0)
            cell->set_material_id(1);
          else
            cell->set_material_id(2);
          // }
          // else
          // {
          //   if (p[1] > 1.25)
          //     cell->set_material_id(1);
          //   else
          //     cell->set_material_id(2);
          // }
        }

        for (const auto &cell : triangulation.active_cell_iterators())
        {
          for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
          {
            // if (std::abs(cell->vertex(v)[0] - 0.25) < 0.00001 && std::abs(cell->vertex(v)[1] < 1))
            // {
            //   cell->vertex(v)[0] = 1.2;
            // }

            // if (std::abs(cell->vertex(v)[0] - 0.75) < 0.00001 && std::abs(cell->vertex(v)[1] > 1))
            // {
            //   cell->vertex(v)[0] = 1.2;
            // }

            //  if (std::abs(cell->vertex(v)[0] - 0.5) < 0.00001 && std::abs(cell->vertex(v)[1] - 1) < 0.00001)
            //     {
            //       cell->vertex(v)[0] = 0.4;
            //     }

            //   if (std::abs(cell->vertex(v)[0] - 0.5) < 0.00001 && std::abs(cell->vertex(v)[1] - 1.25) < 0.00001)
            //     {
            //       cell->vertex(v)[0] = 0.6;
            //     }
          }
        }
      }

      if (par.geometry == "peeling")
      {
        Point<dim> P1(0.0, -3), P2(10, 2);
        std::vector<unsigned int> repetitions{8, 5};

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

      if (par.geometry == "hole")
      {
        GridIn<2> gridin;
        gridin.attach_triangulation(triangulation);
        std::ifstream f("hole.msh");
        // std::ifstream f("fine_hole.msh");

        gridin.read_msh(f);

        triangulation.refine_global(par.initial_refinement);
      }

      if (par.geometry == "half_matrix")
      {
        GridIn<2> gridin;
        gridin.attach_triangulation(triangulation);
        //   std::ifstream f("half_matrix.msh");
        std::ifstream f("half_matrix_fine.msh");

        gridin.read_msh(f);

        triangulation.refine_global(par.initial_refinement);
      }

      if (par.geometry == "full_matrix")
      {
        GridIn<2> gridin;
        gridin.attach_triangulation(triangulation);
        std::ifstream f("matrix.msh");

        gridin.read_msh(f);

        triangulation.refine_global(par.initial_refinement);
      }

      std::ofstream out("grid-1.vtu");
      GridOut grid_out;
      grid_out.write_vtu(triangulation, out);
      pcout << " written to "
            << "grid-1.vtu" << std::endl
            << std::endl;

      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int face : cell->face_indices())
          if ((cell->face(face)->at_boundary() == false) && (cell->material_id() != cell->neighbor(face)->material_id()))
            cell->face(face)->set_user_flag();
          else
            cell->face(face)->clear_user_flag();

      if (par.is_distorted)
      {
        GridTools::transform(&grid_y_transform<dim>, triangulation);

        // GridTools::distort_random(0.3, triangulation, true, 2882);
        //   GridTools::rotate(-M_PI/2.8, triangulation);
      }

      setup_quadrature_point_history();
      double error;
      unsigned int non_lin = 0;
      std::ofstream forces("forces" + std::to_string(refinement_level) + ".txt");
      std::ofstream residual("residual" + std::to_string(refinement_level) + ".txt");
      unsigned int max_nonlin = 0;
      unsigned int at_cycle = 0;
      unsigned int max_nonlin_iter = 150;
      unsigned int internal_iterations;
      // 1 means do one more extra innner iteration
     

      for (unsigned int cycle = 0; cycle < par.cycles; ++cycle)
      {
        error = 1;
        non_lin = 0;
        

        if (cycle < par.unloading || cycle > par.reloading)
        {
          disp[0] += par.displacementx;
          disp[1] += par.displacementy;
        }
        else
        {
          disp[0] -= par.displacementx;
          disp[1] -= par.displacementy;
        }

        pcout << " ####### cycle = " << cycle << " and displacement = " << disp[1] << " ###### \n";

        
        while ((error > par.error && non_lin < max_nonlin_iter))
        {

          if (cycle == 0 && non_lin == 0)
            setup_system(true);
          else
            setup_system(false);

          max_nonlin = std::max(max_nonlin, non_lin);

          if (max_nonlin == non_lin)
            at_cycle = cycle;

          if (max_nonlin == max_nonlin_iter - 1)
            throw std::runtime_error("max non_lin iterations reached b4 convergence");

          non_lin++;

          assemble_system(cycle, non_lin, error);
          error = system_rhs.l2_norm();

          if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          {
            residual.open("residual" + std::to_string(refinement_level) + ".txt", std::ios_base::app);
            if (residual.is_open() == false)
            {
              pcout << "File didn't open\n";
              exit(1);
            }

            else
              residual << error << "\n";

            residual.close();
          }

          //  pcout << "  Residual: " << compute_residual(0) << std::endl;
          pcout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;

          if (error > par.error)
          {
            
            pcout << "  non lin = " << non_lin << std::endl;
            solve(cycle);
            pcout << "  update l2_norm " << distributed_newton_update.l2_norm() << std::endl;
          }

          pcout << "cycle " << cycle << "\n------- \n ";
          pcout << "max nonlin iterations " << max_nonlin << "  at_cycle " << at_cycle << "\n------- \n ";

          if (error > 1e+20)
            throw std::runtime_error("Divergence in the solution");

          if (cycle % par.output_frequency == 0)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results(cycle, non_lin, refinement_level);
          }
        }
        q_point_PPV(cycle, non_lin);
        interface_jump_and_traction_at_qpoints(cycle);
        Tensor<1, dim> ten;
        Tensor<1, dim> jump;
        Tensor<1, dim> inter;

        reaction(ten, 3);
        interface_jump_and_traction(jump, inter);
        pcout << "   reaction force = " << ten[0] << "\t" << ten[1] << "\n";

        pcout << std::endl;
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          forces.open("forces" + std::to_string(refinement_level) + ".txt", std::ios_base::app);
          if (forces.is_open() == false)
          {
            pcout << "File didn't open\n";
            exit(1);
          }

          if (cycle == 1)
            forces << 0 << "\t" << 0 << "\t"
                   << 0 << "\t" << 0 << "\t"
                   << 0 << "\t" << 0 << "\t"
                   << 0 << "\t" << 0 << "\n";

          forces << disp[0] << "\t" << disp[1] << "\t"
                 << ten[0] << "\t" << ten[1] << "\t"
                 << jump[0] << "\t" << jump[1] << "\t"
                 << inter[0] << "\t" << inter[1] << "\n";

          forces.close();
        }
      }
      computing_timer.print_summary();
      computing_timer.reset();
    }
  }
} // closing namespace

#endif
