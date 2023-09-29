#ifndef FRICTION_ADAPTIVITY_H
#define FRICTION_ADAPTIVITY_H

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
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
} // namespace LA

namespace Friction_adaptivity
{

  double heaviside(double x)
  {
    return (x < 0) ? 0.0 : 1.0;
  }

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
      add_parameter("sig_c", sig_c, " ", this->prm, Patterns::Double());
      add_parameter("sig_ci", sig_ci, " ", this->prm, Patterns::Double());
      add_parameter("delta_c", delta_c, " ", this->prm, Patterns::Double());
      add_parameter("delta_ci", delta_ci, " ", this->prm, Patterns::Double());
      add_parameter("ku", ku, " ", this->prm, Patterns::Double());
      add_parameter("kui", kui, " ", this->prm, Patterns::Double());
      add_parameter("friction_coff", friction_coff, " ", this->prm, Patterns::Double());
      add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
      add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
      add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
      add_parameter("quadrature", quadrature, " ", this->prm, Patterns::Selection("trap|gauss|lobatto|gauss-reduced|lobatto-over"));
      add_parameter("symmetry", symmetry, " ", this->prm, Patterns::Integer());
      add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
      add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
      add_parameter("error", error, " ", this->prm, Patterns::Double());
      add_parameter("is_distorted", is_distorted, " ", this->prm, Patterns::Bool());
      add_parameter("adaptivity", adaptivity, " ", this->prm, Patterns::Bool());
      add_parameter("penalty_factor", penalty_factor, " ", this->prm, Patterns::Double());
      add_parameter("unloading", unloading, " ", this->prm, Patterns::Integer());
      add_parameter("reloading", reloading, " ", this->prm, Patterns::Integer());
      add_parameter("output_frequency", output_frequency, " ", this->prm, Patterns::Integer());
    }

    int cycles = 0, initial_refinement = 0, degree = 0, symmetry = 0, unloading = 0, reloading = 0, output_frequency = 0;
    double E = 0, nu = 0, Ep = 0, nup = 0, sig_c = 0, sig_ci = 0,
           delta_c = 0, delta_ci = 0, ku = 0, kui = 0, friction_coff,
           displacementx = 0, displacementy = 0, error = 0, penalty_factor = 0;
    bool is_distorted = 0, adaptivity = 0;
    std::string quadrature = "trap";
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

  template <int dim>
  struct PointHistory
  {
    bool is_damaged;
    bool is_reunload;
    double max_damage;
    double Freddi_g;

    Tensor<1, dim> traction_init;
    Tensor<1, dim> law_g;

    // output
    Tensor<1, dim> cords_out;

    void set_to_zero()
    {
      is_damaged = 0;
      is_reunload = 0;
      max_damage = 0;
      Freddi_g = 0;
      traction_init = 0;
      law_g = 0;
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
  class Friction_adaptivity
  {
  public:
    Friction_adaptivity(const unsigned int degree, const ProblemParameters<dim> &par);
    void run();

  private:
    void setup_system(const bool initial_step);
    void assemble_system(const unsigned int cycle, const unsigned int non_lin);
    void solve();
    void refine_mesh(const unsigned int cycle);
    void output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const;
    void q_point_PPV(const unsigned int cycle, const unsigned int non_lin);

    void reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id);
    void interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress);

    const ProblemParameters<dim> &par;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    const unsigned int degree;
    const QGauss<dim> quadrature;
    const MappingQ1<dim> mapping;

    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;

    std::unique_ptr<const Quadrature<dim - 1>> face_quadrature;

    LA::MPI::SparseMatrix system_matrix;

    LA::MPI::Vector ghosted_solution;
    LA::MPI::Vector distributed_solution;

    LA::MPI::Vector ghosted_newton_update;
    LA::MPI::Vector distributed_newton_update;

    LA::MPI::Vector system_rhs;

    AffineConstraints<double> constraints;
    ConvergenceTable convergence_table;

    Tensor<1, dim> disp;
    std::vector<PointHistory<dim>> quadrature_point_history;
    void setup_quadrature_point_history();

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    ConditionalOStream pcout;
    TimerOutput computing_timer;
  };

  template <int dim>
  Friction_adaptivity<dim>::Friction_adaptivity(const unsigned int degree, const ProblemParameters<dim> &par)
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
  void Friction_adaptivity<dim>::setup_quadrature_point_history()
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
  void Friction_adaptivity<dim>::setup_system(const bool initial_step)
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
  void Friction_adaptivity<dim>::assemble_system(const unsigned int cycle, const unsigned int non_lin)
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

      double penalty = (par.penalty_factor) * elasticity / (cell->measure());

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
      const FEFaceValuesBase<dim> &fe_fv1 = fe_iv.get_fe_face_values(0);
      const FEFaceValuesBase<dim> &fe_fv2 = fe_iv.get_fe_face_values(1);
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

      double penalty = (par.penalty_factor) * ((elasticity < nelasticity) ? nelasticity : elasticity) / (cell->measure());

      std::vector<Tensor<2, dim>> grads_1(n_q_points);
      fe_fv1[displacements].get_function_gradients(ghosted_solution, grads_1);
      std::vector<Tensor<2, dim>> grads_2(n_q_points);
      fe_fv2[displacements].get_function_gradients(ghosted_solution, grads_2);

      Tensor<1, dim> tr, tr_eff;
      Tensor<2, dim> Rot;
      Tensor<1, dim> tangential;
      Tensor<1, dim> TCZ_res;
      Tensor<2, dim> TCZ;
      Tensor<1, dim> law_g, law_p, law_q;

      for (unsigned int point = 0; point < n_q_points; ++point)
      {

        // if (q_points[point][0] < 0.25)
        //   if (cell->material_id() != ncell->material_id())
        //   {
        //     pcout << "History data before  \n";
        //     pcout << "At q point " << q_points[point][0] << "\n";
        //     pcout << "history_is_damaged " << quadrature_points_history[point].is_damaged << "\n";
        //     pcout << "history_is_reunaload " << quadrature_points_history[point].is_reunload << "\n";
        //     pcout << "qp_history_max_damage " << quadrature_points_history[point].max_damage << "\n";
        //     pcout << "qp_history_traction_init " << quadrature_points_history[point].traction_init << "\n";
        //     pcout << "qp_history_law_g " << quadrature_points_history[point].law_g << "\n";
        //     pcout << "qp_history_Freddi_g " << quadrature_points_history[point].Freddi_g << "\n \n";
        //   }
        // get sig_c delta_c and ku from aparameter file
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

        Tensor<2, dim> K;
        K[0][0] = sig_c / delta_c;
        K[1][1] = sig_c / delta_c;

        // Calculation for stress across interfaces to check for damage
        Tensor<2, dim> strain1 = 0.5 * (grads_1[point] + transpose(grads_1[point]));
        Tensor<2, dim> strain2 = 0.5 * (grads_2[point] + transpose(grads_2[point]));
        Tensor<2, dim> stress1 = (lambda * trace(strain1) * Identity + 2 * mu * strain1);
        Tensor<2, dim> stress2 = (nlambda * trace(strain2) * Identity + 2 * nmu * strain2);
        Tensor<2, dim> old_stress_avg = stress1 / 2 + stress2 / 2;

        tangential = cross_product_2d(normals[point]);
        tr[0] = tangential * (old_stress_avg * normals[point]);
        tr[1] = normals[point] * (old_stress_avg * normals[point]);
        tr_eff[0] = tr[0];
        tr_eff[1] = heaviside(tr[1]) * tr[1];

        if (tr_eff.norm() > sig_c)
          quadrature_points_history[point].is_damaged = true;

        // save traction as long as q point hasn't failed
        // This will be used to calculate the shifting
        if (quadrature_points_history[point].is_damaged == false)
          quadrature_points_history[point].traction_init = tr_eff;

        // calculation of rotation
        const double normaly = normals[point][1];
        const int sign_normalx = ((normals[point][0] > 0) ? 1 : ((normals[point][0] < 0) ? -1 : 0));
        Rot[0][0] = normaly;
        Rot[0][1] = -sign_normalx * std::sqrt(1 - normaly * normaly);
        Rot[1][0] = sign_normalx * std::sqrt(1 - normaly * normaly);
        Rot[1][1] = normaly;
        Rot = -Rot;

        law_g = Rot * old_solution_jumps[point];
        // Calculation of shift
        Tensor<1, dim> shift = invert(K) * quadrature_points_history[point].traction_init;
        law_g += shift;

        // calculate p and q
        law_p[0] = law_g[0];
        if (law_g[1] > 0)
          law_q[1] = heaviside(law_g[1]) * law_g[1];
        else
          law_p[0] = quadrature_points_history[point].Freddi_g;

        // calculate damge
        Tensor<1, dim> sig_u = K * (law_g);
        double damage = ku * (sig_c - sig_u.norm()) / ((sig_u.norm()) * (sig_c - ku));
        damage = std::max(0.0, std::min(1.0, damage));

        if (non_lin == 1)
          quadrature_points_history[point].max_damage = std::max(damage, quadrature_points_history[point].max_damage);

        // check for unloading and reloading
        // if (non_lin == 1)
        if ((quadrature_points_history[point].max_damage - damage) > 1e-8)
          quadrature_points_history[point].is_reunload = true;
        else
          quadrature_points_history[point].is_reunload = false;

        // we need this value for unlaoding-relaoding
        if (non_lin == 1)
          if (quadrature_points_history[point].is_reunload == false)
            quadrature_points_history[point].law_g = law_g;

        // calculation of sig trail
        Tensor<1, dim> subs_g;
        subs_g[1] = heaviside(-law_g[1]) * law_g[1];
        Tensor<1, dim> sig_trial = K * (law_g - (law_p + law_q)) + penalty * subs_g;

        // Columb and columb criteria
        double coff = par.friction_coff;
        double Coulomb = std::abs(sig_trial[0]) + coff * sig_trial[1];

        if (law_g[1] < 0)
        {
          if (Coulomb <= 0)
          {
            // do nothing to p
          }
          else
          {

            quadrature_points_history[point].Freddi_g += (1 / K[0][0]) * (coff * sig_trial[1] + std::abs(sig_trial[0])) * (std::abs(sig_trial[0]) / sig_trial[0]);
            law_p[0] = quadrature_points_history[point].Freddi_g;
          }
        }

        // TCZ_res law . Penlaty term is added in the end
        TCZ_res = K * (law_g - damage * (law_p + law_q));

        // TCZ calculation

        TCZ = 0;
        TCZ = K;

        // dTCZ /ddamage x ddamage/du
        if (damage < 1 - 1e-10 && quadrature_points_history[point].is_reunload == false)
        {
          Tensor<1, dim> ddamage_du;
          ddamage_du[0] = -K[0][0] * K[0][0] * law_g[0];
          ddamage_du[1] = -K[1][1] * K[1][1] * law_g[1];
          ddamage_du *= (1 / sig_u.norm() + 1 / (sig_c - sig_u.norm())) * damage / sig_u.norm();
          TCZ += outer_product(-K * (law_p + law_q), ddamage_du);
        }

        // dTCZ/dp x dp/du
        if (law_g[1] > 0)
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
            dp_du[0][0] = 1.0; //! where did this 1 come form
            dp_du[0][1] = coff * (K[1][1] / K[0][0]) * std::abs(law_g[0] - quadrature_points_history[point].Freddi_g) / (law_g[0] - quadrature_points_history[point].Freddi_g);
            TCZ += -K * damage * dp_du;
          }
        }

        //  dTCZ/dq * dq/du
        Tensor<2, dim> dq_du;
        dq_du[1][1] = heaviside(law_g[1]);
        TCZ += -K * damage * dq_du;

        if (quadrature_points_history[point].is_reunload)
        {

          damage = quadrature_points_history[point].max_damage;
          law_g -= shift;

          TCZ_res = K * (1 - damage) * law_g;
          TCZ_res[0] *= 1 + shift[0] / (quadrature_points_history[point].law_g[0] - shift[0]);
          TCZ_res[1] *= 1 + shift[1] / (quadrature_points_history[point].law_g[1] - shift[1]);

          TCZ = K * (1 - damage);
          TCZ[0][0] *= 1 + shift[0] / (quadrature_points_history[point].law_g[0] - shift[0]);
          TCZ[1][1] *= 1 + shift[1] / (quadrature_points_history[point].law_g[1] - shift[1]);
        }

        // If we are in compression
        if (law_g[1] < 0)
        {
          TCZ_res[1] = 0;
          TCZ[1][1] = 0;
          TCZ[0][1] = 0;
          TCZ[1][0] = 0;

          TCZ_res[1] += penalty * law_g[1];
          TCZ[1][1] += penalty;
        }

        // if (q_points[point][0] < 0.25)
        //   if (cell->material_id() != ncell->material_id())
        //   {
        //     pcout << "History data after  \n";
        //     pcout << "At q point " << q_points[point][0] << "\n";
        //     pcout << "history_is_damaged " << quadrature_points_history[point].is_damaged << "\n";
        //     pcout << "history_is_reunaload " << quadrature_points_history[point].is_reunload << "\n";
        //     pcout << "qp_history_max_damage " << quadrature_points_history[point].max_damage << "\n";
        //     pcout << "qp_history_traction_init " << quadrature_points_history[point].traction_init << "\n";
        //     pcout << "qp_history_law_g " << quadrature_points_history[point].law_g << "\n";
        //     pcout << "qp_history_Freddi_g " << quadrature_points_history[point].Freddi_g << "\n \n";
        //   }

        // if (q_points[point][0] < 0.014)
        //   if (cell->material_id() != ncell->material_id())
        //   {
        //     pcout << "data  \n";
        //     pcout << "At q point " << q_points[point][0] << "\n";
        //     pcout << "damage " << damage << "\n";
        //     pcout << "law_g " << law_g << "\n";
        //     pcout << "shift " << shift << "\n";
        //     pcout << "law_g shifted " << law_g - shift << "\n";
        //     pcout << "law_p " << law_p << "\n";
        //     pcout << "tr_eff " << tr_eff << "\n";
        //     pcout << "TCZ_res " << TCZ_res << "\n";
        //     pcout << "TCZ " << TCZ << "\n";

        //     pcout << "\n";
        //   }

        quadrature_points_history[point].cords_out = q_points[point];

        if (cell->material_id() == ncell->material_id())
          quadrature_points_history[point].is_damaged = false;

        for (unsigned int i = 0; i < n_dofs_face; ++i)
        {
          Tensor<2, dim> stress_i;
          if (i < n_dofs_face / 2)
          {
            Tensor<2, dim> grads_i_avg = 0.5 * fe_fv1[displacements].gradient(i, point);
            Tensor<2, dim> strain_i = 0.5 * (grads_i_avg + transpose(grads_i_avg));
            stress_i = (lambda * trace(strain_i) * Identity + 2 * mu * strain_i);
          }
          else
          {
            Tensor<2, dim> grads_i_avg = 0.5 * fe_fv2[displacements].gradient(i - n_dofs_face / 2, point);
            Tensor<2, dim> strain_i = 0.5 * (grads_i_avg + transpose(grads_i_avg));
            stress_i = (nlambda * trace(strain_i) * Identity + 2 * nmu * strain_i);
          }

          for (unsigned int j = 0; j < n_dofs_face; ++j)
          {

            Tensor<2, dim> stress_j;
            if (j < n_dofs_face / 2)
            {
              Tensor<2, dim> grads_j_avg = 0.5 * fe_fv1[displacements].gradient(j, point);
              Tensor<2, dim> strain_j = 0.5 * (grads_j_avg + transpose(grads_j_avg));
              stress_j = (lambda * trace(strain_j) * Identity + 2 * mu * strain_j);
            }
            else
            {
              Tensor<2, dim> grads_j_avg = 0.5 * fe_fv2[displacements].gradient(j - n_dofs_face / 2, point);
              Tensor<2, dim> strain_j = 0.5 * (grads_j_avg + transpose(grads_j_avg));
              stress_j = (nlambda * trace(strain_j) * Identity + 2 * nmu * strain_j);
            }

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
  void Friction_adaptivity<dim>::solve()
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
  void Friction_adaptivity<dim>::refine_mesh(const unsigned int cycle)
  {

    TimerOutput::Scope t(computing_timer, "Refine");
    Tensor<2, dim> Identity;
    Identity[0][0] = 1;
    Identity[1][1] = 1;

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
    CopyData copy_data;

    const auto face_flagging = [&](const auto &cell,
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

      double elasticity = 0, poisson = 0, nelasticity = 0, npoisson = 0, sig_c = 0;
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

      if (cell->material_id() == ncell->material_id())
        sig_c = par.sig_c;
      else
        sig_c = par.sig_ci;

      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      // PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());

      const FEFaceValuesBase<dim> &fe_fv1 = fe_iv.get_fe_face_values(0);
      const FEFaceValuesBase<dim> &fe_fv2 = fe_iv.get_fe_face_values(1);
      const FEValuesExtractors::Vector displacements(0);
      std::vector<Tensor<2, dim>> grads_1(n_q_points);
      fe_fv1[displacements].get_function_gradients(ghosted_solution, grads_1);
      std::vector<Tensor<2, dim>> grads_2(n_q_points);
      fe_fv2[displacements].get_function_gradients(ghosted_solution, grads_2);

      double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
      double mu = elasticity / (2 * (1 + poisson));
      double nlambda = nelasticity * npoisson / ((1 + npoisson) * (1 - 2 * npoisson));
      double nmu = nelasticity / (2 * (1 + npoisson));

      // if (quadrature_points_history != nullptr)
      if (cell->material_id() != ncell->material_id())
        for (unsigned int point = 0; point < n_q_points; ++point)
        // if (quadrature_points_history[point].traction_init.norm() > 0.8 * sig_c)
        //  if (cycle > 7)
        {
          Tensor<2, dim> strain1 = 0.5 * (grads_1[point] + transpose(grads_1[point]));
          Tensor<2, dim> strain2 = 0.5 * (grads_2[point] + transpose(grads_2[point]));
          Tensor<2, dim> stress1 = (lambda * trace(strain1) * Identity + 2 * mu * strain1);
          Tensor<2, dim> stress2 = (nlambda * trace(strain2) * Identity + 2 * nmu * strain2);
          Tensor<2, dim> old_stress_avg = stress1 / 2 + stress2 / 2;

          Tensor<1, dim> tangential = cross_product_2d(normals[point]);
          Tensor<1, dim> tr, tr_eff;
          tr[0] = tangential * (old_stress_avg * normals[point]);
          tr[1] = normals[point] * (old_stress_avg * normals[point]);
          tr_eff[0] = tr[0];
          tr_eff[1] = heaviside(tr[1]) * tr[1];

          if (tr_eff.norm() > 0.5 * sig_c)
          {
            cell->set_refine_flag();
            ncell->set_refine_flag();
          }

          if (cell->refine_flag_set() && cell->level() == 5)
            cell->clear_refine_flag();

          if (ncell->refine_flag_set() && ncell->level() == 5)
            ncell->clear_refine_flag();
        }
    };

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          nullptr,
                          nullptr,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_interior_faces_both | MeshWorker::assemble_ghost_faces_both,
                          nullptr,
                          face_flagging);

    // for (const auto cell : dof_handler.active_cell_iterators())
    //  // if (cell->is_locally_owned())
    //     for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
    //       if (cell->face(face)->at_boundary() == false)
    //       {
    //         const auto ncell = cell->neighbor(face);

    //         double sig_c = 0;
    //         if (cell->material_id() == ncell->material_id())
    //         {
    //           sig_c = par.sig_c;
    //         }
    //         else
    //         {
    //           sig_c = par.sig_ci;
    //         }

    //         PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(face)->user_pointer());

    //         if (quadrature_points_history != nullptr)
    //           if (cell->material_id() != ncell->material_id())
    //             for (unsigned int point = 0; point < face_quadrature->size(); ++point)
    //               if (quadrature_points_history[point].traction_init.norm() > 0.8 * sig_c)
    //               // if (cycle > 40)
    //               {
    //                 cell->set_refine_flag();
    //               //  ncell->set_refine_flag();
    //               }

    //         if (cell->refine_flag_set() && cell->level() == 6)
    //           cell->clear_refine_flag();

    //         // if (ncell->refine_flag_set() && ncell->level() == 6)
    //         //   ncell->clear_refine_flag();
    //       }

    unsigned int history_index = 0;
    for (auto &face : triangulation.active_face_iterators())
    {

      if (std::abs(face->center()[1] - 1) > 1e-5 || face->at_boundary())
        for (unsigned int q = 0; q < face_quadrature->size(); ++q)
        {
        //  quadrature_point_history[history_index].set_to_zero();
         // quadrature_point_history[history_index + q].set_to_zero();
        }

      history_index += face_quadrature->size();
    }

    triangulation.prepare_coarsening_and_refinement();

    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(distributed_solution);

    FE_DGQ<dim> history_fe(1);
    DoFHandler<dim> history_dof_handler(triangulation);
    history_dof_handler.distribute_dofs(history_fe);

    std::vector<LA::MPI::Vector> distributed_history_max_damage_field(4),
        ghosted_history_max_damage_field(4),
        distributed_history_Freddi_g_field(4),
        ghosted_history_Freddi_g_field(4);

    for (unsigned int j = 0; j < 4; ++j)
    {
      distributed_history_max_damage_field[j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_max_damage_field[j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      distributed_history_Freddi_g_field[j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_Freddi_g_field[j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
    }

    std::vector<std::vector<LA::MPI::Vector>>
        distributed_history_traction_init_field(dim, std::vector<LA::MPI::Vector>(4)),
        ghosted_history_traction_init_field(dim, std::vector<LA::MPI::Vector>(4)),
        distributed_history_law_g_field(dim, std::vector<LA::MPI::Vector>(4)),
        ghosted_history_law_g_field(dim, std::vector<LA::MPI::Vector>(4));

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < 4; ++j)
      {
        distributed_history_traction_init_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
        ghosted_history_traction_init_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);

        distributed_history_law_g_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
        ghosted_history_law_g_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      }

    std::vector<Vector<double>>
        local_history_max_damage_values_at_qpoints(4),
        local_history_max_damage_fe_values(4),
        local_history_Freddi_g_values_at_qpoints(4),
        local_history_Freddi_g_fe_values(4);

    for (unsigned int j = 0; j < 4; ++j)
    {
      local_history_max_damage_values_at_qpoints[j].reinit(face_quadrature->size());
      local_history_max_damage_fe_values[j].reinit(history_fe.dofs_per_cell);
      local_history_Freddi_g_values_at_qpoints[j].reinit(face_quadrature->size());
      local_history_Freddi_g_fe_values[j].reinit(history_fe.dofs_per_cell);
    }

    std::vector<std::vector<Vector<double>>>
        local_history_traction_init_values_at_qpoints(dim, std::vector<Vector<double>>(4)),
        local_history_traction_init_fe_values(dim, std::vector<Vector<double>>(4)),
        local_history_law_g_values_at_qpoints(dim, std::vector<Vector<double>>(4)),
        local_history_law_g_fe_values(dim, std::vector<Vector<double>>(4));

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < 4; ++j)
      {
        local_history_traction_init_values_at_qpoints[i][j].reinit(face_quadrature->size());
        local_history_traction_init_fe_values[i][j].reinit(history_fe.dofs_per_cell);

        local_history_law_g_values_at_qpoints[i][j].reinit(face_quadrature->size());
        local_history_law_g_fe_values[i][j].reinit(history_fe.dofs_per_cell);
      }

    FullMatrix<double> qpoint_to_dof_matrix(history_fe.dofs_per_face,
                                            face_quadrature->size());

    for (const auto cell : history_dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
          FETools::compute_projection_from_face_quadrature_points_matrix(history_fe,
                                                                         *face_quadrature,
                                                                         *face_quadrature,
                                                                         cell,
                                                                         face,
                                                                         qpoint_to_dof_matrix);

          PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(face)->user_pointer());

          if (quadrature_points_history != nullptr)
            for (unsigned int q = 0; q < face_quadrature->size(); ++q)
            {
              local_history_max_damage_values_at_qpoints[face][q] = quadrature_points_history[q].max_damage;
              local_history_Freddi_g_values_at_qpoints[face][q] = quadrature_points_history[q].Freddi_g;
            }

          qpoint_to_dof_matrix.vmult_add(local_history_max_damage_fe_values[face],
                                         local_history_max_damage_values_at_qpoints[face]);

          qpoint_to_dof_matrix.vmult_add(local_history_Freddi_g_fe_values[face],
                                         local_history_Freddi_g_values_at_qpoints[face]);

          cell->set_dof_values(local_history_max_damage_fe_values[face],
                               distributed_history_max_damage_field[face]);

          cell->set_dof_values(local_history_Freddi_g_fe_values[face],
                               distributed_history_Freddi_g_field[face]);

          local_history_max_damage_fe_values[face] = 0;
          local_history_Freddi_g_fe_values[face] = 0;
        }

    for (unsigned int i = 0; i < dim; ++i)
      for (const auto cell : history_dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          {
            FETools::compute_projection_from_face_quadrature_points_matrix(history_fe,
                                                                           *face_quadrature,
                                                                           *face_quadrature,
                                                                           cell,
                                                                           face,
                                                                           qpoint_to_dof_matrix);

            PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(face)->user_pointer());

            if (quadrature_points_history != nullptr)
              for (unsigned int q = 0; q < face_quadrature->size(); ++q)
              {
                local_history_traction_init_values_at_qpoints[i][face][q] = quadrature_points_history[q].traction_init[i];
                local_history_law_g_values_at_qpoints[i][face][q] = quadrature_points_history[q].law_g[i];
              }

            qpoint_to_dof_matrix.vmult_add(local_history_traction_init_fe_values[i][face],
                                           local_history_traction_init_values_at_qpoints[i][face]);

            qpoint_to_dof_matrix.vmult_add(local_history_law_g_fe_values[i][face],
                                           local_history_law_g_values_at_qpoints[i][face]);

            cell->set_dof_values(local_history_traction_init_fe_values[i][face],
                                 distributed_history_traction_init_field[i][face]);

            cell->set_dof_values(local_history_law_g_fe_values[i][face],
                                 distributed_history_law_g_field[i][face]);

            local_history_traction_init_fe_values[i][face] = 0;
            local_history_law_g_fe_values[i][face] = 0;
          }

    for (unsigned int j = 0; j < 4; ++j)
    {
      distributed_history_max_damage_field[j].compress(VectorOperation::insert);
      ghosted_history_max_damage_field[j] = distributed_history_max_damage_field[j];
      distributed_history_Freddi_g_field[j].compress(VectorOperation::insert);
      ghosted_history_Freddi_g_field[j] = distributed_history_Freddi_g_field[j];
    }

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < 4; ++j)
      {
        distributed_history_traction_init_field[i][j].compress(VectorOperation::insert);
        ghosted_history_traction_init_field[i][j] = distributed_history_traction_init_field[i][j];

        distributed_history_law_g_field[i][j].compress(VectorOperation::insert);
        ghosted_history_law_g_field[i][j] = distributed_history_law_g_field[i][j];
      }

    // DataOut<dim> data_out1;
    // data_out1.add_data_vector(history_dof_handler, distributed_history_max_damage_field[2], "max_damage_projected");
    // data_out1.build_patches();
    // data_out1.write_vtu_with_pvtu_record("./output/", "max_damage_projected", cycle, mpi_communicator, 3, 0);

    std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>> history_Damage_transfer(4, history_dof_handler);
    std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>> history_Freddi_g_transfer(4, history_dof_handler);

    for (unsigned int j = 0; j < 4; ++j)
    {
      history_Damage_transfer[j].prepare_for_coarsening_and_refinement(distributed_history_max_damage_field[j]);
      history_Freddi_g_transfer[j].prepare_for_coarsening_and_refinement(distributed_history_Freddi_g_field[j]);
    }

    std::vector<std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>>> history_traction_init_transfer(dim, std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>>(4, history_dof_handler));
    std::vector<std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>>> history_law_g_transfer(dim, std::vector<parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>>(4, history_dof_handler));
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < 4; ++j)
      {
        history_traction_init_transfer[i][j].prepare_for_coarsening_and_refinement(distributed_history_traction_init_field[i][j]);
        history_law_g_transfer[i][j].prepare_for_coarsening_and_refinement(distributed_history_law_g_field[i][j]);
      }

     triangulation.execute_coarsening_and_refinement();

     setup_quadrature_point_history();

    dof_handler.distribute_dofs(fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
    ghosted_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    distributed_solution.reinit(locally_owned_dofs, mpi_communicator);
    solution_transfer.interpolate(distributed_solution);
    ghosted_solution = distributed_solution;

    history_dof_handler.distribute_dofs(history_fe);
    for (unsigned int j = 0; j < 4; ++j)
    {
      distributed_history_max_damage_field[j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_max_damage_field[j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      history_Damage_transfer[j].interpolate(distributed_history_max_damage_field[j]);
      ghosted_history_max_damage_field[j] = distributed_history_max_damage_field[j];
      distributed_history_Freddi_g_field[j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_Freddi_g_field[j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      history_Freddi_g_transfer[j].interpolate(distributed_history_Freddi_g_field[j]);
      ghosted_history_Freddi_g_field[j] = distributed_history_Freddi_g_field[j];
    }

    for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < 4; ++j)
    {
      distributed_history_traction_init_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_traction_init_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      history_traction_init_transfer[i][j].interpolate(distributed_history_traction_init_field[i][j]);
      ghosted_history_traction_init_field[i][j] = distributed_history_traction_init_field[i][j];

      distributed_history_law_g_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), mpi_communicator);
      ghosted_history_law_g_field[i][j].reinit(history_dof_handler.locally_owned_dofs(), DoFTools::extract_locally_relevant_dofs(history_dof_handler), mpi_communicator);
      history_law_g_transfer[i][j].interpolate(distributed_history_law_g_field[i][j]);
      ghosted_history_law_g_field[i][j] = distributed_history_law_g_field[i][j];
    }

    // // DataOut<dim> data_out2;
    // // data_out2.add_data_vector(history_dof_handler, ghosted_history_traction_init_field[1], "tr1_transfered");
    // // data_out2.build_patches();
    // // data_out2.write_vtu_with_pvtu_record("./output/", "tr1_transfered", cycle, mpi_communicator, 3, 0);

    FEFaceValues<dim> fe_fv_damage(history_fe,
                                   *face_quadrature,
                                   update_values | update_quadrature_points);

    std::vector<double> max_damage_values(face_quadrature->size());
    std::vector<double> Freddi_g_values(face_quadrature->size());

    std::vector<std::vector<double>> traction_init_values(dim, std::vector<double>(face_quadrature->size()));
    std::vector<std::vector<double>> law_g_values(dim, std::vector<double>(face_quadrature->size()));

    for (const auto cell : history_dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        //     // pcout << "\n";
        //     // pcout << "cell index " << cell->index() << "\n";
        // pcout << "cell center " << cell->center() << "\n";

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
          // pcout << "face index " << cell->face(face)->index() << "\n";
          // pcout << "face center " << cell->face(face)->center() << "\n";

          fe_fv_damage.reinit(cell, face);
          const auto &q_points = fe_fv_damage.get_quadrature_points();

          fe_fv_damage.get_function_values(ghosted_history_max_damage_field[face], max_damage_values);
          fe_fv_damage.get_function_values(ghosted_history_Freddi_g_field[face], Freddi_g_values);
          for (unsigned int i = 0; i < dim; ++i)
          {
            fe_fv_damage.get_function_values(ghosted_history_traction_init_field[i][face], traction_init_values[i]);
            fe_fv_damage.get_function_values(ghosted_history_law_g_field[i][face], law_g_values[i]);
          }

          PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(face)->user_pointer());

          if (quadrature_points_history != nullptr)
            for (unsigned int q = 0; q < face_quadrature->size(); ++q)
            {
              // pcout << "qpoint " << q_points[q] << "\n";

              quadrature_points_history[q].max_damage = max_damage_values[q];
              quadrature_points_history[q].Freddi_g = Freddi_g_values[q];

              if (quadrature_points_history[q].max_damage > 0)
                quadrature_points_history[q].is_damaged = true;

              for (unsigned int i = 0; i < dim; ++i)
              {
                quadrature_points_history[q].traction_init[i] = traction_init_values[i][q];
                quadrature_points_history[q].law_g[i] = law_g_values[i][q];
              }
            }
        }
      }

    // const auto quadrature_reset = [&](const auto &cell,
    //                                   const unsigned int &f,
    //                                   const unsigned int &sf,
    //                                   const auto &ncell,
    //                                   const unsigned int &nf,
    //                                   const unsigned int &nsf,
    //                                   ScratchData<dim> &scratch_data,
    //                                   CopyData &copy_data)
    // {

    // };

    // MeshWorker::mesh_loop(dof_handler.begin_active(),
    //                       dof_handler.end(),
    //                       nullptr,
    //                       nullptr,
    //                       scratch_data,
    //                       copy_data,
    //                       MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
    //                       nullptr,
    //                       quadrature_reset);
  }

  template <int dim>
  void Friction_adaptivity<dim>::output_results(const unsigned int cycle, const unsigned int non_lin, const unsigned int refinement_level) const
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

  template <int dim>
  void Friction_adaptivity<dim>::q_point_PPV(const unsigned int cycle, const unsigned int non_lin)
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
      for (unsigned int point = 0; point < n_q_points; ++point)
      {

        std::ofstream q_point;

        q_point.open("output/cycle_" + std::to_string(cycle) + "_" + std::to_string(non_lin) + ".txt", std::ios_base::app);
        q_point << quadrature_points_history[point].cords_out[0] << "," << quadrature_points_history[point].cords_out[1] << ","
                << quadrature_points_history[point].max_damage
                << "\n";
        q_point.close();
      }
    };

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
    CopyData copy_data;

    std::ofstream q_point;
    if (cycle == 0)
    {
      // create files
      for (int i = 1; i < par.cycles; i++)
        for (int j = 1; j < 25; j++)
        {
          q_point.open("output/cycle_" + std::to_string(i) + "_" + std::to_string(j) + ".txt");
          q_point << "x"
                  << ","
                  << "y"
                  << ","
                  << "max_damage"
                  << "\n";
          q_point.close();
        }
    }
    else
      MeshWorker::mesh_loop(dof_handler.begin_active(),
                            dof_handler.end(),
                            nullptr,
                            nullptr,
                            scratch_data,
                            copy_data,
                            MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                            nullptr,
                            face_worker);
  }

  template <int dim>
  void Friction_adaptivity<dim>::reaction(Tensor<1, dim> &reaction_stress, const types::boundary_id &boundary_id)
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
  void Friction_adaptivity<dim>::interface_jump_and_traction(Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
  {
    TimerOutput::Scope t(computing_timer, "other");

    interface_stress = 0;
    interface_jump = 0;
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

      if (cell->material_id() != ncell->material_id())
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

    ScratchData<dim> scratch_data(mapping, fe, quadrature, *face_quadrature);
    CopyData copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          nullptr,
                          nullptr,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                          nullptr,
                          face_worker);

    counter = Utilities::MPI::sum(counter, mpi_communicator);
    interface_jump = Utilities::MPI::sum(interface_jump, mpi_communicator);
    interface_stress = Utilities::MPI::sum(interface_stress, mpi_communicator);
    interface_jump = interface_jump / counter;
  }

  template <int dim>
  void Friction_adaptivity<dim>::run()
  {

    for (unsigned int refinement_level = 0; refinement_level < 1; refinement_level++)
    {
      triangulation.clear();

      Point<dim> P1, P2(1, 2);
      std::vector<unsigned int> repetitions{1, 2};
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

      std::ofstream out("grid-1.vtu");
      GridOut grid_out;
      grid_out.write_vtu(triangulation, out);
      pcout << " written to "
            << "grid-1.vtu" << std::endl
            << std::endl;

      if (par.is_distorted)
      {
        // GridTools::transform(&grid_y_transform<dim>, triangulation);
        GridTools::distort_random(0.3, triangulation, true, 2882);
        //   GridTools::rotate(-M_PI/2, triangulation);
      }

      setup_quadrature_point_history();
      double error;
      unsigned int non_lin = 0;
      std::ofstream forces("forces" + std::to_string(refinement_level) + ".txt");
      unsigned int max_nonlin = 0;
      unsigned int at_cycle = 0;
      unsigned int max_nonlin_iter = 25;

      for (int cycle = 0; cycle < par.cycles; ++cycle)
      {
        error = 1;
        non_lin = 0;

        if (cycle < par.unloading || cycle > par.reloading)
        {
          disp[0] += par.displacementx;
          disp[1] = par.displacementy;
        }
        else
        {
          disp[0] -= par.displacementx;
          disp[1] = par.displacementy;
        }

        pcout << " ####### cycle = " << cycle << " and displacement = " << disp[1] << " ###### \n";

        if (par.adaptivity == true)
          if (cycle != 0)
            refine_mesh(cycle);

        while (error > par.error && non_lin < max_nonlin_iter)
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

          assemble_system(cycle, non_lin);
          error = system_rhs.l2_norm();

          pcout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;

          if (error > par.error)
          {
            pcout << "  non lin = " << non_lin << std::endl;
            solve();
            pcout << "  update l2_norm " << distributed_newton_update.l2_norm() << std::endl;
          }

          pcout << "max nonlin iterations " << max_nonlin << "  at_cycle " << at_cycle << "\n------- \n ";

          if (error > 1e+20)
            throw std::runtime_error("Divergence in the solution");

          if (cycle % par.output_frequency == 0)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results(cycle, non_lin, refinement_level);
          }
          // q_point_PPV(cycle, non_lin);
        }

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
}

#endif
