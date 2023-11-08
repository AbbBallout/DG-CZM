#if !defined Friction_adaptivity_everywhere_H
#define Friction_adaptivity_everywhere_H

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

////////////// WARNING ////////////////////
//!!!!!!!!!!!!!!!!!!!!
//////// Don't run in paralle with mesh adaptivity //////////
///////////////////////////////////////
//////////////////////////////////////

namespace LA
{
    using namespace dealii::LinearAlgebraPETSc;
} // namespace LA

template <int dim>
dealii::Point<dim> grid_y_transform(const dealii::Point<dim> &pt_in)
{

    const double max_slope = 1.0;
    double slope = (pt_in[1] - 0.5 * pt_in[1] * pt_in[1]) * max_slope;
    dealii::Point<dim> pt_out = pt_in;

    if (std::abs(pt_in[1] - 1) < 0.99)
        pt_out[1] = pt_in[1] + slope * (pt_in[0] - 0.5);

    return pt_out;
}

namespace Friction_adaptivity_everywhere
{
    using namespace dealii;

    double heaviside(double x)
    {
        return (x < 0) ? 0.0 : 1.0;
    }

    double sign(double x)
    {
        return (x < 0) ? -1.0 : 1.0;
    }

    template <int dim>
    class ProblemParameters : public ParameterAcceptor
    {
    public:
        ProblemParameters() : ParameterAcceptor("main")
        {
            add_parameter("testing_I", testing_I, " ", this->prm, Patterns::Integer());
            add_parameter("testing_F", testing_F, " ", this->prm, Patterns::Double());

            add_parameter("degree", degree, " ", this->prm, Patterns::Integer());
            add_parameter("quadrature", quadrature, " ", this->prm, Patterns::Selection("trap|gauss|lobatto|gauss-reduced|lobatto-over"));
            add_parameter("symmetry", symmetry, " ", this->prm, Patterns::Integer());
            add_parameter("penalty_factor", penalty_factor, " ", this->prm, Patterns::Double());

            enter_subsection("geometry");
            {
                add_parameter("geometry", geometry, " ", this->prm, Patterns::Selection("reactangle|half_matrix|matrix|hole|ENF|EN_10"));
                add_parameter("E1", E1, " ", this->prm, Patterns::Double());
                add_parameter("nu1", nu1, " ", this->prm, Patterns::Double());
                add_parameter("E2", E2, " ", this->prm, Patterns::Double());
                add_parameter("nu2", nu2, " ", this->prm, Patterns::Double());
                add_parameter("initial_refinement", initial_refinement, " ", this->prm, Patterns::Integer());
                add_parameter("is_distorted", is_distorted, " ", this->prm, Patterns::Bool());
            }
            leave_subsection();

            enter_subsection("CZM");
            {
                add_parameter("is_everywhere", is_everywhere, " ", this->prm, Patterns::Bool());
                add_parameter("type", type, " ", this->prm, Patterns::Selection("extrinsic|intrinsic"));

                add_parameter("penetration_penalty", penetration_penalty, " ", this->prm, Patterns::Double());

                add_parameter("sig_ics", sig_ics, " ", this->prm, Patterns::Double());
                add_parameter("delta_ics", delta_ics, " ", this->prm, Patterns::Double());
                add_parameter("G_ics", G_ics, " ", this->prm, Patterns::Double());
                add_parameter("sig_icn", sig_icn, " ", this->prm, Patterns::Double());
                add_parameter("delta_icn", delta_icn, " ", this->prm, Patterns::Double());
                add_parameter("G_icn", G_icn, " ", this->prm, Patterns::Double());

                add_parameter("sig_cs", sig_cs, " ", this->prm, Patterns::Double());
                add_parameter("delta_cs", delta_cs, " ", this->prm, Patterns::Double());
                add_parameter("G_cs", G_cs, " ", this->prm, Patterns::Double());
                add_parameter("sig_cn", sig_cn, " ", this->prm, Patterns::Double());
                add_parameter("delta_cn", delta_cn, " ", this->prm, Patterns::Double());
                add_parameter("G_cn", G_cn, " ", this->prm, Patterns::Double());

                add_parameter("friction_coff", friction_coff, " ", this->prm, Patterns::Double());
            }
            leave_subsection();

            enter_subsection("robustness");
            {
                add_parameter("always_check_for_damage_and_unloading", always_check_for_damage_and_unloading, " ", this->prm, Patterns::Bool());
                add_parameter("update_damage_once", update_damage_once, " ", this->prm, Patterns::Bool());
                add_parameter("newton_relaxation", newton_relaxation, " ", this->prm, Patterns::Double());
                add_parameter("with_adaptive_relaxation", with_adaptive_relaxation, " ", this->prm, Patterns::Bool());
                add_parameter("max_external_iterations", max_external_iterations, " ", this->prm, Patterns::Integer());
                add_parameter("damage_error", damage_error, " ", this->prm, Patterns::Double());
                add_parameter("external_iterations_error", external_iterations_error, " ", this->prm, Patterns::Double());
                add_parameter("ignore_non_convergence", ignore_non_convergence, " ", this->prm, Patterns::Bool());
                add_parameter("regularization1", regularization1, " ", this->prm, Patterns::Double());
                add_parameter("regularization2", regularization2, " ", this->prm, Patterns::Double());
                add_parameter("regularization3", regularization3, " ", this->prm, Patterns::Double());
            }
            leave_subsection();

            enter_subsection("adaptivity");
            {
                add_parameter("with_adaptivity", with_adaptivity, " ", this->prm, Patterns::Bool());
                add_parameter("threshold_stress_factor", threshold_stress_factor, " ", this->prm, Patterns::Double());
                add_parameter("max_cell_level", max_cell_level, " ", this->prm, Patterns::Integer());
            }

            leave_subsection();

            enter_subsection("solver_control");
            {
                add_parameter("cycles", cycles, " ", this->prm, Patterns::Integer());
                add_parameter("error", error, " ", this->prm, Patterns::Double());
                add_parameter("max_nonlin_iter", max_nonlin_iter, " ", this->prm, Patterns::Integer());
                add_parameter("displacementx", displacementx, " ", this->prm, Patterns::Double());
                add_parameter("displacementy", displacementy, " ", this->prm, Patterns::Double());
                add_parameter("unloading", unloading, " ", this->prm, Patterns::Integer());
                add_parameter("reloading", reloading, " ", this->prm, Patterns::Integer());
                add_parameter("output_directory", output_directory, " ", this->prm, Patterns::Anything());
                add_parameter("file_name", file_name, " ", this->prm, Patterns::Anything());
                add_parameter("output_frequency", output_frequency, " ", this->prm, Patterns::Integer());
            }
            leave_subsection();
        }

        int testing_I = 0, degree = 0, initial_refinement = 0, cycles = 0, max_nonlin_iter = 0, unloading = 0, reloading = 0,
            output_frequency = 0, symmetry = 0, max_cell_level = 0, max_external_iterations = 0;
        double testing_F = 0, E1 = 0, E2 = 0, nu1 = 0, nu2 = 0, penalty_factor = 0, error = 0, displacementx = 0, displacementy = 0,
               sig_ics = 0, delta_ics = 0, k_ius = 0, sig_icn = 0, G_ics = 0, delta_icn = 0, k_iun = 0, G_icn = 0,
               sig_cs = 0, delta_cs = 0, k_us = 0, sig_cn = 0, G_cs = 0, delta_cn = 0, k_un = 0, G_cn = 0,
               CZM_penalty_factor = 0, friction_coff = 0, threshold_stress_factor = 0, damage_error = 0, external_iterations_error = 0,
               newton_relaxation = 0, penetration_penalty = 0, regularization1 = 0, regularization2 = 0, regularization3 = 0;
        std::string quadrature = "gauss", type = "extrinsic", geometry = "reactangle", output_directory = "output", file_name = "forces";
        bool is_distorted = 0, with_adaptivity = 0, is_everywhere = 0, update_damage_once = 0,
             always_check_for_damage_and_unloading = 0, ignore_non_convergence = 0, with_adaptive_relaxation = 0;
    };

    template <int dim>
    class StressPostprocessor : public DataPostprocessorTensor<dim>
    {
    public:
        StressPostprocessor(const ProblemParameters<dim> &par)
            : DataPostprocessorTensor<dim>("stress",
                                           update_gradients),
              par(par)
        {
        }

        virtual void evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                           std::vector<Vector<double>> &computed_quantities) const override
        {
            AssertDimension(input_data.solution_gradients.size(),
                            computed_quantities.size());

            Tensor<2, dim> strain;
            Tensor<2, dim> stress;

            for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
            {
                AssertDimension(computed_quantities[p].size(),
                                (Tensor<2, dim>::n_independent_components));
                Tensor<2, dim> Identity;

                Identity[0][0] = 1;
                Identity[1][1] = 1;

                const typename DoFHandler<dim>::cell_iterator current_cell =
                    input_data.template get_cell<dim>();

                double elasticity = 0, poisson = 0;

                if (current_cell->material_id() == 1)
                {
                    elasticity = par.E1;
                    poisson = par.nu1;
                }
                else
                {
                    elasticity = par.E2;
                    poisson = par.nu2;
                }

                double lambda = elasticity * poisson / ((1 + poisson) * (1 - 2 * poisson));
                double mu = elasticity / (2 * (1 + poisson));

                for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                        strain[d][e] = (input_data.solution_gradients[p][d][e] + input_data.solution_gradients[p][e][d]) / 2;

                stress = lambda * trace(strain) * Identity + 2 * mu * strain;

                for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                        computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = stress[d][e];
            }

            // Tensor<2, dim> str_avg;
            // for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
            // {

            //     for (unsigned int d = 0; d < dim; ++d)
            //         for (unsigned int e = 0; e < dim; ++e)
            //             str_avg[d][e] += computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))];
            // }

            // for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
            // {

            //     for (unsigned int d = 0; d < dim; ++d)
            //         for (unsigned int e = 0; e < dim; ++e)
            //             computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))] = str_avg[d][e] / (input_data.solution_gradients.size());
            // }
        }

    private:
        const ProblemParameters<dim> &par;
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
        bool is_damaged;
        bool is_fully_damaged;
        bool is_reunload;

        double max_damage;
        double Freddi_g;

        double damage;
        double damage_init;

        Tensor<1, dim> traction_init;
        Tensor<1, dim> law_g;
    };

    template <int dim>
    class Friction_adaptivity_everywhere
    {
    public:
        Friction_adaptivity_everywhere(const ProblemParameters<dim> &par);
        void run();

    private:
        void setup_system(const bool initial_step);
        void assemble_system(const unsigned int cycle, const unsigned int non_lin);
        void solve(double prev_error);
        double calculate_damage_error();
        void refine_mesh();
        void output_results(const unsigned int cycle) const;

        void reaction_and_traction(const types::boundary_id &boundary_id, Tensor<1, dim> &reaction_stress, Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress);

        void ENLM(const unsigned int material_ID, double &E, double &nu, double &lambda, double &mu);
        void setup_quadrature_point_history();

        const ProblemParameters<dim> &par;
        MPI_Comm mpi_communicator;
        parallel::distributed::Triangulation<dim> triangulation;
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

        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        ConditionalOStream pcout;
        TimerOutput computing_timer;

        Tensor<2, dim> Identity;
        double newton_relaxation;
    };

    template <int dim>
    Friction_adaptivity_everywhere<dim>::Friction_adaptivity_everywhere(const ProblemParameters<dim> &par)
        : par(par),
          mpi_communicator(MPI_COMM_WORLD),
          triangulation(mpi_communicator),
          quadrature(par.degree + 1),
          mapping(),
          fe(FE_DGQ<dim>(par.degree), dim),
          dof_handler(triangulation),
          pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
          computing_timer(mpi_communicator,
                          pcout,
                          TimerOutput::never,
                          TimerOutput::wall_times)
    {
        TimerOutput::Scope t(computing_timer, "other");
        disp = 0.0;
        newton_relaxation = par.newton_relaxation;

        Identity[0][0] = 1;
        Identity[1][1] = 1;

        if (par.quadrature == "gauss")
            face_quadrature = std::make_unique<const QGauss<dim - 1>>(par.degree + 1);
        if (par.quadrature == "lobatto")
            face_quadrature = std::make_unique<const QGaussLobatto<dim - 1>>(par.degree + 1);
        if (par.quadrature == "gauss-reduced")
            face_quadrature = std::make_unique<const QGauss<dim - 1>>(par.degree);
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::ENLM(const unsigned int material_ID, double &E, double &nu, double &lambda, double &mu)
    {
        if (material_ID == 1)
        {
            E = par.E1;
            nu = par.nu1;
        }
        else
        {
            E = par.E2;
            nu = par.nu2;
        }

        lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        mu = E / (2 * (1 + nu));
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::setup_quadrature_point_history()
    {
        TimerOutput::Scope t(computing_timer, "other");

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
    void Friction_adaptivity_everywhere<dim>::setup_system(const bool initial_step)
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
    void Friction_adaptivity_everywhere<dim>::assemble_system(const unsigned int cycle, const unsigned int non_lin)
    {
        TimerOutput::Scope t(computing_timer, "Assemble");

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

            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

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

            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

            double penalty = (par.penalty_factor) * elasticity / (cell->diameter());

            if (par.geometry == "reactangle" || par.geometry == "matrix" || par.geometry == "hole")
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

            if (par.geometry == "half_matrix")
            {
                if ((cell->face(face_no)->boundary_id() == 0) || (cell->face(face_no)->boundary_id() == 2))
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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      (0.0) * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * (0.0) *
                                      JxW[point]; // dx
                        }
                    }

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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      (disp * normals[point]) * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * (disp * normals[point]) *
                                      JxW[point]; // dx
                        }
                    }
            }

            if (par.geometry == "ENF")
            {
                double rad = 2;

                if (cell->face(face_no)->boundary_id() == 2 && std::abs(cell->face(face_no)->center()[0] - 0) < rad)
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
                                      (0.0) * disp * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * (0.0) * disp *
                                      JxW[point]; // dx
                        }
                    }

                if (cell->face(face_no)->boundary_id() == 2 && std::abs(cell->face(face_no)->center()[0] - 140) < rad)
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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      (0.0) * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * (0.0) *
                                      JxW[point]; // dx
                        }
                    }

                if (cell->face(face_no)->boundary_id() == 3 && (std::abs(cell->face(face_no)->center()[0] - 70) < rad / 2))
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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      disp * normals[point] * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * disp * normals[point] *
                                      JxW[point]; // dx
                        }
                    }
            }

            if (par.geometry == "EN_10")
            {
                double h = 55.0, l = 240.0, b = 10.0;

                Tensor<1, dim> compressive_force;
                compressive_force[1] = -0.6;

                if (cell->face(face_no)->boundary_id() == 0 && std::abs(cell->face(face_no)->center()[1]) < h / 2 && cell->material_id() == 1)
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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      disp * normals[point] * JxW[point] // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * disp * normals[point] *
                                      JxW[point]; // dx
                        }
                    }

                if ((cell->face(face_no)->boundary_id() == 1 && std::abs(cell->face(face_no)->center()[1]) > h / 2 + b / 2 && cell->material_id() == 1) || (cell->face(face_no)->boundary_id() == 2))
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
                                    ((-fe_fv[displacements].value(i, point) * normals[point]) *
                                         (((lambda * trace(strainj) * Identity + 2 * mu * strainj) * normals[point]) * normals[point])

                                     + par.symmetry *
                                           (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                           (fe_fv[displacements].value(j, point) * normals[point]) // Symetry term

                                     + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (fe_fv[displacements].value(j, point) * normals[point])) *
                                    JxW[point];
                            }

                            copy_data.cell_rhs(i) +=

                                -(
                                    -(fe_fv[displacements].value(i, point) * normals[point]) *
                                        (((lambda * trace(strain) * Identity + 2 * mu * strain) * normals[point]) * normals[point])

                                    + par.symmetry *
                                          (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                          (old_solution_values[point] * normals[point]) // Symetry term

                                    + penalty * (fe_fv[displacements].value(i, point) * normals[point]) * (old_solution_values[point] * normals[point])) *
                                    JxW[point]

                                + par.symmetry * (((lambda * trace(straini) * Identity + 2 * mu * straini) * normals[point]) * normals[point]) *
                                      disp * normals[point] * JxW[point] * 0.0 // Symetry term

                                + penalty *
                                      fe_fv[displacements].value(i, point) * normals[point] * disp * normals[point] * 0.0 *
                                      JxW[point]; // dx
                        }
                    }

                if (cell->face(face_no)->boundary_id() == 3)
                    for (unsigned int point = 0; point < n_q_points; ++point)
                    {

                        for (unsigned int i = 0; i < n_facet_dofs; ++i)
                            copy_data.cell_rhs(i) += fe_fv[displacements].value(i, point) * compressive_force * JxW[point]; // dx
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

            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

            double nelasticity = 0, npoisson = 0, nlambda = 0, nmu = 0;
            ENLM(ncell->material_id(), nelasticity, npoisson, nlambda, nmu);

            const FEValuesExtractors::Vector displacements(0);

            std::vector<Tensor<1, dim>> old_solution_jumps(n_q_points);
            fe_iv[displacements].get_jump_in_function_values(ghosted_solution, old_solution_jumps);

            std::vector<Tensor<2, dim>> gradu_avg(n_q_points);
            fe_iv[displacements].get_average_of_function_gradients(ghosted_solution, gradu_avg);

            std::vector<Tensor<2, dim>> grads_1(n_q_points);
            fe_fv1[displacements].get_function_gradients(ghosted_solution, grads_1);
            std::vector<Tensor<2, dim>> grads_2(n_q_points);
            fe_fv2[displacements].get_function_gradients(ghosted_solution, grads_2);

            PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());

            double penalty = (par.penalty_factor) * (elasticity / 2 + nelasticity / 2) / std::min(cell->diameter(), ncell->diameter());
            Tensor<1, dim> sig_c, delta_c, G_c;
            if (cell->material_id() == ncell->material_id())
            {
                sig_c[0] = par.sig_cs;
                delta_c[0] = par.delta_cs;
                G_c[0] = par.G_cs;
                sig_c[1] = par.sig_cn;
                delta_c[1] = par.delta_cn;
                G_c[1] = par.G_cn;
            }
            else
            {
                sig_c[0] = par.sig_ics;
                delta_c[0] = par.delta_ics;
                G_c[0] = par.G_ics;
                sig_c[1] = par.sig_icn;
                delta_c[1] = par.delta_icn;
                G_c[1] = par.G_icn;
            }

            Tensor<2, dim> slope;
            slope[0][0] = sig_c[0] / delta_c[0];
            slope[1][1] = sig_c[1] / delta_c[1];

            Tensor<1, dim> tangential;
            Tensor<1, dim> traction, traction_eff;
            Tensor<2, dim> Rot;
            Tensor<1, dim> law_g, law_p, law_q;
            Tensor<1, dim> TCZ_res;
            Tensor<2, dim> TCZ;

            for (unsigned int point = 0; point < n_q_points; ++point)
            {
                Tensor<2, dim> strain1 = 0.5 * (grads_1[point] + transpose(grads_1[point]));
                Tensor<2, dim> strain2 = 0.5 * (grads_2[point] + transpose(grads_2[point]));
                Tensor<2, dim> stress1 = (lambda * trace(strain1) * Identity + 2 * mu * strain1);
                Tensor<2, dim> stress2 = (nlambda * trace(strain2) * Identity + 2 * nmu * strain2);
                Tensor<2, dim> old_stress_avg = stress1 / 2 + stress2 / 2;

                tangential = cross_product_2d(normals[point]);

                traction[0] = tangential * (old_stress_avg * normals[point]);
                traction[1] = normals[point] * (old_stress_avg * normals[point]);
                traction_eff[0] = traction[0];
                traction_eff[1] = heaviside(traction[1]) * traction[1];

                const double normaly = normals[point][1];
                const int sign_normalx = ((normals[point][0] > 0) ? 1 : ((normals[point][0] < 0) ? -1 : 0));
                Rot[0][0] = normaly;
                Rot[0][1] = -sign_normalx * std::sqrt(1 - normaly * normaly);
                Rot[1][0] = sign_normalx * std::sqrt(1 - normaly * normaly);
                Rot[1][1] = normaly;
                Rot = -Rot;

                law_g = Rot * old_solution_jumps[point];

                Tensor<1, dim> shift = invert(slope) * quadrature_points_history[point].traction_init;
                if (par.type == "intrinsic")
                    shift = 0;

                Tensor<1, dim> law_g_unshifted = law_g;
                law_g += shift;

                double coff = par.friction_coff;
                // for robustness
                if (non_lin <= 1 && coff > 0 + 1e-6)
                    coff = 10;

                if (coff > 1e-6)
                    law_p[0] = heaviside(law_g[1]) * law_g[0] + heaviside(-law_g[1]) * quadrature_points_history[point].Freddi_g;
                else
                    law_p[0] = law_g[0];
                law_q[1] = law_g[1];

                // Columb and columb criteria
                Tensor<1, dim> tau_trial = slope * (law_g - (law_p + law_q));
                tau_trial[1] += heaviside(-law_g_unshifted[1]) * par.penetration_penalty * law_g_unshifted[1];

                double Coulomb = std::abs(tau_trial[0]) + coff * tau_trial[1];

                if (Coulomb > 0 && law_g[1] < 0 && coff > 1e-6)
                {
                    quadrature_points_history[point].Freddi_g += (1 / slope[0][0]) * Coulomb * sign(tau_trial[0]);
                    law_p[0] = quadrature_points_history[point].Freddi_g;
                }

                // caclulation of damage
                double damage = 0, delta_ec = 0, delta_eu = 0;

                delta_ec = std::pow(slope[0][0] * (law_g[0] / law_g.norm()) / sig_c[0], 2);
                delta_ec += std::pow(slope[1][1] * (law_g[1] / law_g.norm()) / sig_c[1], 2);
                delta_ec = 1 / sqrt(delta_ec);

                delta_eu = slope[0][0] * delta_ec * std::pow(law_g[0] / law_g.norm(), 2) / (2 * G_c[0]);
                delta_eu += slope[1][1] * delta_ec * std::pow(law_g[1] / law_g.norm(), 2) / (2 * G_c[1]);
                delta_eu = 1 / delta_eu;

                damage = delta_eu * (law_g.norm() - delta_ec) / (law_g.norm() * (delta_eu - delta_ec));

                if (cycle == 0 /*&& non_lin == 1 */) // initially everything is nan
                    damage = 0;
                damage = std::max(0.0, std::min(1.0, damage));

                quadrature_points_history[point].damage = damage;

                if (non_lin == 1)
                    quadrature_points_history[point].damage_init = damage;

                if (par.update_damage_once == true)
                    damage = quadrature_points_history[point].damage_init;
                ////// #######  CONTROLS ######## ///////

                if ((non_lin == 1 && par.always_check_for_damage_and_unloading == false) || par.always_check_for_damage_and_unloading == true)
                    if (par.type == "extrinsic")
                        if (std::pow(traction_eff[0] / sig_c[0], 2) + std::pow(traction_eff[1] / sig_c[1], 2) > 1 || damage > 0 + 1e-6)
                            quadrature_points_history[point].is_damaged = true;

                if (par.type == "intrinsic")
                    quadrature_points_history[point].is_damaged = true;

                if (quadrature_points_history[point].is_damaged == false)
                    quadrature_points_history[point].traction_init = traction_eff;

                if (non_lin == 1)
                    quadrature_points_history[point].max_damage = std::max(damage, quadrature_points_history[point].max_damage);

                // check for unloading and reloading
                // if (non_lin == 1)
                if ((non_lin == 1 && par.always_check_for_damage_and_unloading == false) || par.always_check_for_damage_and_unloading == true)
                    if (quadrature_points_history[point].is_fully_damaged == false)
                    {
                        if ((quadrature_points_history[point].max_damage - damage) > 1e-8)
                            quadrature_points_history[point].is_reunload = true;
                        else
                            quadrature_points_history[point].is_reunload = false;
                    }

                // we need this value for unlaoding-relaoding
                if (non_lin == 1)
                    if (quadrature_points_history[point].is_reunload == false)
                        quadrature_points_history[point].law_g = law_g;

                // if (non_lin == 1)
                if (quadrature_points_history[point].is_damaged && damage > 1 - 1e-12)
                    quadrature_points_history[point].is_fully_damaged = true;

                if (par.geometry == "ENF" && cell->material_id() != ncell->material_id() && cell->face(f)->center()[0] > 100)
                {
                    quadrature_points_history[point].is_fully_damaged = true;
                    quadrature_points_history[point].is_damaged = true;
                    quadrature_points_history[point].max_damage = 1.0;
                }

                if (par.geometry == "EN_10")
                {
                    if (par.is_everywhere == false)
                        if (cell->face(f)->center()[1] > 55 / 2)
                            quadrature_points_history[point].is_damaged = false;

                    if (par.is_everywhere == true)
                        if (cell->material_id() == 1 && ncell->material_id() == 1)
                            quadrature_points_history[point].is_damaged = false;
                }

                if (quadrature_points_history[point].is_fully_damaged)
                    damage = 1.0;

                //  END CONTROLS ### ///////////////

                // TCZ && TCZ_res calculation
                if (quadrature_points_history[point].is_reunload)
                {
                    damage = quadrature_points_history[point].max_damage;

                    // shift is zero for an intrinsic law so this is consitent

                    TCZ_res = slope * (law_g - shift - damage * (law_p + law_q - shift));

                    TCZ = 0;
                    // dTCZ/du
                    TCZ[0][0] = 1.0;
                    TCZ[1][1] = 1.0;

                    // dTCZ/dp * dp/du
                    if (law_g[1] > 0 || coff < 1e-6)
                    {
                        Tensor<2, dim> dp_du;
                        dp_du[0][0] = 1.0;
                        TCZ += -damage * dp_du;
                    }
                    else
                    {
                        if (Coulomb > 0)
                        {
                            Tensor<2, dim> dp_du;
                            dp_du[0][0] = 1.0;
                            dp_du[0][1] = coff * ((slope[1][1]) / slope[0][0]) * sign(law_g[0] - quadrature_points_history[point].Freddi_g);

                            dp_du[0][1] += coff * ((par.penetration_penalty) / slope[0][0]) * sign(law_g[0] - quadrature_points_history[point].Freddi_g);

                            TCZ += -damage * dp_du;
                        }
                    }

                    // dTCZ/dq * dq/du
                    Tensor<2, dim> dq_du;
                    dq_du[1][1] = 1.0;
                    TCZ += -damage * dq_du;

                    TCZ = slope * TCZ;

                    double small = 1e-10; // for regularization
                    TCZ_res[0] *= 1 + shift[0] / (quadrature_points_history[point].law_g[0] - shift[0] + small);
                    TCZ_res[1] *= 1 + shift[1] / (quadrature_points_history[point].law_g[1] - shift[1] + small);
                    TCZ[0] *= 1 + shift[0] / (quadrature_points_history[point].law_g[0] - shift[0] + small);
                    TCZ[1] *= 1 + shift[1] / (quadrature_points_history[point].law_g[1] - shift[1] + small);

                    TCZ_res[1] += heaviside(-law_g_unshifted[1]) * par.penetration_penalty * law_g_unshifted[1];
                    TCZ[1][1] += heaviside(-law_g_unshifted[1]) * par.penetration_penalty;
                }
                else
                {
                    TCZ_res = slope * (law_g - damage * (law_p + law_q));

                    TCZ = 0;

                    // dTCZ/du
                    TCZ[0][0] = 1.0;
                    TCZ[1][1] = 1.0;

                    // dTCZ/Ds * dDs/du
                    if (par.update_damage_once == false || non_lin == 1)
                        if (damage > 0 && quadrature_points_history[point].is_fully_damaged == false)
                        {
                            Tensor<1, dim> ddamage_du;

                            double ddelta_dus = law_g[0] / law_g.norm();
                            double dec_dus = law_g[0] * delta_ec * (1 - std::pow(slope[0][0] * delta_ec / sig_c[0], 2)) / (law_g.norm() * law_g.norm());
                            double duc_dus = 2 * law_g[0] * delta_eu * (1 - slope[0][0] * delta_ec * delta_eu / (2 * G_c[0])) - (delta_eu / delta_ec) * dec_dus;

                            double ddelta_dun = law_g[1] / law_g.norm();
                            double dec_dun = law_g[1] * delta_ec * (1 - std::pow(slope[1][1] * delta_ec / sig_c[1], 2)) / (law_g.norm() * law_g.norm());
                            double duc_dun = 2 * law_g[1] * delta_eu * (1 - slope[1][1] * delta_ec * delta_eu / (2 * G_c[1])) - (delta_eu / delta_ec) * dec_dun;

                            ddamage_du[0] = duc_dus * (law_g.norm() - delta_ec) * law_g.norm() * (delta_eu - delta_ec);
                            ddamage_du[0] += (ddelta_dus - dec_dus) * delta_eu * law_g.norm() * (delta_eu - delta_ec);
                            ddamage_du[0] += -ddelta_dus * (delta_eu - delta_ec) * delta_eu * (law_g.norm() - delta_ec);
                            ddamage_du[0] += -(duc_dus - dec_dus) * law_g.norm() * delta_eu * (law_g.norm() - delta_ec);

                            ddamage_du[1] = duc_dun * (law_g.norm() - delta_ec) * law_g.norm() * (delta_eu - delta_ec);
                            ddamage_du[1] += (ddelta_dun - dec_dun) * delta_eu * law_g.norm() * (delta_eu - delta_ec);
                            ddamage_du[1] += -ddelta_dun * (delta_eu - delta_ec) * delta_eu * (law_g.norm() - delta_ec);
                            ddamage_du[1] += -(duc_dun - dec_dun) * law_g.norm() * delta_eu * (law_g.norm() - delta_ec);

                            ddamage_du *= 1 / std::pow(law_g.norm() * (delta_eu - delta_ec), 2);

                            TCZ += outer_product(-law_p - law_q, ddamage_du);
                        }

                    // dTCZ/dp * dp/du
                    if (law_g[1] > 0 || coff < 1e-6)
                    {
                        Tensor<2, dim> dp_du;
                        dp_du[0][0] = 1.0;
                        TCZ += -damage * dp_du;
                    }
                    else
                    {

                        if (Coulomb > 0)
                        {
                            Tensor<2, dim> dp_du;
                            dp_du[0][0] = 1.0;
                            dp_du[0][1] = coff * ((slope[1][1]) / slope[0][0]) * sign(law_g[0] - quadrature_points_history[point].Freddi_g);

                            dp_du[0][1] += coff * ((par.penetration_penalty) / slope[0][0]) * sign(law_g[0] - quadrature_points_history[point].Freddi_g);

                            TCZ += -damage * dp_du;
                        }
                    }

                    // dTCZ/dq * dq/du
                    Tensor<2, dim> dq_du;
                    dq_du[1][1] = 1.0;
                    TCZ += -damage * dq_du;

                    TCZ = slope * TCZ;

                    TCZ_res[1] += heaviside(-law_g_unshifted[1]) * par.penetration_penalty * law_g_unshifted[1];
                    TCZ[1][1] += heaviside(-law_g_unshifted[1]) * par.penetration_penalty;
                }

                // Regularize
                if (law_g[1] < 0)
                    if (damage > 0)
                    {
                        TCZ[0][0] += par.regularization1;

                        if (law_g[1] > -1e-12)
                            TCZ[1][0] += par.regularization2;

                        TCZ[0][0] += par.regularization3;
                        TCZ[1][0] += par.regularization3;
                    }

                if (par.is_everywhere == false)
                    if (cell->material_id() == ncell->material_id())
                    {
                        quadrature_points_history[point].is_damaged = false;
                        quadrature_points_history[point].max_damage = 0;
                        quadrature_points_history[point].damage_init = 0;
                    }

                if (par.geometry == "matrix")
                {
                    if (cell->material_id() != 1 && ncell->material_id() != 1)
                        quadrature_points_history[point].is_damaged = false;

                    if (std::abs(cell->center()[1] - 0.5) > 0.48)
                        quadrature_points_history[point].is_damaged = false;
                }
                // if (cell->material_id() != ncell->material_id())
                // if (q_points[point][1] > 0.91)
                // if (quadrature_points_history[point].is_reunload)
                if (quadrature_points_history[point].is_damaged)
                {
                    //   pcout << "At q point " << q_points[point] << "\n";
                    //  pcout << "is damaged " << quadrature_points_history[point].is_damaged << "\n";
                    // pcout << "is reloaded " << quadrature_points_history[point].is_reunload << "\n";
                    //  pcout << "is fully damaged " << quadrature_points_history[point].is_fully_damaged << "\n";
                    // pcout << "traction " << traction << "\n";

                    // pcout << "damage " << damage << "\n";
                    //   pcout << "history.max_damage " << quadrature_points_history[point].max_damage << "\n";
                    // pcout << "shift " << shift << "\n";
                    // pcout << "law_g shifted " << law_g << "\n";
                    // pcout << "law_g unshifted " << law_g_unshifted << "\n";
                    //  pcout << "history law_g " << quadrature_points_history[point].law_g << "\n";
                    //  pcout << "history traction_init " << quadrature_points_history[point].traction_init << "\n";

                    //  pcout << "law_p " << law_p << "\n";
                    // pcout << "law_q " << law_q << "\n";
                    // pcout << "tau_trial " << tau_trial << "\n";
                    // pcout << "tau_u " << tau_u << "\n";
                    // pcout << "Coulomb " << Coulomb << "\n";
                    //   pcout << "Freddi_g " << quadrature_points_history[point].Freddi_g << "\n";

                    //  pcout << "TCZ_res " << TCZ_res << "\n";
                    //  pcout << "TCZ " << TCZ << "\n";

                    //   pcout << "\n";
                }

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
    void Friction_adaptivity_everywhere<dim>::solve(double prev_error)
    {

        TimerOutput::Scope t(computing_timer, "solve");
        pcout << "solving" << std::endl;
        SolverControl solver_control;
        PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
        solver.solve(system_matrix, distributed_newton_update, system_rhs);
        ghosted_newton_update = distributed_newton_update;

        if (par.with_adaptive_relaxation)
        {

            if (system_rhs.l2_norm() > prev_error)
                newton_relaxation = newton_relaxation * 0.75;
            else if (system_rhs.l2_norm() < 100)
                newton_relaxation = newton_relaxation * 1.2;

            if (newton_relaxation < 0.05)
                newton_relaxation = 0.05;

            if (newton_relaxation > 1.0)
                newton_relaxation = 1.0;
        }

        pcout << "newton_relaxation = " << newton_relaxation << "\n";

        distributed_solution.add(newton_relaxation, distributed_newton_update);
        ghosted_solution = distributed_solution;
    }

    template <int dim>
    double Friction_adaptivity_everywhere<dim>::calculate_damage_error()
    {
        double damage_error = 0;
        unsigned int counter = 0;

        const auto face_damage_error = [&](const auto &cell,
                                           const unsigned int &f,
                                           const unsigned int &sf,
                                           const auto &ncell,
                                           const unsigned int &nf,
                                           const unsigned int &nsf,
                                           ScratchData<dim> &scratch_data,
                                           CopyData & /* copy_data */)
        {
            FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
            fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

            const auto &q_points = fe_iv.get_quadrature_points();
            const unsigned int n_q_points = q_points.size();

            PointHistory<dim> *quadrature_points_history = reinterpret_cast<PointHistory<dim> *>(cell->face(f)->user_pointer());

            if ((par.is_everywhere == false && cell->material_id() != ncell->material_id()) || par.is_everywhere == true)
                for (unsigned int point = 0; point < n_q_points; ++point)
                {
                    counter++;
                    damage_error += std::abs(quadrature_points_history[point].damage_init - quadrature_points_history[point].damage);
                }
        };

        MeshWorker::mesh_loop(dof_handler.begin_active(),
                              dof_handler.end(),
                              nullptr,
                              nullptr,
                              ScratchData<dim>(mapping, fe, quadrature, *face_quadrature),
                              CopyData(),
                              MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                              nullptr,
                              face_damage_error);
        counter = Utilities::MPI::sum(counter, mpi_communicator);
        damage_error = Utilities::MPI::sum(damage_error, mpi_communicator);
        damage_error = damage_error / counter;

        return damage_error;
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::refine_mesh()
    {
        TimerOutput::Scope t(computing_timer, "Refine");
        double counter = 0;
        const auto face_flagging = [&](const auto &cell,
                                       const unsigned int &f,
                                       const unsigned int &sf,
                                       const auto &ncell,
                                       const unsigned int &nf,
                                       const unsigned int &nsf,
                                       ScratchData<dim> &scratch_data,
                                       CopyData & /* copy_data */)
        {
            FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
            fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

            double nelasticity = 0, npoisson = 0, nlambda = 0, nmu = 0;
            ENLM(ncell->material_id(), nelasticity, npoisson, nlambda, nmu);

            const auto &q_points = fe_iv.get_quadrature_points();
            const unsigned int n_q_points = q_points.size();
            const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

            const FEFaceValuesBase<dim> &fe_fv1 = fe_iv.get_fe_face_values(0);
            const FEFaceValuesBase<dim> &fe_fv2 = fe_iv.get_fe_face_values(1);
            const FEValuesExtractors::Vector displacements(0);
            std::vector<Tensor<2, dim>> grads_1(n_q_points);
            fe_fv1[displacements].get_function_gradients(ghosted_solution, grads_1);
            std::vector<Tensor<2, dim>> grads_2(n_q_points);
            fe_fv2[displacements].get_function_gradients(ghosted_solution, grads_2);

            Tensor<1, dim> sig_c;
            if (cell->material_id() == ncell->material_id())
            {
                sig_c[0] = par.sig_cs;
                sig_c[1] = par.sig_cn;
            }
            else
            {
                sig_c[0] = par.sig_ics;
                sig_c[1] = par.sig_icn;
            }

            Tensor<1, dim> tr, traction_eff;
            if ((par.is_everywhere == false && cell->material_id() != ncell->material_id()) || par.is_everywhere == true)
                for (unsigned int point = 0; point < n_q_points; ++point)
                {
                    Tensor<2, dim> strain1 = 0.5 * (grads_1[point] + transpose(grads_1[point]));
                    Tensor<2, dim> strain2 = 0.5 * (grads_2[point] + transpose(grads_2[point]));
                    Tensor<2, dim> stress1 = (lambda * trace(strain1) * Identity + 2 * mu * strain1);
                    Tensor<2, dim> stress2 = (nlambda * trace(strain2) * Identity + 2 * nmu * strain2);
                    Tensor<2, dim> old_stress_avg = stress1 / 2 + stress2 / 2;

                    Tensor<1, dim> tangential = cross_product_2d(normals[point]);

                    tr[0] = tangential * (old_stress_avg * normals[point]);
                    tr[1] = normals[point] * (old_stress_avg * normals[point]);
                    traction_eff[0] = tr[0];
                    traction_eff[1] = heaviside(tr[1]) * tr[1];

                    if (std::pow(traction_eff[0] / sig_c[0], 2) + std::pow(traction_eff[1] / sig_c[1], 2) > std::pow(par.threshold_stress_factor, 2))
                    {
                        cell->set_refine_flag();
                        ncell->set_refine_flag();
                    }

                    if (cell->refine_flag_set() && cell->level() == par.max_cell_level)
                        cell->clear_refine_flag();

                    if (ncell->refine_flag_set() && ncell->level() == par.max_cell_level)
                        ncell->clear_refine_flag();

                    if (par.is_everywhere == false)
                        if (cell->material_id() == ncell->material_id())
                        {
                            cell->clear_refine_flag();
                            ncell->clear_refine_flag();
                        }

                    if (par.geometry == "EN_10")
                        if (par.is_everywhere == false)
                            if (cell->face(f)->center()[1] > 55 / 2)
                            {
                                cell->clear_refine_flag();
                                ncell->clear_refine_flag();
                            }
                }
            if (cell->refine_flag_set() || ncell->refine_flag_set())
                counter++;
        };

        MeshWorker::mesh_loop(dof_handler.begin_active(),
                              dof_handler.end(),
                              nullptr,
                              nullptr,
                              ScratchData<dim>(mapping, fe, quadrature, *face_quadrature),
                              CopyData(),
                              MeshWorker::assemble_own_interior_faces_both | MeshWorker::assemble_ghost_faces_both,
                              nullptr,
                              face_flagging);

        counter = Utilities::MPI::sum(counter, mpi_communicator);

        if (counter != 0)
        {
            triangulation.prepare_coarsening_and_refinement();
            parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(dof_handler);
            solution_transfer.prepare_for_coarsening_and_refinement(distributed_solution);

            // create and init history vectors for trasnfer
            FE_DGQ<dim> history_fe(1);
            DoFHandler<dim> history_dof_handler(triangulation);
            history_dof_handler.distribute_dofs(history_fe);

            // $ you might not need all of these vecots
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

            // create and init more history stuff for trasnfer and projection
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

            // project from local data on q_points to the history vectors defined above
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

            // create solution transfer objects
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

            // reset q point data on the new mesh.
            setup_quadrature_point_history();

            // transfer solution  new mesh
            dof_handler.distribute_dofs(fe);
            locally_owned_dofs = dof_handler.locally_owned_dofs();
            locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
            ghosted_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
            distributed_solution.reinit(locally_owned_dofs, mpi_communicator);
            solution_transfer.interpolate(distributed_solution);
            ghosted_solution = distributed_solution;

            // transfer history data to new mesh
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

            // interpolate form solution vectors to q point.

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

                    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                    {

                        fe_fv_damage.reinit(cell, face);
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
        }
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::output_results(const unsigned int cycle) const
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

            data_out.write_vtu_with_pvtu_record(par.output_directory + "/", "solution", cycle, mpi_communicator, 3, 0);
        }
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::reaction_and_traction(const types::boundary_id &boundary_id, Tensor<1, dim> &reaction_stress,
                                                                    Tensor<1, dim> &interface_jump, Tensor<1, dim> &interface_stress)
    {
        TimerOutput::Scope t(computing_timer, "other");
        reaction_stress = 0.0;
        interface_stress = 0;
        interface_jump = 0;
        unsigned int counter = 0;
        const auto boundary_worker = [&](const auto &cell,
                                         const unsigned int &face_no,
                                         ScratchData<dim> &scratch_data,
                                         CopyData &copy_data)
        {
            (void)copy_data;
            scratch_data.fe_interface_values.reinit(cell, face_no);
            const FEFaceValuesBase<dim> &fe_fv =
                scratch_data.fe_interface_values.get_fe_face_values(0);

            const auto &q_points = fe_fv.get_quadrature_points();
            const unsigned int n_q_points = q_points.size();
            const std::vector<double> &JxW = fe_fv.get_JxW_values();
            const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();
            const FEValuesExtractors::Vector displacements(0);

            std::vector<Tensor<2, dim>> gradu(n_q_points);
            fe_fv[displacements].get_function_gradients(ghosted_solution, gradu);

            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

            Tensor<1, dim> tangential;
            Tensor<1, dim> local_reaction;

            if ((cell->face(face_no)->boundary_id() == boundary_id))
                //  if (cell->face(face_no)->boundary_id() == boundary_id && (std::abs(cell->face(face_no)->center()[0] - 70) < rad / 2))
                // if (cell->face(face_no)->boundary_id() == 0 && std::abs(cell->face(face_no)->center()[1]) < 55 / 2 && cell->material_id() == 1)
                for (unsigned int point = 0; point < n_q_points; ++point)
                {

                    Tensor<2, dim> strain = 0.5 * (gradu[point] + transpose(gradu[point]));
                    Tensor<2, dim> stress = (lambda * trace(strain) * Identity + 2 * mu * strain);

                    tangential = cross_product_2d(normals[point]);
                    reaction_stress[0] += tangential * (stress * normals[point]) * JxW[point];
                    reaction_stress[1] += normals[point] * (stress * normals[point]) * JxW[point];
                }
        };

        const auto face_worker = [&](const auto &cell,
                                     const unsigned int &f,
                                     const unsigned int &sf,
                                     const auto &ncell,
                                     const unsigned int &nf,
                                     const unsigned int &nsf,
                                     ScratchData<dim> &scratch_data,
                                     CopyData &copy_data)
        {
            (void)copy_data;
            FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
            fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

            const FEFaceValuesBase<dim> &fe_fv1 = fe_iv.get_fe_face_values(0);
            const FEFaceValuesBase<dim> &fe_fv2 = fe_iv.get_fe_face_values(1);

            const auto &q_points = fe_iv.get_quadrature_points();
            const unsigned int n_q_points = q_points.size();

            const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
            const std::vector<double> &JxW = fe_iv.get_JxW_values();
            double elasticity = 0, poisson = 0, lambda = 0, mu = 0;
            ENLM(cell->material_id(), elasticity, poisson, lambda, mu);

            double nelasticity = 0, npoisson = 0, nlambda = 0, nmu = 0;
            ENLM(ncell->material_id(), nelasticity, npoisson, nlambda, nmu);

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

                    interface_stress[0] += tangential * (stress * normals[point]) * JxW[point];
                    interface_stress[1] += normals[point] * (stress * normals[point]) * JxW[point];
                    interface_jump[0] += jumpu[point] * normals[point];
                    interface_jump[1] += jumpu[point] * tangential;
                }
        };

        ScratchData<dim> scratch_data(mapping, fe, quadrature, QGauss<dim - 1>(5));
        CopyData copy_data;

        MeshWorker::mesh_loop(dof_handler.begin_active(),
                              dof_handler.end(),
                              nullptr,
                              nullptr,
                              scratch_data,
                              copy_data,
                              MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once | MeshWorker::assemble_ghost_faces_once,
                              boundary_worker,
                              face_worker);

        counter = Utilities::MPI::sum(counter, mpi_communicator);
        reaction_stress = Utilities::MPI::sum(reaction_stress, mpi_communicator);
        interface_jump = Utilities::MPI::sum(interface_jump, mpi_communicator);
        interface_stress = Utilities::MPI::sum(interface_stress, mpi_communicator);
        interface_jump = interface_jump / counter;
    }

    template <int dim>
    void Friction_adaptivity_everywhere<dim>::run()
    {
        if (par.geometry == "reactangle")
        {
            Point<dim> P1, P2(1, 2);
            std::vector<unsigned int> repetitions{1, 2};
            GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
            triangulation.refine_global(par.initial_refinement);

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
        else if (par.geometry == "half_matrix" || par.geometry == "matrix" || par.geometry == "hole")
        {
            GridIn<2> gridin;
            gridin.attach_triangulation(triangulation);
            std::ifstream f;
            if (par.geometry == "half_matrix")
                f.open("half_matrix.msh");
            else if (par.geometry == "matrix")
                f.open("matrix1.msh");
            else if (par.geometry == "hole")
                f.open("hole.msh");
            gridin.read_msh(f);
            triangulation.refine_global(par.initial_refinement);
        }
        else if (par.geometry == "ENF")
        {
            Point<dim> P1(0, 0), P2(140, 6);
            std::vector<unsigned int> repetitions{70, 3};
            GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
            triangulation.refine_global(par.initial_refinement);

            Point<dim> p;
            for (const auto &cell : triangulation.active_cell_iterators())
            {
                p = cell->center();
                if (p[1] > 3)
                    cell->set_material_id(1);
                else
                    cell->set_material_id(2);
            }
        }
        else if (par.geometry == "EN_10")
        {
            double h = 55.0, l = 240.0, b = 10.0;
            Point<dim> P1(0, 0), P2(l, h + h / 2 + b);
            std::vector<unsigned int> repetitions{2, 1};
            GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, P1, P2, true);
            triangulation.refine_global(par.initial_refinement);

            Point<dim> p;
            for (const auto &cell : triangulation.active_cell_iterators())
            {
                p = cell->center();
                if (std::abs(p[1] - h / 2 - b / 2) > b)
                    cell->set_material_id(1);
                else
                {
                    cell->set_material_id(2);
                    //  cell->set_refine_flag();
                }
            }
            // triangulation.execute_coarsening_and_refinement();
        }
        else
            throw std::runtime_error("no input mesh");

        std::ofstream out(par.output_directory + "/grid-1.vtu");
        GridOut grid_out;
        grid_out.write_vtu(triangulation, out);
        pcout << " written to "
              << "grid-1.vtu" << std::endl
              << std::endl;

        if (par.is_distorted)
        {
            GridTools::transform(&grid_y_transform<dim>, triangulation);

            // GridTools::distort_random(0.1, triangulation, true, 2882);
            // triangulation.refine_global(1);
            // GridTools::rotate(-M_PI/2, triangulation);
        }

        setup_quadrature_point_history();

        double error;
        unsigned int non_lin = 0;
        unsigned int max_nonlin_iter = par.max_nonlin_iter;
        unsigned int max_nonlin_output = 0;
        unsigned int max_non_lin_at_cycle = 0;
        unsigned int max_external_iterations_output = 0;
        unsigned int max_external_iterations_at_cycle = 0;

        std::ofstream forces(par.output_directory + "/" + par.file_name + ".txt");

        for (int cycle = 0; cycle < par.cycles; ++cycle)
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

            if (par.with_adaptivity == true)
                if (cycle != 0)
                    refine_mesh();

            unsigned int max_external_iterations = par.max_external_iterations;
            unsigned int external_iterations = 0;

            bool doing_external_iterations = external_iterations < max_external_iterations;

            while ((error > par.error && non_lin < max_nonlin_iter) || doing_external_iterations)
            {
                if (external_iterations >= max_external_iterations)
                    doing_external_iterations = false;

                double prev_error = system_rhs.l2_norm();
                if (non_lin == 0)
                    prev_error = 1e+20;

                if (cycle == 0 && non_lin == 0)
                    setup_system(true);
                else
                    setup_system(false);

                max_nonlin_output = std::max(max_nonlin_output, non_lin);
                max_external_iterations_output = std::max(max_external_iterations_output, external_iterations);

                if (max_nonlin_output == non_lin)
                    max_non_lin_at_cycle = cycle;

                if (max_external_iterations_output == external_iterations)
                    max_external_iterations_at_cycle = cycle;

                if (par.ignore_non_convergence == false)
                    if (max_nonlin_output == max_nonlin_iter - 1)
                        throw std::runtime_error("max non_lin iterations reached b4 convergence");

                non_lin++;
                assemble_system(cycle, non_lin);
                error = system_rhs.l2_norm();

                pcout << "  system_rhs l2_norm " << system_rhs.l2_norm() << std::endl;
                pcout << "  damage_error " << calculate_damage_error() << "\n";

                if ((error > par.error && !doing_external_iterations) || (error > par.external_iterations_error && doing_external_iterations))
                {
                    pcout << "  external_iterations = " << external_iterations << std::endl;
                    pcout << "  non lin = " << non_lin << std::endl;

                    solve(prev_error);
                    pcout << "  update l2_norm " << distributed_newton_update.l2_norm() << std::endl;
                }
                else if (doing_external_iterations)
                {
                    if (calculate_damage_error() < par.damage_error)
                    {
                        doing_external_iterations = false;
                    }
                    else
                    {
                        external_iterations++;
                        non_lin = 0;
                    }
                }
                pcout << "max external_iterationses " << max_external_iterations_output << " at_cycle " << max_external_iterations_at_cycle << "\n";
                pcout << "max nonlin iterations " << max_nonlin_output << "  at_cycle " << max_non_lin_at_cycle << "\n------- \n ";

                if (error > 1e+30)
                    throw std::runtime_error("Divergence in the solution");

                if (cycle % par.output_frequency == 0 || cycle == 0)
                {
                    TimerOutput::Scope t(computing_timer, "output");
                    output_results(cycle);
                }
            }

            Tensor<1, dim> reaction, interface_jump, interface_traction;
            reaction_and_traction(3, reaction, interface_jump, interface_traction);
            pcout << "   reaction force = " << reaction[0] << "\t" << reaction[1] << "\n";
            pcout << std::endl;

            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            {
                forces.open(par.output_directory + "/" + par.file_name + ".txt", std::ios_base::app);
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
                       << reaction[0] << "\t" << reaction[1] << "\t"
                       << interface_jump[0] << "\t" << interface_jump[1] << "\t"
                       << interface_traction[0] << "\t" << interface_traction[1] << "\n";

                forces.close();
            }
        }
        computing_timer.print_summary();
        computing_timer.reset();
    }
} // namespace namespace Friction_adaptivity_everywhere

#endif //  Friction_adaptivity_everywhere_H