/* ---------------------------------------------------------------------
 *
 *   CODE FOR PHASE FIELD SIMULATION
 *   
*  TRANSFER THE SOLUTION See step 15
 * SNES FOR DAMAGE See step 15
 * SAVE THE POST PROCESSING (deformations, tensions)
 https://groups.google.com/g/dealii/c/uME5kwaJIEU/m/M-hMBIUHCAAJ 
https://groups.google.com/g/dealii/c/ymDfJIumYw4/m/hITK6VhXCgAJ

In the future, build a ConstitutiveLaw class such as
https://www.dealii.org/developer/doxygen/deal.II/code_gallery_goal_oriented_elastoplasticity.html

// NOTE
Both fields share the same quadrature formula

 * ---------------------------------------------------------------------
 *
 */
 
#ifndef NRPHASE
#define NRPHASE


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>    
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/physics/transformations.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <chrono>



namespace phasefieldNR
{   

   const bool plane_strain = true;
   const unsigned int model_type = 1; // 1 for full elastic energy degradation // 2 for masonry like formulation // 3 Deviatoric // 4 full-deviatoric
   const unsigned int damage_type = 2; // 1 for quadratic // 2 for linear



    using namespace dealii;

    template <int dim>
    struct PointHistory_stresses
    {
        SymmetricTensor<2, dim> old_stress;
    };

    template <int dim>
    SymmetricTensor<4, dim> copy_tensor_to_symmetric_tensor(const  Tensor<4, dim> &to_be_copied)
    {
    SymmetricTensor<4, dim> tmp;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            tmp[i][j][k][l] = to_be_copied[i][j][k][l];
    return tmp;
    }

    template <int dim>
    SymmetricTensor<2, dim> copy_tensor_to_symmetric_tensor(const  Tensor<2, dim> &to_be_copied)
    {
    SymmetricTensor<2, dim> tmp;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
            tmp[i][j] = to_be_copied[i][j];
    return tmp;
    }
  

    template <int dim>
    struct PointHistory
    {
        SymmetricTensor<2, dim>  old_stress;
        SymmetricTensor<2, dim>  old_strain;
        double old_phase_field; 
    };

    template <int dim>
    SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda, const double mu)
    {
        SymmetricTensor<4, dim> tmp;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                        tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
                            ((i == l) && (j == k) ? mu : 0.0) +
                            ((i == j) && (k == l) ? lambda : 0.0));
        return tmp;
    }


    template <int dim>
    SymmetricTensor<4, dim> get_stress_strain_tensor_dev(const double mu)
    {
        SymmetricTensor<4, dim> tmp;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                        for (unsigned int l = 0; l < dim; ++l)
                            tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
                              ((i == l) && (j == k) ? mu : 0.0) +
                              ((i == j) && (k == l) ? -2./3.*mu : 0.0));
        return tmp;
    }

    template <int dim>
    SymmetricTensor<4, dim> get_stress_strain_tensor_sph(const double lambda,
                                                  const double mu)
    {
        double bulk = lambda + 2.*mu/3.;
        SymmetricTensor<4, dim> tmp;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                        tmp[i][j][k][l] = ((i == j) && (k == l) ? bulk : 0.0);
        return tmp;
    }


    template <int dim>
    SymmetricTensor<4,dim> Identity ()
    {
        SymmetricTensor<4,dim> tmp;
        for (unsigned int i=0;i<dim;++i)
            for (unsigned int j=0;j<dim;++j)
                for (unsigned int k=0;k<dim;++k)
                    for (unsigned int l=0;l<dim;++l)
                        {
                        double a=0,b=0;
                        if (i==k & j==l)
                            a=0.5;
                        if (i==l & j==k)
                            b=0.5;
                        tmp[i][j][k][l]=a+b;
                        }
        return tmp;
    }

    template <int dim>
    inline
    Point<dim> get_tractions (const SymmetricTensor<2,dim> &stress, const Point<dim> &normal)
    {
        Assert (stress.size() == dim, ExcInternalError());

        Point<dim> traction;

        for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
	            traction[i] += stress[i][j] * normal[j];

        return traction;
  } 
 
    template <int dim>
    inline SymmetricTensor<2, dim> get_strain(const FEValues<dim>& fe_values,
        const unsigned int   shape_func,
        const unsigned int   q_point)
    {
        SymmetricTensor<2, dim> tmp;
        for (unsigned int i = 0; i < dim; ++i)
            tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = i + 1; j < dim; ++j)
                tmp[i][j] =
                (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
                    fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
                2;
        return tmp;
    }


    template <int dim>
    inline SymmetricTensor<2, dim>
        get_strain(const std::vector<Tensor<1, dim>>& grad)
    {
        Assert(grad.size() == dim, ExcInternalError());
        SymmetricTensor<2, dim> strain;
        for (unsigned int i = 0; i < dim; ++i)
            strain[i][i] = grad[i][i];
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = i + 1; j < dim; ++j)
                strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
        return strain;
    }
    Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>>& grad_u)
    {
        const double curl = (grad_u[1][0] - grad_u[0][1]);
        const double angle = std::atan(curl);
        return Physics::Transformations::Rotations::rotation_matrix_2d(-angle);
    }
    Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>>& grad_u)
    {
        const Point<3> curl(grad_u[2][1] - grad_u[1][2],
            grad_u[0][2] - grad_u[2][0],
            grad_u[1][0] - grad_u[0][1]);
        const double tan_angle = std::sqrt(curl * curl);
        const double angle = std::atan(tan_angle);
        if (std::abs(angle) < 1e-9)
        {
            static const double rotation[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
            static const Tensor<2, 3> rot(rotation);
            return rot;
        }
        const Point<3> axis = curl / tan_angle;
        return Physics::Transformations::Rotations::rotation_matrix_3d(axis,
            -angle);
    }
 
    // ------------------------------------------------------------------------------------------
    // DAMAGE FUNCTIONS
    // ------------------------------------------------------------------------------------------
    inline double g_alpha (double &alpha)
    {
        return (1-alpha)*(1-alpha);
    }

    inline double g_prime_alpha (double &alpha)
    {
        return -2.+2.*alpha;
    }

    inline double g_second_alpha ()
    {
        return 2.;
    }

    inline double w_alpha (double &alpha)
    {   
        double tmp;
        if (damage_type == 1)
            tmp =  alpha*alpha;
        else if (damage_type == 2)
            tmp = alpha;    
        return tmp;
    }

    inline double w_prime_alpha (double &alpha)
    {   
        double tmp;
        if (damage_type == 1)
            tmp = 2.*alpha;
        else if (damage_type == 2)
            tmp = 1.;    
        return tmp;
    }

    inline double w_second_alpha ()
    {   
        double tmp;
        if (damage_type == 1)
            tmp =  2.;
        else if (damage_type == 2)
            tmp = 0.;    
        return tmp;
    }

    inline double heaviside_neg (double &x)
    {
        double result = 0;
        if (x < 0) result = 1;
        return result;
    }

    // ------------------------------------------------------------------------------------------
    // Class definition
    // ------------------------------------------------------------------------------------------
    template <int dim>
    class TopLevel
    {
    public:
        TopLevel();
        ~TopLevel();
        void run();

    private:
        void create_coarse_grid();
        void refine_initial_grid();
        void set_data ();
        void set_initial_values_alpha (double &a);

        void evaluate_eps_c (SymmetricTensor<2,dim> &eps_pos, SymmetricTensor<2,dim> &eps_neg, const SymmetricTensor<2,dim> &eps, const double nu);
        void evaluate_stress (SymmetricTensor<2,dim> &sigma_pos, SymmetricTensor<2,dim> &sigma_neg, const SymmetricTensor<2,dim> &eps, const SymmetricTensor<2,dim> &sigma, const double EE, const double nu);
        void evaluate_stiffness_matrix (SymmetricTensor<4,dim> &C_plus, SymmetricTensor<4,dim> &C_minus, const SymmetricTensor<2,dim> &eps, const SymmetricTensor<4,dim> &C, const double EE, const double nu);   
        void stiffness_matrix(SymmetricTensor<4,dim> &C_plus, SymmetricTensor<4,dim> &C_minus, const SymmetricTensor<2,dim> &eps, const SymmetricTensor<4,dim> &C, const double EE, const double nu);         

        void setup_system_elas();
        void assemble_system_elas(const bool flag_iter);
        void solve_linear_problem(const bool flag_elastic_iter);

        void setup_system_alpha();
        void assemble_system_alpha(const bool flag_iter_alpha);
        unsigned int solve_alpha(const bool flag_iter_alpha);
        void check_alpha(Vector<double> &alpha);
        
        void output_results_elas(const unsigned int cycle) const;
        void output_results_alpha(const unsigned int cycle) const;     
        void output_stress(const unsigned int cycle);   
        void reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id);
        void energy (double &bulk_energy, double &surface_energy);

        void solve_timestep();
        void do_initial_timestep();
        void do_timestep();

        parallel::shared::Triangulation<dim> triangulation;
        FESystem<dim> fe;
        FESystem<dim> fe_vec;
        DoFHandler<dim> dof_handler;
        DoFHandler<dim> dof_handler_vec;

        AffineConstraints<double> hanging_node_constraints;
        AffineConstraints<double> hanging_node_constraints_vec;  

        const QGauss<dim> quadrature_formula;
        std::vector<PointHistory<dim>> quadrature_point_history;
        std::vector<PointHistory_stresses<dim>> quadrature_point_history_stress;

        // Elasticity Matrices - Tensor - Vectors
        PETScWrappers::MPI::SparseMatrix system_matrix_elas;
        PETScWrappers::MPI::Vector system_rhs_elas;

        static const SymmetricTensor<4,dim> stress_strain_tensor;
        static const SymmetricTensor<4,dim> null_s_fourth_order_tensor;
        static const SymmetricTensor<4,dim> stress_strain_tensor_dev;
        static const SymmetricTensor<4,dim> stress_strain_tensor_sph;
        SymmetricTensor<4, dim> C_Pos, C_Neg;

        Vector<double> solution_u;
        Vector<double> newton_update_u;

        // Damage Matricies - Tensor - Vectors
        PETScWrappers::MPI::SparseMatrix system_matrix_alpha;
        PETScWrappers::MPI::Vector system_rhs_alpha;

        Vector<double> solution_alpha;
        Vector<double> solution_alpha_previous_step;
        Vector<double> newton_update_alpha;

        PETScWrappers::MPI::Vector alpha_vec;

        // Timestep parameters
        double  real_time;
        double  real_time_step;
        double  real_time_final;
        unsigned int real_timestep_no = 0;

        // MPI communicators
        MPI_Comm mpi_communicator;
        const unsigned int n_mpi_processes;
        const unsigned int this_mpi_process;


        // MPI Output on shell 
        ConditionalOStream pcout;

        // Indexset
        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        IndexSet locally_owned_dofs_vec;
        IndexSet locally_relevant_dofs_vec;

        // Elasticity / Damage parameters
        double Poisson;
        double Young;

        double k_res; 
        double alpha0;
        double ell;
        double Gc;
        double c_w;
        double gamma_penalty;

        double error_tol;
    };

    template <int dim>
    void TopLevel<dim>::set_data ()
    {
    // Inizializzazione parametri.
        Poisson = 0.2;
        Young = 70000;

        k_res=1.e-6; 
        alpha0=0.; 
        ell=6;
        Gc=0.5;
        gamma_penalty = 1.e4;
        
        error_tol=1.e-5;
        
        if (damage_type == 1)
            c_w = 2.;
        else if (damage_type == 2)    
            c_w = 8./3.;

        real_time = 0.;
        real_time_final = 30;
        real_time_step = 1; 
    }

    // ------------------------------------------------------------------------------------------
    // SETUP INITIAL / BC FUNCTIONS
    // ------------------------------------------------------------------------------------------
    // Alpha setup value
    template <int dim>
    void TopLevel<dim>::set_initial_values_alpha (double &alpha0)
    {  
        solution_alpha_previous_step=alpha0;     
        solution_alpha=alpha0;        
    }

    // Elasticity - Body force class
    template <int dim>
    class BodyForce : public Function<dim>
    {
    public:
        BodyForce();
        virtual void vector_value(const Point<dim>& p,
            Vector<double>& values) const override;
        virtual void
            vector_value_list(const std::vector<Point<dim>>& points,
                std::vector<Vector<double>>& value_list) const override;
    };

    template <int dim>
    BodyForce<dim>::BodyForce()
        : Function<dim>(dim)
    {}

    template <int dim>
    inline void BodyForce<dim>::vector_value(const Point<dim>& /*p*/,
        Vector<double>& values) const
    {
        Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
        const double g = 9.81;
        const double rho = 7700;
        values = 0;
        values(dim - 1) = -rho * g * 0;
    }

    template <int dim>
    void BodyForce<dim>::vector_value_list(
        const std::vector<Point<dim>>& points,
        std::vector<Vector<double>>& value_list) const
    {
        const unsigned int n_points = points.size();
        Assert(value_list.size() == n_points,
            ExcDimensionMismatch(value_list.size(), n_points));
        for (unsigned int p = 0; p < n_points; ++p)
            BodyForce<dim>::vector_value(points[p], value_list[p]);
    }

    template <int dim>
    class IncrementalBoundaryValues : public Function<dim>
    {
    public:
        IncrementalBoundaryValues(const double present_time,
            const double present_timestep);
        virtual void vector_value(const Point<dim>& p,
            Vector<double>& values) const override;
        virtual void
            vector_value_list(const std::vector<Point<dim>>& points,
                std::vector<Vector<double>>& value_list) const override;
    private:
        const double velocity;
        const double present_time;
        const double present_timestep;
    };

    template <int dim>
    IncrementalBoundaryValues<dim>::IncrementalBoundaryValues(
        const double present_time,
        const double present_timestep)
        : Function<dim>(dim)
        , velocity(0.01)
        , present_time(present_time)
        , present_timestep(present_timestep)
    {}

    template <int dim>
    void
        IncrementalBoundaryValues<dim>::vector_value(const Point<dim>& /*p*/,
            Vector<double>& values) const
    {
        Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
        values = 0;
        values(1) = + present_time * velocity;
    }

    template <int dim>
    void IncrementalBoundaryValues<dim>::vector_value_list(
        const std::vector<Point<dim>>& points,
        std::vector<Vector<double>>& value_list) const
    {
        const unsigned int n_points = points.size();
        Assert(value_list.size() == n_points,
            ExcDimensionMismatch(value_list.size(), n_points));
        for (unsigned int p = 0; p < n_points; ++p)
            IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]);
    }

    // ------------------------------------------------------------------------------------------
    // ELASTICITY FUNCTIONS
    // ------------------------------------------------------------------------------------------
    // Tensors evaluations
    template <int dim>
    const SymmetricTensor<4,dim>
    TopLevel<dim>::stress_strain_tensor = get_stress_strain_tensor<dim> (/*lambda = */ (70000*0.2)/((1+0.2)*(1-2*0.2)),
                                                                         /*mu     = */ 70000/(2.*(1+0.2)));

    template <int dim>
    const SymmetricTensor<4,dim>
    TopLevel<dim>::stress_strain_tensor_dev = get_stress_strain_tensor_dev<dim> (/*mu  = */ 70000/(2.*(1+0.2)));

    template <int dim>
    const SymmetricTensor<4,dim>
    TopLevel<dim>::stress_strain_tensor_sph = get_stress_strain_tensor_sph<dim> (/*lambda = */ (70000*0.2)/((1+0.2)*(1-2*0.2)),
                                                                                 /*mu     = */ 70000/(2.*(1+0.2)));

    template <int dim>
    const SymmetricTensor<4,dim>
    TopLevel<dim>::null_s_fourth_order_tensor = get_stress_strain_tensor<dim> (0.,0.);

    template <int dim>
    void TopLevel<dim>::evaluate_eps_c(SymmetricTensor<2,dim> &eps_pos, SymmetricTensor<2,dim> &eps_neg, const SymmetricTensor<2,dim> &eps, const double nu)
    {
        double aa;
        double aa_cap;

        if (plane_strain)
            aa = nu / (1-2*nu);
        else
            aa = nu / (1-nu);
        
        aa_cap = aa/(1+aa);

        double eig_1, eig_2;
        Tensor<1,dim> n1; 
        Tensor<1,dim> n2; 

        const std::array<std::pair<double, Tensor<1, dim, double>>, dim>  eigen_eps = eigenvectors(eps);
        eig_1 = eigen_eps[0].first;
        eig_2 = eigen_eps[1].first;
        n1 = eigen_eps[0].second;
        n2 = eigen_eps[1].second;

        if (eig_2 >= 0)  // Prima condizione
        {
            eps_pos = eps;
            eps_neg = 0.;
        }
        else if ( (1+aa)*eig_1+aa*eig_2 >= 0)    //Seconda condizione
        {
            Tensor<2,dim> M1, M2;
            M1 = outer_product(n1,n1);
            M2 = outer_product(n2,n2);
            eps_pos = symmetrize((eig_1 + aa_cap*eig_2)*M1);    
            eps_neg = eps-eps_pos;
        }
        else      // Terza condizione
        {
            eps_pos = 0.;
            eps_neg = eps;
        }
    }    

    template <int dim>
    void TopLevel<dim>::evaluate_stress (SymmetricTensor<2,dim> &sigma_pos, SymmetricTensor<2,dim> &sigma_neg, const SymmetricTensor<2,dim> &eps, const SymmetricTensor<2,dim> &sigma, const double EE, const double nu)
    {
        double aa;
        double aa_cap;

        double lmbda = (EE*nu)/((1+nu)*(1-2*nu));
        double mu = EE/(2*(1+nu));

        if (plane_strain)
            aa = nu / (1-2*nu);
        else
            aa = nu / (1-nu);
        
        aa_cap = aa/(1+aa);

        double eig_1, eig_2;
        Tensor<1,dim> n1; 
        Tensor<1,dim> n2; 

        const std::array<std::pair<double, Tensor<1, dim, double>>, dim> eigen_eps = eigenvectors(eps);
        eig_1 = eigen_eps[0].first;
        eig_2 = eigen_eps[1].first;
        n1 = eigen_eps[0].second;
        n2 = eigen_eps[1].second;

        if (eig_2 >= 0)  // Prima condizione
        {
            sigma_pos = sigma;
            sigma_neg = 0.;
        }
        else if ( (1+aa)*eig_1+aa*eig_2 >= 0)    //Seconda condizione
        {
            Tensor<2,dim> M1, M2;
            M1 = outer_product(n1,n1);
            M2 = outer_product(n2,n2);
            sigma_pos = symmetrize((lmbda +2*mu)*(eig_1 + aa_cap*eig_2)*(M1 + aa_cap*M2));
            sigma_neg = symmetrize(((2.0*mu*(1.0 + 2.0*aa*(1.0 + aa)) + lmbda)/((1.0 + aa)*(1.0 + aa)))*eig_2*M2);   
        }
        else      // Terza condizione
        {
            sigma_pos = 0.;
            sigma_neg = sigma;
        }
    } 

    template <int dim>
    void TopLevel<dim>::evaluate_stiffness_matrix(SymmetricTensor<4,dim> &C_plus, SymmetricTensor<4,dim> &C_minus, const SymmetricTensor<2,dim> &eps, const SymmetricTensor<4,dim> &C, const double EE, const double nu) 
    {
    // Switch parametro alpha per plane strain/stress
        double aa;
        double aa_cap;

        double lmbda = (EE*nu)/((1+nu)*(1-2*nu));
        double mu = EE/(2*(1+nu));

        if (plane_strain)
            aa = nu / (1-2*nu);
        else
            aa = nu / (1-nu);

        aa_cap = aa/(1+aa);

        double eig_1, eig_2;
        Tensor<1,dim> n1; 
        Tensor<1,dim> n2; 

        const std::array<std::pair<double, Tensor<1, dim, double>>, dim>  eigen_eps = eigenvectors(eps);
        eig_1 = eigen_eps[0].first;
        eig_2 = eigen_eps[1].first;
        n1 = eigen_eps[0].second;
        n2 = eigen_eps[1].second;    

        if (eig_2 >= 0)    // Prima condizione
        {
            C_plus = C;
            C_minus = null_s_fourth_order_tensor;
        }  
        else if ( ((1+aa)*eig_1+aa*eig_2) >= 0)    //Seconda condizione
        {
            Tensor<4, dim>  C_plus_temp, C_minus_temp;
            Tensor<2,dim> M1, M2, M12, M21;
            M1 = outer_product(n1,n1);
            M2 = outer_product(n2,n2);
            M12 = outer_product(n1, n2);
            M21 = outer_product(n2, n1);

            Tensor<4,dim> hh, S1212, S1221, S2112, S2121, TT1, TT2;
            S1212 = outer_product(M12, M12);
            S1221 = outer_product(M12, M21);
            S2112 = outer_product(M21, M12);
            S2121 = outer_product(M21, M21);

            TT1 = (1/(2*(eig_1-eig_2)))*(S1212 + S1221 + S2112 + S2121);
            TT2 = -TT1;
        
            double coeff_pos;
            coeff_pos = lmbda + 2*mu;
            Tensor<4, dim> H = outer_product(M1 + aa_cap*M2, M1 + aa_cap*M2);
            C_plus_temp =  coeff_pos * (H + (eig_1+aa_cap*eig_2)*(TT1 + aa_cap*TT2));

            double coeff_neg;
            coeff_neg = ((2.0*mu*(1.0 + 2.0*aa*(1.0 + aa)) + lmbda)/((1.0 + aa)*(1.0 + aa)));
            C_minus_temp = coeff_neg*(outer_product(M2,M2) + eig_2*TT2);

            C_plus = copy_tensor_to_symmetric_tensor(C_plus_temp);
            C_minus = copy_tensor_to_symmetric_tensor(C_minus_temp);
        }
        else // Terza condizione
        {
            C_plus = null_s_fourth_order_tensor;
            C_minus = C;
        }
    }  

    template <int dim>
    void TopLevel<dim>::stiffness_matrix(SymmetricTensor<4,dim> &C_plus, SymmetricTensor<4,dim> &C_minus, const SymmetricTensor<2,dim> &eps_u, const SymmetricTensor<4,dim> &C, const double EE, const double nu) 
    {
        if (model_type == 1)
        {
            C_plus = C;
            C_minus = null_s_fourth_order_tensor;
        }
        else if (model_type == 2)
        {
            evaluate_stiffness_matrix(C_plus, C_minus, eps_u, C, EE, nu);
        }
        else if (model_type == 3)
        {    
            C_plus = stress_strain_tensor_dev;
            C_minus = stress_strain_tensor_sph;
        }
        else if (model_type == 4)
        {    
            if (trace(eps_u)>0)
            {
                C_plus = C;
                C_minus = null_s_fourth_order_tensor;            
            }
            else 
            {
                C_plus = stress_strain_tensor_dev;
                C_minus = stress_strain_tensor_sph;
            } 
        }
    }

    // ------------------------------------------------------------------------------------------
    // DAMAGE FUNCTIONS
    // ------------------------------------------------------------------------------------------
    template <int dim>
    void TopLevel<dim>::check_alpha (Vector<double> &alpha)
    {
        for (unsigned int i=0; i<alpha.size(); ++i)
        {
            if (alpha(i) < 0.) 
                alpha(i) = 0.;
            else if (alpha(i) > 1.0) 
                alpha(i) = 1.e0;       
        }         

    }

    // ------------------------------------------------------------------------------------------
    // CONSTRUCTOR
    // ------------------------------------------------------------------------------------------
    template <int dim>
    TopLevel<dim>::TopLevel()
        : triangulation(MPI_COMM_WORLD)
        , fe(FE_Q<dim>(1), 1)
        , fe_vec(FE_Q<dim>(1), dim)
        , dof_handler(triangulation)
        , dof_handler_vec(triangulation)
        , quadrature_formula(fe.degree + 1)
        , mpi_communicator(MPI_COMM_WORLD)
        , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
        , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
        , pcout(std::cout, this_mpi_process == 0)
    {}

    // Destructor
    template <int dim>
    TopLevel<dim>::~TopLevel()
    {
        dof_handler_vec.clear();
    }

    // ------------------------------------------------------------------------------------------
    // CREATE MESH
    // ------------------------------------------------------------------------------------------
    template <int dim>
    void TopLevel<dim>::create_coarse_grid()
    {
        std::vector<unsigned int> sudd(dim);
        unsigned int sudd_x=15, sudd_y=45;
        sudd[0]=sudd_x;
        sudd[1]=sudd_y;
 
        GridGenerator::subdivided_hyper_rectangle(triangulation, sudd,
                                                 Point<dim>(0,0),
        	                                    Point<dim>(100,300),
                                                           false);  
        
        for (const auto& cell : triangulation.active_cell_iterators())
            for (const auto& face : cell->face_iterators())
                if (face->at_boundary())
                {
                    const Point<dim> face_center = face->center();
                    if (face_center[1] == 0)
                        face->set_boundary_id(1);
                    else if (face_center[1] == 300)
                        face->set_boundary_id(3);
                }
        triangulation.refine_global(0);
    }

    // ------------------------------------------------------------------------------------------
    // ELASTICITY - Setup / Assemble / Solve
    // ------------------------------------------------------------------------------------------
    // Elasticity - Initialize vectors and matrix
    template <int dim>
    void TopLevel<dim>::setup_system_elas()
    {
        dof_handler_vec.distribute_dofs(fe_vec);
        locally_owned_dofs_vec = dof_handler_vec.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler_vec, locally_relevant_dofs_vec);

        hanging_node_constraints_vec.clear();
        DoFTools::make_hanging_node_constraints(dof_handler_vec, hanging_node_constraints_vec);
        hanging_node_constraints_vec.close();

        DynamicSparsityPattern sparsity_pattern_vec(locally_relevant_dofs_vec);
        DoFTools::make_sparsity_pattern(dof_handler_vec,
                                        sparsity_pattern_vec,
                                        hanging_node_constraints_vec,
                                        /*keep constrained dofs*/ false);
        SparsityTools::distribute_sparsity_pattern(sparsity_pattern_vec,
                                                   locally_owned_dofs_vec,
                                                   mpi_communicator,
                                                   locally_relevant_dofs_vec);
        
        system_matrix_elas.reinit(locally_owned_dofs_vec,
                                  locally_owned_dofs_vec,
                                  sparsity_pattern_vec,
                                  mpi_communicator);
        system_rhs_elas.reinit(locally_owned_dofs_vec, mpi_communicator);
        solution_u.reinit(dof_handler_vec.n_dofs());
        newton_update_u.reinit(dof_handler_vec.n_dofs());
    }

    // Assemble system
    template <int dim>
    void TopLevel<dim>::assemble_system_elas(const bool flag_iter)
    {
        system_rhs_elas = 0;
        system_matrix_elas = 0;
        
        FEValues<dim> fe_values(fe_vec,
                                quadrature_formula,
                                update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

        FEValues<dim> fe_values_alpha(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                      update_quadrature_points | update_JxW_values);  
        
        const unsigned int dofs_per_cell = fe_vec.dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();
       
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     cell_rhs(dofs_per_cell);
        
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        std::vector<std::vector<Tensor<1, dim>>> previous_gradient(n_q_points,std::vector<Tensor<1,dim>>(dim));
        std::vector<double> local_solution_alpha(n_q_points);  

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_vec.begin_active(),
						                               endc = dof_handler_vec.end();

        typename DoFHandler<dim>::active_cell_iterator cell_alpha = dof_handler.begin_active();
        
        BodyForce<dim>              body_force; // Sono gia presenti nella formulazione ma sono settate pari a 0
        std::vector<Vector<double>> body_force_values(n_q_points, Vector<double>(dim));
        
        for (unsigned int index = 0; cell!=endc; ++cell, ++cell_alpha, ++index)
        {
            if (cell->is_locally_owned())
            {
                cell_matrix = 0;
                cell_rhs = 0;

                fe_values.reinit(cell);
                fe_values_alpha.reinit(cell_alpha);        

                fe_values.get_function_gradients(solution_u, previous_gradient);
                fe_values_alpha.get_function_values(solution_alpha, local_solution_alpha);
                
                for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                {
                    
                    const double g_alpha_gauss = g_alpha(local_solution_alpha[q_point]);
                    const SymmetricTensor<2,dim> eps_u = get_strain (previous_gradient[q_point]);
                    SymmetricTensor<4, dim>  C_Pos, C_Neg;

                    stiffness_matrix(C_Pos, C_Neg, eps_u, stress_strain_tensor, Young, Poisson);
                    
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                        const SymmetricTensor<2, dim>   eps_phi_i = get_strain(fe_values, i, q_point);
                        const SymmetricTensor<2, dim>   sigma_i =  (g_alpha_gauss + k_res) * C_Pos * eps_phi_i + C_Neg * eps_phi_i;    

                        for (unsigned int j = 0; j < dofs_per_cell; ++j)         
                        {
                            const SymmetricTensor<2, dim>  eps_phi_j = get_strain(fe_values, j, q_point);
                            cell_matrix(i, j) +=    sigma_i*eps_phi_j * fe_values.JxW(q_point);
                        }

                        if (flag_iter == false) // RHS is populated only by the second iter, the first iter shooting and imbalance generation
                        {
                            const SymmetricTensor<2, dim>  previous_sigma_u =  (g_alpha_gauss + k_res) *  C_Pos * eps_u +  C_Neg * eps_u; 
                            cell_rhs(i) -=previous_sigma_u* eps_phi_i * fe_values.JxW(q_point);
                        }
                    }
                }

                cell->get_dof_indices(local_dof_indices);
                hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                                    cell_rhs,
                                                                    local_dof_indices,
                                                                    system_matrix_elas,
                                                                    system_rhs_elas);
            }
        }
        system_matrix_elas.compress(VectorOperation::add);
        system_rhs_elas.compress(VectorOperation::add);

        std::vector<bool> BC_x = {true, false};
        std::vector<bool> BC_y = {false, true};

        std::map<types::global_dof_index, double> boundary_values;

        Point<dim> point_one, point_two, point_three;
        point_one(0) = 10000.;
        point_one(1) = 0.;
        point_two(0) = 0.;
        point_two(1) = 300.;
        point_three(0) = 1000.;
        point_three(1) = 1000.;

        cell = dof_handler_vec.begin_active(),
        endc = dof_handler_vec.end();

        bool evaluation_point_found = false;
        for (; (cell!=endc) && !evaluation_point_found; ++cell)
        // for (;cell!=endc;++cell)
        {
            if (cell->is_locally_owned())
            {
                for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
                {
                    if (cell->vertex(vertex).distance (point_one) < cell->diameter() * 1e-12)
                    {
                        boundary_values[cell->vertex_dof_index(vertex,0)]=0.;
                        // boundary_values[cell->vertex_dof_index(vertex,1)]=-1.*what_time();
                        evaluation_point_found = true;
                        // break;
                    };
                    if (cell->vertex(vertex).distance (point_two) <	cell->diameter() * 1e-12)
                    {
                        boundary_values[cell->vertex_dof_index(vertex,0)]=0.;
                        evaluation_point_found = true;
                        // break;
                    };
                    if (cell->vertex(vertex).distance (point_three) <	cell->diameter() * 1e-12)
                    {
                	    // boundary_values[cell->vertex_dof_index(vertex,0)]=0;
                        boundary_values[cell->vertex_dof_index(vertex,1)]=0.;
                        evaluation_point_found = true;
                	    // break;
                    };
                };
            }
        }

        VectorTools::interpolate_boundary_values(dof_handler_vec,
                                                 1,
                                                 Functions::ZeroFunction<dim>(dim),
                                                 boundary_values,
                                                 BC_y);

        if (flag_iter == true)
        {
            VectorTools::interpolate_boundary_values(dof_handler_vec,
                                                 3,
                                                 IncrementalBoundaryValues<dim>(real_time, real_time_step),
                                                 boundary_values,
                                                 BC_y);
        }
        else
        {
            VectorTools::interpolate_boundary_values(dof_handler_vec,
                                                 3,
                                                 Functions::ZeroFunction<dim>(dim),
                                                 boundary_values,
                                                 BC_y);
        } 

        PETScWrappers::MPI::Vector tmp(locally_owned_dofs_vec, mpi_communicator);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_elas, tmp, system_rhs_elas, false);
        if (flag_iter == true)
            solution_u = tmp;   
    }

    // Solve elasticity problem
    template <int dim>
    void TopLevel<dim>::solve_linear_problem(const bool flag_elastic_iter)
    {
        PETScWrappers::MPI::Vector distributed_solution_u(locally_owned_dofs_vec, mpi_communicator);
        SolverControl cn;
        PETScWrappers::SparseDirectMUMPS mumps(cn, mpi_communicator);
        // mumps.set_symmetric_mode(true);

        if (flag_elastic_iter == true)
        {
            distributed_solution_u = solution_u;
            mumps.solve(system_matrix_elas, distributed_solution_u, system_rhs_elas);
            solution_u = distributed_solution_u;
            hanging_node_constraints.distribute(solution_u);
        }
        else
        {
            distributed_solution_u = newton_update_u;
            mumps.solve(system_matrix_elas, distributed_solution_u, system_rhs_elas);
            newton_update_u = distributed_solution_u;
            hanging_node_constraints.distribute(newton_update_u);
        }
    }

    // ------------------------------------------------------------------------------------------
    // DAMAGE - Setup / Assemble / Solve
    // ------------------------------------------------------------------------------------------
    // Setup system damage
    template <int dim>
    void TopLevel<dim>::setup_system_alpha()
    {
        dof_handler.distribute_dofs(fe);
        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

        hanging_node_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        hanging_node_constraints.close();
        
        DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs);
        DoFTools::make_sparsity_pattern(dof_handler,
                                        sparsity_pattern,
                                        hanging_node_constraints,
                                        false);
        SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
                                                   locally_owned_dofs,
                                                   mpi_communicator,
                                                   locally_relevant_dofs);

        system_matrix_alpha.reinit(locally_owned_dofs,
                                   locally_owned_dofs,
                                   sparsity_pattern,
                                   mpi_communicator);
        system_rhs_alpha.reinit(locally_owned_dofs, mpi_communicator);

        solution_alpha.reinit(dof_handler.n_dofs());
        solution_alpha_previous_step.reinit(dof_handler.n_dofs());
        newton_update_alpha.reinit(dof_handler.n_dofs());

        alpha_vec.reinit(locally_owned_dofs, mpi_communicator);
    }

    // Assemble system damage
    template <int dim>
    void TopLevel<dim>::assemble_system_alpha(const bool flag_iter_alpha)
    {
        system_rhs_alpha = 0;
        system_matrix_alpha = 0;
        
        FEValues<dim> fe_values_alpha(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                      update_quadrature_points | update_JxW_values);
        
        FEValues<dim> fe_values(fe_vec,
                                quadrature_formula,
                                update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();
        
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     cell_rhs(dofs_per_cell);
        
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        std::vector<double> local_solution_alpha(n_q_points);
        std::vector<double> previous_step_local_solution_alpha(n_q_points);
        std::vector<Tensor<1, dim>> previous_grad_local_solution_alpha (n_q_points);
        std::vector< std::vector<Tensor<1, dim>>> previous_gradient (n_q_points, std::vector<Tensor<1,dim>>(dim));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_vec.begin_active();

        typename DoFHandler<dim>::active_cell_iterator cell_alpha = dof_handler.begin_active(),
					    	                           endc_alpha = dof_handler.end();

        for (unsigned int index = 0; cell_alpha!=endc_alpha; ++cell, ++cell_alpha, ++index)
        {
            if (cell_alpha->is_locally_owned())
            {
                cell_matrix = 0;
                cell_rhs = 0;
            
                fe_values_alpha.reinit(cell_alpha);
                fe_values.reinit (cell);

                fe_values.get_function_gradients(solution_u,previous_gradient);            
            
                fe_values_alpha.get_function_values(solution_alpha,local_solution_alpha);
                fe_values_alpha.get_function_values(solution_alpha_previous_step,previous_step_local_solution_alpha);
                fe_values_alpha.get_function_gradients (solution_alpha,previous_grad_local_solution_alpha);

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                {
                    const SymmetricTensor<2,dim> eps_u = get_strain (previous_gradient[q_point]);
                    double alpha_diff=local_solution_alpha[q_point]-previous_step_local_solution_alpha[q_point];

                    SymmetricTensor<4, dim>  C_Pos, C_Neg;
                    stiffness_matrix(C_Pos, C_Neg, eps_u, stress_strain_tensor, Young, Poisson);                           
              
                    double elastic_source=0.;
                    elastic_source = 0.5* eps_u*C_Pos*eps_u;    

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            cell_matrix(i,j)+= (fe_values_alpha.shape_grad(i,q_point) *
                                                fe_values_alpha.shape_grad(j,q_point) *
                                                            (2.*Gc*ell/c_w)           
                                                +
                                                fe_values_alpha.shape_value(i,q_point)  *
                                                fe_values_alpha.shape_value(j,q_point)  *
                                                (Gc/(ell*c_w)*w_second_alpha()+
                                                g_second_alpha()*elastic_source+
                                                gamma_penalty*heaviside_neg(alpha_diff))     
                                                ) * fe_values_alpha.JxW(q_point);
                        }
                    }
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        cell_rhs(i) -=  (g_prime_alpha(local_solution_alpha[q_point])*elastic_source+
                                            w_prime_alpha(local_solution_alpha[q_point])*Gc/(ell*c_w)+
                                            gamma_penalty*std::min(alpha_diff,0.))  *   
                                        fe_values_alpha.shape_value(i,q_point) *
                                        fe_values_alpha.JxW(q_point);

                        cell_rhs(i) -=   (previous_grad_local_solution_alpha[q_point]*2.*Gc*ell/c_w) *
                                        fe_values_alpha.shape_grad(i,q_point) *
                                        fe_values_alpha.JxW(q_point);
                    }                                                     
                }
                cell_alpha->get_dof_indices(local_dof_indices);
                hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                                    cell_rhs,
                                                                    local_dof_indices,
                                                                    system_matrix_alpha,
                                                                    system_rhs_alpha);
            }
        }
        system_matrix_alpha.compress(VectorOperation::add);
        system_rhs_alpha.compress(VectorOperation::add);
        std::map<types::global_dof_index, double> boundary_values;

        if (flag_iter_alpha == true)
        {
            VectorTools::interpolate_boundary_values(dof_handler,
                                                    1,
                                                    Functions::ZeroFunction<dim>(),
                                                    boundary_values);

            VectorTools::interpolate_boundary_values(dof_handler,
                                                    3,
                                                    Functions::ZeroFunction<dim>(),
                                                    boundary_values);
        }
        else
        {
            VectorTools::interpolate_boundary_values(dof_handler,
                                                    1,
                                                    Functions::ZeroFunction<dim>(),
                                                    boundary_values);

            VectorTools::interpolate_boundary_values(dof_handler,
                                                    3,
                                                    Functions::ZeroFunction<dim>(),
                                                    boundary_values);

        }

        PETScWrappers::MPI::Vector tmp_alpha(locally_owned_dofs, mpi_communicator);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_alpha, tmp_alpha, system_rhs_alpha, false);
        if (flag_iter_alpha == true) solution_alpha = tmp_alpha;
    }

    // Solve system alpha 
    template <int dim>
    unsigned int TopLevel<dim>::solve_alpha(const bool flag_iter_alpha)
    {
        PETScWrappers::MPI::Vector distributed_alpha(locally_owned_dofs, mpi_communicator);
        SolverControl solver_control(dof_handler.n_dofs(), 1e-12* system_rhs_alpha.l2_norm());
        PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
        // PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix_alpha);
        PETScWrappers::PreconditionSOR preconditioner(system_matrix_alpha,1.2);
  
         if (flag_iter_alpha == true)
        {
            distributed_alpha = solution_alpha;
            cg.solve(system_matrix_alpha, distributed_alpha, system_rhs_alpha, preconditioner);
            solution_alpha = distributed_alpha;
            hanging_node_constraints.distribute(solution_alpha);
        }
        else
        {
            distributed_alpha =  newton_update_alpha;
            cg.solve(system_matrix_alpha, distributed_alpha, system_rhs_alpha, preconditioner);
            newton_update_alpha = distributed_alpha;
            hanging_node_constraints.distribute(newton_update_alpha);
        } 
        return solver_control.last_step();
    }

    // ------------------------------------------------------------------------------------------
    // OUTPUT RESULTS
    // ------------------------------------------------------------------------------------------
    // Elasticity
    template <int dim>
    void TopLevel<dim>::output_results_elas(const unsigned int cycle) const
    {
        const Vector<double> localized_solution(solution_u);
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler_vec);

        std::vector<std::string> solution_names(dim, "u");
        std::vector<DataComponentInterpretation::DataComponentInterpretation> 
                    interpretation(dim,
                                   DataComponentInterpretation::component_is_part_of_vector);
   
        data_out.add_data_vector(localized_solution,
                                   solution_names,
                                   DataOut<dim>::type_dof_data,
                                   interpretation);
        data_out.build_patches();

        const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record("./", "solution_u", cycle, mpi_communicator, 4);
        if (this_mpi_process == 0)
        {
            static std::vector<std::pair<double, std::string>> times_and_names;
            times_and_names.push_back(std::pair<double, std::string>(cycle, pvtu_filename));
            std::ofstream pvd_output("./solution_u.pvd");
            DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }
    }
 
    // Damage
    template <int dim>
    void TopLevel<dim>::output_results_alpha(const unsigned int cycle) const
    {
        const Vector<double> localized_solution(solution_alpha);
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);

        std::vector<std::string> solution_names;
        solution_names.push_back ("alpha");
        data_out.add_data_vector(localized_solution,
                                solution_names);
        data_out.build_patches();

        const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record("./", "solution_alpha", cycle, mpi_communicator, 4);
        if (this_mpi_process == 0)
        {
            static std::vector<std::pair<double, std::string>> times_and_names;
            times_and_names.push_back(std::pair<double, std::string>(cycle, pvtu_filename));
            std::ofstream pvd_output("./solution_alpha.pvd");
            DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }
    }

    // Stresses
    template <int dim>
    void TopLevel<dim>::output_stress(const unsigned int cycle)
    {
        Vector<double> s_xx (triangulation.n_active_cells());
        Vector<double> s_xy (triangulation.n_active_cells());
        Vector<double> s_yy (triangulation.n_active_cells());
        Vector<double> source_plot (triangulation.n_active_cells());
        Vector<double> source_plot_neg (triangulation.n_active_cells());
        
        FEValues<dim> fe_values(fe,
                                quadrature_formula,
                                update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);
        
        FEValues<dim> fe_values_vec(fe_vec,
                                    quadrature_formula,
                                    update_values | update_gradients |
                                    update_quadrature_points | update_JxW_values);
        
        const unsigned int n_q_points = quadrature_formula.size();
        
        std::vector<std::vector<Tensor< 1, dim>>> previous_gradient (n_q_points,
                                                                    std::vector<Tensor<1,dim>>(dim));
        std::vector<double> local_solution_alpha(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_vec.begin_active(),
                                                       endc = dof_handler_vec.end();
        typename DoFHandler<dim>::active_cell_iterator cell_alpha = dof_handler.begin_active();

        for (unsigned int index = 0; cell!=endc; ++cell, ++cell_alpha, ++index)
        {
            if (cell->is_locally_owned())
            {
                fe_values_vec.reinit(cell);
                fe_values.reinit (cell_alpha);

                fe_values_vec.get_function_gradients(solution_u, previous_gradient);
                fe_values.get_function_values(solution_alpha, local_solution_alpha);

                double c_s_xx = 0;
                double c_s_xy = 0;
                double c_s_yy = 0;
                double el_cell = 0;
                double el_cell_neg = 0;

                for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
                {
                    const double g_alpha_gauss = g_alpha(local_solution_alpha[q_point]);
                    
                    const SymmetricTensor<2,dim> eps_u = get_strain (previous_gradient[q_point]);
                    
                    SymmetricTensor<4, dim>  C_pos, C_neg;
                    stiffness_matrix(C_pos, C_neg, eps_u, stress_strain_tensor, Young, Poisson);

                    el_cell += 0.5 * eps_u * C_pos * eps_u + 1.e-16;
                    el_cell_neg += 0.5 * (eps_u) * C_neg * (eps_u) + 1.e-16;
                    const SymmetricTensor<2, dim> stresses = (g_alpha_gauss + k_res) * C_pos * eps_u + C_neg * eps_u;

                    c_s_xx += stresses[0][0];
                    c_s_xy += stresses[0][1];
                    c_s_yy += stresses[1][1];
                }

                s_xx(index) = (c_s_xx / quadrature_formula.size());
                s_xy(index) = (c_s_xy / quadrature_formula.size());
                s_yy(index) = (c_s_yy / quadrature_formula.size());
                source_plot(index) =(el_cell / quadrature_formula.size());
                source_plot_neg(index) =(el_cell_neg / quadrature_formula.size());
            }
        }

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        switch (dim)
        {
        case 1:
            data_out.add_data_vector(s_xx, "s_xx");
            break;
        case 2:
            data_out.add_data_vector(s_xx, "s_xx");
            data_out.add_data_vector(s_xy, "s_xy");
            data_out.add_data_vector(s_yy, "s_yy");
            data_out.add_data_vector(source_plot, "source");
            data_out.add_data_vector(source_plot_neg, "source_neg");
            break; 
        default:
            Assert (false, ExcInternalError());
        }
        data_out.build_patches();

        const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record("./", "stresses", cycle, mpi_communicator, 4);
        if (this_mpi_process == 0)
        {
            static std::vector<std::pair<double, std::string>> times_and_names;
            times_and_names.push_back(std::pair<double, std::string>(cycle, pvtu_filename));
            std::ofstream pvd_output("./stresses.pvd");
            DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }
    }

    // ------------------------------------------------------------------------------------------
    // RUN
    // ------------------------------------------------------------------------------------------
    template <int dim>
    void TopLevel<dim>::run()
    {
        set_data ();
        do_initial_timestep();

        std::ofstream output_text("output.txt", std::ios::out);
        output_text << " il dato riportato e' il seguente:  bulk - surface - total energy" << std::endl;

        while (real_time < real_time_final)
        {
            do_timestep();

            // Output
            output_results_elas(real_timestep_no);
            output_results_alpha(real_timestep_no);            
            
            double bulk_energy=0.;
            double surface_energy=0.;      
            energy (bulk_energy, surface_energy);
            
            Tensor<1,dim> reaction_stress_id_3;
            reaction(reaction_stress_id_3, 3);   
            
            if (this_mpi_process == 0)
            {
            output_text << real_timestep_no << "  " << bulk_energy  << "  " << surface_energy << "  " << bulk_energy+ surface_energy 
                                              << "  " << reaction_stress_id_3[0] << "  " << reaction_stress_id_3[1] << std::endl;
            }
        }
    }

    // -------------------------------------------------------------------------------------------
    // Timestep solve - Elasticity / Damage
    // -------------------------------------------------------------------------------------------
    template <int dim>
    void TopLevel<dim>::do_initial_timestep()
    {
        real_time += real_time_step;
        ++real_timestep_no;
        pcout << "Timestep " << real_timestep_no << " at time " << real_time << std::endl;

        create_coarse_grid();

        pcout << "    Number of active cells:       "
              << triangulation.n_active_cells() << " (by partition:";
              for (unsigned int p = 0; p < n_mpi_processes; ++p)
                    pcout << (p == 0 ? ' ' : '+') << (GridTools::count_cells_with_subdomain_association(triangulation, p));
        pcout << ")" << std::endl;

        setup_system_elas();
        setup_system_alpha();
        set_initial_values_alpha(alpha0);

        pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs() << " (by partition:";
        for (unsigned int p = 0; p < n_mpi_processes; ++p)
            pcout << (p == 0 ? ' ' : '+') << (DoFTools::count_dofs_with_subdomain_association(dof_handler, p));
        pcout << ")" << std::endl;

        solve_timestep();
        
        real_time -= real_time_step;
        --real_timestep_no;
    }

    // Timestep function
    template <int dim>
    void TopLevel<dim>::do_timestep()
    {
        real_time += real_time_step;
        ++real_timestep_no;

        pcout << "Timestep " << real_timestep_no << " at time " << real_time << std::endl;

        solve_timestep();
    }

    template <int dim>
    void TopLevel<dim>::solve_timestep()
    {
        // Control variables declaration
        double error_alternate = 1.0;
        double  error_elastic;
        double  error_alpha; 
        const double error_toll =1.e-4;
        const double error_toll_elastic=1.e-6;
        const double error_toll_alpha=1.e-6;        
        unsigned int iter_counter_am = 0;
        unsigned int iter_elastic;
        unsigned int iter_alpha; 
        const unsigned int max_iteration = 2000;
        const unsigned int max_iteration_elastic = 500;
        const unsigned int max_iteration_alpha = 500;    
        bool solve_step;
        bool solve_step_alpha;

                while (error_alternate > error_toll && iter_counter_am <max_iteration)
        {
            
            pcout << "      Solving Elastic problem..." << std::endl;
            
            iter_elastic = 0;            
            solve_step = true;

            assemble_system_elas(solve_step);
            solve_linear_problem(solve_step);

            do
            {
                ++iter_elastic;
                solve_step = false;

                assemble_system_elas(solve_step);
                solve_linear_problem(solve_step);
                
                solution_u.add(1.0,newton_update_u);
                error_elastic = newton_update_u.l2_norm();

            }while (error_elastic > error_toll_elastic && iter_elastic <max_iteration_elastic);

            pcout << "          Iterations: " << iter_elastic << std::endl;
            pcout << "          --- Error_on_Newton_update_u: "   << error_elastic << std::endl; 
            pcout << "          --- rhs_error_elastic: "   <<  system_rhs_elas.l2_norm()  <<  std::endl; 
            
            iter_alpha = 0;  
            solve_step_alpha = true;
            Vector<double> temp_alpha;
            temp_alpha = solution_alpha;
            pcout << "      Solving Damage problem..." << std::endl;  

            assemble_system_alpha(solve_step_alpha);
            solve_alpha(solve_step_alpha);

            do
            {
                ++iter_alpha;
                solve_step_alpha = false;

                assemble_system_alpha(solve_step_alpha);
                solve_alpha(solve_step_alpha);
                
                solution_alpha.add(1.0,newton_update_alpha);
                error_alpha = newton_update_alpha.l2_norm();

                alpha_vec = solution_alpha;
                
                pcout << "          Total Damage Iteration: " << iter_alpha << std::endl;
                pcout << "          --- Error_on_Newton_update_alpha: " << error_alpha  << std::endl;
                pcout << "          --- rhs_error_alpha: "   <<  system_rhs_alpha.l2_norm()  <<  std::endl; 

              //  pcout << "          --- Min_alpha: "   <<  alpha_vec.min()  <<  std::endl;
              //  pcout << "          --- Max_alpha: "   <<  alpha_vec.max()  <<  std::endl;
                

            } while (error_alpha > error_toll_alpha && iter_alpha < max_iteration_alpha);

            temp_alpha -=solution_alpha;
            error_alternate =  temp_alpha.linfty_norm();
            iter_counter_am++;
            
            pcout << " Number of iteration: " << iter_counter_am << std::endl;
            pcout << " Error_on_alpha:  " << error_alternate << std::endl;
            pcout << " Alpha_max:  " << solution_alpha.linfty_norm() << std::endl;
        }                        
        check_alpha(solution_alpha);
        solution_alpha_previous_step = solution_alpha;
    }

    // ------------------------------------------------------------------------------------------
    // REACTION - ENERGY EVALUATION
    // ------------------------------------------------------------------------------------------
    // Reaction
    template <int dim>
    void TopLevel<dim>::reaction (Tensor<1,dim> &reaction_stress, const types::boundary_id &boundary_id)
    {
        QGauss<dim-1> face_quadrature_formula(5);

        FEFaceValues<dim> fe_face_values (fe_vec, 
                                          face_quadrature_formula,
                                          UpdateFlags(update_values             | 				
                                                      update_gradients          | 
                                                      update_quadrature_points  |
                                                      update_normal_vectors     |
                                                      update_JxW_values));

        FEFaceValues<dim> fe_face_values_scalar (fe, 
                                                face_quadrature_formula,
                                                UpdateFlags(update_values            | 
                                                            update_quadrature_points |
                                                            update_normal_vectors    |
                                                            update_JxW_values));

        const unsigned int n_face_q_points = face_quadrature_formula.size();

        std::vector<double> local_solution_alpha(n_face_q_points);
        std::vector< std::vector< Tensor< 1, dim > > > previous_gradient (n_face_q_points,std::vector<Tensor<1,dim> >(dim));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_vec.begin_active(),
						                               endc = dof_handler_vec.end();

        typename DoFHandler<dim>::active_cell_iterator cell_alpha = dof_handler.begin_active();
  
        Tensor<1,dim> reaction_mpi;

        for (; cell!=endc; ++cell, ++cell_alpha)
            if (cell->is_locally_owned())
            {    
                for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    	        if (cell->face(face)->boundary_id() == boundary_id)
	            {
                    fe_face_values.reinit (cell, face);
                    fe_face_values_scalar.reinit (cell_alpha, face);        
                    fe_face_values_scalar.get_function_values(solution_alpha,local_solution_alpha);     
                    fe_face_values.get_function_gradients(solution_u,previous_gradient);                        
       
	                for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                    {
                        const SymmetricTensor<2,dim> eps_u = get_strain (previous_gradient[q_point]);
                        SymmetricTensor<4, dim>  C_Pos, C_Neg;        
                        const Tensor<1,dim> &N = fe_face_values.normal_vector(q_point);
                        const double JxW_f = fe_face_values.JxW(q_point);                        
                        stiffness_matrix(C_Pos, C_Neg, eps_u, stress_strain_tensor, Young, Poisson);
                        
                        const SymmetricTensor<2,dim> sigma = (g_alpha(local_solution_alpha[q_point])+k_res)*C_Pos*eps_u + C_Neg*eps_u;    
                        reaction_mpi +=sigma * N * JxW_f;;
                    }
	            }
            }
        reaction_stress = Utilities::MPI::sum(reaction_mpi,mpi_communicator);    
    }

    template <int dim>
    void TopLevel<dim>::energy (double &bulk_energy, double &surface_energy)
    {
        FEValues<dim> fe_values(fe_vec,
                                quadrature_formula,
                                update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

        FEValues<dim> fe_values_alpha(fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                      update_quadrature_points | update_JxW_values);
    
        const unsigned int n_q_points = quadrature_formula.size();
       
        std::vector<std::vector<Tensor< 1, dim >>> previous_gradient(n_q_points,std::vector<Tensor<1,dim> >(dim));
        std::vector<double> local_solution_alpha(n_q_points);
        std::vector<Tensor< 1, dim > > local_solution_grad_alpha (n_q_points);

        double el_en, ph_el;          

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_vec.begin_active(),
						                               endc = dof_handler_vec.end();

        typename DoFHandler<dim>::active_cell_iterator cell_alpha = dof_handler.begin_active();
        
        for (; cell!=endc; ++cell, ++cell_alpha)
            if (cell->is_locally_owned())
            {
                fe_values.reinit(cell);
                fe_values_alpha.reinit (cell_alpha);        
                fe_values_alpha.get_function_values(solution_alpha,local_solution_alpha);
                fe_values_alpha.get_function_gradients (solution_alpha,local_solution_grad_alpha);          
                fe_values.get_function_gradients(solution_u,previous_gradient);                         
                
                for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                {
                    SymmetricTensor<4, dim>  C_Pos, C_Neg;            
                    const SymmetricTensor<2,dim> eps_u = get_strain (previous_gradient[q_point]);                                
                    double grad_alpha_square=local_solution_grad_alpha[q_point]*local_solution_grad_alpha[q_point];
                    stiffness_matrix(C_Pos, C_Neg, eps_u, stress_strain_tensor, Young, Poisson);
                   
                    double elastic_energy_density =  0.5*(g_alpha(local_solution_alpha[q_point])+k_res)*eps_u*C_Pos*eps_u + 0.5*eps_u*C_Neg*eps_u; 
                    el_en +=elastic_energy_density*fe_values.JxW(q_point); 
                    ph_el += (Gc/c_w)*(ell*grad_alpha_square+w_alpha(local_solution_alpha[q_point])/ell)*fe_values.JxW(q_point);                        

                }              
            }
        bulk_energy = Utilities::MPI::sum(el_en,mpi_communicator);
        surface_energy = Utilities::MPI::sum(ph_el,mpi_communicator);    
    }
} // namespace phasefield

#endif 
