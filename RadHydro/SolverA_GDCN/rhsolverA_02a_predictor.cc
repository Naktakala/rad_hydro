#include "rhsolverA.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/chi_math_banded_solvers.h"

#include "ChiTimer/chi_timer.h"

void chi_radhydro::SolverA_GDCN::
Predictor(const std::map<uint64_t, BCSetting>& bc_setttings,
          const std::vector<double>&      kappa_a_n,
          const std::vector<double>&      kappa_t_n,
          double dt,
          const std::vector<UVector>&     U_n,
          const std::vector<GradUTensor>& grad_U_n,
          std::vector<UVector>&           U_nph,

          const std::vector<double>&            rad_E_n,
          const std::vector<chi_mesh::Vector3>& grad_rad_E_n,
          std::vector<double>&                  rad_E_nph)
{
  std::vector<UVector> U_n_star = U_n;
  std::vector<double>  rad_E_n_star = rad_E_n;

  const size_t num_local_nodes = rad_E_n.size();

  //####################################################### PREDICTOR
  const double tau = 1/(0.5*dt);

  //=================================== Advection of U and rad_E
  MHM_HydroRadEPredictor(
    *grid, *fv, m_gamma,     //Stuff
    tau,                     //tau input
    U_n, grad_U_n,           //Hydro inputs
    rad_E_n, grad_rad_E_n,   //RadE inputs
    U_n_star, rad_E_n_star   //Outputs
    );

  //=================================== Update density and momentum to nph
  DensityMomentumUpdateWithRadMom(
    *grid, *fv, bc_setttings, kappa_t_n,    //Stuff
    tau,                                    //tau input
    U_n, U_n_star,                          //Hydro inputs
    rad_E_n,                                //RadE input
    m_kappa_s_function, m_kappa_a_function, //Opacity functions
    U_nph                                   //Output
  );

  //=================================== Internal energy and radiation energy
  {
    std::vector<double> k5,k6;
    MatDbl A(num_local_nodes, VecDbl(num_local_nodes,0.0));
    VecDbl b(num_local_nodes, 0.0);

    AssembleGeneralEnergySystem(
      *grid, fv, bc_setttings,
      kappa_a_n  , kappa_t_n,
      kappa_a_n  , kappa_t_n, m_Cv,                   //Stuff
      tau, /*theta1=*/1.0, /*theta2=*/0.0,          //tau, theta input
      U_n, U_n, U_n_star, U_nph, grad_U_n,          //Hydro inputs
      rad_E_n, rad_E_n, rad_E_n_star, grad_rad_E_n, //RadE inputs
      k5, k6, A, b);                                //Outputs

    rad_E_nph = chi_math::TDMA(A,b);

    //============================ Back sub for internal energy
    for (const auto& cell_c : grid->local_cells)
    {
      const uint64_t c     = cell_c.local_id;
      double rho_c_nph     = U_nph[c][RHO];
      Vec3   u_c_nph       = VelocityFromCellU(U_nph[c]);
      double u_abs_sqr_nph = u_c_nph.NormSquare();
      double e_c_nph       = k5[c]*rad_E_nph[c] + k6[c];

      double& E_c_nph = U_nph[c](MAT_E);

      E_c_nph = rho_c_nph*(0.5 * u_abs_sqr_nph + e_c_nph);
    }//for cell c

//    //================================= Compute energy conservation
//    auto system_init_energy = ComputeSysEnergyChange(dt,U_n,rad_E_n,
//                                                    bc_settings);
//    auto new_system_energy = ComputeSysEnergyChange(dt,U_nph,rad_E_nph,
//                                                    bc_settings);
//    auto dEmat = new_system_energy.Emat - system_init_energy.Emat;
//    auto dErad = new_system_energy.Erad - system_init_energy.Erad;
//    chi::log.Log()
//      << "me dEmat: " <<std::setw(10)<<std::scientific<<std::setprecision(8)<< dEmat << "\n"
//      << "me dErad: " <<std::setw(10)<<std::scientific<<std::setprecision(8)<< dErad << "\n";

  }//scope internal e and rad_E
}

