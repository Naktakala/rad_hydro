#include "rhsolverA.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/chi_math_banded_solvers.h"

void chi_radhydro::SolverA_GDCN::
  Corrector(const std::map<uint64_t, BCSetting>& bc_setttings,
            const std::vector<double>&      kappa_a_n,
            const std::vector<double>&      kappa_t_n,
            const std::vector<double>&      kappa_a_nph,
            const std::vector<double>&      kappa_t_nph,
            double dt,
            const std::vector<UVector> &U_n,
            const std::vector<UVector> &U_nph,
            const std::vector<GradUTensor> &grad_U_nph,
            std::vector<UVector> &U_np1,

            const std::vector<double> &rad_E_n,
            const std::vector<double> &rad_E_nph,
            const std::vector<chi_mesh::Vector3> &grad_rad_E_nph,
            std::vector<double> &rad_E_np1)
{
  chi::log.Log0Verbose1() << "Executing corrector";

  std::vector<UVector> U_nph_star = U_n;
  std::vector<double>  rad_E_nph_star = rad_E_n;

  auto& Cv = scalar_fields.at("Cv");
  auto& gamma = scalar_fields.at("gamma");

  const size_t num_local_nodes = rad_E_n.size();

  //####################################################### PREDICTOR
  const double   tau     = 1/dt;

  //=================================== Advection of U and rad_E
  MHM_HydroRadECorrector(
    *grid, *fv, bc_setttings, gamma,    //Stuff
    tau,                                //tau input
    U_n, U_nph, grad_U_nph,             //Hydro inputs
    rad_E_n, rad_E_nph, grad_rad_E_nph, //RadE inputs
    U_nph_star, rad_E_nph_star          //Outputs
  );

  //=================================== Update density and momentum to np1
  DensityMomentumUpdateWithRadMom(
    *grid, *fv, bc_setttings, kappa_t_nph,  //Stuff
    tau,                                    //tau input
    U_nph, U_nph_star,                      //Hydro inputs
    rad_E_nph,                              //RadE input
    U_np1                                   //Output
  );

  //=================================== Internal energy and radiation energy
  {
    std::vector<double> k5,k6;
    MatDbl A(num_local_nodes, VecDbl(num_local_nodes,0.0));
    VecDbl b(num_local_nodes, 0.0);

    AssembleGeneralEnergySystem(
      *grid, fv, bc_setttings,
      kappa_a_n  , kappa_t_n,
      kappa_a_nph, kappa_t_nph, Cv,                       //Stuff
      tau, /*theta1=*/0.5, /*Crank-Nicolson*/             //tau, theta input
      U_n, U_nph, U_nph_star, U_np1, grad_U_nph,          //Hydro inputs
      rad_E_n, rad_E_nph, rad_E_nph_star, grad_rad_E_nph, //RadE inputs
      k5, k6, A, b);                                      //Outputs

    rad_E_np1 = chi_math::TDMA(A,b);

    for (const auto& cell_c : grid->local_cells)
    {
      const uint64_t c = cell_c.local_id;
      double rho_c_np1 = U_np1[c][0];
      Vec3   u_c_np1   = chi_radhydro::VelocityFromCellU(U_np1[c]);
      double u_abs_sqr_np1 = u_c_np1.NormSquare();
      double e_c_np1 = k5[c]*rad_E_np1[c] + k6[c];

      double& E_c_np1 = U_np1[c](4);

      E_c_np1 = rho_c_np1*(0.5 * u_abs_sqr_np1 + e_c_np1);
    }//for cell c
  }//scope internal e and rad_E

  chi::log.Log0Verbose1() << "Done executing corrector";
}