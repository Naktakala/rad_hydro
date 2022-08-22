#include "rhsolverC.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

//###################################################################
/**Generic update of density and momentum from lagged
 * radiation momentum deposition.*/
void chi_radhydro::SolverC_SNCN_MFEM::
  DensityMomentumUpdateWithRadF0(
    SimRefs&                              sim_refs,
    const std::vector<double>&            kappa_t,
    double                                tau,
    const std::vector<UVector>&           U_old,
    const std::vector<UVector>&           U_int,
    const std::vector<Vec3>&              rad_F0_old,
    std::vector<UVector>&                 U_new
  )
{
  U_new = U_int;

  const auto& Cv          = sim_refs.Cv;
  const auto& bc_settings = sim_refs.bc_settings;
  const auto& grid        = sim_refs.grid;
  auto& fv                = sim_refs.fv;
  auto& pwld              = *m_pwld;

  const auto& kappa_s_function = sim_refs.kappa_s_function;
  const auto& kappa_a_function = sim_refs.kappa_a_function;

  for (const auto& cell : sim_refs.grid.local_cells)
  {
    const uint64_t c        = cell.local_id;
    const auto&    fem_data = pwld.GetUnitIntegrals(cell);
    const size_t   Nn       = fem_data.NumNodes();

    const UVector& U_c_old     = U_old[c];
    const double   rho_c_old   = U_old[c][RHO];

    const double   sigma_t_c = rho_c_old * kappa_t[c];

    Vec3 rad_F0_c_old;
    for (size_t i=0; i<Nn; ++i)
    {
      const auto imap = pwld.MapDOFLocal(cell, i);
      rad_F0_c_old += rad_F0_old[imap];
    }
    rad_F0_c_old /= static_cast<double>(Nn);

    UVector U_c_new_01 = U_new[c];

    if (std::fabs(kappa_t[c]) < 1.0e-10) continue;

    const double sigma_t_div_c = sigma_t_c/speed_of_light_cmpsh;

    U_c_new_01(1) += (1 / tau) * sigma_t_div_c * rad_F0_c_old.x;
    U_c_new_01(2) += (1 / tau) * sigma_t_div_c * rad_F0_c_old.y;
    U_c_new_01(3) += (1 / tau) * sigma_t_div_c * rad_F0_c_old.z;

    U_new[c] = U_c_new_01;
  }//for cell
}