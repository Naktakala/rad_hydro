#include "radtranDO.h"

#include "ChiMath/Quadratures/product_quadrature.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"

void chi_radtran::RadTranDO::Initialize()
{
  //=================================== Process options
  options.num_angles       = basic_options("num_angles").IntegerValue();
  options.scattering_order = basic_options("scattering_order").IntegerValue();
  options.num_groups       = basic_options("num_groups").IntegerValue();

  ChiInvalidArgument(options.num_angles%2,
    "Angular-quadrature number of angles must be divisible by 2")

  //=================================== Create unknown manager
  num_moments = options.scattering_order+1;
  for (size_t m=0; m<num_moments; ++m)
    uk_man_phi.AddUnknown(chi_math::UnknownType::VECTOR_N, options.num_groups);

  //=================================== Create spatial discretization
  grid = chi_mesh::GetCurrentHandler().GetGrid();
  sdm = SDM_PWLD::New(grid);

  num_local_phiDOFs = sdm->GetNumLocalDOFs(uk_man_phi);

  //=================================== Initialize discrete vectors
  phi_old.assign(num_local_phiDOFs, 0.0);

  //=================================== Create angular quadrature
  auto product_quad = std::make_shared<chi_math::ProductQuadrature>();
  product_quad->InitializeWithGL(static_cast<int>(options.num_angles/2), false);
  angular_quadrature = product_quad;

  //=================================== Sort angles into angle-sets
  if (grid->Attributes() | chi_mesh::DIMENSION_1)
  {
    std::vector<AngularQuadraturePointInfo> top_hemi;
    std::vector<AngularQuadraturePointInfo> bot_hemi;

    const size_t num_angles = angular_quadrature->abscissae.size();
    for (size_t n=0; n<num_angles; ++n)
    {
      const auto& qp = angular_quadrature->abscissae[n];
      const auto& omega = angular_quadrature->omegas[n];
      const auto& weight = angular_quadrature->weights[n];

      if (omega.z > 0.0) top_hemi.emplace_back(qp, omega, weight, n);
      else               bot_hemi.emplace_back(qp, omega, weight, n);
    }

    if (not top_hemi.empty())
    {
      auto top_hemi_spds =
        std::make_shared<SweepStructure>(*grid, top_hemi.front().omega);
      angle_sets.emplace_back(top_hemi, top_hemi_spds);
    }
    if (not bot_hemi.empty())
    {
      auto bot_hemi_spds =
        std::make_shared<SweepStructure>(*grid, bot_hemi.front().omega);
      angle_sets.emplace_back(bot_hemi, bot_hemi_spds);
    }

  }
  else
    ChiLogicalError(true, "Only 1D meshes are currently supported.")

  //======================================== Sweep scheduler
  sweep_scheduler = std::make_unique<SweepSchedulerFIFO>(angle_sets);

  sweep_scheduler->Initialize();
}