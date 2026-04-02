#include "NavierStokesSolver.H"
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void NavierStokesSolver::predictVelocity()
{
    BL_PROFILE("Predict Velocity");
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
        // amrex::Print() << "DEBUG: predictVelocity starting for level " << lev 
        //               << ", time = " << curr_time << ", dt = " << dt 
        //               << ", xvel_max = " << xvel[lev].max(0) 
        //               << ", yvel_max = " << yvel[lev].max(0) 
        //               << ", xvel_min = " << xvel[lev].min(0) 
        //               << ", yvel_min = " << yvel[lev].min(0) << std::endl;

        // Check inputs before computation
        // AMREX_ALWAYS_ASSERT(xvel_old[lev].ok());
        // AMREX_ALWAYS_ASSERT(yvel_old[lev].ok());
        
        // AMREX_ALWAYS_ASSERT(!xvel_old[lev].contains_nan());
        // AMREX_ALWAYS_ASSERT(!yvel_old[lev].contains_nan());

        
        // if (xvel_old[lev].contains_nan()) {
        //     amrex::Print() << "ERROR: NaN in xvel_old at level " << lev << " before prediction!" << std::endl;
        // }
        // if (yvel_old[lev].contains_nan()) {
        //     amrex::Print() << "ERROR: NaN in yvel_old at level " << lev << " before prediction!" << std::endl;
        // }
        // xvel_old[lev].FillBoundary(geom[lev].periodicity());
        // yvel_old[lev].FillBoundary(geom[lev].periodicity());

        const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
        const amrex::GpuArray<Real, AMREX_SPACEDIM> plo = geom[lev].ProbLoArray();
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dxinv = geom[lev].InvCellSizeArray();

        Real dx2inv = 1.0 / (dx[0] * dx[0]);
        Real dy2inv = 1.0 / (dx[1] * dx[1]);
        
        // Capture scalar values as plain Real (POD types are safe)
        Real Mu_val = Mu;
        Real dt_val = dt;

        // xvel update
        for (MFIter mfi(xvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box &bx = mfi.tilebox(IntVect::TheDimensionVector(0)); // nodal in x
            // Get device arrays - these are safe
            const auto &u = xvel[lev][mfi].array();
            const auto &uo = xvel_old[lev][mfi].array();
            const auto &vo = yvel_old[lev][mfi].array();

            // Capture everything by value
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Real ue = 0.5 * (uo(i + 1, j, k) + uo(i, j, k));
                Real uw = 0.5 * (uo(i - 1, j, k) + uo(i, j, k));
                Real un = 0.5 * (uo(i, j + 1, k) + uo(i, j, k));
                Real us = 0.5 * (uo(i, j - 1, k) + uo(i, j, k));
                Real vs = 0.5 * (vo(i - 1, j, k) + vo(i, j, k));
                Real vn = 0.5 * (vo(i - 1, j + 1, k) + vo(i, j + 1, k));

                Real diff = dx2inv * (uo(i + 1, j, k) - 2.0 * uo(i, j, k) + uo(i - 1, j, k)) 
                          + dy2inv * (uo(i, j + 1, k) - 2.0 * uo(i, j, k) + uo(i, j - 1, k));
                Real adv = dxinv[0] * (ue * ue - uw * uw) + dxinv[1] * (vn * un - vs * us);
                u(i, j, k) = uo(i, j, k) + dt_val * (Mu_val * diff - adv);
            });
        }

        //FillPhysBCs_xnd(xvel[lev], dx, plo, geom[lev].Domain(), u_bc, u_bc_func);
      
      
       // AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());

        // yvel update
        for (MFIter mfi(yvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box &bx = mfi.tilebox(IntVect::TheDimensionVector(1)); // nodal in y
            // Get device arrays - these are safe
            const auto &v = yvel[lev][mfi].array();
            const auto &uo = xvel_old[lev][mfi].array();
            const auto &vo = yvel_old[lev][mfi].array();
            
            // Capture everything by value
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Real ue = 0.5 * (uo(i + 1, j, k) + uo(i+1, j-1, k));
                Real uw = 0.5 * (uo(i, j, k) + uo(i, j-1, k));
                Real vn = 0.5 * (vo(i, j, k) + vo(i, j+1, k));
                Real vs = 0.5 * (vo(i, j - 1, k) + vo(i, j, k));
                Real ve = 0.5 * (vo(i + 1, j, k) + vo(i, j, k));
                Real vw = 0.5 * (vo(i - 1, j , k) + vo(i, j, k));

                Real diff = dx2inv * (vo(i + 1, j, k) - 2.0 * vo(i, j, k) + vo(i - 1, j, k)) 
                          + dy2inv * (vo(i, j + 1, k) - 2.0 * vo(i, j, k) + vo(i, j - 1, k));
                Real adv = dxinv[0] * (ue * ve - uw * vw) + dxinv[1] * (vn * vn - vs * vs);
                v(i, j, k) = vo(i, j, k) + dt_val * (Mu_val * diff - adv);
            });
        }

        //FillPhysBCs_ynd(yvel[lev], dx, plo, geom[lev].Domain(), v_bc, v_bc_func);

        // // Validate yvel after computation
        // if (yvel[lev].contains_nan()) {
        //     amrex::Print() << "ERROR: NaN in yvel after computation at level " << lev << std::endl;
        //     amrex::Print() << "DEBUG: Mu = " << Mu << ", dt = " << dt
        //                   << ", dx = " << dx[0] << "," << dx[1] << std::endl;
            
        //     // Detailed debugging for yvel
        //     for (MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi) {
        //         const auto& fab = yvel[lev][mfi];
        //         if (fab.contains_nan()) {
        //             amrex::Print() << "  NaN in yvel fab " << mfi.index() 
        //                           << ", box = " << mfi.validbox() << std::endl;
                    
        //             // Check specific values in the problematic fab
        //             const auto& arr = fab.array();
        //             const Box& box = mfi.validbox();
        //             IntVect lo = box.smallEnd();
        //             IntVect hi = box.bigEnd();
                    
        //             // Print a few values to see what's happening
        //             for (int j = lo[1]; j <= std::min(lo[1] + 5, hi[1]); ++j) {
        //                 for (int i = lo[0]; i <= std::min(lo[0] + 5, hi[0]); ++i) {
        //                     amrex::Print() << "  yvel(" << i << "," << j << ",0) = " << arr(i,j,0) << std::endl;
        //                 }
        //             }
        //             break;
        //         }
        //     }
        // }

        // AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
        // Allocate flux register if needed
        //if (lev < max_level) {
          //  if (flux_register.size() <= lev) {
            //    flux_register.resize(max_level + 1);
            //}

            //if (!flux_register[lev]) {
              //  flux_register[lev] = std::make_unique<amrex::FluxRegister>(
                //    grids[lev], dmap[lev], IntVect(ref_ratio[lev], ref_ratio[lev]), lev, xvel[lev].nComp());
            //}
        //}

        // Flux register updates
        //if (lev < finest_level && flux_register[lev]) {
          //  flux_register[lev]->FineAdd(xvel[lev], 0, 0, 0, xvel[lev].nComp(), 1.0);  // dir = 0
          //  flux_register[lev]->FineAdd(yvel[lev], 1, 0, 0, yvel[lev].nComp(), 1.0);  // dir = 1
        //}

        //if (lev > 0 && flux_register[lev - 1]) {
          //  flux_register[lev - 1]->CrseAdd(xvel[lev], 0, 0, 0, xvel[lev].nComp(), 1.0, Geom(lev - 1));  // dir = 0
          //  flux_register[lev - 1]->CrseAdd(yvel[lev], 1, 0, 0, yvel[lev].nComp(), 1.0, Geom(lev - 1));  // dir = 1
        //}
    }
}