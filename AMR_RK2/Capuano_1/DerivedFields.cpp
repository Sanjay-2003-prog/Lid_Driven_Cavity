#include <NavierStokesSolver.H>

#include <AMReX_MultiFabUtil.H>

void NavierStokesSolver::computeCellCenteredVel()
{
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for (MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &vel = U[lev][mfi].array();
            const auto &u = xvel[lev][mfi].const_array();
            const auto &v = yvel[lev][mfi].const_array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                vel(i, j, k, 0) = 0.5 * (u(i,j, k) + u(i+1, j, k));
                vel(i, j, k, 1) = 0.5 * (v(i,j, k) + v(i, j + 1, k)); });
        }
    }
}

void NavierStokesSolver::computeDivergence()
{
    Real divU_max = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const Real *dxinv = geom[lev].InvCellSize();
        for (MFIter mfi(divU[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &u = xvel[lev][mfi].const_array();
            const auto &v = yvel[lev][mfi].const_array();
            const auto &rhs = divU[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        { rhs(i, j, k) = dxinv[0] * (u(i + 1, j, k) - u(i, j, k)) + dxinv[1] * (v(i, j + 1, k) - v(i, j, k)); });
        }
        if (lev < finest_level)
        {
            int coarse_value = 1;
            int fine_value = 0;
            iMultiFab mask = makeFineMask(divU[lev], grids[lev + 1], ref_ratio[lev],
                                          coarse_value, fine_value);
            Real lev_max_norm = divU[lev].norminf(mask);

            divU_max = std::max(divU_max, lev_max_norm);
        }
        else
        {
            Real lev_max_norm = divU[lev].norminf();
            divU_max = std::max(divU_max, lev_max_norm);
        }
    }
    ParallelDescriptor::ReduceRealMax(divU_max);

    PrintToFile("log") << "div max norm : " << divU_max << "\n";
}

void NavierStokesSolver::computeMagnitude()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const BoxArray &ba = U[lev].boxArray();
        const DistributionMapping &dm = U[lev].DistributionMap();

        if (U_magnitude.size() <= lev)
            U_magnitude.resize(lev + 1);
        U_magnitude[lev].define(ba, dm, 1, 0);

        for (MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &u = U[lev][mfi].const_array();
            const auto &mag = U_magnitude[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                Real sum = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    sum += u(i,j,k,d) * u(i,j,k,d);
                }
                mag(i,j,k) = std::sqrt(sum); });
        }
    }
}

void NavierStokesSolver::computeVorticity(int lev)
{
    // for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real *dxinv = geom[lev].InvCellSize();

        for (MFIter mfi(vorticity[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &u = xvel[lev][mfi].const_array();
            const auto &v = yvel[lev][mfi].const_array();
            const auto &omega = vorticity[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                Real vw = 0.25 * (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k));
                Real ve = 0.25 * (v(i,j,k) + v(i+1,j,k) + v(i,j+1,k) + v(i+1,j+1,k));
                Real un = 0.25 * (u(i,j,k) + u(i+1,j,k) + u(i,j+1,k) + u(i+1,j+1,k));
                Real us = 0.25 * (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k));

                omega(i,j,k) = dxinv[0] * (ve - vw) + dxinv[1] * (un - us); });
        }

        // vorticity[lev].FillBoundary(geom[lev].periodicity());
    }
}