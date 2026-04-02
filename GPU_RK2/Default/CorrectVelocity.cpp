#include <NavierStokesSolver.H>
#include <FillBC.H>
// #include <AMReX_FluxRegister.H>

void NavierStokesSolver::correctVelocity()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        BL_PROFILE("Correct5 Velocity");
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dxinv = geom[lev].InvCellSizeArray();
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
        const amrex::GpuArray<Real, AMREX_SPACEDIM> plo = geom[lev].ProbLoArray();

        // xvel
        for (MFIter mfi(xvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &u = xvel[lev][mfi].array();
            const auto &p = Pressure[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        { u(i, j, k) -= dxinv[0] * (p(i, j, k) - p(i - 1, j, k)); });
        }

        // yvel
        for (MFIter mfi(yvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
            const auto &v = yvel[lev][mfi].array();
            const auto &p = Pressure[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        { v(i, j, k) -= dxinv[1] * (p(i, j, k) - p(i, j - 1, k)); });
        }

        // // Apply physical boundary conditions
        // FillPhysBCs_xnd(xvel[lev], dx, plo, geom[lev].Domain(), u_bc, u_bc_funcs);
        // FillPhysBCs_ynd(yvel[lev], dx, plo, geom[lev].Domain(), v_bc, v_bc_funcs);
    }
}