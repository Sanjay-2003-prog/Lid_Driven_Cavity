#include <NavierStokesSolver.H>
#include <FillBC.H>
//#include <AMReX_FluxRegister.H>

void NavierStokesSolver::correctVelocity()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real *dxinv = geom[lev].InvCellSize();
        const Real *dx = geom[lev].CellSize();
        const Real *plo = geom[lev].ProbLo();          

        // xvel
        for (MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi) 
        {
            const Box &bx = mfi.validbox();
            const auto &u = xvel[lev][mfi].array();
            const auto &p = Pressure[lev][mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                u(i,j,k) -= dxinv[0] * (p(i,j,k)-p(i-1,j,k));
            });
        }

        // yvel
        for (MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi) 
        {
            const Box &bx = mfi.validbox();
            const auto &v = yvel[lev][mfi].array();
            const auto &p = Pressure[lev][mfi].array();
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                v(i, j, k) -= dxinv[1] * (p(i,j,k)-p(i,j-1,k));
            });
        }


        // // Apply physical boundary conditions
        // FillPhysBCs_xnd(xvel[lev], dx, plo, geom[lev].Domain(), u_bc, u_bc_func);
        // FillPhysBCs_ynd(yvel[lev], dx, plo, geom[lev].Domain(), v_bc, v_bc_func);
        
    }
}