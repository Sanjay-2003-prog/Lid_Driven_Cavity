#include <NavierStokesSolver.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>
#include <FillBC.H>

void NavierStokesSolver::pressurePoisson()
{
    // AMREX_ALWAYS_ASSERT(Pressure.size() > finest_level);
    // AMREX_ALWAYS_ASSERT(divU.size() > finest_level);
    // AMREX_ALWAYS_ASSERT(dt.size() > finest_level);
    // AMREX_ALWAYS_ASSERT(geom.size() > finest_level);
    // AMREX_ALWAYS_ASSERT(grids.size() > finest_level);
    // AMREX_ALWAYS_ASSERT(dmap.size() > finest_level);

    computeDivergence();

    // Fill ghost cells before solve
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        Pressure[lev].FillBoundary(geom[lev].periodicity());
        divU[lev].FillBoundary(geom[lev].periodicity());
    }

    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMaxCoarseningLevel(30);

    MLPoisson mlpoisson(geom, grids, dmap, info);
    mlpoisson.setMaxOrder(2);
    mlpoisson.setDomainBC({p_bc[0], p_bc[1]}, {p_bc[2], p_bc[3]});

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        AMREX_ALWAYS_ASSERT(Pressure[lev].ok());
        AMREX_ALWAYS_ASSERT(Pressure[lev].nComp() == 1);
        AMREX_ALWAYS_ASSERT(Pressure[lev].nGrow() >= 1);
        mlpoisson.setLevelBC(lev, &Pressure[lev]);
    }

    MLMG mlmg(mlpoisson);
    mlmg.setMaxIter(1000);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(0);

    mlmg.solve(GetVecOfPtrs(Pressure), GetVecOfConstPtrs(divU), 1.e-12, 0.0);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real* dx = geom[lev].CellSize();
        const Real* plo = geom[lev].ProbLo();
        FillPhysBCs_pressure(Pressure[lev], geom[lev], p_bc, p_bc_func);

        // Final ghost sync after physical BCs
        Pressure[lev].FillBoundary(geom[lev].periodicity());
        AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());
    }
}