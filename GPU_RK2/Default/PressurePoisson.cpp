#include <NavierStokesSolver.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>
#include <AMReX_Hypre.H>
#include <AMReX_ParmParse.H>
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
    info.setMaxCoarseningLevel(5); // Make it zero for Hypre

    MLPoisson mlpoisson(geom, grids, dmap, info);
    mlpoisson.setMaxOrder(2);
    mlpoisson.setDomainBC({p_bc[0], p_bc[1]}, {p_bc[2], p_bc[3]});

    for (int lev = 0; lev <= finest_level; ++lev)
    {
//        AMREX_ALWAYS_ASSERT(Pressure[lev].ok());
//        AMREX_ALWAYS_ASSERT(Pressure[lev].nComp() == 1);
//        AMREX_ALWAYS_ASSERT(Pressure[lev].nGrow() >= 1);
        mlpoisson.setLevelBC(lev, &Pressure[lev]);
    }

    MLMG mlmg(mlpoisson);
    mlmg.setMaxIter(1000);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(0);
    mlmg.setBottomVerbose(0); // For Hypre
#ifdef AMREX_USE_HYPRE
    //amrex::Print()<<"Hypre Enabled"<< std::endl;
    mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
    
    ParmParse pp_hypre("hypre");
    
    // Choose Hypre solver type
    pp_hypre.add("solver", "pcg");           // Options: pcg, gmres, flexgmres, lgmr, etc.
    
    // Choose preconditioner
    pp_hypre.add("preconditioner", "boomeramg"); // Options: boomeramg, pilut, parasails, euclid, etc.
    
    // BoomerAMG specific settings (if used)
    pp_hypre.add("boomeramg_max_levels", 25);     // Max AMG levels
    pp_hypre.add("boomeramg_coarsen_type", 3);    // 6 = HMIS coarsening
    pp_hypre.add("boomeramg_interp_type", 6);      // 6 = extended+i interpolation
    pp_hypre.add("boomeramg_relax_type", 3);       // 3 = hybrid symmetric Gauss-Seidel
    pp_hypre.add("boomeramg_num_sweeps", 2);       // Number of relaxation sweeps
    
    pp_hypre.add("verbose", 0);                    // Hypre's own verbosity
#endif
#ifdef AMREX_USE_PETSC
    mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
#endif

    mlmg.solve(GetVecOfPtrs(Pressure), GetVecOfConstPtrs(divU), 1.e-12, 0.0);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
        // const amrex::GpuArray<Real, AMREX_SPACEDIM> plo = geom[lev].ProbLoArray();
        FillPhysBCs_pressure(Pressure[lev], geom[lev], p_bc, p_bc_funcs);

        // Final ghost sync after physical BCs
        // Pressure[lev].FillBoundary(geom[lev].periodicity());
//        AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());
    }
}
