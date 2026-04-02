#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrCore.H>
#include <AMReX_Reduce.H>

#include "NavierStokesSolver.H"
#include "FillBC.H"
#include "BCFunc.H"

using namespace amrex;

NavierStokesSolver::NavierStokesSolver() : AmrCore()
{
    readParameters();
    // Print()<<"1"<<std::endl;

    int nlevels = max_level + 1;

    xvel.resize(nlevels);
    yvel.resize(nlevels);
    xvel_old.resize(nlevels);
    yvel_old.resize(nlevels);
    xvel_k1.resize(nlevels);
    yvel_k1.resize(nlevels);
    xvel_k2.resize(nlevels);
    yvel_k2.resize(nlevels);
    U.resize(nlevels);
    // U_old.resize(nlevels);
    U_magnitude.resize(nlevels);
    Pressure.resize(nlevels);
    // Pressure_old.resize(nlevels);
    divU.resize(nlevels);
    vorticity.resize(nlevels);
    // vorticity_old.resize(nlevels);
    //  flux_x.resize(nlevels);
    //  flux_y.resize(nlevels);
    //  flux_register.resize(max_level);

    // t_new.resize(nlevels, 0.0);
    // t_old.resize(nlevels, 0.0);
    // dt.resize(nlevels, dt);
    // istep.resize(nlevels, 0);

    // nsubsteps.resize(nlevels, 1);
    // for (int lev = 1; lev < nlevels; ++lev)
    // {
    //     nsubsteps[lev] = ref_ratio[lev - 1];  // e.g., 2 for 2x finer
    // }

    // Print()<<"2"<<std::endl;
    // Set BC types - Order: left, bottom, right, top
// In NavierStokesSolver constructor, update the boundary condition setup:

    // Set BC types - Order: left, bottom, right, top
    u_bc = {amrex::LinOpBCType::Dirichlet,  // left
            amrex::LinOpBCType::Dirichlet,  // bottom
            amrex::LinOpBCType::Dirichlet,  // right
            amrex::LinOpBCType::Dirichlet}; // top

    v_bc = {amrex::LinOpBCType::Dirichlet,  // left
            amrex::LinOpBCType::Dirichlet,  // bottom
            amrex::LinOpBCType::Dirichlet,  // right
            amrex::LinOpBCType::Dirichlet}; // top

    p_bc = {amrex::LinOpBCType::Neumann,    // left
            amrex::LinOpBCType::Neumann,    // bottom
            amrex::LinOpBCType::Neumann,    // right
            amrex::LinOpBCType::Neumann};   // top

    // Initialize boundary condition functors
    // // Order: left, bottom, right, top
    // amrex::Array<MovingWallBCFunc, 4> u_bc_func_local;
    // amrex::Array<VelocityBCFunc, 4> v_bc_func_local;
    // amrex::Array<PressureBCFunc, 4> p_bc_func_local;  

    // for (int i = 0; i < 4; i++)
    // {
    //     if (i == 3) { // top wall
    //         u_bc_funcs[i] = MovingWallBCFunc(1.0);  // Moving wall
    //     } else {
    //         u_bc_funcs[i] = MovingWallBCFunc(0.0);  // Stationary walls
    //     }
    //     v_bc_funcs[i] = VelocityBCFunc();
    //     p_bc_funcs[i] = PressureBCFunc();
    // }


    // // For now, we'll store these in the class as member variables
    // u_bc_funcs = u_bc_func_local;
    // v_bc_funcs = v_bc_func_local;
    // p_bc_funcs = p_bc_func_local;

    // Print()<<"5"<<std::endl;
}

void NavierStokesSolver::readParameters()
{

    BL_PROFILE("Read Parameters");
    ParmParse pp;

    // Time stepping
    pp.get("dt", dt);
    pp.get("Re", Re);

    Mu = (1.0 * 1.0 * 1.0) / Re;

    pp.get("maxSteps", maxSteps);
    pp.get("finalTime", finalTime);

    // Optional parameters
    pp.query("SaveEvery", SaveEvery);
    pp.query("steadyTol", steadyTol);
    pp.query("cfl_adv", cfl_adv);
    pp.query("cfl_diff", cfl_diff);

    // pp.query("do_subcycle", do_subcycle);
    // pp.query("do_reflux", do_reflux);

    // Refinement sensitivity
    pp.query("refine_threshold", refine_threshold);

    // AMR settings
    ParmParse pp_amr("amr");
    pp_amr.query("max_level", max_level);
    // Read refinement ratio - single set applies to all levels
    if (max_level > 0)
    {
        ref_ratio.resize(max_level);

        // Read the 2D/3D refinement ratio
        Vector<int> ref_comps(AMREX_SPACEDIM);
        pp_amr.getarr("ref_ratio", ref_comps, 0, AMREX_SPACEDIM);

        // Create IntVect and use for all levels
        IntVect ref_iv(ref_comps.data());
        for (int lev = 0; lev < max_level; ++lev)
        {
            ref_ratio[lev] = ref_iv;
        }
    }
    pp_amr.query("regrid_int", regrid_int);
    pp_amr.query("n_error_buf", n_error_buf);
    pp_amr.query("max_grid_size", max_grid_size);
    pp_amr.query("blocking_factor", blocking_factor);
}

NavierStokesSolver::~NavierStokesSolver() {}

void NavierStokesSolver::Evolve()
{
    // For constant refinement ratio across all levels
    // for (int lev = 0; lev <= max_level; ++lev)
    // {
    //     ref_ratio[lev] = amrex::IntVect(ref_ratio_vals[lev], ref_ratio_vals[lev]);
    // }
    Real residual;
    // readParameters();
    // Print()<<"begin"<<std::endl;
    while (continueRun())
    {
        // Print()<<"1"<<std::endl;

        // computeVorticity();

        // if (do_subcycle)
        //     timeStepWithSubcycling(0, curr_time, iteration);
        // else
        //     timeStepNoSubcycling(curr_time, iteration);

        old_time = curr_time;
        curr_time += dt;
        curr_iter++;

        // Print()<<"1"<<std::endl;

        if (max_level > 0 && regrid_int > 0) // We may need to regrid
        {
            if (curr_iter % regrid_int == 0)
            {
                regrid(0, curr_time);
            }
        }


        // copy new to old
        copyNewToOld();

        //RK Step 1

        // prediction
        predictVelocity_k1();

        for (int lev = 0; lev <= finest_level; lev++)
        {
            FillPatchVelocity(lev, 0, 1);
        }

        AverageDownVelocity();

        //RK Step 2

        // prediction
        predictVelocity_k2();

        for (int lev = 0; lev <= finest_level; lev++)
        {
            FillPatchVelocity(lev, 0, 1);
        }

        AverageDownVelocity();

        // poisson
        pressurePoisson();
        for (int lev = 0; lev <= finest_level; lev++)
        {
            FillPatchPressure(lev, 0, 1);
        }

        // correction
        correctVelocity();
        for (int lev = 0; lev <= finest_level; lev++)
        {
            FillPatchVelocity(lev, 0, 1);
        }
        AverageDownVelocity();
        AverageDownPressure();


        residual = CheckSteadyState();

        // Print()<<"1"<<std::endl;

        amrex::PrintToFile("log") << "Iteration: " << curr_iter
                                  << ", Residual: " << residual
                                  << ", Finest_level: " << finest_level
                                  << ", Time: " << curr_time
                                  << ", dt: " << dt
                                  << ", dt_adv: " << dt_adv
                                  << ", dt_diff: " << dt_diff
                                  << std::endl;

        if (curr_iter % SaveEvery == 0)
            WritePlotFile();

        if (curr_iter > 1 && residual < steadyTol)
        {
            amrex::PrintToFile("log") << "Steady state reached at iteration " << curr_iter << std::endl;
            WritePlotFile();
            break;
        }
        // Print()<<"1"<<std::endl;

        ComputeTimeStep();

        // Print()<<"1"<<std::endl;

        // amrex::Real stop_time = amrex::second() - start_time;
        // const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
        // amrex::ParallelDescriptor::ReduceRealMax(stop_time, IOProc);
        // amrex::PrintToFile("log") << "Iteration curr_time = " << stop_time << "\n\n";
    }
}

// void NavierStokesSolver::timeStepWithSubcycling(int lev, Real curr_time, int iteration)
// {

//     AMREX_ALWAYS_ASSERT(!xvel_old[lev].contains_nan());

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());

//     // Regridding logic
//     if (regrid_int > 0)
//     {
//         static Vector<int> last_regrid_step(max_level+1, 0);

//         if (lev < max_level && istep[lev] > last_regrid_step[lev])
//         {
//             if (istep[lev] % regrid_int == 0)
//             {
//                 int old_finest = finest_level;

//                 amrex::Print() << "DEBUG: Before regrid - finest_level = " << finest_level
//                             << ", old_finest = " << old_finest << std::endl;

//                 regrid(lev, curr_time);

//                 for (int k = lev; k <= finest_level; ++k) {
//                     last_regrid_step[k] = istep[k];
//                 }

//                 for (int k = old_finest+1; k <= finest_level; ++k) {
//                     dt[k] = dt[k-1] / MaxRefRatio(k-1);
//                 }

//                 for (int new_lev = old_finest + 1; new_lev <= finest_level; ++new_lev)
//                 {
//                     amrex::Print() << "DEBUG: Validating new level " << new_lev << std::endl;

//                     // === Allocate and initialize xvel_old and yvel_old ===
//                     int nghost = xvel[new_lev].nGrow(); // match ghost width
//                     xvel_old[new_lev].define(grids[new_lev], dmap[new_lev], 1, nghost);
//                     yvel_old[new_lev].define(grids[new_lev], dmap[new_lev], 1, nghost);

//                     MultiFab::Copy(xvel_old[new_lev], xvel[new_lev], 0, 0, 1, nghost);
//                     MultiFab::Copy(yvel_old[new_lev], yvel[new_lev], 0, 0, 1, nghost);

//                     xvel_old[new_lev].FillBoundary(geom[new_lev].periodicity());
//                     yvel_old[new_lev].FillBoundary(geom[new_lev].periodicity());

//                     const Real* dx_new = geom[new_lev].CellSize();
//                     const Real* plo_new = geom[new_lev].ProbLo();

//                     FillPhysBCs_xnd(xvel_old[new_lev], dx_new, plo_new, geom[new_lev].Domain(), u_bc, u_bc_funcs);
//                     FillPhysBCs_ynd(yvel_old[new_lev], dx_new, plo_new, geom[new_lev].Domain(), v_bc, v_bc_funcs);

//                     // === Validate ===
//                     AMREX_ALWAYS_ASSERT(!xvel[new_lev].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!yvel[new_lev].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!xvel_old[new_lev].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!yvel_old[new_lev].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!Pressure[new_lev].contains_nan());

//                     amrex::Print() << "DEBUG: New level " << new_lev << " validation passed" << std::endl;
//                 }

//                 for (int k = old_finest+1; k <= finest_level; ++k) {
//                     dt[k] = dt[k-1] / MaxRefRatio(k-1);
//                     if (dt[k] <= 0.0 || std::isnan(dt[k]) || dt[k] > 1e6)
//                         dt[k] = 1e-6;
//                 }

//                 for (int k = 0; k <= finest_level; ++k)
//                 {
//                     const Real* dx_k = geom[k].CellSize();
//                     const Real* plo_k = geom[k].ProbLo();

//                     FillPhysBCs_xnd(xvel[k], dx_k, plo_k, geom[k].Domain(), u_bc, u_bc_funcs);
//                     FillPhysBCs_ynd(yvel[k], dx_k, plo_k, geom[k].Domain(), v_bc, v_bc_funcs);

//                     xvel[k].FillBoundary(geom[k].periodicity());
//                     yvel[k].FillBoundary(geom[k].periodicity());

//                     AMREX_ALWAYS_ASSERT(!xvel[k].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!yvel[k].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!xvel_old[k].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!yvel_old[k].contains_nan());
//                     AMREX_ALWAYS_ASSERT(!Pressure[k].contains_nan());
//                 }
//             }
//         }
//     }

//     // Time update
//     t_old[lev] = t_new[lev];
//     t_new[lev] += dt[lev];
//     Real t_nph = t_old[lev] + 0.5 * dt[lev];

//     copyNewToOld();

//     predictVelocity();

//     Pressure[lev].FillBoundary(geom[lev].periodicity());
//     xvel[lev].FillBoundary(geom[lev].periodicity());
//     yvel[lev].FillBoundary(geom[lev].periodicity());

//     pressurePoisson();

//     Pressure[lev].FillBoundary(geom[lev].periodicity());
//     xvel[lev].FillBoundary(geom[lev].periodicity());
//     yvel[lev].FillBoundary(geom[lev].periodicity());

//     correctVelocity();

//     Pressure[lev].FillBoundary(geom[lev].periodicity());
//     xvel[lev].FillBoundary(geom[lev].periodicity());
//     yvel[lev].FillBoundary(geom[lev].periodicity());

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());

//     ++istep[lev];

// if (lev < finest_level)
// {
//     for (int i = 1; i <= nsubsteps[lev + 1]; ++i)
//     {
//         timeStepWithSubcycling(lev + 1, curr_time + (i - 1) * dt[lev + 1], i);
//     }

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());

//     if (flux_register[lev + 1])
//     {
//         flux_register[lev + 1]->FineAdd(xvel[lev + 1], 0, 0, 0, 1, 1.0);
//         flux_register[lev + 1]->FineAdd(yvel[lev + 1], 1, 0, 0, 1, 1.0);
//     }

//     // Only CrseAdd to the current level's register
//     if (flux_register[lev])
//     {
//         flux_register[lev]->CrseAdd(xvel[lev + 1], 0, 0, 0, 1, 1.0, geom[lev]);
//         flux_register[lev]->CrseAdd(yvel[lev + 1], 1, 0, 0, 1, 1.0, geom[lev]);
//     }

//     xvel[lev].FillBoundary(geom[lev].periodicity());
//     yvel[lev].FillBoundary(geom[lev].periodicity());

//     if (do_reflux && flux_register[lev + 1] && flux_register[lev])
//     {
//         const int ncomp = xvel[lev].nComp();

//         if (flux_x.size() <= lev) flux_x.resize(lev + 1);
//         if (flux_y.size() <= lev) flux_y.resize(lev + 1);

//         flux_x[lev].define(xvel[lev].boxArray(), xvel[lev].DistributionMap(), ncomp, 0);
//         flux_y[lev].define(yvel[lev].boxArray(), yvel[lev].DistributionMap(), ncomp, 0);

//         for (MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.validbox();
//             auto const& u = xvel[lev][mfi].const_array();
//             auto const& fx = flux_x[lev][mfi].array();

//             ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
//             {
//                 fx(i,j,k) = u(i,j,k);
//             });
//         }

//         for (MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.validbox();
//             auto const& v = yvel[lev][mfi].const_array();
//             auto const& fy = flux_y[lev][mfi].array();

//             ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
//             {
//                 fy(i,j,k) = v(i,j,k);
//             });
//         }

//         flux_x[lev].FillBoundary(geom[lev].periodicity());
//         flux_y[lev].FillBoundary(geom[lev].periodicity());

//         AMREX_ALWAYS_ASSERT(!flux_x[lev].contains_nan());
//         AMREX_ALWAYS_ASSERT(!flux_y[lev].contains_nan());

//         flux_register[lev + 1]->CrseInit(flux_x[lev], 0, 0, ncomp, 0);
//         flux_register[lev + 1]->CrseInit(flux_y[lev], 0, 0, ncomp, 1);

//         flux_register[lev + 1]->Reflux(xvel[lev], 0, 1.0, 0, ncomp, Geom(lev));
//         flux_register[lev + 1]->Reflux(yvel[lev], 0, 1.0, 0, ncomp, Geom(lev));
//         // Apply boundary conditions after refluxing
//         FillPhysBCs_xnd(xvel[lev], geom[lev].CellSize(), geom[lev].ProbLo(),
//                         geom[lev].Domain(), u_bc, u_bc_funcs);
//         FillPhysBCs_ynd(yvel[lev], geom[lev].CellSize(), geom[lev].ProbLo(),
//                         geom[lev].Domain(), v_bc, v_bc_funcs);
//     }

//     Pressure[lev].FillBoundary(geom[lev].periodicity());
//     xvel[lev].FillBoundary(geom[lev].periodicity());
//     yvel[lev].FillBoundary(geom[lev].periodicity());

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());

//     AverageDownTo(lev);

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());
// }

//     AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
//     AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());

// }

// void NavierStokesSolver::timeStepNoSubcycling(Real curr_time, int iteration)
// {
//     if (max_level > 0 && regrid_int > 0)
//     {
//         if (istep[0] % regrid_int == 0)
//         {
//             regrid(0, curr_time);
//         }
//     }

//     if (Verbose()) {
//         for (int lev = 0; lev <= finest_level; ++lev) {
//             amrex::Print() << "[Level " << lev << " step " << istep[lev] + 1 << "] ";
//             amrex::Print() << "ADVANCE with curr_time = " << t_new[lev]
//                            << " dt = " << dt[0] << std::endl;
//         }
//     }

//     for (int lev = 0; lev <= finest_level; ++lev) {
//         t_old[lev] = t_new[lev];
//         t_new[lev] += dt[0];
//     }

//     copyNewToOld();
//     predictVelocity();
//     pressurePoisson();
//     correctVelocity();

//     for (int lev = finest_level - 1; lev >= 0; --lev)
//     {
//         for (int dir = 0; dir <= 1; ++dir)
//         {
//             if (flux_register[lev + 1])
//             {
//                 if (dir == 0)
//                 {
//                     flux_register[lev + 1]->FineAdd(xvel[lev + 1], dir, 0, 0, xvel[lev + 1].nComp(), 1.0);
//                     flux_register[lev]->CrseAdd(xvel[lev + 1], dir, 0, 0, xvel[lev + 1].nComp(), 1.0, Geom(lev));
//                 }
//                 else
//                 {
//                     flux_register[lev + 1]->FineAdd(yvel[lev + 1], dir, 0, 0, yvel[lev + 1].nComp(), 1.0);
//                     flux_register[lev]->CrseAdd(yvel[lev + 1], dir, 0, 0, yvel[lev + 1].nComp(), 1.0, Geom(lev));
//                 }
//             }
//         }
//     }

//     // Apply flux corrections from fine to coarse
//     for (int lev = finest_level - 1; lev >= 0; --lev)
//     {
//         if (do_reflux && flux_register[lev])
//         {
//             flux_register[lev]->Reflux(xvel[lev], 0, 1.0, 0, xvel[lev].nComp(), Geom(lev));
//             flux_register[lev]->Reflux(yvel[lev], 0, 1.0, 0, yvel[lev].nComp(), Geom(lev));
//         }
//     }

//     AverageDown();

//     for (int lev = 0; lev <= finest_level; ++lev)
//         ++istep[lev];

// }

bool NavierStokesSolver::continueRun()
{
    // Print()<<maxSteps<<std::endl;
    if (maxSteps < 0)
    {
        // Print()<<maxSteps<<std::endl;
        return curr_time < finalTime;
    }
    else
    {
        // Print()<<curr_time<<std::endl;
        return (curr_iter < maxSteps && curr_time < finalTime);
    }
}

void NavierStokesSolver::InitialData()
{
    BL_PROFILE("Inial Data");
    // InitAmrCore();
    // Print()<<"1"<<std::endl;

    InitFromScratch(curr_time);

    // Print()<<"1"<<std::endl;

    // AverageDown();

    // Print()<<"1"<<std::endl;

    // Display(xvel[0], "xvel");
    // Display(yvel[0], "yvel");

    // Print()<<"1"<<std::endl;

    WritePlotFile();
    //Print()<<"1"<<std::endl;
}

void NavierStokesSolver::memoryallocation(int lev, const BoxArray &ba,
                                          const DistributionMapping &dm)
{
    BL_PROFILE("Memory Allocation");
    // Print()<<"1"<<std::endl;

    // Memory Allocation (Important to run on GPU)

    #ifdef AMREX_USE_GPU
        amrex::Print() << "GPU Memory is ENABLED" << std::endl;
        mf_info.SetArena(amrex::The_Device_Arena()); // Allocates memory on GPU
    #endif

    // Pressure (cell-centered)
    Pressure[lev].define(ba, dm, 1, 1, mf_info);
    // Pressure_old[lev].define(ba, dm, 1, 1);
    U[lev].define(ba, dm, AMREX_SPACEDIM, 0, mf_info);
    // U_old[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    U_magnitude[lev].define(ba, dm, 1, 0, mf_info);

    // u velocity (x-face centered)
    BoxArray ba_u = ba;
    ba_u.convert(IntVect::TheDimensionVector(0)); // staggered in x
    xvel[lev].define(ba_u, dm, 1, 1, mf_info);
    xvel_old[lev].define(ba_u, dm, 1, 1, mf_info);
    xvel_k1[lev].define(ba_u, dm, 1, 1, mf_info);
    xvel_k2[lev].define(ba_u, dm, 1, 1, mf_info);
    // flux_x[lev].define(xvel[lev].boxArray(), xvel[lev].DistributionMap(), xvel[lev].nComp(), 1);

    // v velocity (y-face centered)
    BoxArray ba_v = ba;
    ba_v.convert(IntVect::TheDimensionVector(1)); // staggered in y
    yvel[lev].define(ba_v, dm, 1, 1, mf_info);
    yvel_old[lev].define(ba_v, dm, 1, 1, mf_info);
    yvel_k1[lev].define(ba_v, dm, 1, 1, mf_info);
    yvel_k2[lev].define(ba_v, dm, 1, 1, mf_info);
    // flux_y[lev].define(yvel[lev].boxArray(), yvel[lev].DistributionMap(), yvel[lev].nComp(), 1);

    // Derived quantities
    divU[lev].define(ba, dm, 1, 0, mf_info);
    vorticity[lev].define(ba, dm, 1, 0, mf_info);
    // vorticity_old[lev].define(ba, dm, 1, 0);
    // Print()<<"1"<<std::endl;
}

void NavierStokesSolver::InitializeField(int lev, Real time)
{
    BL_PROFILE("Inialize Field");

    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
    const amrex::GpuArray<Real, AMREX_SPACEDIM> problo = geom[lev].ProbLoArray();

    // Initialize pressure (cell-centered)
    for (MFIter mfi(Pressure[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.validbox();
        auto arr = Pressure[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            // Always initialize to 0 for lid-driven cavity
            arr(i, j, k) = 0.0;
        });
    }

    // Initialize u-velocity (x-face centered)
    for (MFIter mfi(xvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.validbox();
        auto arr = xvel[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = problo[0] + i * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            // For lid-driven cavity, set u=1 only at top boundary
            // But during initialization, we can set all to 0
            arr(i, j, k) = 0.0;
        });
    }

    // Initialize v-velocity (y-face centered)
    for (MFIter mfi(yvel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.validbox();
        auto arr = yvel[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + j * dx[1];
            arr(i, j, k) = 0.0;
        });
    }

    // Initialize other fields
    U[lev].setVal(0.0);
    vorticity[lev].setVal(0.0);
    U_magnitude[lev].setVal(0.0);
    xvel_k1[lev].setVal(0.0);
    yvel_k1[lev].setVal(0.0);
    xvel_k2[lev].setVal(0.0);
    yvel_k2[lev].setVal(0.0);
    

    // Apply boundary conditions
    // Note: Boundary order is: 0=left, 1=bottom, 2=right, 3=top

    // In InitializeField function:
    FillPhysBCs_pressure(Pressure[lev], geom[lev], p_bc, p_bc_funcs);
    FillPhysBCs_xvel(xvel[lev], geom[lev], u_bc, u_bc_funcs);
    FillPhysBCs_yvel(yvel[lev], geom[lev], v_bc, v_bc_funcs);

    // if (lev == 0) {
    //     Display(xvel[0], "xvel at level 0");
    // }

    // Copy current fields to old
    copyNewToOld();

    // Fill ghost cells (uncomment if needed)
    // xvel[lev].FillBoundary(geom[lev].periodicity());
    // yvel[lev].FillBoundary(geom[lev].periodicity());
    // Pressure[lev].FillBoundary(geom[lev].periodicity());

    // Sanity checks (uncomment if needed)
    // AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
    // AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
    // AMREX_ALWAYS_ASSERT(!Pressure[lev].contains_nan());
}

void NavierStokesSolver::copyNewToOld()
{
    BL_PROFILE("Copying");

    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy(xvel_old[lev], xvel[lev], 0, 0, xvel[lev].nComp(), xvel[lev].nGrow());
        MultiFab::Copy(yvel_old[lev], yvel[lev], 0, 0, yvel[lev].nComp(), yvel[lev].nGrow());
        // MultiFab::Copy(U_old[lev], U[lev], 0, 0, U[lev].nComp(), U[lev].nGrow());
        // MultiFab::Copy(Pressure_old[lev], Pressure[lev], 0, 0, Pressure[lev].nComp(), Pressure[lev].nGrow());
        // MultiFab::Copy(vorticity_old[lev], vorticity[lev], 0, 0, vorticity[lev].nComp(), vorticity[lev].nGrow());
    }
}

void NavierStokesSolver::WritePlotFile()
{
    BL_PROFILE("WritePlotFile");

    computeCellCenteredVel();
    computeDivergence();
    computeMagnitude();
    for (int lev = 0; lev <= finest_level; lev++)
    {
        computeVorticity(lev);
    }

    int ncomp = AMREX_SPACEDIM + 4; // Ux, Uy, p, divU, (vorticity), (magnitude)

    Vector<std::string> varnames(ncomp);
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        varnames[d] = "U" + std::to_string(d);
    }
    varnames[AMREX_SPACEDIM] = "p";
    varnames[AMREX_SPACEDIM + 1] = "divU";
    varnames[AMREX_SPACEDIM + 2] = "vorticity";
    varnames[AMREX_SPACEDIM + 3] = "magnitude";

    Vector<MultiFab> plotmf(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // AMREX_ALWAYS_ASSERT(U[lev].nComp() >= AMREX_SPACEDIM);
        // AMREX_ALWAYS_ASSERT(Pressure[lev].nComp() >= 1);
        // AMREX_ALWAYS_ASSERT(divU[lev].nComp() >= 1);
        // AMREX_ALWAYS_ASSERT(U_magnitude[lev].nComp() >= 1);
        // AMREX_ALWAYS_ASSERT(vorticity[lev].nComp() >= 1);

        plotmf[lev].define(grids[lev], dmap[lev], ncomp, 0);

        MultiFab::Copy(plotmf[lev], U[lev], 0, 0, AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf[lev], Pressure[lev], 0, AMREX_SPACEDIM, 1, 0);
        MultiFab::Copy(plotmf[lev], divU[lev], 0, AMREX_SPACEDIM + 1, 1, 0);
        MultiFab::Copy(plotmf[lev], vorticity[lev], 0, AMREX_SPACEDIM + 2, 1, 0);
        MultiFab::Copy(plotmf[lev], U_magnitude[lev], 0, AMREX_SPACEDIM + 3, 1, 0);
    }

    std::string filename = Concatenate("Output/plt", curr_iter);
    amrex::Print() << "Writing plotfile to " << filename << " with finest_level = " << finest_level << "\n";

    WriteMultiLevelPlotfile(filename, finest_level + 1,
                            GetVecOfConstPtrs(plotmf), varnames,
                            Geom(), curr_time,
                            Vector<int>(1, curr_iter),
                            Vector<IntVect>(1, IntVect(1)));
}

void NavierStokesSolver::ComputeTimeStep()
{
    BL_PROFILE("Compute Time Step");
    dt_adv = std::numeric_limits<Real>::max();
    dt_diff = std::numeric_limits<Real>::max();
    Real global_dt_min = std::numeric_limits<Real>::max();

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();

        // Compute max velocities using built-in reductions (GPU-compatible)
        Real max_u_lev = xvel[lev].norm0(0, 0, true); 
        Real max_v_lev = yvel[lev].norm0(0, 0, true); 

        Real vel_scale = (max_u_lev / dx[0] + max_v_lev / dx[1] + 1e-15);
        Real dt_adv_lev = cfl_adv / vel_scale;

        Real deltax = std::min(dx[0], dx[1]);
        Real dt_diff_lev = cfl_diff * (deltax * deltax) / (4.0 * Mu);

        dt_adv = std::min(dt_adv, dt_adv_lev);
        dt_diff = std::min(dt_diff, dt_diff_lev);

        Real dt_lev = std::min(dt_adv_lev, dt_diff_lev);
        global_dt_min = std::min(global_dt_min, dt_lev);
    }

    // Parallel reduction across MPI ranks
    ParallelDescriptor::ReduceRealMin(dt_adv);
    ParallelDescriptor::ReduceRealMin(dt_diff);
    amrex::ParallelDescriptor::ReduceRealMin(global_dt_min);

    dt = global_dt_min;
}



Real NavierStokesSolver::CheckSteadyState()
{
    BL_PROFILE("Check Steady State");
    
    Real max_diff = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // // Fill ghost cells before stencil access
        // xvel[lev].FillBoundary(geom[lev].periodicity());
        // yvel[lev].FillBoundary(geom[lev].periodicity());
        // xvel_old[lev].FillBoundary(geom[lev].periodicity());
        // yvel_old[lev].FillBoundary(geom[lev].periodicity());

        // AMREX_ALWAYS_ASSERT(!xvel[lev].contains_nan());
        // AMREX_ALWAYS_ASSERT(!yvel[lev].contains_nan());
        // AMREX_ALWAYS_ASSERT(!xvel_old[lev].contains_nan());
        // AMREX_ALWAYS_ASSERT(!yvel_old[lev].contains_nan());
        
        for (MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();

            auto const &un = xvel[lev].const_array(mfi);
            auto const &uo = xvel_old[lev].const_array(mfi);
            auto const &vn = yvel[lev].const_array(mfi);
            auto const &vo = yvel_old[lev].const_array(mfi);

            // Use a Gpu::DeviceScalar for thread-safe reduction
            amrex::Gpu::DeviceScalar<Real> ds_local_max_diff(0.0);
            Real* p_local_max_diff = ds_local_max_diff.dataPtr();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real u_cc = 0.5 * (un(i,j,k) + un(i+1,j,k));
                Real uo_cc = 0.5 * (uo(i,j,k) + uo(i+1,j,k));
                Real v_cc = 0.5 * (vn(i,j,k) + vn(i,j+1,k));
                Real vo_cc = 0.5 * (vo(i,j,k) + vo(i,j+1,k));

                Real vel_mag_n = std::sqrt(u_cc*u_cc + v_cc*v_cc);
                Real vel_mag_o = std::sqrt(uo_cc*uo_cc + vo_cc*vo_cc);

                Real diff = (vel_mag_n > vel_mag_o) ? (vel_mag_n - vel_mag_o) : (vel_mag_o - vel_mag_n);  

                amrex::Gpu::Atomic::Max(p_local_max_diff, diff);
            });


            Real local_max_diff = ds_local_max_diff.dataValue();
            
            max_diff = std::max(max_diff, local_max_diff);
        }
    }

    // Parallel reduce max_diff across all MPI ranks
    amrex::ParallelDescriptor::ReduceRealMax(max_diff);

    // Uncomment for debugging
    // amrex::Print() << "Steady state check max_diff = " << max_diff << std::endl;

    return max_diff;
}

void NavierStokesSolver::cleanup()
{
    for (int lev = 0; lev < xvel.size(); ++lev) {
        xvel[lev].clear();
        yvel[lev].clear();
        xvel_old[lev].clear();
        yvel_old[lev].clear();
        U[lev].clear();
        U_magnitude[lev].clear();
        Pressure[lev].clear();
        divU[lev].clear();
        vorticity[lev].clear();
    }
    
    xvel.clear();
    yvel.clear();
    xvel_old.clear();
    yvel_old.clear();
    U.clear();
    U_magnitude.clear();
    Pressure.clear();
    divU.clear();
    vorticity.clear();
    
#ifdef AMREX_USE_GPU
    amrex::Gpu::synchronize();
#endif
}

// void NavierStokesSolver::Display(const MultiFab &mf, const std::string &name)
// {
//     amrex::Print() << "\n--- Displaying MultiFab: " << name << " ---\n";

//     // Create a host copy of the MultiFab for printing
//     MultiFab mf_host;
//     mf_host.define(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);
//     amrex::ParallelDescriptor::Barrier();
    
//     // Copy from device to host
//     mf_host.copy(mf, 0, 0, mf.nComp());
    
//     for (MFIter mfi(mf_host); mfi.isValid(); ++mfi)
//     {
//         const Box &bx = mfi.validbox();
//         const auto &fab = mf_host[mfi];
        
//         amrex::Print() << "Box " << mfi.index() << ": " << bx << "\n";
        
//         // 2D only - no k dimension
//         for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
//         {
//             for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
//             {
//                 IntVect iv(i, j);  // 2D IntVect
//                 // Print all components
//                 for (int n = 0; n < mf_host.nComp(); ++n)
//                 {
//                     amrex::Print() << fab(iv, n) << " ";
//                 }
//             }
//             amrex::Print() << "\n";
//         }
//         amrex::Print() << "\n";
//     }

//     amrex::Print() << "--- End of MultiFab: " << name << " ---\n";
// }