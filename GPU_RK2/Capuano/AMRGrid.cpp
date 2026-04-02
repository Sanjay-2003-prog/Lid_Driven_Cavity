#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Math.H>

#include "NavierStokesSolver.H"
#include "FillBC.H"

using namespace amrex;

// Create a new level (initial or regridded)
void NavierStokesSolver::MakeNewLevelFromScratch(int lev, Real time,
                                                 const BoxArray &ba,
                                                 const DistributionMapping &dm)
{
    memoryallocation(lev, ba, dm);
    InitializeField(lev, time);
}

void NavierStokesSolver::MakeNewLevelFromCoarse(int lev, Real time,
                                                const BoxArray &ba,
                                                const DistributionMapping &dm)
{
    // Allocate memory for all fields at this level
    memoryallocation(lev, ba, dm);

    FillCoarsePatch(lev, time, Pressure[lev], Pressure[lev - 1], 0, Pressure[lev].nComp());
    FillCoarsePatch(lev, time, xvel[lev], xvel[lev - 1], 0, xvel[lev].nComp());
    FillCoarsePatch(lev, time, yvel[lev], yvel[lev - 1], 0, yvel[lev].nComp());
}

void NavierStokesSolver::RemakeLevel(int lev, Real time,
                                     const BoxArray &ba,
                                     const DistributionMapping &dm)
{
    // Create staggered BoxArrays for face-centered data
    BoxArray xba = ba;
    xba.convert(IntVect::TheDimensionVector(0)); // x-face centered
    BoxArray yba = ba;
    yba.convert(IntVect::TheDimensionVector(1)); // y-face centered

    // --- Step 1: Create all temporaries ---

    // Pressure (cell-centered)
    MultiFab tmp_Pressure(ba, dm, Pressure[lev].nComp(), Pressure[lev].nGrow(), mf_info);
    tmp_Pressure.setVal(0.0);

    // Velocity (face-centered)
    MultiFab tmp_xvel(xba, dm, xvel[lev].nComp(), xvel[lev].nGrow(), mf_info);
    MultiFab tmp_yvel(yba, dm, yvel[lev].nComp(), yvel[lev].nGrow(), mf_info);
    tmp_xvel.setVal(0.0);
    tmp_yvel.setVal(0.0);

    // Old velocity (for time integration)
    MultiFab tmp_xvel_old(xba, dm, xvel_old[lev].nComp(), xvel_old[lev].nGrow(), mf_info);
    MultiFab tmp_yvel_old(yba, dm, yvel_old[lev].nComp(), yvel_old[lev].nGrow(), mf_info);
    tmp_xvel_old.setVal(0.0);
    tmp_yvel_old.setVal(0.0);

    // Other fields (cell-centered)
    MultiFab tmp_U(ba, dm, U[lev].nComp(), U[lev].nGrow(), mf_info);
    MultiFab tmp_divU(ba, dm, divU[lev].nComp(), divU[lev].nGrow(), mf_info);
    MultiFab tmp_vorticity(ba, dm, vorticity[lev].nComp(), vorticity[lev].nGrow(), mf_info);

    tmp_U.setVal(0.0);
    tmp_divU.setVal(0.0);
    tmp_vorticity.setVal(0.0);

    // --- Step 2: Fill data using FillPatch ---

    if (lev == 0)
    {
        // Level 0: only fine data for pressure (scalar)
        Vector<MultiFab *> fmf_pressure = {&Pressure[lev]};
        Vector<Real> ft_pressure = {time};
        Vector<MultiFab *> cmf_empty;
        Vector<Real> ct_empty;

        FillPatch(lev, time, tmp_Pressure, cmf_empty, ct_empty, fmf_pressure, ft_pressure, 0, Pressure[lev].nComp());

        // For velocities at level 0, use array-based FillPatch
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};

        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fmf_vel;
        fmf_vel.resize(1);
        fmf_vel[0] = {&xvel[lev], &yvel[lev]};
        Vector<Real> ft_vel = {time};

        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cmf_vel_empty;
        Vector<Real> ct_vel_empty;

        // Use MovingWallBCFunc for u-velocity, VelocityBCFunc for v-velocity
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> u_bf(
            struct_FillPhysBCs<MovingWallBCFunc>{u_bc, u_bc_funcs, amrex::IntVect(AMREX_D_DECL(1, 0, 0))});
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> v_bf(
            struct_FillPhysBCs<VelocityBCFunc>{v_bc, v_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 1, 0))});

        FillPatch(lev, time, tmp_vel, cmf_vel_empty, ct_vel_empty, fmf_vel, ft_vel, 0, xvel[lev].nComp(), u_bf, v_bf);
    }
    else
    {
        // Level > 0: both coarse and fine data for pressure (scalar)
        Vector<MultiFab *> cmf_pressure = {&Pressure[lev - 1]};
        Vector<Real> ct_pressure = {time};
        Vector<MultiFab *> fmf_pressure = {&Pressure[lev]};
        Vector<Real> ft_pressure = {time};

        FillPatch(lev, time, tmp_Pressure, cmf_pressure, ct_pressure, fmf_pressure, ft_pressure, 0, Pressure[lev].nComp());

        // For velocities at level > 0, use array-based FillPatch with coarse and fine data
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};

        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cmf_vel;
        cmf_vel.resize(1);
        cmf_vel[0] = {&xvel[lev - 1], &yvel[lev - 1]};
        Vector<Real> ct_vel = {time};

        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fmf_vel;
        fmf_vel.resize(1);
        fmf_vel[0] = {&xvel[lev], &yvel[lev]};
        Vector<Real> ft_vel = {time};

        // Use MovingWallBCFunc for u-velocity, VelocityBCFunc for v-velocity
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> u_bf(
            struct_FillPhysBCs<MovingWallBCFunc>{u_bc, u_bc_funcs, amrex::IntVect(AMREX_D_DECL(1, 0, 0))});
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> v_bf(
            struct_FillPhysBCs<VelocityBCFunc>{v_bc, v_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 1, 0))});

        FillPatch(lev, time, tmp_vel, cmf_vel, ct_vel, fmf_vel, ft_vel, 0, xvel[lev].nComp(), u_bf, v_bf);
    }

    // --- Step 3: Apply physical boundary conditions ---
    // Use PressureBCFunc for pressure
    FillPhysBCs_pressure(tmp_Pressure, geom[lev], p_bc, p_bc_funcs);
    // Use MovingWallBCFunc for x-velocity
    FillPhysBCs_xvel(tmp_xvel, geom[lev], u_bc, u_bc_funcs);
    // Use VelocityBCFunc for y-velocity
    FillPhysBCs_yvel(tmp_yvel, geom[lev], v_bc, v_bc_funcs);

    // Fill boundary ghost cells for periodicity
    tmp_Pressure.FillBoundary(geom[lev].periodicity());
    tmp_xvel.FillBoundary(geom[lev].periodicity());
    tmp_yvel.FillBoundary(geom[lev].periodicity());

    // --- Step 4: Swap temporaries with existing data ---

    std::swap(tmp_Pressure, Pressure[lev]);
    std::swap(tmp_xvel, xvel[lev]);
    std::swap(tmp_yvel, yvel[lev]);

    std::swap(tmp_xvel_old, xvel_old[lev]);
    std::swap(tmp_yvel_old, yvel_old[lev]);

    std::swap(tmp_U, U[lev]);
    std::swap(tmp_divU, divU[lev]);
    std::swap(tmp_vorticity, vorticity[lev]);
}

// Error estimation for AMR tagging
void NavierStokesSolver::ErrorEst(int lev, TagBoxArray &tags, Real time, int ngrow)
{
    computeVorticity(lev);

    static bool first = true;
    static Vector<Real> vorterr;

    // Only read thresholds once
    if (first)
    {
        first = false;
        ParmParse pp("amr");
        int n = pp.countval("vorterr");
        if (n > 0)
        {
            pp.getarr("vorterr", vorterr, 0, n);
        }
    }

    if (lev >= vorterr.size())
        return;

    const int tagval = TagBox::SET;
    const MultiFab &vort = vorticity[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vort, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        const auto vortfab = vort.array(mfi);
        const auto tagfab = tags.array(mfi);
        Real threshold = vorterr[lev];

        // Threshold-based tagging
        amrex::ParallelFor(bx,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                           {
                               if (amrex::Math::abs(vortfab(i, j, k)) > threshold)
                               {
                                   tagfab(i, j, k) = tagval;
                               }
                           });
    }
}

void NavierStokesSolver::ClearLevel(int lev)
{
    xvel[lev].clear();
    yvel[lev].clear();
    xvel_old[lev].clear();
    yvel_old[lev].clear();
    U[lev].clear();
    Pressure[lev].clear();
    divU[lev].clear();
    vorticity[lev].clear();
}

// Average down from fine to coarse

void NavierStokesSolver::AverageDownVelocity()
{
    for (int lev = finest_level; lev > 0; lev--)
    {
        amrex::Array<const amrex::MultiFab *, AMREX_SPACEDIM> vel = {&xvel[lev], &yvel[lev]};
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> velc = {&xvel[lev - 1], &yvel[lev - 1]};
        amrex::average_down_faces(vel, velc, ref_ratio[lev - 1], geom[lev - 1]);
    }
}

void NavierStokesSolver::AverageDownPressure()
{
    for (int lev = finest_level - 1; lev >= 0; --lev)
    {
        amrex::average_down(Pressure[lev + 1], Pressure[lev],
                            geom[lev + 1], geom[lev],
                            0, Pressure[lev].nComp(), ref_ratio[lev]);
    }
}

void NavierStokesSolver::FillPatchVelocity(int lev, int icomp, int ncomp)
{
    if (lev > 0)
    {
        // Create temporary face-centered velocity fields
        BoxArray ba_u = grids[lev];
        ba_u.convert(IntVect::TheDimensionVector(0)); // staggered in x
        amrex::MultiFab tmp_xvel(ba_u, dmap[lev], ncomp, 1, mf_info);

        BoxArray ba_v = grids[lev];
        ba_v.convert(IntVect::TheDimensionVector(1)); // staggered in y
        amrex::MultiFab tmp_yvel(ba_v, dmap[lev], ncomp, 1, mf_info);

        tmp_xvel.setVal(0.0);
        tmp_yvel.setVal(0.0);

        // Create arrays for FillPatch (using the array-based version)
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};

        // Prepare coarse and fine data arrays
        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cmf_vel;
        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fmf_vel;
        amrex::Vector<Real> ct_vel, ft_vel;

        cmf_vel.resize(1);
        fmf_vel.resize(1);
        ct_vel.resize(1);
        ft_vel.resize(1);

        // Coarse data from lev-1
        cmf_vel[0] = {&xvel[lev - 1], &yvel[lev - 1]};
        ct_vel[0] = curr_time;

        // Fine data from lev (source for interpolation)
        fmf_vel[0] = {&xvel[lev], &yvel[lev]};
        ft_vel[0] = curr_time;

        // Boundary condition functors - different types for u and v
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> u_bf(
            struct_FillPhysBCs<MovingWallBCFunc>{u_bc, u_bc_funcs, amrex::IntVect(AMREX_D_DECL(1, 0, 0))});
        amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> v_bf(
            struct_FillPhysBCs<VelocityBCFunc>{v_bc, v_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 1, 0))});

        // Call the array-based FillPatch for velocities
        FillPatch(lev, curr_time, tmp_vel,
                  cmf_vel, ct_vel,
                  fmf_vel, ft_vel,
                  icomp, ncomp,
                  u_bf, v_bf);

        // Apply physical boundary conditions with appropriate functor types
        FillPhysBCs_xvel(tmp_xvel, geom[lev], u_bc, u_bc_funcs);
        FillPhysBCs_yvel(tmp_yvel, geom[lev], v_bc, v_bc_funcs);

        // Swap temporaries with permanent storage
        std::swap(tmp_xvel, xvel[lev]);
        std::swap(tmp_yvel, yvel[lev]);
    }
    else
    {
        // For lev == 0, just apply boundary conditions
        FillPhysBCs_xvel(xvel[lev], geom[lev], u_bc, u_bc_funcs);
        FillPhysBCs_yvel(yvel[lev], geom[lev], v_bc, v_bc_funcs);
    }
}

void NavierStokesSolver::FillPatchPressure(int lev, int icomp, int ncomp)
{
    if (lev > 0)
    {
        // Create temporary pressure field
        amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], ncomp, 1, mf_info);
        tmp_Pressure.setVal(0.0);

        // Prepare coarse and fine data vectors
        Vector<MultiFab *> cmf_pressure = {&Pressure[lev - 1]};
        Vector<Real> ct_pressure = {curr_time};
        Vector<MultiFab *> fmf_pressure = {&Pressure[lev]};
        Vector<Real> ft_pressure = {curr_time};

        // Call the scalar FillPatch for pressure
        FillPatch(lev, curr_time, tmp_Pressure,
                  cmf_pressure, ct_pressure,
                  fmf_pressure, ft_pressure,
                  icomp, ncomp);

        // Apply physical boundary conditions with PressureBCFunc
        FillPhysBCs_pressure(tmp_Pressure, geom[lev], p_bc, p_bc_funcs);

        // Swap with permanent storage
        std::swap(tmp_Pressure, Pressure[lev]);
    }
    else
    {
        // For lev == 0, just apply boundary conditions
        FillPhysBCs_pressure(Pressure[lev], geom[lev], p_bc, p_bc_funcs);
    }
}

void NavierStokesSolver::FillPatch(
    int lev,
    Real time,
    amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> const &mf,                 // dest mf
    const amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> &cmf, // coarse mf
    const amrex::Vector<amrex::Real> &ct,                                      // coarse time
    const amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> &fmf, // fine mf
    const amrex::Vector<amrex::Real> &ft,                                      // fine time
    int icomp,
    int ncomp,
    amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> phys_bc_fill_u,
    amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> phys_bc_fill_v)
{
    // No need to fill external boundary values
    amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM> bcs;
    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        bcs[n].resize(ncomp);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            for (int ic = 0; ic < ncomp; ++ic)
            {
                bcs[n][ic].setLo(idim, amrex::BCType::int_dir);
                bcs[n][ic].setHi(idim, amrex::BCType::int_dir);
            }
        }
    }

    if (lev == 0)
    {
        amrex::PhysBCFunctNoOp physbc;
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            amrex::Vector<amrex::MultiFab *> fmf_tmp(fmf.size());
            for (size_t i = 0; i < fmf_tmp.size(); i++)
            {
                fmf_tmp[i] = fmf[i][dir];
            }
            amrex::FillPatchSingleLevel(*mf[dir], time, fmf_tmp, ft, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        amrex::Interpolater *mapper = &amrex::face_divfree_interp;

        amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM> bcs1;

        for (int n = 0; n < AMREX_SPACEDIM; ++n)
        {
            bcs1[n].resize(ncomp);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                for (int ic = 0; ic < ncomp; ++ic)
                {
                    bcs1[n][ic].setLo(idim, amrex::BCType::ext_dir); // dummy
                    bcs1[n][ic].setHi(idim, amrex::BCType::ext_dir); // dummy
                }
            }
        }

        // Use different types for u and v
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>>
            cphysbc_u(geom[lev - 1], bcs1[0], phys_bc_fill_u);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>>
            cphysbc_v(geom[lev - 1], bcs1[1], phys_bc_fill_v);

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>>
            fphysbc_u(geom[lev], bcs1[0], phys_bc_fill_u);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>>
            fphysbc_v(geom[lev], bcs1[1], phys_bc_fill_v);

        // Create arrays of the correct types
        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>>, AMREX_SPACEDIM>
            cphysbc_arr_u = {cphysbc_u, cphysbc_u}; // Note: This might need adjustment
        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>>, AMREX_SPACEDIM>
            cphysbc_arr_v = {cphysbc_v, cphysbc_v};

        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>>, AMREX_SPACEDIM>
            fphysbc_arr_u = {fphysbc_u, fphysbc_u};
        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>>, AMREX_SPACEDIM>
            fphysbc_arr_v = {fphysbc_v, fphysbc_v};

        // This part needs to be handled carefully - you might need to call FillPatchTwoLevels separately for u and v
        // For simplicity, consider using the BCFunc approach instead of separate types
    }
}

void NavierStokesSolver::FillPatch(
    int lev,
    Real time,
    MultiFab &mf,                  // dest mf
    const Vector<MultiFab *> &cmf, // coarse mf
    const Vector<Real> &ct,        // coarse time
    const Vector<MultiFab *> &fmf, // fine mf
    const Vector<Real> &ft,        // fine time
    int icomp,
    int ncomp)
{
    const IndexType &ixt = mf.boxArray().ixType();

    // Select appropriate interpolator based on field type
    Interpolater *mapper = nullptr;
    if (ixt.cellCentered())
    {
        mapper = &cell_cons_interp;
    }
    else if (ixt.nodeCentered(0) || ixt.nodeCentered(1))
    {
        mapper = &face_linear_interp;
    }
    else
    {
        amrex::Abort("FillPatch: unsupported index type");
    }

    // Boundary conditions - all internal dirichlet for interior filling
    Vector<BCRec> bcs(ncomp);
    for (int n = 0; n < ncomp; ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[n].setLo(idim, BCType::int_dir);
            bcs[n].setHi(idim, BCType::int_dir);
        }
    }

    if (lev == 0)
    {
        // Single level fill - use fine data only
        // Use PressureBCFunc for pressure fields
        GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>> dummy_bf(
            struct_FillPhysBCs<PressureBCFunc>{p_bc, p_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 0, 0))});
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> physbc(geom[lev], bcs, dummy_bf);

        amrex::FillPatchSingleLevel(mf, time, fmf, ft,
                                    0, icomp, ncomp, geom[lev],
                                    physbc, 0);
    }
    else
    {
        // Two-level fill with interpolation
        // Use PressureBCFunc for pressure fields
        GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>> dummy_bf(
            struct_FillPhysBCs<PressureBCFunc>{p_bc, p_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 0, 0))});
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> cphysbc(geom[lev - 1], bcs, dummy_bf);
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> fphysbc(geom[lev], bcs, dummy_bf);

        amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                  0, icomp, ncomp,
                                  geom[lev - 1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  ref_ratio[lev - 1], mapper,
                                  bcs, 0);
    }
}

void NavierStokesSolver::FillCoarsePatch(
    int lev,
    Real time,
    MultiFab &mf,        // dest mf
    const MultiFab &cmf, // coarse mf
    int icomp,
    int ncomp)
{
    BL_ASSERT(lev > 0);

    const IndexType &ixt = mf.boxArray().ixType();

    // Select appropriate interpolator based on field type
    Interpolater *mapper = nullptr;
    if (ixt.cellCentered())
    {
        mapper = &cell_cons_interp;
    }
    else if (ixt.nodeCentered(0) || ixt.nodeCentered(1))
    {
        // For face-centered data, use face_linear_interp for FillCoarsePatch
        mapper = &face_linear_interp;
    }
    else
    {
        amrex::Abort("FillCoarsePatch: unsupported index type");
    }

    // Set boundary conditions - use ext_dir as dummy for interpolation
    Vector<BCRec> bcs(ncomp);
    for (int n = 0; n < ncomp; ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[n].setLo(idim, BCType::ext_dir);
            bcs[n].setHi(idim, BCType::ext_dir);
        }
    }

    // Perform interpolation with appropriate boundary condition functors
    if (&mf == &xvel[lev])
    {
        // Use MovingWallBCFunc for x-velocity
        GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> u_bf(
            struct_FillPhysBCs<MovingWallBCFunc>{u_bc, u_bc_funcs, amrex::IntVect(AMREX_D_DECL(1, 0, 0))});
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>> cphysbc(geom[lev - 1], bcs, u_bf);
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>> fphysbc(geom[lev], bcs, u_bf);

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp,
                                     geom[lev - 1], geom[lev],
                                     cphysbc, 0, fphysbc, 0,
                                     ref_ratio[lev - 1], mapper, bcs, 0);
    }
    else if (&mf == &yvel[lev])
    {
        // Use VelocityBCFunc for y-velocity
        GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> v_bf(
            struct_FillPhysBCs<VelocityBCFunc>{v_bc, v_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 1, 0))});
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>> cphysbc(geom[lev - 1], bcs, v_bf);
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>> fphysbc(geom[lev], bcs, v_bf);

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp,
                                     geom[lev - 1], geom[lev],
                                     cphysbc, 0, fphysbc, 0,
                                     ref_ratio[lev - 1], mapper, bcs, 0);
    }
    else
    {
        // For non-velocity fields (pressure, etc.) use PressureBCFunc
        GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>> dummy_bf(
            struct_FillPhysBCs<PressureBCFunc>{p_bc, p_bc_funcs, amrex::IntVect(AMREX_D_DECL(0, 0, 0))});
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> cphysbc(geom[lev - 1], bcs, dummy_bf);
        PhysBCFunct<GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> fphysbc(geom[lev], bcs, dummy_bf);

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp,
                                     geom[lev - 1], geom[lev],
                                     cphysbc, 0, fphysbc, 0,
                                     ref_ratio[lev - 1], mapper, bcs, 0);
    }
}

std::string NavierStokesSolver::getFieldName(MultiFab &mf, int lev)
{
    if (&mf == &xvel[lev])
        return "xvel";
    if (&mf == &xvel_old[lev])
        return "xvel_old";
    if (&mf == &yvel[lev])
        return "yvel";
    if (&mf == &yvel_old[lev])
        return "yvel_old";
    if (&mf == &U[lev])
        return "U";
    if (&mf == &Pressure[lev])
        return "Pressure";
    if (&mf == &divU[lev])
        return "divU";
    if (&mf == &vorticity[lev])
        return "vorticity";

    return "unknown";
}