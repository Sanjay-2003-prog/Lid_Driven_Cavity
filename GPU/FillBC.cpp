#include "FillBC.H"
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Print.H>

void FillPhysBCs_pressure(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<PressureBCFunc, 4> &bc_funcs)  // Using PressureBCFunc
{
    mf.FillBoundary();
    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;
    bcs1.resize(1);
    
    // Set dummy boundary conditions (these aren't actually used)
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir);
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    // Create boundary function with nodality (0,0) for cell-centered pressure
    struct_FillPhysBCs<PressureBCFunc> bf(bctype, bc_funcs, amrex::IntVect(0, 0));
    amrex::GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>> gpu_bf(bf);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<PressureBCFunc>>> 
        physbcf(geom, bcs1, gpu_bf);
    
    // Apply boundary conditions
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}

void FillPhysBCs_xvel(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<MovingWallBCFunc, 4> &bc_funcs)
{
    mf.FillBoundary();
    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;
    bcs1.resize(1);
    
    // Set dummy boundary conditions
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir);
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    // Create boundary function with nodality (1,0) for x-velocity (staggered in x)
    struct_FillPhysBCs<MovingWallBCFunc> bf(bctype, bc_funcs, amrex::IntVect(1, 0));
    amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>> gpu_bf(bf);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<MovingWallBCFunc>>> 
        physbcf(geom, bcs1, gpu_bf);
    
    // Apply boundary conditions
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}

void FillPhysBCs_yvel(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<VelocityBCFunc, 4> &bc_funcs)
{
    mf.FillBoundary();
    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;
    bcs1.resize(1);
    
    // Set dummy boundary conditions
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir);
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    // Create boundary function with nodality (0,1) for y-velocity (staggered in y)
    struct_FillPhysBCs<VelocityBCFunc> bf(bctype, bc_funcs, amrex::IntVect(0, 1));
    amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>> gpu_bf(bf);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs<VelocityBCFunc>>> 
        physbcf(geom, bcs1, gpu_bf);
    
    // Apply boundary conditions
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}