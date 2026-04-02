#ifndef _fillbc_H_
#define _fillbc_H_

#include <AMReX_MultiFab.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_BCRec.H>
#include "defs.H"
#include "BCFunc.H"

// Template functor that holds boundary condition functors
template <typename BC_FUNC>
struct struct_FillPhysBCs
{
    // Store boundary condition functors
    BC_FUNC bc_func[4];
    amrex::LinOpBCType bc_type[4];
    amrex::IntVect nodality;

    // Constructor - takes array of functors
    AMREX_GPU_HOST
    struct_FillPhysBCs(
        const amrex::Array<amrex::LinOpBCType, 4>& bctype,
        const amrex::Array<BC_FUNC, 4>& bc_funcs,
        const amrex::IntVect &nd)
        : nodality(nd)
    {
        // Copy boundary types and functors
        for (int i = 0; i < 4; i++) {
            bc_type[i] = bctype[i];
            bc_func[i] = bc_funcs[i];
        }
    }

    // Default constructor needed for GPU
    AMREX_GPU_HOST_DEVICE
    struct_FillPhysBCs() = default;

    // Copy constructor for GPU
    AMREX_GPU_HOST_DEVICE
    struct_FillPhysBCs(const struct_FillPhysBCs& other) = default;

    // Main operator() called by AMReX for each ghost cell
    AMREX_GPU_DEVICE
    void operator()(const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
                    const int dcomp, const int numcomp,
                    amrex::GeometryData const &geom, const amrex::Real time,
                    const amrex::BCRec *bcr, const int bcomp,
                    const int orig_comp) const
    {
        // Ensure we're only writing to the first component
        if (dcomp != 0) return;
        
        // Call the appropriate fill function based on nodality
        if (nodality == amrex::IntVect(0, 0)) {
            fill_cc(iv, f, geom, time);
        }
        else if (nodality == amrex::IntVect(1, 0)) {
            fill_xnd(iv, f, geom, time);
        }
        else if (nodality == amrex::IntVect(0, 1)) {
            fill_ynd(iv, f, geom, time);
        }
    }

    //(AMREX_D_DECL(x, y, 0.0), time);

    // Fill x-nodal (staggered in x) - for u-velocity
    AMREX_GPU_DEVICE
    void fill_xnd(const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
                  amrex::GeometryData const &geom, const amrex::Real time) const
    {
        const amrex::Box& domain = geom.Domain();
        const amrex::Real* dx = geom.CellSize();
        const amrex::Real* plo = geom.ProbLo();
        const amrex::Real* phi = geom.ProbHi();
        
        int i = iv[0];
        int j = iv[1];
        
        // Left boundary (x-direction)
        if (i <= domain.smallEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = plo[0];
            
            // Call the stored functor
            amrex::Real bval = bc_func[0](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
            {
                if (i == domain.smallEnd(0))
                {
                    f(i,j,0) = bval;
                } else {
                    f(i,j,0) = 2.0 * bval - f(i+2, j, 0);
                }
            }
            else if (bc_type[0] == amrex::LinOpBCType::Neumann)
            {
                if (i == domain.smallEnd(0))
                {
                    
                } 
                else 
                {
                    f(i,j,0) = f(i+2, j, 0) - 2.0* bval * dx[0];
                }
            }
        }
        
        // Right boundary (x-direction)
        if (i > domain.bigEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = phi[0];
            
            amrex::Real bval = bc_func[2](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
            {
                if (i == domain.bigEnd(0) + 1) {
                    f(i,j,0) = bval;
                } else {
                    f(i,j,0) = 2.0 * bval - f(i-2, j, 0);
                }
            }
            else if (bc_type[2] == amrex::LinOpBCType::Neumann)
            {
                if (i == domain.bigEnd(0) + 1)
                {
                    
                } else {
                    f(i,j,0) = f(i-2, j, 0) + 2.0 * bval * dx[0];
                }
            }
        }
        
        // Bottom boundary (y-direction) - for x-nodal, this is cell-centered in y
        if (j < domain.smallEnd(1))
        {
            amrex::Real y = plo[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[1](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i, j+1, 0);
            }
            else if (bc_type[1] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i, j+1, 0) - bval * dx[1];
            }
        }
        
        // Top boundary (y-direction)
        if (j > domain.bigEnd(1))
        {
            amrex::Real y = phi[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[3](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i, j-1, 0);
            }
            else if (bc_type[3] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i, j-1, 0) + bval * dx[1];
            }
        }
    }

    // Fill y-nodal (staggered in y) - for v-velocity
    AMREX_GPU_DEVICE
    void fill_ynd(const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
                  amrex::GeometryData const &geom, const amrex::Real time) const
    {
        const amrex::Box& domain = geom.Domain();
        const amrex::Real* dx = geom.CellSize();
        const amrex::Real* plo = geom.ProbLo();
        const amrex::Real* phi = geom.ProbHi();
        
        int i = iv[0];
        int j = iv[1];
        
        // Left boundary (x-direction) - for y-nodal, this is cell-centered in x
        if (i < domain.smallEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = plo[0];
            
            amrex::Real bval = bc_func[0](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i+1, j, 0);
            }
            else if (bc_type[0] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i+1, j, 0) - bval * dx[0];
            }
        }
        
        // Bottom boundary (y-direction)
        if (j <= domain.smallEnd(1))
        {
            amrex::Real y = plo[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[1](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
            {
                if (j == domain.smallEnd(1)) {
                    f(i,j,0) = bval;
                } else {
                    f(i,j,0) = 2.0 * bval - f(i, j+2, 0);
                }
            }
            else if (bc_type[1] == amrex::LinOpBCType::Neumann)
            {
                if (j == domain.smallEnd(1)) {
                    f(i,j,0) = f(i, j+1, 0) - bval * dx[1];
                } else {
                    f(i,j,0) = f(i, j+2, 0) - 2.0 * bval * dx[1];
                }
            }
        }
        
        // Right boundary (x-direction)
        if (i > domain.bigEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = phi[0];
            
            amrex::Real bval = bc_func[2](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i-1, j, 0);
            }
            else if (bc_type[2] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i-1, j, 0) + bval * dx[0];
            }
        }
        
        // Top boundary (y-direction)
        if (j > domain.bigEnd(1))
        {
            amrex::Real y = phi[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[3](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
            {
                if (j == domain.bigEnd(1)+1)
                {
                    f(i,j,0) = bval;
                } else {
                    f(i,j,0) = 2.0 * bval - f(i, j-2, 0);
                }
            }
            else if (bc_type[3] == amrex::LinOpBCType::Neumann)
            {
                if (j == domain.bigEnd(1)+1)
                {
                    
                } else {
                    f(i,j,0) = f(i, j-2, 0) + 2.0 * bval * dx[1];
                }
            }
        }
    }

    // Fill cell-centered - for pressure
    AMREX_GPU_DEVICE
    void fill_cc(const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
                 amrex::GeometryData const &geom, const amrex::Real time) const
    {
        const amrex::Box& domain = geom.Domain();
        const amrex::Real* dx = geom.CellSize();
        const amrex::Real* plo = geom.ProbLo();
        const amrex::Real* phi = geom.ProbHi();
        
        int i = iv[0];
        int j = iv[1];
        
        // Left boundary
        if (i < domain.smallEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = plo[0];
            
            amrex::Real bval = bc_func[0](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i+1, j, 0);
            }
            else if (bc_type[0] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i+1, j, 0) - bval * dx[0];
            }
        }
        
        // Bottom boundary
        if (j < domain.smallEnd(1))
        {
            amrex::Real y = plo[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[1](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i, j+1, 0);
            }
            else if (bc_type[1] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i, j+1, 0) - bval * dx[1];
            }
        }
        
        // Right boundary
        if (i > domain.bigEnd(0))
        {
            amrex::Real y = plo[1] + (j + 0.5) * dx[1];
            amrex::Real x = phi[0];
            
            amrex::Real bval = bc_func[2](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i-1, j, 0);
            }
            else if (bc_type[2] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i-1, j, 0) + bval * dx[0];
            }
        }
        
        // Top boundary
        if (j > domain.bigEnd(1))
        {
            amrex::Real y = phi[1];
            amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            
            amrex::Real bval = bc_func[3](AMREX_D_DECL(x, y, 0.0), time);
            
            if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
            {
                f(i,j,0) = 2.0 * bval - f(i, j-1, 0);
            }
            else if (bc_type[3] == amrex::LinOpBCType::Neumann)
            {
                f(i,j,0) = f(i, j-1, 0) + bval * dx[1];
            }
        }
    }
};

// Function declarations for filling boundary conditions
// These will be implemented in FillBC.cpp
void FillPhysBCs_pressure(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<PressureBCFunc, 4> &bc_funcs);  

void FillPhysBCs_xvel(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<MovingWallBCFunc, 4> &bc_funcs);

void FillPhysBCs_yvel(
    amrex::MultiFab &mf,
    const amrex::Geometry &geom,
    const amrex::Array<amrex::LinOpBCType, 4> &bctype,
    const amrex::Array<VelocityBCFunc, 4> &bc_funcs);

#endif