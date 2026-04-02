#ifndef BCFUNC_H
#define BCFUNC_H

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>



// Functor for velocity boundary conditions (no-slip walls)
struct VelocityBCFunc {
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    amrex::Real operator()(AMREX_D_DECL(double x, double y, double z), double time) const
    {
        // No-slip walls - return 0 for all boundaries
        return 0.0;
    }
};

// Functor for pressure boundary conditions (Neumann)
struct PressureBCFunc {
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    amrex::Real operator()(AMREX_D_DECL(double x, double y, double z), double time) const
    {
        // Neumann BC for pressure - value not used directly
        return 0.0;
    }
};

// Functor for moving wall (lid-driven cavity top wall)
struct MovingWallBCFunc {
    // Data member (state)
    amrex::Real m_wall_velocity;
    
    // Default constructor (required by AMReX)
    AMREX_GPU_HOST_DEVICE
    MovingWallBCFunc() : m_wall_velocity(0.0) {}
    
    // Constructor with parameter
    AMREX_GPU_HOST_DEVICE
    MovingWallBCFunc(double vel) : m_wall_velocity(vel) {}
    
    // Call operator uses the state
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    amrex::Real operator()(AMREX_D_DECL(double x, double y, double z), double time) const
    {
        // Return the stored velocity
        return m_wall_velocity;
    }
};

// struct MovingWallBCFunc {
//     double m_wall_velocity;
    
//     // Default constructor (needed for arrays)
//     AMREX_GPU_HOST_DEVICE
//     MovingWallBCFunc() : m_wall_velocity(0.0) {}
    
//     // Constructor with parameter
//     AMREX_GPU_HOST_DEVICE
//     MovingWallBCFunc(double vel) : m_wall_velocity(vel) {}
    
//     AMREX_GPU_HOST_DEVICE
//     AMREX_INLINE
//     double operator()(AMREX_D_DECL(double x, double y, double z), double time) const
//     {
//         return m_wall_velocity;
//     }
// };

// For backward compatibility, we can still define a type alias
// But this won't be used in the new implementation
typedef double (*LegacyBCFuncPtr)(AMREX_D_DECL(double, double, double), double);

#endif 