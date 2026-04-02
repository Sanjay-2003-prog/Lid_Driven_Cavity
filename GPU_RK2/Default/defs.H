#ifndef _defs_H_
#define _defs_H_

#include <AMReX_LO_BCTYPES.H>

#ifndef AMREX_USE_GPU
#include <functional>
#endif

#ifdef AMREX_USE_GPU
typedef double (*FieldFunc)(AMREX_D_DECL(double, double, double), double);
#else
typedef std::function<double(AMREX_D_DECL(double, double, double), double)> FieldFunc;
#endif

typedef amrex::Array<amrex::LinOpBCType, 2 * AMREX_SPACEDIM> BCArray;
typedef amrex::Array<FieldFunc, 2 * AMREX_SPACEDIM> BCValuesArray;

#endif