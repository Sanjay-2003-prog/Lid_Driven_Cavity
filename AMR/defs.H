#ifndef _defs_H_
#define _defs_H_

#include <AMReX_LO_BCTYPES.H>

typedef std::function<double(AMREX_D_DECL(double, double, double), double)> FieldFunc;
typedef amrex::Array<amrex::LinOpBCType, 2 * AMREX_SPACEDIM> BCArray;
typedef amrex::Array<FieldFunc, 2 * AMREX_SPACEDIM> BCValuesArray;

#endif
