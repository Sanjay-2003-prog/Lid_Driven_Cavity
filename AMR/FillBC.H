#ifndef _fillbc_H_
#define _fillbc_H_

#include <AMReX_MultiFab.H>
#include <AMReX_LO_BCTYPES.H>   
#include <defs.H>

#include <AMReX_BCRec.H>


struct struct_FillPhysBCs
{
    struct_FillPhysBCs
    (
        const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
        const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values,
        const amrex::IntVect& nd 
    )
    :
    bc_type(bctype),
    bc_func(bc_values),
    nodality(nd)
    {}

    BCArray bc_type;
    BCValuesArray bc_func;
    amrex::IntVect nodality;

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const;


    AMREX_GPU_DEVICE
    void fill_xnd(const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const;                         

    AMREX_GPU_DEVICE
    void fill_ynd(const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const;   

    AMREX_GPU_DEVICE
    void fill_cc(const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const;  

};

void FillPhysBCs_pressure
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
);

void FillPhysBCs_xvel
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
);

void FillPhysBCs_yvel
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
);

#endif
