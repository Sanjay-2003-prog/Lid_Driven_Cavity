#include "FillBC.H"

#include <AMReX_PhysBCFunct.H>

AMREX_GPU_DEVICE
void struct_FillPhysBCs::operator() 
(
    const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
    const int dcomp, const int numcomp,
    amrex::GeometryData const& geom, const amrex::Real time,
    const amrex::BCRec* bcr, const int bcomp,
    const int orig_comp
) const
{
    if (nodality == amrex::IntVect(0,0))
    {
        fill_cc(iv, dest, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp);
    }
    if (nodality == amrex::IntVect(1,0))
    {
        fill_xnd(iv, dest, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp);
    }
    if (nodality == amrex::IntVect(0,1))
    {
        fill_ynd(iv, dest, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp);
    }
}


/************************************************************************************/

AMREX_GPU_DEVICE
void struct_FillPhysBCs::fill_xnd 
(
    const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& f,
    const int dcomp, const int numcomp,
    amrex::GeometryData const& geom, const amrex::Real time,
    const amrex::BCRec* bcr, const int bcomp,
    const int orig_comp
) const
{
    amrex::Box domain = geom.Domain();
    const amrex::Real *dx = geom.CellSize();
    const amrex::Real *plo = geom.ProbLo();
    const amrex::Real *phi = geom.ProbHi();

    amrex::Dim3 dlo = lbound(domain);
    amrex::Dim3 dhi = ubound(domain);

    // left boundary
    if (iv[0] <= dlo.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = plo[0];

        amrex::Real bval = bc_func[0](x, y, time);

        if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
        {
            if (iv[0] == dlo.x)
            {
                f(iv[0], iv[1], 0) =  bval ;
            }
            else
            {
                f(iv[0], iv[1], 0) = 2.0* bval - f(iv[0] + 2, iv[1], 0);                
            }
        }
        else if (bc_type[0] == amrex::LinOpBCType::Neumann)
        {
            if (iv[0] == dlo.x)
            {

            }
            else
            {
                f(iv[0], iv[1], 0) = f(iv[0] + 2, iv[1], 0) -  2.0* bval * dx[0];
            }
        }
    }

    // bottom
    if (iv[1] < dlo.y)
    {
        amrex::Real y = plo[1];
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];

        amrex::Real bval = bc_func[1](x, y, time);
        if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0] , iv[1]+1, 0);
        }
        else if (bc_type[1] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0], iv[1] + 1, 0) -  bval * dx[1];
        }
    }
    // right
    if (iv[0] > dhi.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = phi[0];

        amrex::Real bval = bc_func[2](x, y, time);

        if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
        {
            if (iv[0] == dhi.x+1)
            {
                f(iv[0], iv[1], 0) = bval;
            }
            else
            {
                f(iv[0], iv[1], 0) = 2.0* bval - f(iv[0] - 2, iv[1], 0);
                
            }
        }
        else if (bc_type[2] == amrex::LinOpBCType::Neumann)
        {
            if (iv[0] == dhi.x+1)
            {
            }
            else
            {
                f(iv[0], iv[1], 0) = f(iv[0] - 2, iv[1], 0) + 2.0* bval * dx[0];
            }
        }
    }

    // top
    if (iv[1] > dhi.y)
    {
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];
        amrex::Real y = plo[1];

        amrex::Real bval = bc_func[3](x, y, time);
        if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0], iv[1] - 1, 0);
        }
        else if (bc_type[3] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0], iv[1] - 1, 0) + bval * dx[1];
        }
    }
}

AMREX_GPU_DEVICE
void struct_FillPhysBCs::fill_ynd(
    const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
    const int dcomp, const int numcomp,
    amrex::GeometryData const &geom, const amrex::Real time,
    const amrex::BCRec *bcr, const int bcomp,
    const int orig_comp) const
{
    amrex::Box domain = geom.Domain();
    const amrex::Real *dx = geom.CellSize();
    const amrex::Real *plo = geom.ProbLo();
    const amrex::Real *phi = geom.ProbHi();

    amrex::Dim3 dlo = lbound(domain);
    amrex::Dim3 dhi = ubound(domain);

    // left boundary
    if (iv[0] < dlo.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = plo[0];

        amrex::Real bval = bc_func[0](x, y, time);

        if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0] + 1, iv[1], 0);
        }
        else if (bc_type[0] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0] + 1, iv[1], 0) - bval * dx[0];
        }
    }

    // bottom
    if (iv[1] <= dlo.y)
    {
        amrex::Real y = plo[1];
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];

        amrex::Real bval = bc_func[1](x, y, time);
        if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
        {
            if (iv[1] == dlo.y)
            {
                f(iv[0], iv[1], 0) = bval;
            }
            else
            {
                f(iv[0], iv[1], 0) = 2 * bval - f(iv[0], iv[1] + 2, 0);
            }
        }
        else if (bc_type[1] == amrex::LinOpBCType::Neumann)
        {
            if (iv[1] == dlo.y)
            {

            }
            else
            {
                f(iv[0], iv[1], 0) = f(iv[0], iv[1] + 2, 0) - 2.0 * bval * dx[1];
            }
        }
    }
    // right
    if (iv[0] > dhi.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = phi[0];

        amrex::Real bval = bc_func[2](x, y, time);

        if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0] - 1, iv[1], 0);
        }
        else if (bc_type[2] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0] - 1, iv[1], 0) + bval * dx[0];
        }
    }

    // top
    if (iv[1] > dhi.y)
    {
        amrex::Real y = phi[1];
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];

        amrex::Real bval = bc_func[3](x, y, time);
        if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
        {
            if (iv[1] == dhi.y + 1)
            {
                f(iv[0], iv[1], 0) = bval;
            }
            else
            {
                f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0], iv[1] - 2, 0);
            }
        }
        else if (bc_type[3] == amrex::LinOpBCType::Neumann)
        {
            if (iv[1] == dhi.y + 1)
            {

            }
            else
            {
                f(iv[0], iv[1], 0) = f(iv[0], iv[1] - 2, 0) + 2.0 * bval * dx[1];
            }
        }
    }
}

AMREX_GPU_DEVICE
void struct_FillPhysBCs::fill_cc(
    const amrex::IntVect &iv, amrex::Array4<amrex::Real> const &f,
    const int dcomp, const int numcomp,
    amrex::GeometryData const &geom, const amrex::Real time,
    const amrex::BCRec *bcr, const int bcomp,
    const int orig_comp) const
{
    amrex::Box domain = geom.Domain();
    const amrex::Real *dx = geom.CellSize();
    const amrex::Real *plo = geom.ProbLo();
    const amrex::Real *phi = geom.ProbHi();

    amrex::Dim3 dlo = lbound(domain);
    amrex::Dim3 dhi = ubound(domain);

    // left boundary
    if (iv[0] < dlo.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = plo[0];

        amrex::Real bval = bc_func[0](x, y, time);

        if (bc_type[0] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0] + 1, iv[1], 0);
        }
        else if (bc_type[0] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0] + 1, iv[1], 0) - bval * dx[0];
        }
    }

    // bottom
    if (iv[1] < dlo.y)
    {
        amrex::Real y = plo[1];
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];

        amrex::Real bval = bc_func[1](x, y, time);
        if (bc_type[1] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0], iv[1] + 1, 0);
        }
        else if (bc_type[1] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0], iv[1] + 1, 0) - bval * dx[1];
        }
    }

    // right
    if (iv[0] > dhi.x)
    {
        amrex::Real y = plo[1] + (iv[1] + 0.5) * dx[1];
        amrex::Real x = phi[0];

        amrex::Real bval = bc_func[2](x, y, time);

        if (bc_type[2] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0] - 1, iv[1], 0);
        }
        else if (bc_type[2] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0] - 1, iv[1], 0) + bval * dx[0];
        }
    }

    // top
    if (iv[1] > dhi.y)
    {
        amrex::Real y = phi[1];
        amrex::Real x = plo[0] + (iv[0] + 0.5) * dx[0];

        amrex::Real bval = bc_func[3](x, y, time);
        if (bc_type[3] == amrex::LinOpBCType::Dirichlet)
        {
            f(iv[0], iv[1], 0) = 2.0 * bval - f(iv[0], iv[1] - 1, 0);
        }
        else if (bc_type[3] == amrex::LinOpBCType::Neumann)
        {
            f(iv[0], iv[1], 0) = f(iv[0], iv[1] - 1, 0) + bval * dx[1];
        }
    }
}


void FillPhysBCs_pressure
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
)
{
    mf.FillBoundary();
    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;

    bcs1.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir); // these are only dummy, not used to fill actual bcs
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    amrex::GpuBndryFuncFab<struct_FillPhysBCs> bf(struct_FillPhysBCs{bctype, bc_values, amrex::IntVect(0,0)});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs> > physbcf(geom, bcs1, bf);
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}

void FillPhysBCs_xvel
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
)
{
    mf.FillBoundary();

    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;

    bcs1.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir); // these are only dummy, not used to fill actual bcs
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    amrex::GpuBndryFuncFab<struct_FillPhysBCs> bf(struct_FillPhysBCs{bctype, bc_values, amrex::IntVect(1,0)});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs> > physbcf(geom, bcs1, bf);
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}

void FillPhysBCs_yvel
(
    amrex::MultiFab& mf, 
    const amrex::Geometry& geom,
    const amrex::Array<amrex::LinOpBCType, 2*AMREX_SPACEDIM>& bctype,
    const amrex::Array<FieldFunc, 2*AMREX_SPACEDIM>& bc_values 
)
{
    mf.FillBoundary();

    amrex::Real time = 0.0;

    amrex::Vector<amrex::BCRec> bcs1;

    bcs1.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs1[0].setLo(idim, amrex::BCType::ext_dir); // these are only dummy, not used to fill actual bcs
        bcs1[0].setHi(idim, amrex::BCType::ext_dir);
    }

    amrex::GpuBndryFuncFab<struct_FillPhysBCs> bf(struct_FillPhysBCs{bctype, bc_values, amrex::IntVect(0,1)});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<struct_FillPhysBCs> > physbcf(geom, bcs1, bf);
    physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), time, 0);
}


