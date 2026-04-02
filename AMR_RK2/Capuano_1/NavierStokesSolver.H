#ifndef NAVIERSTOKES_SOLVER_H_
#define NAVIERSTOKES_SOLVER_H_

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FillPatchUtil.H>
//#include <AMReX_FluxRegister.H>

#include "FillBC.H"
#include <limits>
#include <defs.H>

using namespace amrex;
using namespace std;

class NavierStokesSolver : public AmrCore
{
public:
    NavierStokesSolver();
    virtual ~NavierStokesSolver();

    void InitialData();     // Initializes all levels via InitFromScratch
    void Evolve();          // Main time-stepping loop
    void WritePlotFile();   // Output plotfile

    // AMR lifecycle overrides
    virtual void MakeNewLevelFromScratch(int lev, Real time,
                                         const BoxArray& ba,
                                         const DistributionMapping& dm) override;

    virtual void MakeNewLevelFromCoarse(int lev, Real time,
                                        const BoxArray& ba,
                                        const DistributionMapping& dm) override;

    virtual void RemakeLevel(int lev, Real time,
                             const BoxArray& ba,
                             const DistributionMapping& dm) override;

    virtual void ClearLevel(int lev) override;

    virtual void ErrorEst(int lev, TagBoxArray& tags,
                          Real time, int ngrow) override;

    // AMR data management

    // In NavierStokesSolver class declaration
    void FillPatchPressure( int lev, int icomp, int ncomp);

    void FillPatchVelocity( int lev, int icomp, int ncomp);

    void FillPatch( int lev, Real time, 
                    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> const& mf,                  // dest mf
                    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& cmf,  // coarse mf
                    const amrex::Vector<amrex::Real>& ct,                                      // coarse time
                    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& fmf,  // fine mf
                    const amrex::Vector<amrex::Real>& ft,                                      // fine time
                    int icomp, int ncomp,
                    amrex::GpuBndryFuncFab<struct_FillPhysBCs> phys_bc_fill_u,
                    amrex::GpuBndryFuncFab<struct_FillPhysBCs> phys_bc_fill_v);
    
    void FillPatch( int lev, Real time, 
                    MultiFab& mf,                        // dest mf
                    const Vector<MultiFab*>& cmf,         // coarse mf
                    const Vector<Real>& ct,               // coarse time
                    const Vector<MultiFab*>& fmf,         // fine mf
                    const Vector<Real>& ft,               // fine time
                    int icomp, int ncomp);  
    void FillCoarsePatch( int lev, Real time, 
                        MultiFab& mf,                        // dest mf
                        const MultiFab& cmf,                  // coarse mf
                        int icomp, int ncomp);
    
    void AverageDownVelocity();
    void AverageDownPressure();
    //void AverageDown();                     // Average all levels
    //void AverageDownTo(int crse_lev);       // Average fine to coarse level
    

    // Time stepping
    //void timeStepWithSubcycling(int lev, Real time, int iteration);
    // void timeStepNoSubcycling(Real time, int iteration);

    // User-defined initial conditions
    void SetUFunction(FieldFunc f) { u_func = f; }
    void SetVFunction(FieldFunc f) { v_func = f; }
    void SetPFunction(FieldFunc f) { p_func = f; }

private:
    // Initialization
    void readParameters();
    void memoryallocation(int lev, const BoxArray& ba, const DistributionMapping& dm);
    void InitializeField(int lev, Real time);

    // Time stepping internals
    void predictVelocity_k1();
    void predictVelocity_k2();
    void pressurePoisson();
    void correctVelocity();
    void copyNewToOld();
    void ComputeTimeStep();
    Real CheckSteadyState();
    bool continueRun();

    // Data access helpers
    void GetPressureData(int lev, Real time, Vector<MultiFab*>& mf, Vector<Real>& t, bool single);
    void GetFaceVelocityData(int lev, Real time, Vector<MultiFab*>& xmf, Vector<MultiFab*>& ymf, Vector<Real>& t, bool single);
    void GetVorticityData(int lev, Real time, Vector<MultiFab*>& mf, Vector<Real>& t, bool single);

    std::string getFieldName(amrex::MultiFab& mf, int lev);


    // Diagnostics
    void computeCellCenteredVel();
    void computeDivergence();
    void computeMagnitude();
    void computeVorticity(int lev);
    void Display(const MultiFab& mf, const std::string& name);

    // Simulation control
    int curr_iter = 0;
    int maxSteps = -1;
    int SaveEvery = 100;
    Real curr_time = 0.0;
    //Real dt = 1.0e-6;
    Real finalTime = 10.0;
    int Re = 100;
    Real Mu;
    // Vector<Real> t_new;  // Current time at each level
    // Vector<Real> t_old;  // Previous time at each level
    Real dt = 1.0e-6; // inital time step
    // Vector<Real> dt;     // Time step at each level
    // Vector<int>  istep;  // Number of steps taken at each level
    // Vector<int> nsubsteps;   // Number of substeps per level (used in subcycling)
    // bool do_subcycle = true; // Enable/disable subcycling
    // bool do_reflux = true;   // Enable/disable flux correction



    Real dt_adv = 0.0, dt_diff = 0.0;
    Real cfl_adv = 0.5, cfl_diff = 0.25;
    Real steadyTol = 1e-6;

    Real refine_threshold = 0.01;

    // Physical domain bounds
    Real prob_lo[AMREX_SPACEDIM];
    Real prob_hi[AMREX_SPACEDIM];

    // Field data
    Vector<MultiFab> xvel, xvel_old, xvel_k1, xvel_k2;
    Vector<MultiFab> yvel, yvel_old, yvel_k1, yvel_k2;
    Vector<MultiFab> U, U_magnitude; //U_old 
    Vector<MultiFab> Pressure; //Pressure_old;
    Vector<MultiFab> divU;
    Vector<MultiFab> vorticity; //vorticity_old;
    // Vector<MultiFab> flux_x;
    // Vector<MultiFab> flux_y;


    // Boundary conditions
    FieldFunc u_func, v_func, p_func;
    BCArray u_bc, v_bc, p_bc;
    BCValuesArray u_bc_func, v_bc_func, p_bc_func;

    // Flux register for AMR refluxing
    // Vector<std::unique_ptr<FluxRegister>> flux_register;

    // AMR settings (parsed from inputs)
    int max_level = 0;
    Vector<int> ref_ratio_vals;
    Vector<amrex::IntVect> ref_ratio; 

    int regrid_int = 0;
    Real old_time = 0.0;
    int n_error_buf = 0;
    int max_grid_size = 32;
    int blocking_factor = 8;
};

#endif
