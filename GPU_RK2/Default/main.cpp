#include <NavierStokesSolver.H>
// #include <AMReX_ParmParse.H>
// #include <AMReX_BLProfiler.H>
using namespace amrex;

int main(int argc, char** argv)
{
    Initialize(argc, argv);
    {
        auto strt_time = ParallelDescriptor::second();

        #ifdef AMREX_TINY_PROFILING
            amrex::Print() << "Tiny Profiling is ENABLED" << std::endl;
            // amrex::Print() << "Profile output will be written to: ";
            
            // // You can optionally read the output filename to confirm
            // ParmParse pp("tiny_profiler");
            // std::string output_file;
            // if (pp.query("output_file", output_file)) {
            //     amrex::Print() << output_file << std::endl;
            // } else {
            //     amrex::Print() << "default location" << std::endl;
            // }
        #else
            amrex::Print() << "Tiny Profiling is DISABLED" << std::endl;
        #endif

        //Print()<<"1"<<std::endl;
        //Print()<<"2"<<std::endl;
        NavierStokesSolver solver;
        //Print()<<"1"<<std::endl;
        //Print()<<"2"<<std::endl;
        solver.InitialData();
        //Print()<<"1"<<std::endl;
        solver.Evolve();
        solver.cleanup();
        auto total_time = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(total_time,IOProc);

        // Tell the I/O Processor to write out the "run time"
        amrex::PrintToFile("log") << "Run time = " << total_time << " secs " << std::endl;

        //BL_PROFILE_SET_FILENAME("GPU_Timing.txt");
        //amrex::Print() << "Run time = " << total_time << std::endl;

        // #ifdef AMREX_USE_GPU
        //     // Force GPU cleanup before AMReX finalize
        //     amrex::Gpu::Device::reset();
        // #endif
    }
    Finalize();
    return 0;
}
