#include <NavierStokesSolver.H>
using namespace amrex;

int main(int argc, char** argv)
{
    Initialize(argc, argv);
    auto strt_time = ParallelDescriptor::second();
    //Print()<<"1"<<std::endl;
    //Print()<<"2"<<std::endl;
    NavierStokesSolver solver;
    //Print()<<"1"<<std::endl;
    //Print()<<"2"<<std::endl;
    solver.InitialData();
    //Print()<<"1"<<std::endl;
    solver.Evolve();
    auto total_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(total_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::PrintToFile("log") << "Run time = " << total_time << " secs " << std::endl;
    //amrex::Print() << "Run time = " << total_time << std::endl;
    Finalize();
    return 0;
}
