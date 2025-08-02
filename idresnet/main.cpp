#include <iostream>
#include "rnet.hpp"
#include <chrono>
#include <boost/program_options.hpp>
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <sys/types.h>
#include <unistd.h>
namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv) {
    srand(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    int L = 100,mode,type,verbosity,pc_precision;
    int r = 10;
    int create_images=0;
    int profileOn=0;
    int dynamics=0;
    int write_dot=0;
    int stripAttack=0;
    int random_net=0;
    int error_check=0;
    int triangular=0;
    float q = 1;
    double p = 0.65;
    double p_step,kavg,lambda;
    int steps;
    long unsigned int time_stamp = time(0);
    string config_file;

        try {

        po::options_description desc("Supported options",1024);
        desc.add_options()
            ("help", "produce help message")
	    ("L,L",po::value<int>(&L)->default_value(100), "L")
            ("r,r",po::value<int>(&r)->default_value(10), "r")
	    ("p,p",po::value<double>(&p)->default_value(0.65),"p (initial)")
            ("dpc",po::value<int>(&pc_precision)->default_value(10), "delta pc - pc_precision.  pc can be found to within a single node-removal but it is very costly to get that close.  dragons live in dpc<2")
	    ("kavg,k",po::value<double>(&kavg)->default_value(2.54),"kavg (initial)")
	    ("lambda,l",po::value<double>(&lambda)->default_value(-1),"lambda for connectivity distribution r~ lambda*exp(-lambda*r) if lambda<1 it reverts to the lattices with r")
	    ("q,q",po::value<float>(&q)->default_value(1),"q")
	    ("type,t",po::value<int>(&type)->default_value(1),"dependency type: no dependency(0), no voltage(1), geometric dependency(2), dynamic dependency(3)")
	    ("mode,m",po::value<int>(&mode)->default_value(1),"simulation mode: find_pc then sample above it (1), create matrices to test solvers on(2), sample fixed p values only(3)")
	    //("find,f",po::value<int>(&find_flag)->default_value(1),"find pc: if set to zero, the p variable will be taken as pc with no attempt to find the actual pc (mode 3 only)")
	    ("steps,s",po::value<int>(&steps)->default_value(100),"number of steps for prange (spacing is set by dp) (mode 3 only)")
	    ("dp,dp",po::value<double>(&p_step)->default_value(0.001), "p_step")
	    //("zero,z",po::value<double>(&zero_size)->default_value(0.05),"zero size (the size at which the network is considered equal to zero)")
	    ("triangular","use triangular lattice instead of square lattice (affects lambda networks only)")
	    ("history,h","record_lattice_history")
            ("dynamics","output cascade dynamics for all p (in json)")
            ("profile","output time profiling info to file and terminal if verbosity>0")            
            ("strip","remove strip of width r")
            ("verbosity,v",po::value<int>(&verbosity)->default_value(0),"verbose output 0:silent,1:profiling info")
            ("dot","write graph_viz file Graph[0|1].dot")
            ("random","use random networks -- no lambda, no r, for now ER networks with kavg only")
            ("error","Run error checking code.  This should never be run while trying to obtain results as it adds lots of calculations that should not be necessary if everything is working.")
	    ("timestamp",po::value<long unsigned int>(&time_stamp),"time_stamp (long int used to differentiate the files, if not supplied will taken from current clock)")
	    ("config",po::value<string>(&config_file)->default_value("interdep.cfg"),"path to config file");

        po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
	ifstream ifs(config_file.c_str());
        if (ifs){
	 cout<<"Config file read ("<<config_file<<").\n";
	 po::store(po::parse_config_file(ifs,desc),vm);
	 po::notify(vm);
	}

        if ((vm.count("help"))){
            cout << desc << "\n";
            return 1;
        }
        if (vm.count("history")) 
 	    create_images=1;
        if (vm.count("profile"))
            profileOn=1;
        if (vm.count("dynamics"))
            dynamics=1;
        if (vm.count("strip"))
            stripAttack=1;
        if (vm.count("dot"))
            write_dot=1;
        if (vm.count("error"))
            error_check=1;
	if (vm.count("triangular"))
	    triangular=1;
        if (vm.count("random") || (lambda<0 && r<0))
            random_net=1;
	if (!vm.count("timestamp"))
	    time_stamp = time(0);

    }catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
	return 1;
    }
    int N = L*L;
    RNet rn;
    /** Set various options**/
    if (stripAttack)
        rn.set_strip(1);
    if (profileOn)
        rn.turnProfileOn();

    if (error_check)
	rn.enableErrorCheck();

    rn.setDependencyType(type);

    if(random_net){
        rn.Initialize(L,kavg,q,time_stamp);
    }else{
        rn.Initialize(L,r,q,kavg,lambda,triangular,time_stamp);
    }
    if (write_dot)
        rn.write_graphviz();
    rn.setPcPrecision(pc_precision);
    rn.setVerbosity(verbosity);
    //cout << "RNet initialized" << endl;

    FloatVec pVec;
    double pc=0.5927;
    switch(mode){
	case 1:{
            if(dynamics)
                rn.enableWriteDynamics();
	    pc = rn.find_pc(p,0.01,create_images);
            if(pc<1){
                if(steps<0){
                        steps =  ((1 - pc) / p_step);
                }
                    for(int i=0; i<steps; ++i){
                        double p_ = pc+(steps-i)*p_step;
                        if (p_< 1)
                            pVec.push_back(p_);
                        else{
                            pVec.push_back(1);
                            break;
                        }
                    }
                
                rn.tracked_cascade_multiple_p(pVec,0);
            }
	    rn.write_summary(); 
            if (dynamics)
                rn.closeDynamicsFile();
	    break;
	}
	case 2:
	    rn.setDependencyType(0);
	    rn.enableWriteLaplacian();
	    rn.tracked_cascade_single_p(p,0);
	    break;
    case 3:
        for (int i=0; i<steps; ++i)
        {
            double p_ = p+(steps-i)*p_step;

            if (p_ <1)
                pVec.push_back(p_);
        }

        if (dynamics)
            rn.enableWriteDynamics();
        rn.tracked_cascade_multiple_p(pVec,create_images);
        rn.write_summary();
        if (dynamics)
            rn.closeDynamicsFile();
        break;

	
    }
  //  rn.reset();
    return 0;

}
