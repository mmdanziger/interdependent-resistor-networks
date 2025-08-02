#include <iostream>
#include "rnet.hpp"
#include <chrono>
#include <boost/program_options.hpp>
#include <omp.h>
#include <sstream>
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv) {
    srand(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    int L = 500,mode,verbosity,pc_precision;
    int r = -1;
    int type=1;
    int create_images=0;
    int profileOn=0;
    int dynamics=0;
    int write_dot=0;
    int stripAttack=0;
    int random_net=0;
    int triangular=0;
    int error_check=0;
    int num_samps=30;
    int thread_count=0;
    double min_exp=-1;
    float q = 1;
    double p = 0.65;
    double p_step,kavg=4,lambda;
    int steps;
    long unsigned int time_stamp = time(0);
    string config_file;

        try {

        po::options_description desc("Supported options",1024);
        desc.add_options()
        ("help", "produce help message")
        ("L,L",po::value<int>(&L)->default_value(100), "L")
        ("kavg,k",po::value<double>(&kavg)->default_value(4),"kavg (initial)")
        ("lambda,l",po::value<double>(&lambda)->default_value(-1),"lambda for connectivity distribution r~ lambda*exp(-lambda*r) if lambda<1 it reverts to the lattices with r")
        ("mode,m",po::value<int>(&mode)->default_value(1),"simulation mode: find_pc then sample above it (1), create matrices to test solvers on(2), sample fixed p values only(3)")
        ("samples,s",po::value<int>(&num_samps)->default_value(30),"number of samples")
        ("threads,t",po::value<int>(&thread_count)->default_value(1),"number of threads")
        ("zetamin,z",po::value<double>(&min_exp)->default_value(-1),"minimum power to sample (ie, zeta_min = 10^zetamin")
	("triangular","use triangular lattice instead of square lattice (affects lambda networks only)")
	("dynamics","output cascade dynamics for all p (in json)")
        ("profile","output time profiling info to file and terminal if verbosity>0")            
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
/*        if (vm.count("profile"))
            profileOn=1;
        if (vm.count("dynamics"))
            dynamics=1;
        if (vm.count("strip"))
            stripAttack=1;

        if (vm.count("error"))
            error_check=1;
        */
	if (vm.count("triangular"))
	  triangular=1;
        if (vm.count("dot"))
            write_dot=1;
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

    omp_set_num_threads(thread_count);
    //cout << "RNet initialized" << endl;
    ofstream cc_ofile,radius_ofile,diameter_ofile,all_ofile;
    std::stringstream cc_fname,radius_fname,diameter_fname,all_fname;
    if(triangular==0)
      all_fname << "ClusteringRadiusDiameter_L" << L << "_" <<time(0) <<".json";
    else
      all_fname << "ClusteringRadiusDiameterTriangular_L" << L << "_" <<time(0) <<".json";
   // cc_fname << "Clustering_L" << L << "_" <<time(0) <<".json";
   // radius_fname << "Radius_L" << L << "_" <<time(0) <<".json";
   // diameter_fname << "Diameter_L" << L << "_" <<time(0) <<".json";
    
    
    all_ofile.open(all_fname.str().c_str()); //rvalue black magic, i think
   // length_ofile.open("Lengths.json");
    
    vector<double> zeta_vec;
    //int num_samps=30;
    double max_exp = log10(L/2);
    double step = (max_exp - min_exp) / (num_samps +1);
    for(int i=0; i<num_samps; i++)
        zeta_vec.push_back(pow(10.0, min_exp + i * step));
    //map<double,long> overlaps;
    map<double, double> clustering_coefficients;
    map<double, long> radiuses;
    map<double, long> diameters;
//     RNet dummyr;
//     dummyr.Initialize(L,r,q,kavg,lambda,triangular,time_stamp);
//     for(int samples=5; samples<N; samples*=2){
// 	auto rdpair = dummyr.getRadiusDiameter(samples);
// 	std::cout << samples << "\t" << rdpair.first << "\t" <<rdpair.second << std::endl;
// 	std::cout.flush();
//     }
//     
//     return 0;
    
    //map<double,vector<double>> lengths;
//#pragma omp parallel for
    std::cout << "zeta\t\tcc\tr\td\n";
    for(unsigned i=0; i< zeta_vec.size(); ++i){
        double zeta = zeta_vec[i];
        RNet rn;    
	if(write_dot)
	  rn.setDependencyType(0);
	else
	  rn.setDependencyType(1);
        rn.setVerbosity(verbosity);
        rn.Initialize(L,r,q,kavg,1 / zeta,triangular,time_stamp);
	std::cout << zeta << "\t";
        //auto thisoverlap = rn.getOverlap();
        //overlaps[zeta] = thisoverlap;
        auto thiscc = rn.getClusteringCoefficient();
	if(write_dot){
	  rn.write_graphviz();
	  return 0;
	}
	clustering_coefficients[zeta] = thiscc;
	std::cout << thiscc << "\t";
	auto thisrd = rn.getRadiusDiameter(min(350,N));
	radiuses[zeta] = thisrd.first;
	diameters[zeta] = thisrd.second;
	std::cout << thisrd.first << "\t";
	std::cout << thisrd.second << "\n";
	//lengths[zeta] = rn.getLinkLengths(0);
    }
    
    all_ofile << " { \"clustering\" : " ;
    jsonMap(clustering_coefficients,all_ofile);
    all_ofile << " ,\n \"radius\" : " ;
    jsonMap(radiuses,all_ofile);
    all_ofile << " ,\n \"diameter\" : " ;
    jsonMap(diameters,all_ofile);
    all_ofile << " }\n " ;
    all_ofile.close();
    
    //jsonMapOfArrays(lengths,length_ofile);
}