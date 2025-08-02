//#define NDEBUG 
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <algorithm>
//#include <chrono>
#include <assert.h>
#include <cmath>
#include <iomanip>
#include <boost/bimap.hpp>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/array.hpp>
#include <boost/format.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/exterior_property.hpp>
//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/grid_graph.hpp>
//#include <boost/graph/breadth_first_search.hpp>
//#include <boost/graph/filtered_graph.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/property_map/property_map.hpp>
//#include "custom_components.hpp"
#include "utilities.hpp"
#include "jsonutils.hpp"
#include "trianglevector2d.hpp"
//#include "breadth_first_search.hpp"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>
//#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/LU>
// #define INIT_TIME 0
// #define BFS_TIME 1
// #define LAP_CREATE_TIME 2
// #define LAP_SOLVE_TIME 3
// #define VBFS_TIME 4
// #define PRE_BINSEARCH_TIME 5
// #define BINSEARCH_TIME 6
// #define TIME_SECTIONS 7
static const int MAX_ITERATIONS = 5000;
static const int INIT_TIME = 0;
static const int BFS_TIME =1;
static const int LAP_CREATE_TIME =2;
static const int LAP_SOLVE_TIME =3;
static const int VBFS_TIME =4;
static const int PRE_BINSEARCH_TIME =5;
static const int BINSEARCH_TIME =6;
static const int TIME_SECTIONS =7;
static const float MIN_ALLOWED_P = 0.f;
static const float MAX_ALLOWED_P = 1.f;
std::vector<float> profilingVector(TIME_SECTIONS);
boost::array<string, TIME_SECTIONS > profilingNames = { {"init", "bfs", "laplacian_create", "laplacian_solve", "voltage_bfs", "prebinsearch", "binsearch" } };

static const int MINIMUM_PC_PRECISION = 2;

typedef boost::property<vertex_alive_t, boost::array<float, 2>,  boost::property<boost::vertex_color_t, boost::default_color_type> >  VertexPropertiesWithAlive;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexPropertiesWithAlive> GraphWithAlive;
typedef boost::property<boost::vertex_color_t, boost::default_color_type>   VertexProperties;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperties> Graph;
typedef boost::property_map<Graph, boost::vertex_color_t>::type ColorMap;
typedef boost::property_map<GraphWithAlive, vertex_alive_t>::type AliveMap;
//typedef boost::property_map< ::Graph, boost::vertex_index_t>::type IndexMap;
typedef boost::queue<int> queue_t;
// typedef components_recorder< vector<int> > comp_recorder_t;
typedef boost::property_traits<ColorMap>::value_type ColorValue;
typedef boost::color_traits<ColorValue> Color;
typedef std::pair<int,std::pair<int,int> > DistancePair;
typedef boost::bimap< int, int> IntBimap;
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> EigenSparse;
using namespace std;
boost::array<bool, 2> wrapped = { { true, true } };

   template <class Graph, class ComponentMap>
  inline typename boost::property_traits<ComponentMap>::value_type
 custom_connected_components(const Graph& g, ComponentMap c, ColorMap cmap)
  {
    if (boost::num_vertices(g) == 0) return 0;

    int c_count=0;
    components_recorder<ComponentMap> vis(c, c_count);
    for(unsigned int v =0; v<boost::num_vertices(g); ++v){
	if (cmap[v] == Color::white()){
	    ++c_count;
	    boost::breadth_first_visit(g,v,color_map(cmap).visitor(vis));
	}
    } 
    return c_count + 1;
  }


struct run_summary {
  run_summary(){}
  run_summary(int dep_type, double p, double xi, int hole, int Ng, int iter) : m_p(p), m_xi(xi), m_Ng(Ng),m_iter(iter),m_hole(hole),m_type(dep_type) {}
  run_summary(int dep_type, double p, int Ng, int iter) : m_p(p), m_Ng(Ng), m_iter(iter), m_type(dep_type) { m_hole=0; }
  run_summary(int dep_type, double p, int Ng, int Nb, double S) : m_p(p),m_S(S),m_Ng(Ng), m_Nb(Nb) ,m_type(dep_type){ 
      if (m_type>0){
          m_iter=m_Nb;
          m_Nb=0;}
  }
  string getHeader(){
      switch(m_type){
          case 0:
              return string("p\tNg\tNb\tS");
              break;
          case 1:
              return string("p\tNg\tNOI");
              break;
          case 2:
          case 3:
              return string("p\tNg\tNOI\tS");
              break;
      }
      return string("");
  }
  
  double m_p=0,m_xi=0,m_S=0;
  int m_Ng=0,m_iter=0,m_hole=0,m_Nb=0,m_type=0;
  
};

ostream& operator<< (ostream& stream, run_summary& rs){
  switch (rs.m_type){
    case -1://temporarily(?) deprecated
      stream << setprecision (10) << rs.m_p << "\t" << rs.m_Ng << "\t" << rs.m_iter << "\t" << rs.m_xi << "\t" << rs.m_hole;
      break;
    case 0:
      stream << setprecision (10) << rs.m_p << "\t" << rs.m_Ng << "\t" << rs.m_Nb <<  "\t" << rs.m_S ;
      break;
    case 1:
      stream << setprecision (10) << rs.m_p << "\t" << rs.m_Ng << "\t" << rs.m_iter ;
      break;
    case 2:
    case 3:
      stream << setprecision (10) << rs.m_p << "\t" << rs.m_Ng << "\t" << rs.m_iter <<  "\t" << rs.m_S ;
      break;
   }
  return stream;
}

class RNet
{
  protected:
    string dummy;
    int L,N;
    int no_neighbor_count,stop,iteration_number,last_cascade_vertex,firstP,firstIter,pc_precision,r;
    std::array<int,2> Ng;
    float q;
    int verbosity=0;
    int isRandom=0;
    int isGrid=0;
    int isLambda=0;
    int strip=0;
    int double_attack_flag=1;
    int writeLaplacianOn=0;
    int profileOn=0;
    int writeDynamics=0;
    int errorCheck=0;
    int depType=0;
    int triangular = 0;
    boost::random::exponential_distribution<double> exp;
    double kavg,lambda;
    double real_pc=0;
    long unsigned int time_stamp;
    std::array<Graph,2> G;
    IntBimap interlinks;
    boost::array<vector<int>,2>cascade_order;
    vector<int> component_vec;
    chrono::high_resolution_clock HRClock;
    ColorMap color;
    vector<bool> multiplex_dependency_vec;
    //AliveMap alive;    
    std::vector< std::array<float,2> > alive;
    map<double,run_summary> summary_map;
    void add_triplet(int node1, int node2, int value);
    virtual void kill_all();
    std::vector<T> tripletList;
    std::vector<DistancePair> pencils;
    std::vector<TriDistancePair> tripencils;


    boost::array<vector<int> ,2> node2matrix;

    std::vector< int > v0;
    boost::array<std::vector< double >, 2> voltage;
    double eps = 1e-9;
    float maxRelativeError = 1e-3f;
    double total_current=0;
    int sources,sinks;
    int V_fixed=100;//arbitrary number to make the numeric zeros easier to spot
    ofstream laplacianFile;
    ofstream V0File;
    ofstream dynamicsFile;
    bool periodic_bc;
    
/** METHODS **/    
  public:
    RNet();
    void Initialize(int L_, int r_, float q_,float kavg_, float lambda_, int triangular_,long unsigned int time_stamp_);
    void Initialize( int L_, float kavg_, float q_, long unsigned int time_stamp_);
    void create_interlinks();
    void make_connectivity_links(int lattice_index);
    void create_lattices();
    void restore_lattices();
    void calculate_distances();
    int draw_triangle_target(int sNode);
    void calculate_triangle_distances();
    void start_cascade(double p);
    virtual int trim_to_gc();
    virtual int trim_to_st_clusters();
    int get_Ng();
    vector<vector<int>> get_adjacency_list(int lattice_index);
    int kill_dependent();
    vector<int> get_giant_hole();
    void continue_cascade(double p);
    void start_tracked_cascade(double p);
    void start_tracked_dual_cascade(double pA, double pB);
    void remove_strip();
    unordered_set<int> boost_giant(bool voltage_calc = false);
    void write_history(int lattice_index);
    void set_strip(int s);
    void set_double_attack(int d);
/** Voltage calc methods **/
    bool hasSpanningCluster();
    int findVoltage();
    void enableWriteLaplacian(int number_to_write = 4);
    void writeLaplacian();
    virtual void enableWriteDynamics();
    void enableErrorCheck();
    void disableErrorCheck();
    void closeDynamicsFile();
    void solveLaplacian();
    int trim_to_backbone();
    double findTotalCurrent();
    int getFixedVoltage();
    int draw_k();
    int draw_target(int sNode);
    void setPcPrecision(int pcPrecision);
    void setVerbosity(int verbosity_);
    void setDependencyType(int dependency_type);
    void turnProfileOn();
/** Graph calculations**/    
    map<int,int> getDegreeDistribution(int lattice_index);
    vector<double> getLinkLengths(int lattice_index);
    long getOverlap();
    double getClusteringCoefficient();
    pair<long,long> getRadiusDiameter(int number_of_samples);
/** Methods and objects needed for custom BFS **/
    template <class ComponentMap>
//     inline typename boost::property_traits<ComponentMap>::value_type 
    int voltage_connected_components(ComponentMap c);
    template <class ComponentMap>
    void voltage_breadth_first_visit(queue_t & Q, int v,  components_recorder<ComponentMap> vis);
/** 'main' type methods **/
    virtual void tracked_cascade_single_p(double p,  int create_images=0);
    void continued_cascade(double p0, double p_step, double zero_size, int create_images);
    void tracked_cascade_multiple_p(FloatVec p, int create_images);
    double find_pc(double p0, double step, int create_images);
    void clear_summary();
    void write_summary();
    void write_graphviz();
    void writeVoltage();

};
RNet::RNet() { }

void RNet::Initialize(int L_, int r_, float q_,float kavg_, float lambda_, int triangular_, long unsigned int time_stamp_) 
{
    L = L_;
    r = r_;
    N = L*L;
    q = q_;
    kavg = kavg_;
    lambda = lambda_;
    triangular = triangular_;
    time_stamp = time_stamp_;

    
    auto startTime = HRClock.now();
    component_vec = vector<int>(N,-1);
    stop=0;
    iteration_number = 1;
    if (errorCheck)
        cout << "Constructing lattice-like network - Will do lots of extra calcs, enjoy ;)\n";
    G[0]=Graph(N);
    G[1]=Graph(N);
    periodic_bc = (depType == 1)? true : false;
    if ( lambda>0 ){
        isLambda=1;
        exp.param(lambda);
	if(triangular)
	  calculate_triangle_distances();
        else 
	  calculate_distances();
    //     cout << "calculate_distances();" <<endl;
        make_connectivity_links(0);
    //     cout << "make_connectivity_links(0);" <<endl;
        make_connectivity_links(1);
    //     cout << "make_connectivity_links(1);" <<endl;
        multiplex_dependency_vec.resize(N);
        for(int i = 0; i<N;i++) multiplex_dependency_vec[i] = randfloat()<q;
    }
    else{
        isGrid=1;
        create_interlinks();
    //printf("Completed table of interlinks.\n");
        create_lattices();
    }
    //alive = boost::get(::vertex_alive_t(),G[0]);//I think it's ok to just use the vertex map for one because it contains info on both and is indexed with unique ints
    alive.resize(N);
    color = boost::get(boost::vertex_color_t(),G[0]);//same as above
    Ng={{N,N}};

    firstP=1;
    firstIter=1;

    //D.resize(N);
    //for (int i =0; i<N; ++i)
//	D[N].resize(N);
    /*
    printf("Number of nodes in A %i\n",boost::num_vertices(A));
    printf("Number of nodes in B %i\n",boost::num_vertices(B));
    start_cascade(0.6);*/
    auto endTime = HRClock.now();
    profilingVector[INIT_TIME]=millisecond_res_diff(endTime,startTime);
}
/**
 * Constructor for random networks
 */
void RNet::Initialize(int L_, float kavg_, float q_, long unsigned int time_stamp_)
{
    L = L_;
    r = -1;
    N = L*L;
    q = q_;
    kavg = kavg_;
    lambda = -1;
    isRandom=1;
    time_stamp = time_stamp_;
    periodic_bc = (depType == 1)? true : false;    
    auto startTime = HRClock.now();
    component_vec = vector<int>(N,-1);
    stop=0;
    iteration_number = 1;
    if (errorCheck)
        cout << "Constructing ER network - Will do lots of extra calcs, enjoy ;)\n";
    G[0]=Graph(N);
    G[1]=Graph(N);

    make_connectivity_links(0);
    make_connectivity_links(1);
    multiplex_dependency_vec.resize(N);
    for(int i = 0; i<N;i++) multiplex_dependency_vec[i] = randfloat()<q;
    alive.resize(N);
    color = boost::get(boost::vertex_color_t(),G[0]);//same as above
    Ng={{N,N}};

    firstP=1;
    firstIter=1;
    //D.resize(N);
    //for (int i =0; i<N; ++i)
//      D[N].resize(N);
    /*
    printf("Number of nodes in A %i\n",boost::num_vertices(A));
    printf("Number of nodes in B %i\n",boost::num_vertices(B));
    start_cascade(0.6);*/
    auto endTime = HRClock.now();
    profilingVector[INIT_TIME]=millisecond_res_diff(endTime,startTime);
}

void RNet::set_strip(int s)
{
strip=s;
}
void RNet::setPcPrecision(int pcPrecision)
{
    pc_precision = pcPrecision>MINIMUM_PC_PRECISION?pcPrecision: MINIMUM_PC_PRECISION;
}

void RNet::turnProfileOn()
{
    profileOn=1;
}

void RNet::set_double_attack(int d)
{
    double_attack_flag = d;
}

void RNet::enableErrorCheck()
{
    errorCheck=1;
}
void RNet::disableErrorCheck()
{
    errorCheck=0;
}
void RNet::enableWriteLaplacian(int number_to_write){
    writeLaplacianOn = number_to_write;
}
void RNet::setVerbosity(int verbosity_)
{
    verbosity=verbosity_;
}

/** For resistor networks we typically take open boundary conditions
 *
 **/ 

void RNet::create_lattices(){
    
    for(int i =1; i<N; ++i ){
        if (i%L!=L-1){
	  boost::add_edge(i,i+1,G[0]);
          boost::add_edge(i,i+1,G[1]);
	}
	if (i<L*L-L){
	  boost::add_edge(i,i+L,G[0]);
          boost::add_edge(i,i+L,G[1]);
	}
    }
    //printf("Generated two lattices\n");
}

void RNet::calculate_triangle_distances(){
    int uniquePairs = static_cast<int>(L*(L+1)/2) -1; //-1 to disallow (0,0) -> force at least (0,1)
    tripencils.resize( uniquePairs );
   int m=0;
    for (int i = 1; i<L;++i){//start at ONE: no dist=0 links
        for (int j=0; j<=i;++j){
            try{
                if (m >= uniquePairs){
                    cout <<"I'm sorry Dave. I can't let "<<m<<" do that to ("<<i<<","<<j<<")..."<<endl;
                    return;
                }
                TriangleVector2D x(i,j);
           tripencils[m++] = TriDistancePair(x.get_length_squared(), x);
            }catch(exception &e){
                cout << "Pencil construction failed" <<endl;
                cout <<e.what() <<endl;
            }
        }
    }
    std::sort(tripencils.begin(), tripencils.end(), [](const TriDistancePair& dpA, const TriDistancePair& dpB){ return dpA.first < dpB.first;});

}

int RNet::draw_triangle_target(int sNode)
{
    if(isRandom){
        return randint(N);
    }
    double rexp,rexp2;
  do{
     rexp = exp(gen);
     rexp2 = rexp*rexp;}while(rexp>L/2);
  auto it = std::upper_bound(tripencils.cbegin(), tripencils.cend(), rexp2, [](const double A, const TriDistancePair B){return A<B.first;});
  if(it!=tripencils.cbegin() && ((*it).first + (*(it-1)).first > 2*rexp2 )  )
    it--;
  TriangleVector2D this_vec = it->second;
  if(gen()%2 == 0) this_vec.swap_coords();
  this_vec.rotate(gen());
  int tNode_i = sNode/L+this_vec.ux;
  int tNode_j = sNode%L+this_vec.uy;
  if(!periodic_bc){
    int trials=0;
  while( ((tNode_i >= L || tNode_i <0) || (tNode_j >= L || tNode_j <0)) && trials<6){
      this_vec.rotate(gen());
      tNode_i = sNode/L+this_vec.ux;
      tNode_j = sNode%L+this_vec.uy;
      trials++;
  }
  if (tNode_i >= L || tNode_i <0 || tNode_j >= L || tNode_j <0) //if you still can't find something, return failure
      return -1;
  }
  else{
    tNode_i = (tNode_i + L) %L;
    tNode_j = (tNode_j + L) %L;
  }
  return tNode_i*L+tNode_j;
}


void RNet::calculate_distances()
{
    int uniquePairs = static_cast<int>(L*(L+1)/2) -1; //-1 to disallow (0,0) -> force at least (0,1)
    pencils.resize( uniquePairs );
    //cout << "Resized pencils" <<endl;
    int m=0;
    for (int i = 1; i<L;++i){//start at ONE: no dist=0 links
        for (int j=0; j<=i;++j){
            try{
                if (m >= uniquePairs){
                    cout <<"I'm sorry Dave. I can't let "<<m<<" do that to ("<<i<<","<<j<<")..."<<endl;
                    return;
                }
           pencils[m++] = std::make_pair(i*i+j*j, std::make_pair(i,j));
            }catch(exception &e){
                cout << "Pencil construction failed" <<endl;
                cout <<e.what() <<endl;
            }
        }
    }
    std::sort(pencils.begin(), pencils.end(), [](const DistancePair& dpA, const DistancePair& dpB){ return dpA.first < dpB.first;});

}

int RNet::draw_k()
{
    return randfloat() < (kavg/2 - static_cast<int>(kavg/2))? kavg/2+1 : kavg/2;
}


int RNet::draw_target(int sNode)
{
    if(isRandom){
        return randint(N);
    }
    double rexp,rexp2;
  do{
     rexp = exp(gen);
     rexp2 = rexp*rexp;}while(rexp>L/2);
  auto it = std::upper_bound(pencils.cbegin(), pencils.cend(), rexp2, [](const int A, const DistancePair B){return A<B.first;});
  if(it!=pencils.cbegin() && ((*it).first + (*(it-1)).first > 2*rexp2 )  )
    it--;
  int di=randsign()*(*it).second.first;
  int dj=randsign()*(*it).second.second;
  if (randsign()>0){//switch i and j with 50% chance because only one pair appears in the pencils list
    std::swap(di,dj);
    }
  int tNode_i = sNode/L+di;
  int tNode_j = sNode%L+dj;
  if(!periodic_bc){
  if (tNode_i >= L || tNode_i <0)
    tNode_i-=(2*di);
  if (tNode_j >= L || tNode_j <0)
    tNode_j-=(2*dj);
  if (tNode_i >= L || tNode_i <0 || tNode_j >= L || tNode_j <0) //if you still can't find something, return failure
      return -1;
  }
  else{
    tNode_i = (tNode_i + L) %L;
    tNode_j = (tNode_j + L) %L;
}
  return tNode_i*L+tNode_j;
}



void RNet::make_connectivity_links(int lattice_index){
    int give_up_count = N/2;
    for (int i =0; i<kavg*N/2; ++i){
        bool added=false;
        int attempt_count=0;
        while(!added){
        int sNode = randint(N);
        int tNode = triangular? draw_triangle_target(sNode) : draw_target(sNode);
          if (tNode > 0 && !boost::edge(sNode,tNode,G[lattice_index]).second){
            boost::add_edge(sNode,tNode,G[lattice_index]);//should i make sure that sNode<tNode? don't think it matters
            added=true;
          }
          if (++attempt_count > give_up_count){
              if(verbosity>0)
                  cout<<"Unable to find a neighbor for node "<< sNode <<std::endl;
            break;
          }
            
        }
    }
    if(errorCheck){
        map<int,int> dd = getDegreeDistribution(lattice_index);
        double kavg_measured =0;
        for (auto it = dd.begin(); it != dd.end(); it++){
            cout<<it->first<<"\t"<<it->second<<"\n";
            kavg_measured+=it->first*it->second;
        }
        kavg_measured/=N;
        cout<<"k_avg (measured) : "<<kavg_measured<<endl;
        
        auto lengths = getLinkLengths(lattice_index);
        auto medianLength = lengths[lengths.size() / 2];
        auto meanLength = std::accumulate(lengths.begin(), lengths.end(), 0.0) / lengths.size();
        auto maxLength = *(lengths.crbegin());
        double varLength=0;
        for( double &i: lengths){
            varLength+=(i-meanLength)*(i-meanLength);
        }
        varLength/=lengths.size();
        
        ofstream debugFile;
        debugFile.open("/tmp/linklengths.json");
        debugFile << "{ \"mean\"   : " << meanLength << ",\n";
        debugFile << "  \"median\" : " << medianLength << ",\n";
        debugFile << "  \"max\"    : " << maxLength << ",\n";        
        debugFile << "  \"std\"    : " << sqrt(varLength) << ",\n";
        debugFile << " \"lengths\" : \n" ;
        jsonArray(lengths,debugFile);
        debugFile << "}\n";
        debugFile.close();
    }
      
}

map<int,int> RNet::getDegreeDistribution(int lattice_index){
    map<int,int> freq;
    for(int i = 0; i<N; ++i)
        freq[boost::degree(i,G[lattice_index])]++;
    return freq;
}

vector<double> RNet::getLinkLengths(int lattice_index){
    vector<double> lengths;
    auto edges = boost::edges(G[lattice_index]);
    for(auto it = edges.first; it!= edges.second; ++it){
        auto s = boost::source(*it, G[lattice_index]);
        auto t = boost::target(*it, G[lattice_index]);
        lengths.push_back(sqrt(distance_squared(s,t,L)));
    }
    std::sort(lengths.begin(), lengths.end());
    return lengths;
}

long int RNet::getOverlap()
{
    vector<int> neighbors;
    unsigned long overlap_count=0;
   for(int i=0; i<N; ++i){
       auto adj0 = boost::adjacent_vertices(i,G[0]);
       auto adj1 = boost::adjacent_vertices(i,G[1]);
       for(auto it0=adj0.first; it0!=adj0.second; ++it0){
           for(auto it1=adj1.first; it1!=adj1.second; ++it1){
               if (*it0 == *it1){
                   overlap_count++;
               }
           }
       }
   }
   if(overlap_count %2 != 0){
        throw logic_error("Odd number of 2*overlap encountered.");
   }
       
   return overlap_count/2;
}

double RNet::getClusteringCoefficient()
{
    // The clustering property, container, and map define the containment
    // and abstract accessor for the clustering coefficients of vertices.
  
    typedef boost::exterior_vertex_property<Graph, float> ClusteringProperty;
    typedef ClusteringProperty::container_type ClusteringContainer;
    typedef ClusteringProperty::map_type ClusteringMap;
    // Compute the clustering coefficients of each vertex in the graph
    // and the mean clustering coefficient which is returned from the
    // computation.
    ClusteringContainer coefs(boost::num_vertices(G[0]));
    ClusteringMap cm(coefs, G[0]);
    float cc = boost::all_clustering_coefficients(G[0], cm);
    return cc;
  
}

std::pair< long int, long int > RNet::getRadiusDiameter(int number_of_samples)
{

    typedef long VertSize;
    
   // std::vector<std::vector<VertSize> > distance_map(N);
    //std::vector<long> eccentricities(N);
      long r=N, d=0,e=0;
      trim_to_gc();
      std::vector<long> sampling_order(N);
      for(int i=0; i<N; ++i)
	sampling_order[i]=i;
      std::random_shuffle(sampling_order.begin(), sampling_order.end());
#pragma omp parallel for private(e) reduction(min: r) reduction(max: d)
    for(long i=0; i<number_of_samples; i++){
      int sample_idx=sampling_order[i];
      if(alive[sample_idx][0] != STILL_ALIVE_CODE) //do not include non-connected components
	continue;
      std::vector<long> source_d(N,-1);
      boost::breadth_first_search(G[0], sample_idx, 
	boost::visitor(
	  boost::make_bfs_visitor(	  
	    boost::record_distances(&source_d[0], boost::on_tree_edge())
	  )
	)
      );
      //distance_map[i] = source_d;
      e = * std::max_element(source_d.begin(), source_d.end());
      if (e>0){
	r = std::min(r,e);
	d = std::max(d,e);
      }
    }

  return std::make_pair(r,d);
  
}


void RNet::create_interlinks(){
    int neighborhood_size = (2*r+1)*(2*r+1);
    vector<int> candidates(neighborhood_size);
    vector<int> Anodes(N);
    no_neighbor_count = 0;
    unordered_map<int,bool>already_tried;
    
    for(  int i=0; i<N; ++i){ //Initialize node ids in source
        Anodes[i]=i;
	already_tried[i]=false;
        }
    if (neighborhood_size >= N){
	candidates.resize(N);
	candidates = Anodes;
    }
    random_shuffle ( Anodes.begin(), Anodes.end() );
    int node_coord[2] = {0,0};
    for( int& node_index : Anodes){ //Find neighbors for each node
	if (randfloat()>q)
	  continue;
        node_coord[0]=node_index/L;
	node_coord[1]=node_index%L;
        int xmin = node_coord[0]-r, xmax=node_coord[0]+r,ymin=node_coord[1]-r,ymax=node_coord[1]+r;
	/* restrict the neighborhood around the edges
	 * this is liable to increase the numbers of no-neighbor nodes in the corners
	 * but we'll chalk it up as a finite effect
	 */
	xmin = xmin > 0 ? xmin : 0;
	ymin = ymin > 0 ? ymin : 0;
	xmax = xmax < L ? xmax : L-1;
	ymax = ymax < L ? ymax : L-1;
	
	
        if(verbosity>0){
            printf("Node %i: (%i,%i)\n",node_index,node_coord[0],node_coord[1]);
            printf("Will now loop x from %i to %i and y from %i to %i\n",xmin,xmax,ymin,ymax);
        }
        /*TODO:Redo by by instantiating a constant map of (2r+1)x(2r+1) candidates 
            and making a copy for each node_index. Then set the values to the ids.  This
            will decrease the number of dynamic allocation calls to 1 + number of failed 
            candidates as opposed to the current (2r+1)^2 + number of fails.  
            Not critical because it's fast enough now but it might be necessary if we 
            increase L by an order of magnitude.
        */
	 int neighborhood_counter=0;
	if (neighborhood_size >= N){
	    
	}else{
	    for( int x = xmin; x<=xmax; ++x){
		for( int y = ymin; y<=ymax; ++y){
		    int x0 = x<0? x+L : x;
		    int y0 = y<0? y+L : y; //stupid sign of the dividend! http://en.wikipedia.org/wiki/Modulo_operation
		int neighbor_id = (x0%L)*L + (y0%L);
		    if(neighbor_id <0 || neighbor_id >= N)
			printf("I think that (%i,%i) is %i, is that so so crazy?\n",x,y,neighbor_id);
		    candidates[neighborhood_counter] = neighbor_id;
		    ++neighborhood_counter;

		    if(verbosity>1){        
			printf("(%i,%i)=%i,",x,y,neighbor_id);
		    }
		    
		}    
	    }
	   // assert(neighborhood_counter == neighborhood_size);
	}
	 int at_the_end = 0;
	bool success = false;
        //for( int &selected_neighbor : candidates){
	vector<int> small_candidate_vector;
	/* This section has been optimized (hopefully) by splitting the calculations into two regimes:
	 * When there are many candidates, we don't really care if the same one gets called multiple times.
	 * But once we've done a whole sweep and not found anything, we need to limit ourselves to the real candidates.
	 * At this point, though, it's not a big deal to push_back because it is a very small vector generated for only a small percentage of the nodes.
	 */
	
	while(!success){
	    int selected_neighbor=0,small_candidate_index=0;
	    if(at_the_end++ > neighborhood_counter){//we are in the regime of a small number of candidates, we did a full sweep and found nothing
		if (small_candidate_vector.size()==0){//if we haven't yet constructed a small_candidate_vector
		    for(int& potential_neighbor : candidates){
			if (verbosity>0)
			    printf("Checking neighbor %i\n",potential_neighbor);
			if(!already_tried[potential_neighbor]){
			  //  printf("Adding neighbor %i\n",potential_neighbor);
			    small_candidate_vector.push_back(potential_neighbor);//create a small vector with only the remaining candidates
			}
		    }
		}
		    if(small_candidate_vector.size()==0)//if we still didn't find anything, that means we have no candidates left and the loop must break
			break;
		    small_candidate_index = randint(small_candidate_vector.size());
		    selected_neighbor = small_candidate_vector[small_candidate_index];
		    
		}
	    else{
		selected_neighbor = candidates[randint(neighborhood_counter)];
		if( already_tried[selected_neighbor])
		    
		    continue;
	    }
   	    
            //auto link_it = interlinks.right.find(selected_neighbor);//ie, search for a value of selected_neighbor
            //if(link_it == interlinks.right.end()){//if it hasn't been seen
	    if(!already_tried[selected_neighbor]){
                interlinks.insert(IntBimap::value_type(node_index,selected_neighbor));
		//printf("%i\t%i\n", node_index,selected_neighbor);
		success=true;
		small_candidate_vector.clear();
		already_tried[selected_neighbor]=true;
                if(verbosity>0){
                    printf("From %i candidates, node %i selected.\n",neighborhood_counter,selected_neighbor);
                }
                break;   
            }

            if(at_the_end > neighborhood_counter){
		assert(false);
		printf("This line should never be called, right?\n");
		small_candidate_vector.erase(small_candidate_vector.begin() + small_candidate_index);
	    }
            
            
        }
	//already_tried.clear();
	if (!success){
	    if (verbosity>0){
		printf("Warning: No link for node %i.\n",node_index);
		fflush(stdout);
	    }
	    ++no_neighbor_count;
	}
        
    }
    

    int real_bad = N-interlinks.size();
    if (q == 1)
      assert(real_bad == no_neighbor_count);
    printf("Percentage without neighbors %.4f (%i/%i)\n",static_cast<float>(real_bad)/N,real_bad,N);
}
 
void RNet::start_cascade(double p)
{
    if (p>=1)
        return;
    
    if(double_attack_flag == 0){    
	for (int i =0; i<N; ++i) {
	    alive[i][0] = (randfloat()>p) ? iteration_number : STILL_ALIVE_CODE;//record the date of death (in iterations)
	    alive[i][1] = STILL_ALIVE_CODE; //initialize B to all alive
	}
    }else{
	for (int i =0; i<N; ++i) {
	    alive[i][0] = (randfloat()>p)? iteration_number : STILL_ALIVE_CODE;//record the date of death (in iterations)
	    alive[i][1] = (randfloat()>p)? iteration_number : STILL_ALIVE_CODE;//record the date of death (in iterations)
	}
    }    
    
    
}
void RNet::kill_all()
{
    int lattice_index = (iteration_number+1)%2;
    for(int node=0; node<N; ++node){
        if( alive[node][lattice_index] == STILL_ALIVE_CODE)
        {
            alive[node][lattice_index]=iteration_number;
        }
    }
    Ng[lattice_index] = 0;
    
}

void RNet::start_tracked_cascade(double p){
    stop=0;
    iteration_number=1;
    Ng={{N,N}};
    if (cascade_order[0].size() == 0){
	cascade_order[0].resize(N);
	for(int i=0; i<N; ++i)
	    cascade_order[0][i] = i;
	random_shuffle(cascade_order[0].begin(), cascade_order[0].end());
	if(double_attack_flag>0){
	    cascade_order[1].resize(N);
	    for(int i=0; i<N; ++i)
		cascade_order[1][i] = i;
	    random_shuffle(cascade_order[1].begin(), cascade_order[1].end());    
	}
    }
    int go_until = (1-p)*N;
    last_cascade_vertex = go_until;
    if(double_attack_flag>0){
	for(int i=0; i<N; ++i)
	    alive[i]={{STILL_ALIVE_CODE,STILL_ALIVE_CODE}};
	for(int i=0; i<go_until; ++i){
	    alive[cascade_order[0][i]][0] = static_cast<float>(iteration_number);
	    alive[cascade_order[1][i]][1] = static_cast<float>(iteration_number);
	}
	Ng={{N-go_until,N-go_until}};
    
    }else{
	for(auto node_it = cascade_order[0].cbegin(); node_it != cascade_order[0].cbegin()+go_until; ++node_it){
	    alive[*node_it] = {{static_cast<float>(iteration_number),STILL_ALIVE_CODE}};
    //	clear_b_vertex(*node_it,0);
	}
	for(auto node_it = cascade_order[0].cbegin()+go_until; node_it != cascade_order[0].cend(); ++node_it){
	    alive[*node_it] = {{STILL_ALIVE_CODE,STILL_ALIVE_CODE}};
	}
	Ng[0]-=go_until;
    }
}


void RNet::start_tracked_dual_cascade(double pA, double pB)
{
    stop=0;
    iteration_number=1;
    Ng={{N,N}};
    if (cascade_order[0].size() == 0){
        cascade_order[0].resize(N);
        for(int i=0; i<N; ++i)
            cascade_order[0][i] = i;
        random_shuffle(cascade_order[0].begin(), cascade_order[0].end());
        cascade_order[1].resize(N);
        for(int i=0; i<N; ++i)
            cascade_order[1][i] = i;
        random_shuffle(cascade_order[1].begin(), cascade_order[1].end());    
    }
    std::array<int,2>  go_until = {{static_cast<int>((1-pA)*N), static_cast<int>((1-pB)*N)}};
    //TODO:make last_cascade_vertex 2dimensional throughout the code
    last_cascade_vertex = go_until[0];
    
        for(int i=0; i<N; ++i)
            alive[i]={{STILL_ALIVE_CODE,STILL_ALIVE_CODE}};
        for(int i=0; i<go_until[0]; ++i){
            alive[cascade_order[0][i]][0] = static_cast<float>(iteration_number);
        }
        for(int i=0; i<go_until[1]; ++i){
            alive[cascade_order[1][i]][1] = static_cast<float>(iteration_number);
        }
    
        Ng={{N-go_until[0],N-go_until[1]}};
    
}


void RNet::remove_strip(){
    for(int i=0;i<L;++i){
	for(int j=0;j<r;++j){
	    if (double_attack_flag>0)
		alive[i*L+j] = {{static_cast<float>(iteration_number),static_cast<float>(iteration_number)}};
	    else
		alive[i*L+j][0]=static_cast<float>(iteration_number);
	}
    }
}

unordered_set<int> RNet::boost_giant(bool voltage_calc){
    int lattice_index = (iteration_number+1)%2;
    //vertex_is_alive<AliveMap> v_filter(boost::get(::vertex_alive_t(), G),lattice_index);
    //boost::filtered_graph<Graph, boost::keep_all, vertex_is_alive<AliveMap> > fg (G, boost::keep_all(), v_filter); 
    fill(component_vec.begin(),component_vec.end(),-1);
    for (int i=0; i<N; ++i)
	color[i] = alive[i][lattice_index] == STILL_ALIVE_CODE? Color::white() : Color::black();
    int num = voltage_calc? voltage_connected_components(&component_vec[0]) : custom_connected_components(G[lattice_index], &component_vec[0],color);
    map<int,int> comp_size;
    struct comp {
	bool operator()(pair<int,int> const &a, pair<int,int> const &b) {
	    return a.second > b.second;
	}
    };
    set <pair<int,int>, comp> sorted_sizes;
    for(int i=0; i<num; i++)
	comp_size[i]=0;
    for(int i=0; i<N; ++i)
	comp_size[ component_vec[i] ] +=1;
    for(int i=0; i<num; i++)
	sorted_sizes.insert(pair<int,int>(i,comp_size[i]));
    //for(auto  i  = sorted_sizes.cbegin(); i != sorted_sizes.cend(); ++i)
    //	printf( "Cluster %i : %i\n" , (*i).first,(*i).second);
    int max_index = (*(sorted_sizes.cbegin())).first;
    int Ng = (*(sorted_sizes.cbegin())).second;
    //printf("Largest cluster has %i nodes.\n", Ng);
    unordered_set<int> giant (Ng);
    for(int i=0; i<N; ++i)
	if (component_vec[i]==max_index)
	    giant.insert(i);
    return giant;

}
int RNet::get_Ng(){
 return( Ng[(iteration_number+1)%2]);   
}

vector<vector<int> > RNet::get_adjacency_list(int lattice_index){
	  vector<vector<int>> adjacency_list(N);
	  auto epair=boost::edges(G[lattice_index]);
	  for(auto ei = epair.first; ei!=epair.second; ++ei)
	  {
	    int sNode = boost::source(*ei,G[lattice_index]);
	    int tNode = boost::target(*ei,G[lattice_index]);
	    adjacency_list[sNode].push_back(tNode);
	    adjacency_list[tNode].push_back(sNode);
	  }
	return adjacency_list;
}


// void RNet::clear_b_vertex(int i, int lattice_index)
// {//do nothing
// return;
// }

/** If the topology is more random, the largest connected component is not 
 *  the sole determinant of voltage carrying.  We need st clusters, ie 
 *  clusters that have sources and sinks (targets).
 **/
int RNet::trim_to_st_clusters(){
    int lattice_index = (iteration_number+1)%2;
    fill(component_vec.begin(),component_vec.end(),-1);
    for (int i=0; i<N; ++i)
        color[i] = alive[i][lattice_index] == STILL_ALIVE_CODE? Color::white() : Color::black();
    int num = custom_connected_components(G[lattice_index], &component_vec[0],color);
    vector<bool> comp_has_source(num,false);
    vector<bool> comp_has_sink(num,false);
    //Iterating over sources
    for(int i=0; i<L; ++i){
        if (component_vec[i]>=0)
            comp_has_source[component_vec[i]]=true;
    }
    //Iterating over sinks
    for(int i=N-L; i<N; ++i){
        if (component_vec[i]>=0)
            comp_has_sink[component_vec[i]]=true;
    }
    //Nodes which were alive but are currently in clusters that do not have both sources and sinks, are now goners
    int number_removed=0;
    for(int i=0; i<N; ++i)
        if (component_vec[i]>=0 && ( !comp_has_source[ component_vec[i] ] || !comp_has_sink[ component_vec[i] ])){
               alive[i][lattice_index]=iteration_number+0.25;
               number_removed++;
        }
    Ng[lattice_index]-=number_removed;
    return number_removed;
}



int RNet::trim_to_gc(){

    int lattice_index = (iteration_number+1)%2;
    auto giant = boost_giant();
    int Ng0 = Ng[lattice_index];
    if (Ng0 == 0){
        return 0;
    }
    int Ng1 = giant.size();
    Ng[lattice_index] =  Ng1;
    for ( int node = 0; node<N; ++node){
	if( alive[node][lattice_index] != STILL_ALIVE_CODE){
	    continue;
	}
	if ( giant.find(node) != giant.cend()){
	    alive[node][lattice_index]=STILL_ALIVE_CODE;
	}
	else{
	    alive[node][lattice_index]=iteration_number+0.25;
// 	    clear_b_vertex(node,lattice_index);
	}
    }
    return(Ng0 - Ng1);
}

bool RNet::hasSpanningCluster(){
    int lattice_index = (iteration_number+1)%2;
    int hasSource=0,hasSink=0;
    for(int node = 0; node<L; ++node){
        if (alive[node][lattice_index] == STILL_ALIVE_CODE){
            hasSource=1;
            break;
        }
    }
    if (hasSource == 0)
    {
        return(false);
    }
    for(int node=N-L; node<N; ++node){
        if (alive[node][lattice_index] == STILL_ALIVE_CODE){
            hasSink=1;
            break;
        }
    }
    if (hasSink == 0)
    {
        return(false);
    }
    return(true);
}

int RNet::findVoltage(){
  auto start = HRClock.now();
  int lattice_index = (iteration_number+1)%2;
  vector<int> ground;

  node2matrix[lattice_index].resize(N,-1);
  fill( node2matrix[lattice_index].begin() , node2matrix[lattice_index].end(), -1);

  //Code to label initial voltage if necessary
//     for node in self.nnLattice:
//         self.nnLattice.node[node]["voltage"]=-1
//         self.nnLattice.node[node]["label"]=str(node)
    sources=0,sinks=0;
// code to figure out which, if any of the nodes are sources/sinks
    int giant_index=0;
for ( int node = 0; node<N; ++node){
  if (alive[node][lattice_index] == STILL_ALIVE_CODE){

    node2matrix[lattice_index][node]=giant_index;

    
    if (node<L){
     sources++; 
    }
/**ONLY ACCURATE FOR PURE LATTICES!!!!!!!**/    
//     else if (node<2*L){
//        if (alive[node-L][lattice_index] == STILL_ALIVE_CODE){
//            v0[giant_index-sources]=V_fixed;//obtained via -L*v_fixed and then removing the source nodes.  check it out, it works.
//        }
//     }
    else if (node>=N-L)
      sinks++;
  giant_index++;
  }
}

v0.resize(Ng[lattice_index]-sinks-sources);
fill(v0.begin(), v0.end(), 0);

for (int sNode = 0; sNode<L; ++sNode){
    if (alive[sNode][lattice_index] == STILL_ALIVE_CODE ){
        auto epair = boost::out_edges(sNode, G[lattice_index]);
        for (auto ei = epair.first; ei!= epair.second; ++ei){
            int tNode = boost::target(*ei,G[lattice_index]);
            if (alive[tNode][lattice_index] == STILL_ALIVE_CODE && tNode < N-L  && tNode >=L){
                if (errorCheck){
                    int tIndex = node2matrix[lattice_index][tNode] - sources;
                    if( tIndex < 0 ||  tIndex > static_cast<int>(v0.size()) ){
                        cerr << "Source "<< sNode << "is connected to node "<< tNode <<" with invalid index " << tIndex << " skipping for now...";
                        continue;
                    }
                }
                v0[node2matrix[lattice_index][tNode] - sources ] += V_fixed;//if we needed non-unity conductivity, this would be sigma*V_fixed
            }
        }

    }
}
if(errorCheck)
    assert(giant_index == Ng[lattice_index]);
if (sources*sinks == 0){
    if (verbosity>2)
        printf("No spanning cluster!\n");
 return(0); 
  
}

// cout<<"Sources\tSinks\n";
// cout<<sources<<"\t"<<sinks<<endl;
//Generate the Laplacian matrix
tripletList.clear();
tripletList.reserve((Ng[lattice_index]-sinks-sources)*7);

//cout<<"matrix dim should be " << v0.size() << endl;
vector<int> degree_vector(v0.size(),0);
auto epair=boost::edges(G[lattice_index]);
for(auto ei = epair.first; ei!=epair.second; ++ei)
{
  int sNode = boost::source(*ei,G[lattice_index]);
  int tNode = boost::target(*ei,G[lattice_index]);   
  if ( (alive[sNode][lattice_index] != STILL_ALIVE_CODE) || (alive[tNode][lattice_index] != STILL_ALIVE_CODE) ){
     continue;
  }
  bool sNodeInclude = ( sNode >= L && sNode < (N -L) ); // include source and target only if not on edges
  bool tNodeInclude = ( tNode >= L && tNode < (N -L) );
  if(sNodeInclude){
      degree_vector[node2matrix[lattice_index][sNode] - sources ]++; //any included node gets degree++ EVEN IF OTHER END IS NOT INCLUDED <--- absolutely necessary for L to be solvable
  }
  if( tNodeInclude ){
      degree_vector[node2matrix[lattice_index][tNode] - sources ]++;
  }
  if (sNodeInclude && tNodeInclude ){
    add_triplet(sNode,tNode,-1);
  }
}
//Save degrees for each relevant node into matrix.  we've already fixed their indices so we don't need the add_triplet function
for(unsigned int i=0;i<degree_vector.size();++i){
    //cout << i << "\t" << degree_vector[i] <<endl;
    tripletList.push_back(T(i,i,degree_vector[i]));
//     if (writeLaplacian>0){
//         laplacianFile << i+1 << "\t" << i+1 << "\t" << degree_vector[i] << "\n";
//     }
}

auto endTime= HRClock.now();
if (profileOn && verbosity>0){
    printf("Loading Laplacian data to triplet list: %.4f s\n", millisecond_res_diff(endTime,start));
}
// if (writeLaplacian>0){
//     laplacianFile.close();
// }
profilingVector[LAP_CREATE_TIME]+=millisecond_res_diff(endTime,start);
start=HRClock.now();
solveLaplacian();
endTime=HRClock.now();
profilingVector[LAP_SOLVE_TIME]+=millisecond_res_diff(endTime,start);

return(1);  
}


void RNet::writeLaplacian(){

    int lattice_index = (iteration_number+1)%2;
    string lfile = boost::str(  boost::format("Rnet_%i_%i_%lu_SPD.mtx") % lattice_index % writeLaplacianOn % time_stamp );
    string vfile = boost::str(  boost::format("Rnet_%i_%i_%lu_b.mtx") % lattice_index % writeLaplacianOn % time_stamp );

    laplacianFile.open(lfile.c_str());
    laplacianFile << "%%MatrixMarket matrix coordinate real symmetric\n";
    laplacianFile << v0.size() <<"\t"<<v0.size()<<"\t"<<tripletList.size()<<"\n"; 
    for (auto it = tripletList.begin(); it!= tripletList.end(); it++)
        laplacianFile << it->row()+1 << "\t" << it->col()+1 << "\t" << it->value() << "\n";
    laplacianFile.close();
   
    V0File.open(vfile.c_str());
    V0File << "\%\%MatrixMarket matrix array real general\n";
    V0File << v0.size() << "\t1\n";
    for (unsigned int i=0; i<v0.size(); i++)
        V0File<< v0[i] <<"\n";
    V0File.close();
    writeLaplacianOn--;
    
          
}

void RNet::solveLaplacian()
{
    auto start = HRClock.now();
    auto endTime= HRClock.now();
    int dim = v0.size();  
    int lattice_index = (iteration_number+1)%2;
    EigenSparse Laplacian(dim,dim);
    Eigen::VectorXd V(dim),V0(dim);



//     }catch(exception& e){
//         cerr << "Unable to initialize Laplacian object"<< endl;
//         cerr << "error: " << e.what() << endl;
//     }


    for (int i=0; i<dim; i++){
	V0[i] = v0[i];
    }
    

    // cout<<"V0:\t"<<V0<<endl;

    Laplacian.setFromTriplets(tripletList.begin(), tripletList.end());
 
try{
//     printf("Generating Laplacian from triplet list: %.4f s\n",millisecond_res_diff(endTime,start));
    // cout<<Laplacian.toDense()<<endl;
    // mat is ready to go!
    start = HRClock.now();
/** ConjugateGradient solver **/    
    //  
    //  Eigen::ConjugateGradient<EigenSparse, Eigen::Upper> solver;
    //  solver.compute(Laplacian);
    //  V = solver.solve(V0);
    //  endTime= HRClock.now();
    //  printf("ConjugateGradient (N=%i): %.4f s\n",dim,millisecond_res_diff(endTime,start));
/** SimplicialCholesky solver code **/
//         Eigen::SimplicialCholesky<EigenSparse> chol(Laplacian);
//         V = chol.solve(V0);
//         endTime= HRClock.now();
//         printf("SimplicalCholesky (N=%i): %.4f s\n",dim,millisecond_res_diff(endTime,start));
//         start = HRClock.now();
/** CholmodSupernodalLLT solver code **/
//     Eigen::CholmodSupernodalLLT<EigenSparse> cholsup(Laplacian);
//     V = cholsup.solve(V0);
//     endTime= HRClock.now();
       
       //if (dim>850000)
//        try{

/** Eigen LLT -- it seems LDLT is faster but since LDLT doesn't assume SPD, I think the marginal savings will disappear in larger matrices **/
    if(!isRandom && (lambda<=0 || lambda >=0.7) ){ //it can't be ER and has to either be a pure lattice lambda<0 or an almost pure lattice lambda>0.7
	Eigen::SimplicialLLT<EigenSparse, Eigen::Lower> solver;
	solver.compute(Laplacian);
	V = solver.solve(V0);
	endTime = HRClock.now();
	if(verbosity>1){
	    cout<<"SimplicialLLT (dim="<<dim<<"): "<<millisecond_res_diff(endTime,start)<<"s\n";}
    }else if(!isRandom ){//&& lambda > 0.3){ // intermediate systems seem to be optimized with bicgstab
	Eigen::BiCGSTAB<EigenSparse> solver;
	solver.setTolerance(eps); 
	solver.setMaxIterations(MAX_ITERATIONS);
	solver.compute(Laplacian);
	V = solver.solve(V0);
	endTime = HRClock.now();
	if(verbosity>1){
	    cout<<"BiCGSTAB (dim="<<dim<<",solver_iters="<<solver.iterations()<<"): "<<millisecond_res_diff(endTime,start)<<"s\n";}
    }else{ //random systems work best with cg
	Eigen::ConjugateGradient<EigenSparse, Eigen::Lower> solver;
	solver.setTolerance(eps); 
	solver.setMaxIterations(MAX_ITERATIONS);
        solver.compute(Laplacian);
	V = solver.solve(V0);
	endTime = HRClock.now();
        if(verbosity>1){
            cout<<"ConjugateGradient (dim="<<dim<<",solver_iters="<<solver.iterations()<<"): "<<millisecond_res_diff(endTime,start)<<"s\n";}
        

        if (solver.iterations() == MAX_ITERATIONS){
                cerr<<"Wait a minute something's wrong " << solver.error() << endl;
                if (solver.error() > eps || std::isnan(solver.error())){
                    throw logic_error("Unable to solve matrix. Probably unsolvable.\n");
                }
        }
    }
    if (errorCheck){
        double err=0;
        Eigen::VectorXd err_vec = Laplacian*V - V0;
        for(int i =0; i< dim; ++i)
            err+= err_vec[i]*err_vec[i];
        cerr<<"Total error per term in |LV - V0| = "<<sqrt(err)/dim<<endl;        
    }
  
       }catch(exception& e){
           cerr << "Unable to solve Laplacian of size " << dim;
           Eigen::SparseLU<EigenSparse> slusolver;
           slusolver.isSymmetric(true);
           slusolver.compute(Laplacian);
           if (slusolver.info() != Eigen::Success)
                cerr<< endl << slusolver.lastErrorMessage() <<endl;
           else{
            double logabsdet = slusolver.logAbsDeterminant();
            cerr << " det = "<< logabsdet << endl;
            cerr << "error: "<<e.what() <<"\n";
            if ((logabsdet) > 0){
                Eigen::ConjugateGradient<EigenSparse, Eigen::Lower> solver;
                solver.setTolerance(eps); 
                solver.setMaxIterations(MAX_ITERATIONS);
                solver.compute(Laplacian);
                V = solver.solveWithGuess(V0,V);
                if (solver.iterations() == MAX_ITERATIONS){
                    stringstream ss;
                    ss << "CG unable to solve matrix of dim "<< dim <<" with nnz "<< tripletList.size() <<"after " << 2*MAX_ITERATIONS << " error is now: " << solver.error() ;
                    throw logic_error(ss.str());
                }
                
            }
        }

           return;
       }

    voltage[lattice_index].resize(dim);
    for (int i = 0; i <dim; ++i){
	voltage[lattice_index][i] = V[i];
    }
    // cout<<"V:\n"<<V<<endl;
}



void RNet::enableWriteDynamics()
{
    writeDynamics=1;
    string fname = boost::str( boost::format("RNetDynamics-type=%i-L=%i-r-%i-k=%.4f-l=%.4f-q=%.4f-%lu.json") % depType % L % r % kavg % lambda % q % time_stamp);
    dynamicsFile.open(fname.c_str());

    dynamicsFile <<  boost::str( boost::format("{\"type\":%i, \"L\":%i, \"q\":%.6f, \"time_stamp\":%lu, ") % depType % L % q % time_stamp);
    if (isGrid)
        dynamicsFile << "\"topology\": \"grid\", \"r\" : " << r <<",";
    else if (isRandom)
        dynamicsFile << "\"topology\": \"ER\", \"k\" : " << kavg <<",";
    else if (isLambda)
        dynamicsFile << "\"topology\": \"lambda\", \"k\" : " << kavg <<", \"lambda\" : "<< lambda << ",";
    if (depType == 3)
        dynamicsFile << "\"V0\" :" << getFixedVoltage() << ",";
    dynamicsFile << "\"runs\" : [\n";
}

void RNet::closeDynamicsFile()
{
    dynamicsFile << "]";
    if (real_pc)
      dynamicsFile << ", \"pc\" : " <<real_pc << ", \"pc_precision\" : " << pc_precision ;
     dynamicsFile << "}";
    dynamicsFile.close();
}

void RNet::add_triplet(int node1, int node2, int value)
{
    int lattice_index =  (iteration_number+1)%2;
    //do not allow nodes on either edge to be calculated
    //they have fixed voltage
    if (node1 >= N-L || node2 >=N-L || node1 <L || node2<L)
      return;

    int i = node2matrix[lattice_index][node1] - sources;

    if (node1 == node2){
      tripletList.push_back(T(i,i,value));
//       if (writeLaplacian>0){
// 	laplacianFile << i+1 << "\t" << i+1 << "\t" << value << "\n";
//       }
//       
//       cout<<i<<"\t"<<i<<"\t"<<node1<<"\t"<<node2<<"\t"<<value<<endl;
    }else{

    int j = node2matrix[lattice_index][node2] - sources;

/**TODO:For some solvers only lower triangular is needed but need to stay synced  **/
//    if (i>j) //only LOWER triangle needs to be saved -- for some solvers!
//        tripletList.push_back(T(i,j,value));
//      else
//        tripletList.push_back(T(j,i,value));
    tripletList.push_back(T(i,j,value));
    tripletList.push_back(T(j,i,value));

    //       if(writeLaplacian>0)
//       {
// 	laplacianFile << j+1 << "\t" << i+1 << "\t" << value << "\n";
//       }
//       cout<<i<<"\t"<<j<<"\t"<<node1<<"\t"<<node2<<"\t"<<value<<endl;

    }
}
   template <class ComponentMap>
//   inline typename boost::property_traits<ComponentMap>::value_type
 int RNet::voltage_connected_components(ComponentMap c)
  {
    int lattice_index =  (iteration_number+1)%2;
    if (boost::num_vertices(G[lattice_index]) == 0) return 0;

    int c_count=0;
    queue_t Q;
    components_recorder<ComponentMap> vis(c, c_count);
    for(unsigned int v =0; v<boost::num_vertices(G[lattice_index]); ++v){
	if (color[v] == Color::white()){
	    ++c_count;
	    voltage_breadth_first_visit(Q,v,vis);
	}
    } 
    return c_count + 1;
  }

template <class ComponentMap>
void RNet::voltage_breadth_first_visit(queue_t & Q, int v, components_recorder<ComponentMap> vis){
    using namespace boost;
    boost::graph_traits<Graph>::out_edge_iterator ei,ei_end;
    int lattice_index =  (iteration_number+1)%2;
    put(color, v, Color::gray());             vis.discover_vertex(v, G[lattice_index]);
    Q.push(v);
    
    double uVoltage,vVoltage;
    while (! Q.empty()) {
        int u = Q.top(); Q.pop();            vis.examine_vertex(u, G[lattice_index]);

        int uVoltageIndex = node2matrix[lattice_index][u];

        if(uVoltageIndex<0)//node is dead, should not happen
        {
             continue;
        } else if(uVoltageIndex<sources)
        {
             uVoltage = getFixedVoltage();//source
        } else if(uVoltageIndex - sources< static_cast<int>( voltage[lattice_index].size()))
        {             
             uVoltage = voltage[lattice_index][uVoltageIndex-sources];//calculated value
        } else {
             uVoltage = 0;//sink 
        }
        for (boost::tie(ei, ei_end) = out_edges(u, G[lattice_index]); ei != ei_end; ++ei) {
            int v = target(*ei, G[lattice_index]);            vis.examine_edge(*ei, G[lattice_index]);

            int vVoltageIndex = node2matrix[lattice_index][v];

        if(vVoltageIndex<0)
        {
            continue;
           
        }else if(vVoltageIndex<sources)
        {
            vVoltage = getFixedVoltage();
        } else if(vVoltageIndex - sources < static_cast<int>(voltage[lattice_index].size())){
            vVoltage = voltage[lattice_index][vVoltageIndex-sources]; 
        } else {
            vVoltage = 0;
        }
            if (fabs(uVoltage-vVoltage) < eps){
                continue;}
            ColorValue v_color = get(color, v);
            if (v_color == Color::white()) {      vis.tree_edge(*ei, G[lattice_index]);
                put(color, v, Color::gray());       vis.discover_vertex(v, G[lattice_index]);
                Q.push(v);
            } else {                              vis.non_tree_edge(*ei, G[lattice_index]);
            if (v_color == Color::gray())       vis.gray_target(*ei, G[lattice_index]);
            else                                vis.black_target(*ei, G[lattice_index]);
            }
        } // end for
        put(color, u, Color::black());          vis.finish_vertex(u, G[lattice_index]);
    } // end while
} // breadth_first_visit
  
double RNet::findTotalCurrent(){
    int lattice_index = (iteration_number+1)%2;
    total_current=0;
    for(int node = N-L; node<N; ++node){
	if( (alive[node][lattice_index] != STILL_ALIVE_CODE) ){
            continue;
        }
        auto eipair = boost::out_edges(node, G[lattice_index]);
        for (auto ei = eipair.first; ei != eipair.second; ei++){
            int tNode = boost::target(*ei, G[lattice_index]);   
            if (alive[tNode][lattice_index] != STILL_ALIVE_CODE || tNode >= N-L)
                continue;
            int u = node2matrix[lattice_index][node];// - sources;
            int v = node2matrix[lattice_index][tNode];// - sources;
            if (u<0 || v<0){
                printf("I thought these were alive... This should not happen!\n");
                continue;
            }
            total_current += fabs(voltage[lattice_index][v - sources]);
         }
    }
    if(errorCheck){
        double total_out_current = total_current;
        total_current = 0;
//  This is to check if the in and out current are the same, when running simulations it's wasted cycles.
        for(int node = 0; node<L; ++node){
            if( (alive[node][lattice_index] != STILL_ALIVE_CODE) ){
                continue;
            }
            auto eipair = boost::out_edges(node, G[lattice_index]);
            for (auto ei = eipair.first; ei != eipair.second; ei++){
                int tNode = boost::target(*ei, G[lattice_index]);   
                if (alive[tNode][lattice_index] != STILL_ALIVE_CODE || tNode <L)
                    continue;
                int u = node2matrix[lattice_index][node];// - sources;
                int v = node2matrix[lattice_index][tNode];// - sources;
                if (u<0 || v<0){
                    printf("I thought these were alive... This should not happen!\n");
                    continue;
                }
                total_current += fabs(V_fixed - voltage[lattice_index][v - sources] );
            }
        }
        if (fabs (total_out_current - total_current)/ (total_current) >  maxRelativeError ){
            cerr << "In/Out current not equal: " << setprecision(10) << total_current << " / " << total_out_current << "\n";
            cerr << "Error :" << fabs (total_out_current - total_current)/ (total_current) << " > " << maxRelativeError << "\n";
            writeVoltage();
            throw logic_error("BAD_CURRENT_CALC");
        } else {
            cout << "In/Out current equal: " << total_current << " / " << total_out_current << "\n";
            writeVoltage();
        }

    }
    return total_current;

}  
int RNet::getFixedVoltage(){
    return(V_fixed);
}

void RNet::writeVoltage(){
    int lattice_index =  (iteration_number+1)%2;
    string fname = boost::str( boost::format("V%i-L=%i-k=%.4f-l=%.4f-q=%.4f-%lu.csv") % lattice_index % L % kavg % lambda % q % time_stamp);
    ofstream vFile(fname.c_str());
    for (int i = 0; i<L; ++i)
        if (alive[i][lattice_index]==STILL_ALIVE_CODE)
            vFile << i << "\t" << getFixedVoltage() << "\n";
    for (int i = L; i<N-L; ++i)
        if (alive[i][lattice_index]==STILL_ALIVE_CODE)
            vFile << i << "\t" << voltage[lattice_index][node2matrix[lattice_index][i] - sources] << "\n";
    for (int i = N-L; i<N; ++i)
        if (alive[i][lattice_index]==STILL_ALIVE_CODE)
            vFile << i << "\t" << 0 << "\n";
    vFile.close();
    
}
int RNet::trim_to_backbone(){

    int lattice_index = (iteration_number+1)%2;
    auto giant = boost_giant(true);
    int Ng0 = Ng[lattice_index];
    int Ng1 = giant.size();
    Ng[lattice_index] =  Ng1;
    for ( int node = 0; node<N; ++node){
        if( alive[node][lattice_index] != STILL_ALIVE_CODE){
            continue;
        }
        if ( giant.find(node) != giant.cend()){
            alive[node][lattice_index]=STILL_ALIVE_CODE;
        }
        else{
            alive[node][lattice_index]=iteration_number+0.5;
//          clear_b_vertex(node,lattice_index);
        }
    }
    return(Ng0 - Ng1);
}


int RNet::kill_dependent()
{
    int lattice_index_from = (iteration_number+1)%2;
    int lattice_index_to = (iteration_number)%2;
    int number_killed = 0, no_neighbor_found = 0;
    if (q>0){
        if (lambda>0 || isRandom){
            for (int i = 0; i<N; ++i){
                if (alive[i][lattice_index_from] != STILL_ALIVE_CODE && multiplex_dependency_vec[i]){
                    if ( alive[i][lattice_index_to] == STILL_ALIVE_CODE){
                        alive[i][lattice_index_to] = iteration_number;
                        number_killed++;
                        
                    }
                }
            }

        }else{
            //TODO: make this nicer.  it's uuuuuugggllly
            if (lattice_index_to){ //ie, if we are going to B
                for( int i =0; i<N;i++){
                    if (alive[i][lattice_index_from] != STILL_ALIVE_CODE){
                        auto left_it = interlinks.left.find(i);
                        if (left_it == interlinks.left.end()){
                            ++no_neighbor_found;
                            continue;}
                        int dependent = left_it->second;
                        if (alive[dependent][lattice_index_to] == STILL_ALIVE_CODE){
        // 		    clear_b_vertex(dependent,lattice_index_to);
                            alive[dependent][lattice_index_to] = iteration_number;
                            number_killed++;
                        }
                    }
                }
                
            }else{
                for( int i =0; i<N;i++){
                    if (alive[i][lattice_index_from] != STILL_ALIVE_CODE){
                        auto right_it = interlinks.right.find(i);
                        if (right_it == interlinks.right.end()){
                            ++no_neighbor_found;
                            continue;}
                        int dependent = right_it->second;
                        if (alive[dependent][lattice_index_to] == STILL_ALIVE_CODE){
        // 		    clear_b_vertex(dependent,lattice_index_to);
                            alive[dependent][lattice_index_to] = iteration_number;
                            number_killed++;
                        }
                    }
                }
            }
        }
    }
    //printf("%i\t[%i -> %i] Total killed (actual): %i (%i)\n", iteration_number,lattice_index_from, lattice_index_to,number_killed + no_neighbor_found,number_killed);
    ++iteration_number;
    Ng[lattice_index_to]-=number_killed;
    stop = (number_killed == 0) ? stop+1 : 0;
    return(number_killed);
}
    
void RNet::continue_cascade(double p){
    stop=0;   

    iteration_number=1;

    int go_until = (1-p)*N;
    if(double_attack_flag>0){
	for(int i=last_cascade_vertex; i<go_until;++i){
	    alive[cascade_order[0][i]][0]=iteration_number;
	    alive[cascade_order[1][i]][1]=iteration_number;
	}
    }else{
	for(auto node_it = cascade_order[0].begin()+last_cascade_vertex; node_it != cascade_order[0].begin()+go_until; ++node_it){
	    alive[*node_it][0] = iteration_number;
    //	clear_b_vertex(*node_it,0);
	}
    }
}



void RNet::write_history(int lattice_index){
    string double_string = (double_attack_flag>0)? "double" : "single";
    string lattice_string = (lattice_index==0)? "A" : "B";
    double p = last_cascade_vertex>0? 1- static_cast<double>(last_cascade_vertex)/N : 0;
    string fname =boost::str(  boost::format("LatticeHistory-L=%i-q=%.5f-lambda=%.4f-k=%.4f-p=%.9f-%s-%s-%lu.matrix") % L % q % lambda % kavg% p% double_string.c_str() % lattice_string.c_str() % time_stamp );
    ofstream file(fname.c_str());
    file << L << "\t" << q << "\t" << r << "\t" << "\n";
    for (int i=0; i<N; ++i){
	file << alive[i][lattice_index] << "\t";
    }
    file.close();
}

// double RNet::get_system_xi(){
//     int lattice_index = 0;//(iteration_number+1)%2;
//     fill(component_vec.begin(),component_vec.end(),-1);
// 
//     for (int i=0; i<N; ++i)
// 	color[i] = alive[i][lattice_index] == STILL_ALIVE_CODE?  Color::black() : Color::white() ;
//     
//     int num = custom_connected_components(G, &component_vec[0],color);
//     
//     uint64_t total_distance=0;
//     uint64_t denominator=0;
// 
//     for (int i=0; i<N; ++i){
// 	if (component_vec[i] == -1)
// 	    continue;
// 	for (int j=i; j<N; ++j){
// 	    if (component_vec[i] == component_vec[j]){
// 		if (total_distance + L > UINT64_MAX){
// 		    printf ("Overflowing!\n");
// 		    return sqrt( static_cast<double>(total_distance)/denominator );
// 		    }
// 		total_distance+=distance_squared(i,j,L);
// 		denominator++;
// 	    }
// 	}
//     }
//     double xi = sqrt( static_cast<double>(total_distance)/denominator );
//     //hole_history.push_back(xi);
//     return xi;
//     
// }

vector<int> RNet::get_giant_hole()
{
    int lattice_index = (iteration_number+1)%2;
    //vertex_is_alive<AliveMap> v_filter(boost::get(::vertex_alive_t(), G),lattice_index);
    //boost::filtered_graph<Graph, boost::keep_all, vertex_is_alive<AliveMap> > fg (G, boost::keep_all(), v_filter); 
    fill(component_vec.begin(),component_vec.end(),-1);
    for (int i=0; i<N; ++i)
	color[i] = alive[i][lattice_index] == STILL_ALIVE_CODE?  Color::black() : Color::white() ;
    int num = custom_connected_components(G[lattice_index], &component_vec[0],color);
    map<int,int> comp_size;
    struct comp {
	bool operator()(pair<int,int> const &a, pair<int,int> const &b) {
	    return a.second > b.second;
	}
    };
    set <pair<int,int>, comp> sorted_sizes;
    for(int i=0; i<num; i++)
	comp_size[i]=0;
    for(int i=0; i<N; ++i)
	comp_size[ component_vec[i] ] +=1;
    for(int i=0; i<num; i++)
	sorted_sizes.insert(pair<int,int>(i,comp_size[i]));
    //for(auto  i  = sorted_sizes.cbegin(); i != sorted_sizes.cend(); ++i)
    //	printf( "Cluster %i : %i\n" , (*i).first,(*i).second);
    int max_index = (*(sorted_sizes.cbegin())).first;
    int Ng = (*(sorted_sizes.cbegin())).second;
    //printf("Largest cluster has %i nodes.\n", Ng);
    vector<int> giant (Ng);
    int k=0;
    for(int i=0; i<N; ++i)
	if (component_vec[i]==max_index)
	    giant[k++]=i;
    return giant;

}
void write_cluster(vector<int> cluster, string filename){
        ofstream file(filename.c_str());
	for (int& c : cluster){
		file << c << "\n";
	}
	file.close();
    
}



double get_cluster_xi(vector<int> cluster, int L){
    uint64_t numerator=0,denominator=0;
    int N=L*L;
    unsigned int M=cluster.size();
    uint64_t expected=M*(static_cast<uint64_t>(M)+1)/2;//stupid 32bit machines x(
    printf("Calculating xi for hole of size %i:\n",M);
    for(unsigned int i=0; i<M; ++i){
	printf("\r[%.4f%%] complete.",100*static_cast<double>(denominator)/expected);
	for (unsigned int j=i; j<M; ++j){
	    if (numerator + N > UINT64_MAX){
		printf ("Overflowing!\n");
		return sqrt( static_cast<long double>(numerator)/denominator);
	    }
	    
	    numerator += distance_squared(cluster[i],cluster[j],L);
	    denominator++;
	}
    }
    printf("\n");
    double xi = sqrt(static_cast<long double>(numerator)/denominator);
    if(xi*xi>M){
	ofstream file("bad_xi.dbg");
	for(unsigned int i=0; i<M; ++i){
	    for (unsigned int j=i; j<M; ++j){
		file << i <<"\t"<<j<<"\t"<< cluster[i] <<"\t"<< cluster[j] << "\t" << distance_squared(cluster[i],cluster[j],L)<<"\n";
	    }
	}
	file.close();
	assert(false);
    }
    
    if(denominator == expected)
	return (xi);
    else{
	cout<<"Numerator:\t"<<numerator<<"\n";
	cout<<"Denominator:\t"<<denominator<<"\n";
	cout<<"Mass:\t\t"<<M<<"\n";
	cout<<"Expected:\t"<<expected<<"\n";
	assert(false);
    }

}

/** Make sure to pass the actual graph G[0] and not the graph array G
 * */
/*DEPRECATED
int num_live_vertices(Graph&  G, int lattice_index){
    boost::graph_traits<Graph>::vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(G);
    
    AliveMap alive = boost::get(::vertex_alive_t(),G);
    int live_vertices;
    live_vertices=0;
    
    for (next = vi; vi != vi_end; vi = next) {
	++next;
	
	if(alive[*vi][lattice_index] == STILL_ALIVE_CODE){
	    live_vertices++;
	    continue;
	}
	
	
    }
    return live_vertices;
    
 
}
*/
void RNet::setDependencyType(int dependency_type)
{
    depType=dependency_type;
}


/*	End the functions in RNet . 
 * 	These should be in a different file but I can't think of a name so here they are.
 * 	These used to comprise main.cpp.
 */


/** The three dependency types are:
 * 0 - no dependency, just calculate the voltage and conductance
 * 1 - standard dependency, ignore voltage, just look for giant cluster
 * 2 - geometric voltage dependency, ignore current levels, just make sure there is a spanning cluster
 * 3 - dynamic voltage dependency, only follow links with current when constructing the giant cluster (ie backbone only)
 * 
 **/
void RNet::tracked_cascade_single_p(double p, int create_images)
{
    start_tracked_cascade(p);
    int delta1=0,delta2=0;
    double S=0;
    //int lattice_index = (iteration_number+1)%2; //uncomment as needed, 
    if(strip)
        remove_strip();
    if(writeDynamics){
     if (!firstP)
         dynamicsFile << ",";
     dynamicsFile << "{\"p\":"<<setprecision(10)<<p<<",\"history\": [ ";
     firstP=0;
    }
    firstIter=1;
    auto total_time_start = HRClock.now();
    auto startTime = HRClock.now();
    auto endTime= HRClock.now();
    while(stop<2){
        
         //Random networks will have current all over, GCC is not a useful filter here -- suspended for now, makes unsolvable matrices...
        startTime = HRClock.now();
        delta1 = isRandom ? trim_to_st_clusters() : trim_to_gc();
        //delta1 = trim_to_gc();
        endTime= HRClock.now();
        profilingVector[BFS_TIME]+=millisecond_res_diff(endTime,startTime);
        if(profileOn && verbosity>0){
             printf("BFS and Giant trim (d=%i): %.4f s\n",delta1,millisecond_res_diff(endTime,startTime));
        }
        
        switch(depType){
            case 0:
                delta2=N-delta1;
                if(findVoltage()){
                    delta2 = trim_to_backbone();
                    S = findTotalCurrent()/V_fixed;
                }
                summary_map[p] = run_summary(0,p,N-delta1,N-delta1-delta2,S);
                stop=2;
                break;
            case 1:
            case 2:
                if(!hasSpanningCluster()){
             //       printf("Lattice %i has no spanning cluster.\n",(iteration_number+1)%2);
                    kill_all();//in this lattice
                }
                break;
            case 3:
                if(findVoltage()){
                    startTime = HRClock.now();
                    if(writeDynamics)
                    {
                        S = findTotalCurrent()/getFixedVoltage();
                    }
                    trim_to_backbone();
                    endTime= HRClock.now();
                    profilingVector[VBFS_TIME]+=millisecond_res_diff(endTime,startTime);
                    if(profileOn && verbosity>0){
                    printf("Trim to backbone : %.4f s\n",millisecond_res_diff(endTime,startTime));
                    }
                    
                }else{
                    kill_all();
                    S=0;
                }
                break;
        }
        kill_dependent();//this guy advances iteration_number, be advised
//        printf("%.6f\t%i\t%i\t%i\n",p,Ng[0],Ng[1],iteration_number);
        if(writeDynamics)
        {
                if(!firstIter)
                    dynamicsFile << ",\n";
                if (depType == 3 || depType == 0)
                {
                    //This calculation of lattice index is offset by one because the iteration number gets advanced by kill_dependent
                    dynamicsFile << "{\"Ng\":"<<Ng[iteration_number%2]<<",\"S\":"<<S<<"}";
                }
                else
                {
                    dynamicsFile << "{\"Ng\":"<<Ng[iteration_number%2]<<"}";
                }
                firstIter=0;
                dynamicsFile.flush();
        }
    }
    if (depType==1){
        summary_map[p] = run_summary(depType,p,Ng[0],iteration_number);
    } else if (depType >1){
        
        S = findVoltage()? findTotalCurrent()/V_fixed : 0;
        summary_map[p] = run_summary(depType,p,Ng[0],iteration_number,S);
    }
    if(writeDynamics){
        dynamicsFile << "]";
        dynamicsFile << ",\"finalS\":" << S;
        dynamicsFile << ",\"finalN\":" << Ng[0];
        dynamicsFile << ",\"time\":"<< millisecond_res_diff(HRClock.now(),total_time_start);
        dynamicsFile << "}";
        dynamicsFile.flush();
    }

    
    
}

/* fairly useless...
void RNet::single_attack_single_p(double p, int create_images) {
    
    char * fname = new char [200];
    string double_name = double_attack_flag? "double.attack" : "single.attack";
    auto start=HRClock.now();
    printf("p = %.9f:  ",p);
    start_tracked_cascade(p);
    if(strip>0)
	    remove_strip();
    auto endTime=HRClock.now();
    sprintf(fname,"DependentLattices-L=%i-r=%i-p=%.9f-q=%.3f-%s.%lu.output",L,r,p,q,double_name.c_str(),time_stamp);
    FILE *f = fopen(fname,"w");
    fprintf(f,"%i\t%i\t%.5f\n",L,r,p);
    fprintf(f,"%i\t%i\t%i\n",iteration_number,Ng[0],Ng[1]);
   while(stop<=1){
        fprintf(f,"%i\t%i\t%i\n",iteration_number,Ng[0],Ng[1]);
        int delta1 = trim_to_gc();
        //printf("BFS trim to gc: %.4f s\n",millisecond_res_diff(endTime,start));
        if( findVoltage()){
          int delta2 = trim_to_backbone();
         kill_dependent();
          double S = findTotalCurrent()/V_fixed;
          printf("p\tN\tNg\tNb\tsigma\n%.4f\t%i\t%i\t%i\t%.4f\n",p,N,N-delta1,N-delta1-delta2,S);
          summary_map[p] = run_summary(p,N-delta1,N-delta1-delta2,S);
        } else {
          printf("p\tN\tNg\tNb\tsigma\n%.4f\t%i\t%i\t%i\t%.4f\n",p,N,N-delta1,0,0);
          summary_map[p] = run_summary(p,N-delta1,0,0);
          return;
        }

   }
   fclose(f);
   delete[] fname;
   printf("(%i,%i)\n",iteration_number, Ng[0]);
    if(create_images)
	    write_history(0);
    if(strip>0){
	    summary_map[p] = run_summary(p,Ng[0],iteration_number);
	    return;
    }

    
    double xi;
   auto hole = get_giant_hole();
    int hole_size = hole.size();
    if (hole_size < 100000){
	write_cluster(hole,string("giant_hole.list"));
	xi = get_cluster_xi(hole,L);
    }else{
	printf( "%i is too big to measure. Assuming circular hole.\n",hole_size);
	xi = sqrt(hole_size/3.14159);
    }
    xi=get_system_xi();

}**/



void RNet::continued_cascade(double p0, double p_step, double zero_size, int create_images){
    
    auto start=HRClock.now();
    char * fname = new char [200];
    double p = p0;
    string double_name = double_attack_flag? "double.attack" : "single.attack";
    printf("Beginning cascade with p = %.3f...\n",p);
    start_tracked_cascade(p);
    auto endTime=HRClock.now();
    sprintf(fname,"DependentLattices-L=%i-r=%i-p=%.5f-q=%.5f.%s.%lu.continued_cascade.output",L,r,p,q,double_name.c_str(),time_stamp);
    FILE *f = fopen(fname,"w");
    fprintf(f,"%i\t%i\t%.5f\t%.5f\n",L,r,p,q);
    
    start=HRClock.now();
    while(stop<2){
	//int not_in_gc = 
	trim_to_gc();
    //  printf("%i : Ng(trimmed) = %i (%i)\t",iteration_number, get_Ng(),not_in_gc);
	kill_dependent();
	if( stop>1 && Ng[0] > zero_size*N ){
	    endTime=HRClock.now();
	    fprintf(f,"%.6f\t%i\t%i\t%i\n",p,Ng[0],Ng[1],iteration_number);
	    printf("%.6f\t%i\t%i\t%i\t%.3f\n",p,Ng[0],Ng[1],iteration_number,millisecond_res_diff(endTime,start));
	    p-=p_step;
	    start = HRClock.now();
	    printf("Starting new cascade with p=%.5f...\n",p);
	    continue_cascade(p);
	    //printf("Number of live vertices in A: %i\n", num_live_vertices(A));
	}

    }
    fprintf(f,"%.6f\t%i\t%i\t%i\n",p,Ng[0],Ng[1],iteration_number);
    printf("%.6f\t%i\t%i\t%i\n",p,Ng[0],Ng[1],iteration_number);
    fclose(f);
    if (create_images){
	start_tracked_cascade(p);
	stop=0;
	while(stop<2){
	    trim_to_gc();
	    kill_dependent();
	}
	write_history(0);
    }

}


void RNet::tracked_cascade_multiple_p(FloatVec p, int create_images) {
    //int i =0;
    for(double& p_ : p){
        //printf("%i: %.4f\n",++i,p_);
	    tracked_cascade_single_p(p_,  create_images);
    }
}

/**  
 * This takes a p_0 as a guess for p_c and then checks values around it until it finds the most precise value
 * recording the output for all of the intermediate steps along the way.
 *
 **/
double RNet::find_pc(double p0, double step,int create_images){
    auto startTime = HRClock.now();
    double p = p0;
    int create_all_images = true? 0 : 1;//TODO: enable this option 
    tracked_cascade_single_p(p,create_all_images);
    double pmax=MAX_ALLOWED_P,pmin=MIN_ALLOWED_P;
    if (Ng[0] == 0){
	while(Ng[0] == 0){
            if (p>=1){
                real_pc = 1;
                return 1;
            }
	    pmin = p;
            step*=2;
	    p+=step;
            if (p>=1)
            {
                p=1;
                break;
            }
	    tracked_cascade_single_p(p,create_all_images);
        if (errorCheck)
            printf("%.10f\t%i\n",p,Ng[0]);	    

	}
	pmax=p;
    }else{

	while(Ng[0]>0){
	    pmax = p;
	    p-=step;
            step*=2;
           if (p<=MIN_ALLOWED_P)
            {
                p=MIN_ALLOWED_P;
                break;
            }
	    tracked_cascade_single_p(p,create_all_images);
        if (errorCheck)
            printf("%.10f\t%i\n",p,Ng[0]);	    
	}
	pmin = p;
    }
    auto endTime = HRClock.now();
    profilingVector[PRE_BINSEARCH_TIME]+=millisecond_res_diff(endTime,startTime);
    printf("Established interval for p (%.9f,%.9f).\n",pmin,pmax);
    startTime=endTime;
    while((pmax - pmin)*N >= pc_precision){
	p = 0.5*(pmax + pmin);	
	tracked_cascade_single_p(p,create_all_images);
    if (errorCheck)
        printf("%.10f\t%i\n",p,Ng[0]);	    
	if(Ng[0] == 0)
	    pmin = p;
	else
	    pmax = p;
    }
    if(create_images>0)
        write_history(0);
    endTime=HRClock.now();
    profilingVector[BINSEARCH_TIME]+=millisecond_res_diff(endTime,startTime);
    real_pc = pmin;
    return pmin;
    
}


/** This is the hybrid version.  It has a problem that there is a slight lack of ergodicity--once you start 
 * following a given cluster as giant, you get different results than you would from a purely random attack.
 * thus the two parts of the hybrid check do not fit perfectly and the results analysis gets screwy.
 **/
/*
double RNet::find_pc(double p0, double step,int create_images,double zero_size){
    
    double p = p0;
    int create_all_images = true? 0 : 1;//TODO: enable this option 
    printf("Beginning cascade with p = %.3f...\n",p);
    start_tracked_cascade(p);
    double pmax=1,pmin=0.55;
    while(stop<2){
	int not_in_gc = trim_to_gc();
    //  printf("%i : Ng(trimmed) = %i (%i)\t",iteration_number, get_Ng(),not_in_gc);
	kill_dependent();
	if( stop>1 && Ng[0] > zero_size*N ){
	    pmax=p;
	    p-=step;
	    continue_cascade(p);
	    summary_map[p] =  run_summary(p,Ng[0],iteration_number);
	    //printf("Number of live vertices in A: %i\n", num_live_vertices(A));
	}
	pmin=p;

    }
    printf("%.6f\t%i\t%i\t%i\n",p,Ng[0],Ng[1],iteration_number);

    

    printf("Established interval for p (%.9f,%.9f).\n",pmin,pmax);
    
    while((pmax - pmin)*N >= 2){
	p = 0.5*(pmax - pmin)+pmin;	
	tracked_cascade_single_p(p,create_all_images);
	if(Ng[0]<zero_size*N)
	    pmin = p;
	else
	    pmax = p;
    }
    if(create_images>0)
        write_history(0);
    
    return pmin;
    
}
*/
void RNet::clear_summary()
{
    summary_map.clear();    
}


void RNet::write_summary(){
    int type = summary_map.begin()->second.m_type;
    string fname =boost::str(  boost::format("DependentLattices-type=%i-L=%i-lambda=%.4f-k=%.4f-q=%.4f-%lu.summary") % type % L % lambda % kavg % q  % time_stamp );
    ofstream file(fname.c_str());
    file << L << "\t" << q << "\t" << r << "\t" << "\n";
    file << summary_map.begin()->second.getHeader()<<"\n";
    for (auto it = summary_map.begin(); it != summary_map.end(); ++it){
	file << fixed << setprecision(10) << it->second << "\n";
    }
    file.close(); 

if (profileOn){
    
   string fname =boost::str(  boost::format("DependentLattices-type=%i-L=%i-lambda=%.4f-k=%.4f-q=%.4f-%lu.profile") % type % L % lambda % kavg % q  % time_stamp );
   ofstream file(fname.c_str());
   for (int i=0; i<TIME_SECTIONS; i++){
        file << profilingNames[i] << "\t" << profilingVector[i] << std::endl;
   }
   file.close();
}    
    
}

void RNet::write_graphviz(){
    ofstream f0("Graph0.dot");
    ofstream f1("Graph1.dot");
 boost::write_graphviz(f0,G[0]);   
 boost::write_graphviz(f1,G[1]);   
}

/*
void bfs_test(Graph& G){
   //IndexMap index_map = boost::get(boost::vertex_index_t(), G);
   int cluster_number=1;
   map<int,int> seen;
   map<int, vector<int> > clusters;
   cluster_visitor vis(&clusters, &seen, &cluster_number);
   printf("Using custom visitor\n");
   boost::breadth_first_search(G,boost::vertex(1,G),visitor(vis));
}
    MY attempt to do the whole bfs thing.  Ended up being orders of magnitude slower than the boost implementation so retired until further notice
    boost::graph_traits<Graph>::vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(G);
    //IndexMap index_map = boost::get(boost::vertex_index_t(), G);
    multimap<int,int> degree_ordering;
    int cluster_number=1;
    int total_seen=0;
    int max_cluster_size=0;
    int max_cluster_index=0;
    int vertex_index;
    
    //Record which nodes need to be checked and use the same map to record clusterID
    map<int,int> seen;
    map<int, vector<int> > clusters;
    
    for (next = vi; vi != vi_end; vi = next) {
	++next;
	vertex_index = *vi;
	if(alive[vertex_index] == STILL_ALIVE_CODE){
	    seen[vertex_index] = 0;
	    clusters[0].push_back(vertex_index);
	    degree_ordering.insert( pair<int,int>(boost::degree(boost::vertex(vertex_index,G),G),vertex_index));

	}
	else{
	    seen[vertex_index] = -1;//These are the dead nodes, don't bother checking them.
	    clusters[-1].push_back(vertex_index);
	}
    }

    int node_degree=0,node_number=0;
    cluster_visitor vis(&clusters, &seen, &cluster_number);
    endTime=time(0);
    printf("Took %i seconds to prepare for BFS.\nChecking clusters for %i nodes...\n",end-start,clusters[0].size());
    //Try the high degree nodes first, they are more likely to lead to the GC
    for( auto node_pair_it = degree_ordering.end(); node_pair_it != degree_ordering.begin(); --node_pair_it){
	node_number = (*node_pair_it).second;
	node_degree = (*node_pair_it).first;
	assert(node_degree != 0);
	if (node_degree<=1){//Since we are checking degree in descending order, if we hit degree 1 we are describing a pair
	    if(node_pair_it == degree_ordering.end()){//Which can only be the GC if it is first on the list
		continue;
	    }else{
		break;
	    }
	}
	if (seen[node_number]!=0)//if not worth checking -1 if previously dead, >0 if already discovered
	    continue;
	//printf("Beginning BFS #%i...\n",cluster_number);

	boost::breadth_first_search(G,boost::vertex(node_number,G),visitor(vis));
	auto last_size = clusters[cluster_number].size();
	total_seen+=last_size;
	//assert(last_size >1);
	if (last_size > max_cluster_size){
	    max_cluster_index = cluster_number;
	    max_cluster_size = last_size;
	}
	if (total_seen + max_cluster_size > Ng){
	    break;
	}
	
	++cluster_number;
    }
    for (auto it = clusters.begin(); it !=clusters.end(); ++it){
	printf("Cluster key %i points to vec of size %i.\n",it->first,it->second.size());
	if (it->second.size() == 2)
	    printf("%i : %i %i\n",it->first,it->second[0], it->second[1]);
    }
    Ng=max_cluster_size;
    printf("Cluster %i is largest and occupies %.3f (%i/%i) of the network.\n",max_cluster_index,1.0*max_cluster_size/N,max_cluster_size,N);
    for( auto node_pair_it = degree_ordering.end(); node_pair_it != degree_ordering.begin(); --node_pair_it){
	node_number = (*node_pair_it).second;
    if( seen[node_number] != max_cluster_index){
	    auto vertex_var = boost::vertex(node_number,G);
	    boost::clear_vertex(vertex_var,G);
	    alive[vertex_var] = iteration_number;
	}
    }

}

class cluster_visitor : public boost::default_bfs_visitor
{
public:
    IndexMap vertex_index; 
    map<int, vector<int> >* clusters_list;
    int* cluster_number;
    map<int,int>* seen;
    
    cluster_visitor(map< int, vector<int> > * clusters, map<int,int> * seen, int * clusterID):
    clusters_list(clusters),cluster_number(clusterID),seen(seen){}
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) 
    {
	(*clusters_list)[*cluster_number].push_back(u);
	(*seen)[u] = *cluster_number;
    }
};

//Old code for deciding which edges are real, bgl version is much cleaner and probably just as fast
for (int node = 0; node<N; ++node){
  if (alive[node][lattice_index] != STILL_ALIVE_CODE)
     continue;
  degree=0;
  if (node%L!=L-1){ //not on right edge
    if (alive[node+1][lattice_index] == STILL_ALIVE_CODE){
      add_triplet(node,node+1,-1);
      degree++;
    }
  }
  if (node%L!=0){ //not on left edge
      if (alive[node-1][lattice_index] == STILL_ALIVE_CODE){
      //add_triplet(node,node-1,-1);
      degree++;
    }
  }
  if (node>=L){ //not in first row
      if (alive[node-L][lattice_index] == STILL_ALIVE_CODE){
      //add_triplet(node,node-L,-1);
      degree++;
    }
  }
  if (node<(N-L)){ //not in last row
    if (alive[node+L][lattice_index] == STILL_ALIVE_CODE){
      degree++;
      add_triplet(node,node+L,-1);
      
    }
  }
  add_triplet(node,node,degree);
}

*/

