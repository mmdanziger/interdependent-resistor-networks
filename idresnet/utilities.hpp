#include <chrono>
#include <vector>
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdint>
#include <sys/types.h>
#include <unistd.h>


static const  int STILL_ALIVE_CODE = 0;
static const  int DEAD_CODE = 1;
using namespace std;
boost::mt19937 gen(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
double millisecond_res_diff(chrono::high_resolution_clock::time_point end, chrono::high_resolution_clock::time_point start) {
    return( chrono::duration_cast<chrono::milliseconds>(end-start).count() / 1000.0);
}

struct vertex_alive_t {
    typedef boost::vertex_property_tag kind;
};
//



typedef vector<double> FloatVec;

template <typename AliveMap>
struct vertex_is_alive {
    vertex_is_alive() { }
    vertex_is_alive(AliveMap alive, int lattice_index) : m_alive(alive),m_lattice_index(lattice_index) { }
    template <typename Vertex>
    bool operator()(const Vertex& v) const {
        return (boost::get(m_alive,v))[m_lattice_index] == STILL_ALIVE_CODE;
    }
    AliveMap m_alive;
    int m_lattice_index;
};


struct vertex_in_cluster {
    vertex_in_cluster() { }
    vertex_in_cluster(vector<int> cluster_map, int ind ) : m_cluster_map(cluster_map) , m_cluster_index(ind) { }
    template <typename Vertex>
    bool operator()(const Vertex& v) const {
        return m_cluster_map[v] == m_cluster_index;
    }
    vector<int> m_cluster_map;
    int m_cluster_index;
};


template <typename EdgeWeightMap>
struct edge_in_graph {
    edge_in_graph() { }
    edge_in_graph(EdgeWeightMap weight, int lattice_index) :m_lattice_index(lattice_index),m_weight(weight) { }
    template <typename Edge>
    bool operator()(const Edge& e) const {
        return (boost::get(m_weight, e)).test(m_lattice_index);
    }
    EdgeWeightMap m_weight;
    int m_lattice_index;
};


// This visitor is used both in the connected_components algorithm
// and in the kosaraju strong components algorithm during the
// second DFS traversal.
template <class ComponentsMap>
class components_recorder : public boost::default_bfs_visitor
{
    typedef typename boost::property_traits<ComponentsMap>::value_type comp_type;
public:
    components_recorder(ComponentsMap c,
                        comp_type& c_count)
        : m_component(c), m_count(c_count) {}

    /*      template <class Vertex, class Graph>
          void start_vertex(Vertex, Graph&) {
            if (m_count == (std::numeric_limits<comp_type>::max)())
              m_count = 0; // start counting components at zero
            else
              ++m_count;
          }*/
    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, Graph&) {

        boost::put(m_component, u, m_count);
    }
protected:
    ComponentsMap m_component;
    comp_type& m_count;

};



inline int distance_squared(int i,int j,int L) {
    int dx = (i/L - j/L);
    int dy = (i%L - j%L);
    if (dx > L/2)
        dx=L-dx;
    else if (dx < -L/2)
        dx=L+dx;
    if (dy > L/2)
        dy=L-dy;
    else if (dy < -L/2)
        dy = L +dy;
    return (dx*dx + dy*dy);
}


int randint(const int N) {
    if( N == 0)
        return 0;
    boost::uniform_int<> dist(0, N-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniform_rand(gen, dist);
    return uniform_rand();
}
float randfloat(float min=0.0, float max=1.0)
{
    boost::uniform_real<float> u(min, max);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > uniform_float(gen, u);
    return uniform_float();
}

double randdouble(double min=0.0, double max=1.0)
{
    boost::uniform_real<double> u(min, max);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > uniform_double(gen, u);
    return uniform_double();
}

int randsign() {
     return rand() > (RAND_MAX/2 +1)? 1 : -1;   
}
int* id_to_coord(int node_id, const int L) {
    int*  coordinates = new int[2];
    coordinates[0]=node_id/L;
    coordinates[1]=node_id%L;
    return(coordinates);
}

