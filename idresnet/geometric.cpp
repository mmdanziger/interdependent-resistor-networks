#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "utilities.hpp"
// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

int main(void)
{
    typedef bg::model::point<double, 2, bg::cs::cartesian> point;
    typedef bg::model::box<point> box;
    typedef std::pair<point, unsigned> value;

    // create the rtree using default constructor
    bgi::rtree< value, bgi::quadratic<16> > rtree;
    unsigned N = 1000000;
    double r = 0.01;
    // create some values
    for ( unsigned i = 0 ; i < N ; ++i )
    {
        // insert new value
        rtree.insert(std::make_pair(point(randdouble(),randdouble()), i));
    }

    // find values intersecting some area defined by a box
    point root_point(randdouble(),randdouble()); 
    auto pxy=bg::model::d2::point_xy<double>();
    
    box query_box(point(0, 0), point(r, r));
    std::vector<value> result_s;
    rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
    
    // find 5 nearest values to a point
    std::vector<value> result_n;
    rtree.query(bgi::nearest(point(0.5, 0.5), 5), std::back_inserter(result_n));
    
    // display results
    std::cout << "spatial query box:" << std::endl;
    std::cout << bg::wkt<box>(query_box) << std::endl;
    std::cout << "spatial query result:" << std::endl;
    std::cout << "Found "<< result_s.size() << "points"<<std::endl;
    int found=0;
    point neighbor;
    while( !found)
    { 
	int ind = randint(0,result_s.size());
	auto dist = bg::comparable_distance(result_s[ind],root_point);
	if (dist<r){
	    found = 1;
	    neighbor=result_s[ind];
	}
        std::cout << bg::wkt<point>(v.first) << " - " << v.second << std::endl;
     }
    std::cout << "knn query point:" << std::endl;
    std::cout << bg::wkt<point>(point(0.5, 0.5)) << std::endl;
    std::cout << "knn query result:" << std::endl;
    BOOST_FOREACH(value const& v, result_n)
        std::cout << bg::wkt<point>(v.first) << " - " << v.second << std::endl;

    return 0;
}