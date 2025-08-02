#include <cmath>
#include <utility>

#define root3over2  0.8660254037844386
/*
 * type to hold vectors in triangular lattice basis, 
 * to be converted to cartesian coordinates as necessary
 */
struct TriangleVector2D;

typedef std::pair<double, TriangleVector2D> TriDistancePair; 


struct TriangleVector2D{
  long ux;
  long uy;
  double x;
  double y;
  TriangleVector2D(){}
  TriangleVector2D(long ux, long uy):ux(ux),uy(uy){update_cartesian();}
  double get_length_squared(){
    //return (ux + uy*0.5) * (ux + uy*0.5) + (root3over2*uy) * (root3over2*uy);
    return x*x + y*y;
  }
  void swap_coords(){
    std::swap(ux,uy);
    update_cartesian();
  }
  void rotate(int number_of_rotations){
    long oldux =ux;
    long olduy = uy;
    switch(number_of_rotations%6){
      case 0:
	return;
      case 1:
	ux = - olduy;
	uy = oldux + olduy;
	break;
      case 2: 
	ux = -olduy - oldux;
	uy = oldux;
	break;
      case 3:
	ux = -oldux;
	uy = -olduy;
	break;
      case 4:
	ux = olduy;
	uy = -oldux - olduy;
	break;
      case 5:
	ux = oldux + olduy;
	uy = -ux;
	break;
      default:
	break;
  }
  update_cartesian();
  }
  void update_cartesian(){
    x = ux + uy*0.5;
    y = root3over2 * uy;
  }

};