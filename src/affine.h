#include "CMU462.h"

namespace CMU462 {

class Affine {
  double t [6]; // transform 
  double pt [6]; // post-transform
  double color [3];
  std::vector<V*> variations;
  public:
  Affine(std::vector<V*> v, double a, double b, double c, double d, double e, double f, double prob, double c1, double c2, double c3)
  {
    variations = v;
    probability = prob;
    t[0] = a;
    t[1] = b;
    t[2] = c;
    t[3] = d;
    t[4] = e;
    t[5] = f;
    pt[0] = 1;
    pt[1] = 0;
    pt[2] = 0;
    pt[3] = 0;
    pt[4] = 1;
    pt[5] = 0;

    color[0] = c1;
    color[1] = c2;
    color[2] = c3;

  }
  void init_pt(double a, double b, double c, double d, double e, double f)
  {
    pt[0] = a;
    pt[1] = b;
    pt[2] = c;
    pt[3] = d;
    pt[4] = e;
    pt[5] = f;
  }
  ~Affine()
  {
    //std::cout << "deleted " << std::endl;
  }
  double probability;
  Vector2D apply_variations(Vector2D p)
  {
    Vector2D out = Vector2D();
    std::vector<V*>::iterator iter;
    for(iter = variations.begin(); iter != variations.end(); iter++)
    {
        out += (*iter)->eval(p);
    }
    return out;
  }
  Vector2D evaluate(Vector2D v){
    Vector2D v1 = Vector2D(t[0]*v.x + t[1]*v.y + t[2], t[3]*v.x + t[4]*v.y + t[5]);
    Vector2D v2 = apply_variations(v1);
    Vector2D v3 = Vector2D(pt[0]*v2.x + pt[1]*v2.y + pt[2], pt[3]*v2.x + pt[4]*v2.y + pt[5]);
    return v3;
  }

  Color evaluate_color(Color in)
  {
    return Color((in.r + color[0])/2, (in.g + color[1])/2, (in.b + color[2])/2, 1);
  }
  
};
}
