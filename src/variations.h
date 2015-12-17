#include "CMU462.h"
#include <stdio.h>

namespace CMU462 {
class V {
  public:
  double weight;

  Vector2D eval(Vector2D p) {
    if(weight == 0)
      return Vector2D();
    else
      return weight * v(p);
  }

  void load(FILE * pFile)
  {
    fscanf(pFile, "%lf", &weight);
    if(weight != 0)
      l(pFile);
  }

  void weight_to_one()
  {
    weight = 1;
  }


  virtual void l(FILE * pFile) =0;
  virtual Vector2D v(Vector2D p) =0;
};

class V0 : public V {
  // Linear
  public:
  //int x, y;
  void l(FILE * pFile){
    //fscanf(pFile, "%i %i",  &x, &y);
    return;
  }
  Vector2D v(Vector2D p){
    return p;
  }
};

class V1 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V2 : public V {
  //Spherical
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    double r2 = p.x*p.x + p.y*p.y;
    if(r2 == 0)
      return p;
    Vector2D out = 1.0/r2 * p;
    return out;
  }
};
class V3 : public V {
  //Swirl
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    double r2 = p.norm2();
    return Vector2D(p.x * sin(r2) - p.y * cos(r2), p.x * cos(r2) + p.y * sin(r2));
  }
};
class V4 : public V {
  //Horseshoe
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    double r = p.norm();
    if(r == 0)
      return p;
    return ((1/r) * Vector2D((p.x - p.y) * (p.x + p.y), 2 * p.x * p.y));
  }
};
class V5 : public V {
  //Polar
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    if(p.y == 0)
      return p;
    double theta = atan(p.x / p.y);
    double r = p.norm();
    double PI = 3.14159;
    return Vector2D(theta/PI, r - 1);
  }
};

class V6 : public V {
  //Handkerchief
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    double r = p.norm();
    if(p.y == 0)
      return p;
    double theta = atan(p.x / p.y);
    return r * Vector2D(sin(theta + r), cos(theta - r));
  }
};
class V7 : public V {
  //Heart
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    double r = p.norm();
    if(p.y == 0)
    {
      return p;
    }
    double theta = atan(p.x / p.y);
    Vector2D out = r * Vector2D(sin(r * theta), -1 * cos(r * theta));
    return out;
  }
};

class V8 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V9 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V10 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V11 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V12 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V13 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V14 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};
class V15 : public V {
  //Sinusoidal
  public:
  void l(FILE * pFile) {
    return;
  }
  Vector2D v(Vector2D p){
    return Vector2D(sin(p.x), sin(p.y));
  }
};



}
