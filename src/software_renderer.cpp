#include "software_renderer.h"


#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "triangulation.h"
#include "variations.h"
#include "affine.h"

using namespace std;

namespace CMU462 {
  
vector<V*> init_V (FILE * pFile)
{
  vector<V*> out;

  V* v0 = new V0;
  v0->load(pFile);
  out.push_back(v0);

  V* v1 = new V1;
  v1->load(pFile);
  out.push_back(v1);

  v1 = new V2;
  v1->load(pFile);
  out.push_back(v1);
  
  v1 = new V3;
  v1->load(pFile);
  out.push_back(v1);
  
  v1 = new V4;
  v1->load(pFile);
  out.push_back(v1);
  
  v1 = new V5;
  v1->load(pFile);
  out.push_back(v1);

  v1 = new V6;
  v1->load(pFile);
  out.push_back(v1);
  
  v1 = new V7;
  v1->load(pFile);
  out.push_back(v1);

  v1 = new V8;
  v1->load(pFile);
  out.push_back(v1);
 
  v1 = new V9;
  v1->load(pFile);
  out.push_back(v1);

  v1 = new V10;
  v1->load(pFile);
  out.push_back(v1);


  v1 = new V11;
  v1->load(pFile);
  out.push_back(v1);


  v1 = new V12;
  v1->load(pFile);
  out.push_back(v1);


  v1 = new V13;
  v1->load(pFile);
  out.push_back(v1);


  v1 = new V14;
  v1->load(pFile);
  out.push_back(v1);

  v1 = new V15;
  v1->load(pFile);
  out.push_back(v1);


  return out;
}

struct temp {
  Vector2D v;
  Color c;
};

void SoftwareRendererImp::draw_svg( SVG& svg) {


  irradiance_buffer =  new float[4*target_w*target_h];
  for (size_t i = 0; i < 4*target_w*target_h; i++){
    irradiance_buffer[i] = 0.0;
  }

  this->max_irr = 0;
  FILE * pFile;
  pFile = fopen(svg.variations, "r+");
  vector<V*> variations;
  if(pFile == NULL)
  {
    cout << "Variation File not Found" << endl;
     V* v0 = new V0;
     v0->weight_to_one();
    variations.push_back(v0);
  }
  else
  {
    variations = init_V(pFile);
    fclose(pFile);
  }
  FILE * tFile = fopen(svg.file, "r+");
  if(tFile == NULL)
    cout << "AAAAAAAAA" << endl;
  int n;
  fscanf(tFile, "%i", &n);

  std::vector<Affine*> transform;
  std::vector<Affine*>::iterator iter;
  double a, b, c, d, e, f;
  for(int i = 0; i < n; i++)
  {
    double a, b, c, d, e, f, p;
    double c1, c2, c3;
    fscanf(tFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &p, &c1, &c2, &c3); 
    Affine* t = new Affine(variations, a, b, c, d, e, f, p, c1, c2, c3);
    fscanf(tFile, "%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f);
    t->init_pt(a, b, c, d, e, f);
    transform.push_back(t);
  }

  fscanf(tFile, "%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f);

  fclose(tFile);

/*  Color colors[n];

  for (int i = 0; i < n; i ++){
    colors[i] = Color::White;
  }

  colors[1] = Color(0.925, 0.0, 0.18), colors[2](0.992, 0.322, 0.0), colors[2](0.992, 0.322, 0.0)*/;

  Vector2D center = Vector2D(target_w/2, target_h/2);
  Vector2D p = Vector2D(0, 0);
  Color col = Color::White;
  int iterations = svg.iterations;
  double max_x = 0;
  double max_y = 0;
  double min_x = 0;
  double min_y = 0;

  vector<struct temp*> points;
  vector<struct temp*>::iterator point_iter;
 

  for ( size_t i = 0; i < iterations; ++i ) {
    double r =  rand() % 100;
    for(iter = transform.begin(); iter != transform.end(); iter++)
    {
      if(r <= (*iter)->probability * 100)
      {
        p = (*iter)->evaluate(p);
        col = (*iter)->evaluate_color(col);
        break;
      }
      else
        r = r - (*iter)->probability * 100;
    }

    /*
    if(p.x < min_x)
      min_x = p.x;
    if(p.x > max_x)
      max_x = p.x;
    if(p.y < min_y)
      min_y = p.y;
    if(p.y > max_y)
      max_y = p.y;

    struct temp *t1 = new struct temp;
    t1->v = p;
    t1->c = col;
    points.push_back(t1);
    */
    //Vector2D p_trans = p * 500;
    //p_trans += center;
    Vector2D p_trans = Vector2D(a*p.x + b*p.y + c, d*p.x + e*p.y + f);
    //cout << p_trans.x << "   " << p_trans.y << endl;
    rasterize_point(p_trans.x, p_trans.y, col);
  }

  /*
  double mid_x = (min_x + max_x)/2;
  double mid_y = (min_y + max_y)/2;

  for(point_iter = points.begin(); point_iter != points.end(); point_iter++)
  {
    Vector2D p_trans = (*point_iter)->v * target_h * 3.0/4.0;
    p_trans += center - 3.0/4.0*target_h*Vector2D((max_x - min_x)/2, (max_y - min_y)/2);
    rasterize_point(p_trans.x, p_trans.y, (*point_iter)->c);
  }
*/
  //Deallocate everything
  for(int i = 0; i < transform.size(); i++)
  {
    delete transform[i];
  }

  for(int i = 0; i < variations.size(); i++)
  {
    delete variations[i];
  }

  
  // resolve and send to render target
  resolve();

}



void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  this->sample_rate = sample_rate;

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 3: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;

  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 4 (part 1):
  // Modify this to implement the transformation stack
  Matrix3x3 old = transformation;
  transformation =  transformation * element->transform;
  

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }
  transformation = old;
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {
    

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //
// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates


void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

// fill in the nearest pixel
  int s_x = (int)floor(x);
  int s_y = (int)floor(y);

  if ( s_x < 1 || s_x >= target_w - 1 ) return;
  if ( s_y < 1 || s_y >= target_h - 1 ) return;

  float center_x = (float)s_x + 0.5;
  float center_y = (float)s_y + 0.5;

  int u = (x < center_x) ? s_x - 1 : s_x + 1;
  int v = (y < center_y) ? s_y - 1 : s_y + 1;

  float dx = abs(x - center_x);
  float dy = abs(y - center_y);



  irradiance_buffer[4*(s_x + s_y * target_w)] += (1 - dx)*(1 - dy)*color.r;
  irradiance_buffer[4*(s_x + v * target_w)] += (1 - dx)*(dy)*color.r;
  irradiance_buffer[4*(u + s_y * target_w)] += (dx)*(1 - dy)*color.r;
  irradiance_buffer[4*(u + v * target_w)] += (dx)*(dy)*color.r;


  irradiance_buffer[4*(s_x + s_y * target_w) + 1] += (1 - dx)*(1 - dy)*color.g;
  irradiance_buffer[4*(s_x + v * target_w) + 1] += (1 - dx)*(dy)*color.g;
  irradiance_buffer[4*(u + s_y * target_w) + 1] += (dx)*(1 - dy)*color.g;
  irradiance_buffer[4*(u + v * target_w) + 1] += (dx)*(dy)*color.g;


  irradiance_buffer[4*(s_x + s_y * target_w) + 2] += (1 - dx)*(1 - dy)*color.b;
  irradiance_buffer[4*(s_x + v * target_w) + 2] += (1 - dx)*(dy)*color.b;
  irradiance_buffer[4*(u + s_y * target_w) + 2] += (dx)*(1 - dy)*color.b;
  irradiance_buffer[4*(u + v * target_w) + 2] += (dx)*(dy)*color.b;


  irradiance_buffer[4*(s_x + s_y * target_w) + 3] += (1 - dx)*(1 - dy);
  irradiance_buffer[4*(s_x + v * target_w) + 3] += (1 - dx)*(dy);
  irradiance_buffer[4*(u + s_y * target_w) + 3] += (dx)*(1 - dy);
  irradiance_buffer[4*(u + v * target_w) + 3] += (dx)*(dy);


  max_irr = max(max_irr, max(irradiance_buffer[4*(s_x + s_y * target_w) + 3], max(irradiance_buffer[4*(s_x + v * target_w) + 3], 
                max(irradiance_buffer[4*(u + s_y * target_w) + 3], irradiance_buffer[4*(u + v * target_w) + 3]))));



/*  // check bounds
  if ( s_x < 0 || s_x >= target_w ) return;
  if ( s_y < 0 || s_y >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (s_x + s_y * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (s_x + s_y * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (s_x + s_y * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (s_x + s_y * target_w) + 3] = (uint8_t) (color.a * 255);*/
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
}

int sign(float x)
{
    if(x < 0)
        return -1;
    else if (x > 0)
        return 1;
    else
        return 0;
}


void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {
  // Task 3: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 3".
  int distance_buffer[target_w * target_h];

  for (size_t i = 0; i < target_w; i++){
    for (size_t j = 0; j < target_h; j++){
      float irradiance = irradiance_buffer[i + j * target_w];
      if (irradiance != 0.0){
        distance_buffer[i + j * target_w] = 0;
      }
      else{
        distance_buffer[i + j * target_w] = target_w * target_h;
      }
    }
  }
/*
  bool changed = true;
  int max_dist = 0;
  while (changed){
    changed = false;
    for (size_t i = 0; i < target_w; i++){
      for (size_t j = 0; j < target_h; j++){
        int update_dist = distance_buffer[i + j * target_w] + 1;
        if (i > 0 && distance_buffer[i - 1 + j * target_w] > update_dist){
          distance_buffer[i - 1 + j * target_w] = update_dist;
          changed = true; max_dist = max(max_dist, update_dist);
        }
        if (i < target_w - 1 && distance_buffer[i + 1 + j * target_w] > update_dist){
          distance_buffer[i + 1 + j * target_w] = update_dist;
          changed = true; max_dist = max(max_dist, update_dist);
        }
        if (j > 0 && distance_buffer[i + (j - 1) * target_w] > update_dist){
          distance_buffer[i + (j - 1) * target_w] = update_dist;
          changed = true; max_dist = max(max_dist, update_dist);
        }
        if (j < target_h - 1 && distance_buffer[i + (j + 1) * target_w] > update_dist){
          distance_buffer[i + (j + 1) * target_w] = update_dist;
          changed = true; max_dist = max(max_dist, update_dist);
        }
      }
    }
  }

*/
  Color color;
  Color color1;
  color1.r = 0;
  color1.g = 0;
  color1.b = 0;
  color1.a = 1;
  for (size_t i = 0; i < target_w; i++){
    for (size_t j = 0; j < target_h; j++){
      float irr = irradiance_buffer[4*(i + j * target_w) + 3];
      float r = irradiance_buffer[4*(i + j * target_w)];
      float g = irradiance_buffer[4*(i + j * target_w) + 1];
      float b = irradiance_buffer[4*(i + j * target_w) + 2];
      //float dist_f = (float)distance_buffer[i + j * target_w] * 1.0/(float)max_dist;
      //color1.b = (1 - dist_f);
      float alpha = (log2(1 + 1023*irr/ (max_irr)))/10.0;
      color.r = (log2(1 + 1023*r/ (max_irr)))/10.0 * alpha + (1 - alpha) * color1.r;
      color.g = (log2(1 + 1023*g/ (max_irr)))/10.0 * alpha + (1 - alpha) * color1.g;
      color.b = (log2(1 + 1023*b/ (max_irr)))/10.0 * alpha + (1 - alpha) * color1.b;
      color.a = 1;

      render_target[4 * (i + j * target_w)    ] = (uint8_t) (color.r * 255);
      render_target[4 * (i + j * target_w) + 1] = (uint8_t) (color.g * 255);
      render_target[4 * (i + j * target_w) + 2] = (uint8_t) (color.b * 255);
      render_target[4 * (i + j * target_w) + 3] = (uint8_t) (color.a * 255);
    }
  }
  delete [] this->irradiance_buffer;
  return;
}



} // namespace CMU462
