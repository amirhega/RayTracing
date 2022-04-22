/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Amir>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
using namespace std;
#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define MAXT 1e30

unsigned char buffer[HEIGHT][WIDTH][3];
double aspectRatio, tanConst;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

// Customized 3D vector
struct Vec
{
	double x;
	double y;
	double z;

	// default constructor
  Vec() { x = y = z = 0.0; }

	// customized constructor
  Vec(double x, double y, double z)
  {
      this->x = x;
      this->y = y;
      this->z = z;
  }

  // Add two vectors
  Vec add(Vec vec)
  {
      return Vec(this->x + vec.x, this->y + vec.y, this->z + vec.z);
  }

  // Scale this vector
	Vec mult(double scale)
  {
   return Vec(scale * this->x, scale * this->y, scale * this->z);
  }

  // Negate this vector
  Vec neg() { return this->mult(-1.0f); }

  // Subtract a vector from this vector
  Vec minus(Vec vec)
  {
    return this->add(vec.neg());
  }

  // Get the length of this vector
  double magnitude(){
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
  }

    // Normalize this vector
	Vec normalize()
	{
    return this->mult(1.0f / this->magnitude());
	}

	// Dot product
	double dot(Vec v)
	{
		return this->x * v.x + this->y * v.y + this->z * v.z;
	}

	// Cross product
	Vec crossP(Vec p)
	{
		double dx = this->y * p.z - this->z * p.y;
		double dy = - (this->x * p.z - this->z * p.x);
		double dz = this->x * p.y - this->y * p.x;

		return Vec(dx, dy, dz);
	}
  Vec clamp() {
    if(this->x > 1.0) this->x = 1.0;
    if(this->y > 1.0) this->y = 1.0;
    if(this->z > 1.0) this->z = 1.0;
    return Vec(this->x,this->y,this->z); 
  }
};

struct Ray {
  Vec origin, direction;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
bool isRightPixel = false;
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double minimum(double v1, double v2) {
  if(v1 < v2) return v1;
  return v2;
}

Vec triangleNormal(Triangle t) {
  Vec v0,v1,v2;
  v0.x = t.v[0].position[0];
  v0.y = t.v[0].position[1];
  v0.z = t.v[0].position[2];

  v1.x = t.v[1].position[0];
  v1.y = t.v[1].position[1];
  v1.z = t.v[1].position[2];

  v2.x = t.v[2].position[0];
  v2.y = t.v[2].position[1];
  v2.z = t.v[2].position[2];

  Vec t1 = v1.minus(v0);
  Vec t2 = v2.minus(v0);

  Vec tNorm = t1.crossP(t2);
  tNorm = tNorm.normalize();
  return tNorm;

}

enum Plane {XY, YZ, ZX};
enum Hit {NONE, SPHERE, TRIANGLE};

void barycentricCoord(const Vec& p, const Triangle t, double *alpha, double *beta, double *gamma) {
  Vec triangleN = triangleNormal(t);
  Vec NYZ,NZX,NXY;
  NYZ.x = 1.0; NYZ.y=0.0; NYZ.z =0.0;
  NZX.x = 0.0; NZX.y=1.0; NZX.z =0.0;
  NXY.x = 0.0; NXY.y=0.0; NXY.z =1.0;

  Plane plane;
  double xyPlane, yzPlane, zxPlane;
  xyPlane = fabs(triangleN.dot(NXY));
  yzPlane = fabs(triangleN.dot(NYZ));
  zxPlane = fabs(triangleN.dot(NZX));

  if(xyPlane > yzPlane) {
    if(xyPlane > zxPlane) plane = XY;
    else plane = ZX;
  } else {
    if(yzPlane > zxPlane) plane = YZ;
    else plane = ZX;
  }
  
  double ax,ay,bx,by,cx,cy, px,py;
  if(plane == XY) {
    ax = t.v[0].position[0];
    ay = t.v[0].position[1];

    bx = t.v[1].position[0];
    by = t.v[1].position[1];

    cx = t.v[2].position[0];
    cy = t.v[2].position[1];

    px = p.x;
    py = p.y;
  } else if(plane == ZX) {
    ax = t.v[0].position[2];
    ay = t.v[0].position[0];

    bx = t.v[1].position[2];
    by = t.v[1].position[0];

    cx = t.v[2].position[2];
    cy = t.v[2].position[0];

    px = p.z;
    py = p.x;
  } else if(plane == YZ) {
    ax = t.v[0].position[1];
    ay = t.v[0].position[2];

    bx = t.v[1].position[1];
    by = t.v[1].position[2];

    cx = t.v[2].position[1];
    cy = t.v[2].position[2];
    
    px = p.y;
    py = p.z;
  }

  double areaAllV = (0.5 * ((bx - ax)*(cy-ay) - (cx-ax)*(by-ay)));
  double areaCC1C2 = (0.5 * ((bx - px)*(cy-py) - (cx-px)*(by-py)));
  double areaC0CC2 = (0.5 * ((px - ax)*(cy-ay) - (cx-ax)*(py-ay)));
  double areaC0C1C = (0.5 * ((bx - ax)*(py-ay) - (px-ax)*(by-ay)));
  *alpha = areaCC1C2 / areaAllV;
  *beta = areaC0CC2 / areaAllV;
  *gamma = areaC0C1C / areaAllV;
}


bool sphereIntersection(Vec v, double& minT, int& minIdx) {
  double xc,yc,zc,r,b,c, min;

  bool found = false;
  for(int i = 0; i < num_spheres; i++) {
    double t0,t1;
    double eps = 1e-8;
    //get sphere info
    xc = spheres[i].position[0];
    yc = spheres[i].position[1];
    zc = spheres[i].position[2];
    r = spheres[i].radius;

    // cout << "v.x: " << v.x << endl;
    // cout << "v.y: " << v.y << endl;
    // cout << "v.z: " << v.z << endl;
    //  cout << "v.x: " << xc << endl;
    // cout << "v.y: " << yc << endl;
    // cout << "v.z: " << zc << endl;
    b = -2.0 * (v.x * xc+v.y*yc+v.z*zc);
    c =  pow(xc,2)+pow(yc,2)+pow(zc,2)-pow(r,2);
    double disc = pow(b,2) - (4*c);
    // cout << "b " << b << endl;
    // cout << "c " << c << endl;

    //if root is greater than 0, else abort
    if(disc >= 0) {
      t0 = (-b + sqrt(disc)) /2.0;
      t1 = (-b - sqrt(disc)) /2.0;

      if((t0 > eps && t1 > eps)) min = minimum(t0,t1);
      if(t0 > eps) min = t0;
      if(t1 > eps) min = t1;
      if(t0 <=eps && t1 <=eps) continue;

      found = true;

      if(min < minT) {
        minT = min;
        minIdx = i;
      }
    } else continue;
  }
  return found;

}

bool triangleIntersection( Vec v, double &min_t, int &minIdx, Hit &hitSphere, int ind) {
  bool found = false;
  for(int i = 0; i < num_triangles; i++) {
    if(isRightPixel) cout << "CHECK 3: " << hitSphere << endl;
    Triangle t = triangles[i];
    Vec planeNormal = triangleNormal(t);
    double a,b,c, d, denom;
    a = planeNormal.x;
    b = planeNormal.y;
    c = planeNormal.z;

    d = -1.0 * (a*t.v[0].position[0]+b*t.v[0].position[1]+c*t.v[0].position[2]);
    denom = ((a*v.x) +(b*v.y)+(c*v.z));

    if(denom == 0) continue;

    double t_value = -1 * (a+b+c+d)/denom;
    if(t_value <= 0) continue;
    
    //check if i have to multiply t
    if(t_value < min_t) {
      // Vec v2(v.x * t_value,v.y * t_value,v.z * t_value);
      Vec v1 = v.mult(t_value);
      double alpha, beta, gamma;
      barycentricCoord(v1, t, &alpha,&beta,&gamma);
      if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) continue;
      found = true;
      min_t = t_value;
      minIdx = i;
      hitSphere = TRIANGLE;
    }
    
  }
  return found;
}

bool triangleIntersection2(const Ray v, int i, double& tGot) {
  bool found = false;
  // for(int i = 0; i < num_triangles; i++) {
    // if(isRightPixel) cout << "CHECK 3: " << hitSphere << endl;
    Triangle t = triangles[i];
    Vec planeNormal = triangleNormal(t);
    double a,b,c, d, denom;
    a = planeNormal.x;
    b = planeNormal.y;
    c = planeNormal.z;

    d = -1.0 * (a*t.v[0].position[0]+b*t.v[0].position[1]+c*t.v[0].position[2]);
    denom = ((a*v.direction.x) +(b*v.direction.y)+(c*v.direction.z));

    if(denom == 0) return false;

    double t_value = -1 * (a*v.origin.x+b*v.origin.y+c*v.origin.z+d)/denom;
    if(t_value <= 0) return false;
    
    //check if i have to multiply t
   
      Vec v2(v.direction.x * t_value,v.direction.y * t_value,v.direction.z * t_value);
      v2 = v2.add(v.origin);
      double alpha, beta, gamma;
      barycentricCoord(v2, t, &alpha,&beta,&gamma);
      if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) return false;
      found = true;
      tGot = t_value;
    
  // }
  return found;
}

bool checkIntersection(const Vec v, double &min_t, int &minIdx, Hit &hitSphere) {
  // if(isRightPixel) cout << "CHECK: " << hitSphere << endl;
  bool intersect = false, intersect2 = false; 
  if(sphereIntersection(v,min_t, minIdx)) {
    hitSphere = SPHERE;
    intersect = true;
  }
  if(triangleIntersection(v, min_t, minIdx, hitSphere, 0)) {
    intersect2 = true;
  }

  return intersect2 || intersect;
}

bool inShadow(Ray shadowRay, double furthestT) {
  double xc,yc,zc,r,b,c, min;
  double eps = 1e-4, epsilon = 1e-4;  

  bool found = false;
  //in sphere intersection
  for(int i = 0; i < num_spheres; i++) {
    double t0,t1;
    //get sphere info
    xc = spheres[i].position[0];
    yc = spheres[i].position[1];
    zc = spheres[i].position[2];
    r = spheres[i].radius;

    b = 2.0 * (shadowRay.direction.x * (shadowRay.origin.x-xc) + shadowRay.direction.y*(shadowRay.origin.y-yc) 
                + shadowRay.direction.z*(shadowRay.origin.z-zc));
    c =  pow(shadowRay.origin.x-xc,2)+pow(shadowRay.origin.y-yc,2)+pow(shadowRay.origin.z-zc,2)-pow(r,2);
    double disc = pow(b,2) - (4*c);

    //if root is greater than 0, else abort
    if(disc >= 0) {
      t0 = (-b + sqrt(disc)) /2.0;
      t1 = (-b - sqrt(disc)) /2.0;
      Vec temp1 = shadowRay.direction;
      Vec temp2 = shadowRay.direction;

      //greater than valid point (nearly 0) and not behind the light
      if(t0 > eps && temp1.mult(t0).magnitude() < furthestT) return true;
      if(t1 > eps && temp2.mult(t1).magnitude() < furthestT) return true;
    } else continue;
  }
  
  double t;
  for(int i = 0; i < num_triangles; i++) {

    Vec n = triangleNormal(triangles[i]);
       	double nd = n.dot(shadowRay.direction);

       	// ray is parallel to the plane
       	if (-epsilon < nd && nd < epsilon) continue;

       	// Arbitrary point on the triangle in order to calculate d
        // Vec p0 = get_vertex_position(triangle_idx, 0);
        Triangle tr = triangles[i];
        Vec p0(tr.v[0].position[0], tr.v[0].position[1], tr.v[0].position[2]);
        double d = n.dot(p0);

       	double t = (d - n.dot(shadowRay.origin)) / nd;

       	// ray intersects a plane behind the scene
       	if (t <= epsilon) continue;
        // shadowRay.direction = shadowRay.direction.mult(t);

    if(triangleIntersection2(shadowRay,i,t) && t > eps && t <furthestT) { 
    if(isRightPixel)cout << "brett sucks" << endl;
    return true; 
    }

    // Triangle t = triangles[i];
    // Vec planeNormal = triangleNormal(t);
    // double ndotd = planeNormal.dot(shadowRay.direction);
    // //  cout << "barycentric2" << endl;

    // if(ndotd > -eps && ndotd < eps) continue;
    // Vec temp(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
    // double d = planeNormal.dot(temp);
    // double t_value = (d-planeNormal.dot(shadowRay.origin))/ ndotd;



    // double a,b,c, d, denom;
    // a = planeNormal.x;
    // b = planeNormal.y;
    // c = planeNormal.z;

    // d = -1.0 * (a*t.v[0].position[0]+b*t.v[0].position[1]+c*t.v[0].position[2]);
    // denom = ((a*shadowRay.direction.x) +(b*shadowRay.direction.y)+(c*shadowRay.direction.z));

    // if(denom == 0) continue;

    // double t_value = -1 * ((a*shadowRay.origin.x)+(b*shadowRay.origin.y)+(c*shadowRay.origin.z)+d)/denom;
    // if(t_value <= eps) continue;
   
    // //check if i have to multiply t
    //   shadowRay.direction = shadowRay.direction.mult(t_value);

    // if(temp1.mult(t_value).magnitude() < furthestT) {
    //   // double alpha, beta, gamma;
    //   // barycentricCoord(shadowRay.direction.mult(t_value).add(shadowRay.origin), t, &alpha,&beta,&gamma);
    //   // if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) continue;
    //   if(t_value > eps && shadowRay.di.mult(t_value).magnitude() < furthestT) return true;
    // }    
  }
  return false;

}

Vec phongLight(Vec primRay, const double minT, int minIdx, Hit hitSphere) {

  Vec intersectionP = primRay.mult(minT);
  Vec phong;
  // cout << minT << endl;
  
  //hits something other than the background
  Ray shadowRay;
  for(int i = 0; i < num_lights; i++) {
    
    shadowRay.origin = primRay.mult(minT);
    Vec lightPoint;
    lightPoint.x = lights[i].position[0];
    lightPoint.y = lights[i].position[1];
    lightPoint.z = lights[i].position[2];

    shadowRay.direction = lightPoint.minus(shadowRay.origin).normalize();
    double furthestPoint = lightPoint.minus(shadowRay.origin).magnitude();
    if(isRightPixel) cout << "Before shadow HERE: " << hitSphere << endl;
    if(!inShadow(shadowRay, furthestPoint)) {
      if(isRightPixel) cout << "GOT HERE: " << hitSphere << endl;
      if(hitSphere == SPHERE) {
        double xc,yc,zc,r;
        Sphere sphere = spheres[minIdx];
        xc = sphere.position[0];
        yc = sphere.position[1];
        zc = sphere.position[2];
        r = sphere.radius;
        Vec sphereVec;
        sphereVec.x = xc;
        sphereVec.y = yc;
        sphereVec.z = zc;

        //do i normalize
        Vec normal = shadowRay.origin.minus(sphereVec).mult(1/r).normalize();
        Vec L = lightPoint.minus(intersectionP).normalize();
        Vec reflect = (normal.mult(2*L.dot(normal)).minus(L)).normalize();

        double LdotN = L.dot(normal);
        if(LdotN <= 1e-8) LdotN = 0;
        // cout << "LdotN: " << LdotN << endl;

        Vec diffuseColor;
        diffuseColor.x = sphere.color_diffuse[0];
        diffuseColor.y = sphere.color_diffuse[1];
        diffuseColor.z = sphere.color_diffuse[2];
        diffuseColor = diffuseColor.mult(LdotN);
        // cout << "diffuse color: " << diffuseColor.x<<endl;
        // cout << "diffuse color: " << diffuseColor.y<<endl;
        // cout << "diffuse color: " << diffuseColor.z<<endl;

        double RdotV = reflect.dot(primRay.mult(-1.0));
        if(RdotV <= 1e-8) RdotV = 0.0;
        // cout << "RdotV: " << RdotV << endl;

        Vec specColor;
        specColor.x = sphere.color_specular[0];
        specColor.y = sphere.color_specular[1];
        specColor.z = sphere.color_specular[2];
        specColor = specColor.mult(pow(RdotV, sphere.shininess));
        // cout << "specColor color: " << specColor.x<<endl;
        // cout << "specColor color: " << specColor.y<<endl;
        // cout << "specColor color: " << specColor.z<<endl;

        Vec color;
        color.x = lights[i].color[0];
        color.y = lights[i].color[1];
        color.z = lights[i].color[2];

        color.x = color.x * ((diffuseColor.x) +specColor.x);
        color.y = color.y * ((diffuseColor.y) +specColor.y);
        color.z = color.z * ((diffuseColor.z) +specColor.z);

        // cout << "x color: " << color.x<<endl;
        // cout << "y color: " << color.y<<endl;
        // cout << "z color: " << color.z<<endl;
        phong = phong.add(color);
        
      } else if(hitSphere == TRIANGLE) {
        // cout << "got here" <<endl;
        double alpha, beta, gamma;
        Triangle triangle = triangles[minIdx];
        barycentricCoord(intersectionP, triangle, &alpha, &beta, &gamma);
         if(isRightPixel) {
        cout << "aplha " << alpha << " " << beta << " " << gamma<<endl;
        }
        Vec triangleNorm = triangleNormal(triangle);
        Vec normal0(triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2]);
        normal0 = normal0.mult(alpha);
        Vec normal1(triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2]);
        normal1 = normal1.mult(beta);
        Vec normal2(triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2]);
        normal2 = normal2.mult(gamma);

        //do i normalize
        Vec normal = normal1.add(normal0).add(normal2).normalize();
        Vec L = lightPoint.minus(intersectionP).normalize();
        Vec reflect = (normal.mult(2*L.dot(normal)).minus(L)).normalize();

        double LdotN = L.dot(normal);
        if(LdotN <= 1e-8) LdotN = 0;
        if(LdotN > 1) LdotN = 1.0;

        Vec diffuseColor;
        diffuseColor.x = alpha*triangle.v[0].color_diffuse[0] + beta*triangle.v[0].color_diffuse[0] + gamma*triangle.v[0].color_diffuse[0];
        diffuseColor.y = alpha*triangle.v[1].color_diffuse[1] + beta*triangle.v[1].color_diffuse[1] + gamma*triangle.v[1].color_diffuse[1];
        diffuseColor.z = alpha*triangle.v[2].color_diffuse[2] + beta*triangle.v[2].color_diffuse[2] + gamma*triangle.v[2].color_diffuse[2];
        diffuseColor = diffuseColor.mult(LdotN);

        double RdotV = reflect.dot(primRay.mult(-1.0));
        if(RdotV <= 1e-8) RdotV = 0.0; //clamp 
        if(RdotV > 1) RdotV = 1.0;

        Vec specColor;
        specColor.x = alpha*triangle.v[0].color_specular[0] + beta*triangle.v[0].color_specular[0] + gamma*triangle.v[0].color_specular[0];
        specColor.y = alpha*triangle.v[1].color_specular[1] + beta*triangle.v[1].color_specular[1] + gamma*triangle.v[1].color_specular[1];
        specColor.z = alpha*triangle.v[2].color_specular[2] + beta*triangle.v[2].color_specular[2] + gamma*triangle.v[2].color_specular[2];
        specColor = specColor.mult(pow(RdotV, alpha*triangle.v[0].shininess + beta*triangle.v[1].shininess + gamma*triangle.v[2].shininess));

        Vec color;
        color.x = lights[i].color[0];
        color.y = lights[i].color[1];
        color.z = lights[i].color[2];

        //formula for phone shading
        color.x = color.x * ((diffuseColor.x) +specColor.x);
        color.y = color.y * ((diffuseColor.y) +specColor.y);
        color.z = color.z * ((diffuseColor.z) +specColor.z);

        phong = phong.add(color);
      }
     } else if(isRightPixel) cout << "DID NOT GET HERE: " << hitSphere << endl;// else if(hitSphere == TRIANGLE && minIdx == 0) cout << "IN SHADOW BUT TRIANGLE" << endl;
  }

  return phong;

}

Vec directionRay(int x, int y) {
  Vec direction;
  double x_cord = (1.0*x/WIDTH)*2*aspectRatio*tanConst;
  double y_cord = (1.0*y/HEIGHT)*2*tanConst;
  x_cord += (-aspectRatio * tanConst);
  y_cord += -tanConst;

  direction.x = x_cord;
  direction.y = y_cord;
  direction.z = -1.0;
  direction = direction.normalize();

  return direction;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  aspectRatio = (double)WIDTH/(double)HEIGHT;
  tanConst = tan(fov/2.0 *3.14159265/180);
  
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      double minT = MAXT;
      int minIdx = -1;
      double r,g,b;
      Vec ray = directionRay(x,y);
      Hit hitSphere = NONE;
      if(x == 202 && y == 318) isRightPixel = true;
      else isRightPixel = false;
      if(isRightPixel) cout << "SHOULD: " << hitSphere << endl;
      const Vec temp(ray.x,ray.y,ray.z);
      if(checkIntersection(ray, minT, minIdx, hitSphere)) {
        // if(temp.x != ray.x || temp.y != ray.y || temp.z !=ray.z) {
        if(isRightPixel) cout << "WHAT IS CHANGING " << minT << endl;
        // }
        Vec color = phongLight(ray, minT, minIdx, hitSphere);
        color.x +=ambient_light[0];
        color.y +=ambient_light[1];
        color.z +=ambient_light[2];
        color = color.clamp();
        color = color.mult(255);
        r = color.x;
        g = color.y;
        b = color.z;
        if(minIdx == 0 && hitSphere == TRIANGLE) {
          // cout << x << " " << y << endl;
        //  cout << "RAY color: " << r<<endl;
        //   cout << "RAY color: " << g<<endl;
        //   cout << "RAY color: " << b<<endl;
        //   cout << "ambient color: " << ambient_light[0]<<endl;
        }
        
      }
      else  {if(isRightPixel) {
        cout << "RAY color: " << r<<endl;
          cout << "RAY color: " << g<<endl;
          cout << "RAY color: " << b<<endl;}r = 255; g= 255; b =255;}
      plot_pixel(x, y, r,g,b);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}