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

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double minimum(double v1, double v2) {
  if(v1 < v2) return v1;
  return v2;
}

bool sphereIntersection(Vec v, double& minT, int& minIdx) {
  double xc,yc,zc,r,b,c, min;


  for(int i = 0; i < num_spheres; i++) {
    double t0,t1;
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

      if((t0 > 0 && t1 > 0)) min = minimum(t0,t1);
      if(t0 > 0) min = t0;
      if(t1 > 0) min = t1;

      if(min < minT) {
        minT = min;
        minIdx = i;
      }
      // else if(t1 > 0) min = t1;
      // else if(t0 > 0) min = t0;
      // else continue;
    } else continue;
    if(t0 <=0 || t1 <=0) continue;
    

  }
  return true;

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

void barycentricCoord(Vec p, Triangle t, double &alpha, double &beta, double &gamma) {
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
  alpha = areaCC1C2 / areaAllV;
  beta = areaC0CC2 / areaAllV;
  gamma = areaC0C1C / areaAllV;
}

bool triangleIntersection(Vec v, Triangle t) {
  Vec planeNormal = triangleNormal(t);
  double a,b,c, d, denom;
  a = planeNormal.x;
  b = planeNormal.y;
  c = planeNormal.z;

  d = -1.0 * (a*t.v[0].position[0]+b*t.v[0].position[1]+c*t.v[0].position[2]);
  denom = ((a*v.x) +(b*v.y)+(c*v.z));

  if(denom == 0) return false;

  double t_value = -1 * (a+b+c+d)/denom;
  if(t_value <= 1e-8) return false;

  //check if i have to multiply t
  v = v.mult(t_value);
  double alpha, beta, gamma;
  barycentricCoord(v, t, alpha,beta,gamma);
  if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) return false;
  return true;
}

bool checkIntersection(Vec v, &min_t) {
  
  for(int i = 0; i < num_triangles; i++) {
    if(triangleIntersection(v,triangles[i])) return true;
  }
  return (sphereIntersection(v));
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

  //  cout << "x: " << direction.x << endl;
  //   cout << "y: " << direction.y << endl;
  //   cout << "z: " << direction.z << endl;

  direction = direction.normalize();
  //  cout << "x: " << direction.x << endl;
  //   cout << "y: " << direction.y << endl;
  //   cout << "z: " << direction.z << endl;
  //   cout <<endl;
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
      int r,g,b;
      Vec ray = directionRay(x,y);
      double minT;
      if(checkIntersection(ray)) {
        r = 87; g= 0; b =0;
      }
      else  {r = 0; g= 0; b =255;}
      plot_pixel(x, y, r,g,b);
      // plot_pixel(x, y, x % 256, y % 256, (x+y) % 256); // replace with some logic. sphere
      //if sphere intersection plot pixel green 
      //else plot pixel blue
      //
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

