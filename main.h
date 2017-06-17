#ifndef MAIN_H//main header file; contains various constants, structs, and macros; all configurable parameters should be defined at the end of this file
#define MAIN_H

#include<math.h>
#include"tinyexpr.h"//needed for displacement functions

#undef complex//suppresses warning
#ifdef __cplusplus
#include<complex>//treat complex in C++ as std::complex<double>
#undef I
#define I COMPLEX_I
#define cexp(...) std::exp(__VA_ARGS__)
#define cpow(...) std::pow(__VA_ARGS__)
#define cabs(...) std::abs(__VA_ARGS__)
#define csqrt(...) std::sqrt(__VA_ARGS__)
typedef std::complex<double>complex;
static const complex I(0,1);
extern"C"{
#else//treat complex in C as double _Complex
#include<stdbool.h>
#include<complex.h>
typedef double _Complex complex;
#endif//both forms of complex are required by C/C++ standards to have the same memory layout

enum CrystalShape{CUBOID,SPHERE,HEMISPHERE,CYLINDER,CONE,BICONE,PRISM,PYRAMID,BIPYRAMID};
enum BeamType{PLANEWAVE,GAUSSIAN,FZP};
typedef struct{float x,z;}pvector;
typedef struct{double x,y,z;}vector;
typedef struct{double e11,e12,e13,e21,e22,e23,e31,e32,e33;}matrix;
typedef struct{int back,front,num;}boundary;
typedef struct{double scalar,x,y,z;}quaternion;
typedef struct{
	int col;//number of pixels in horizontal direction
	int row;//number of pixels in vertical direction
	double pix_size;//pixel size in um
	double dist;//distance to origin
}detector;
typedef struct{
	size_t col;
	size_t row;
	float* d0;
	float* dh;
	pvector* poynting;
}wavefield;

typedef char char256[256];//typedefs char256 as an array of 256 characters

#define PI 3.1415926535897932384626433832795

//this macro causes a compilation error if its argument is not true during compilation
#define ASSERT(...) extern char ASSERT[1/(int)!!(__VA_ARGS__)];

//conditional code generation macro
#define IF(...)CAT(IF,IS0(__VA_ARGS__))
//usage: IF(condition)(code)
//code will be generated only if condition is non-zero during preprocessing

extern int g_max_cpu_cores;//global maximum number of CPU cores to use in calculation
extern int g_max_gpus;//global maximum number of GPUs to use in calculation

int get_num_procs_omp();
size_t get_num_iterations();

int getDeviceCountCuda();
int getMaxThreadsPerBlockCuda();
void setDeviceCuda(int);
void*mallocCuda(size_t);
void freeCuda(void*);
void memsetCuda(void*,int,size_t);
void kernelCuda(size_t,size_t,vector*,void*,double,double,vector,double,double);
void memcpyHostToDeviceCuda(void*,const void*,size_t);
void memcpyDeviceToHostCuda(void*,const void*,size_t);

vector vec_unit(vector);
vector vec_num(double,vector);
vector vec_plus(vector,vector);
vector vec_minus(vector,vector);
double vec_dist(vector,vector);
double vec_abs(vector);
double vec_dot(vector,vector);
vector vec_cross(vector,vector);
vector rotation(matrix,vector);
matrix rot_m(vector,double);
vector oblique_to_cartesian(double,vector,double,vector,double,vector);
wavefield rec_grid_interpolation(complex*,complex*,vector*,int,int,vector,vector,vector,boundary*,boundary*,int);

vector rot_axis();
vector det_direction();
quaternion direction_to_quaternion(vector);

bool cuboid(vector);
bool sphere(vector);
bool hemisphere(vector);
bool cylinder(vector);
bool cone(vector);
bool bicone(vector);
bool prism(vector);
bool pyramid(vector);
bool bipyramid(vector);

//macros used to support IF(condition)(code) macro
#define IF0(...)
#define IF1(...)__VA_ARGS__
#define CATEX(_,...)_##__VA_ARGS__
#define CAT(_,...)CATEX(_,__VA_ARGS__)
#define GET2EX(x,_,...)_
#define GET2(...)IF1(GET2EX(__VA_ARGS__))
#define _0 ,0
#define IS0(...)GET2(CAT(_,__VA_ARGS__),1,)

#ifdef __cplusplus
}
#endif
#endif//do all of the above only once per file

#ifdef USE_ALL_MACROS
#define CODE(...) __VA_ARGS__
#else//if USE_ALL_MACROS is not defined, process only macros used for declaring parameters
#define CODE(...)
#define LINK(...)
#define MUL(...)
#define GROUP(...)
#define TIP(...)
#define MENU(id,name,var,...) PARAM(id,name,char,var,0)
#endif

#ifdef FOREACH//if FOREACH is defined, apply it to each parameter
#define PARAM FOREACH
#elif __cplusplus//otherwise, declare each of the parameters (standard header functionality, may be done repeatedly)
#define PARAM(id,name,type,var,...) extern"C"{extern type var;}
#else
#define PARAM(id,name,type,var,...) extern type var;
#endif

//declare all configurable global parameters below using the PARAM(parameter) and MENU macros
//parameters declared on the same line will appear on the same line in the GUI
//the variable given to a MENU is set to the index of the currently selected item, beginning at zero
//all enums must be declared in the order they appear in their respective MENU
//follow a PARAM with LINK to link the parameter to the given values of the last added MENU, may be combined with MUL
//unlike MENU, LINK does not care about the order of its arguments
//follow a PARAM with MUL(multiply) to multiply the user input by a conversion value before giving it to the computation (doubles only)
//use the GROUP macro to create a new group of parameters
//the TIP(tooltip) macro adds a helpful tooltip to a parameter
//the CODE(code) macro is used by the GUI implementation to insert irregular widgets into the parameter panels

//GROUP syntax: GROUP("display name of parameter group")
//MENU  syntax:  MENU("variable display name",variable type,variable name,"display names","of each","menu option"...)
//PARAM syntax: PARAM("variable display name",variable type,variable name,variable default value)
//LINK  syntax: PARAM(...)LINK(index0,index1,index2...)
//MUL   syntax: PARAM(...)MUL(conversion constant)
//TIP   syntax: PARAM(...)TIP("tooltip text")

//"hidden" parameters for the d-spacing calculator
PARAM(100,"a(Å)",double,g_a,0) PARAM(101,"b(Å)",double,g_b,0) PARAM(102,"c(Å)",double,g_c,0)
PARAM(103,"α(°)",double,g_A,0) PARAM(104,"β(°)",double,g_B,0) PARAM(105,"γ(°)",double,g_C,0)
PARAM(106,"h",int,g_h,0)       PARAM(107,"k",int,g_k,0)       PARAM(108,"l",int,g_l,0)

GROUP("Computation")
PARAM(0,"Save wave field slices",bool,g_saving_wavefield,false)TIP("Saves and displays the internal wavefield structure, but results in a slower calculation and much larger file size")MENU(39,"\tScan type",g_scan_type,"Spiral","Rectangular")
PARAM(12,"Rotation axis",vector,g_rot_axis,{0,0,0})TIP("An axis value will be calculated in place of all zeroes")
PARAM(15,"Grid cell size (μm)",double,g_ds,.005)                     PARAM(13,"Slice step size (μm)",double,g_dy,.005)
PARAM(16,"Convergence threshold",double,g_tolerance,1e-6)            PARAM(14,"Start deviation angle (°)",double,g_dth_start,0)MUL(PI/180)
PARAM(17,"Angular step (°)",double,g_dth_step,.001*PI/180)MUL(PI/180)PARAM(4,"Angular scan points",int,g_num_angle,1)
PARAM(18,"Spiral scan radius (μm)",double,g_spiral_c,0.1)LINK(0)     PARAM(3,"Spiral scan points",int,g_num_scan,1)LINK(0)
PARAM(40,"x1 (μm)",double,g_scan_x1,0)LINK(1) PARAM(41,"x2 (μm)",double,g_scan_x2,0)LINK(1) PARAM(42,"Horizontal scan points",int,g_scan_nx,1)LINK(1)
PARAM(43,"y1 (μm)",double,g_scan_y1,0)LINK(1) PARAM(44,"y2 (μm)",double,g_scan_y2,0)LINK(1) PARAM(45,"Vertical scan points",int,g_scan_ny,1)LINK(1)

GROUP("Detector")
PARAM(27,"Pixel size (μm)",double,g_pixel_size,55)                   PARAM(29,"Rotation angle (°)",double,g_alpha,0)MUL(PI/180)PARAM(2,"Rows",int,g_row,256)
PARAM(28,"Sample to detector distance (m)",double,g_det_dist,1500000)MUL(1e6)                                                  PARAM(1,"Columns",int,g_column,256)

GROUP("Crystal")
MENU(6,"Shape",g_shape,"Cuboid","Sphere","Hemisphere","Cylinder","Cone","Bicone","Prism","Pyramid","Bipyramid") PARAM(31,"d spacing (Å)",double,g_d_spacing,.000067875)MUL(1e-4) CODE(array[0]=param;group->line->addWidget(spacingButton=new DialogButton("Calculate",spacing,[=]{((DoubleParam*)param)->set(spacing->result());}));)
PARAM(30,"Radius (μm)",double,g_radius,.15)TIP("For shapes with base edges, this is the circumradius of the base")LINK(SPHERE,HEMISPHERE,CYLINDER,CONE,BICONE,PRISM,BIPYRAMID,PYRAMID)   PARAM(35,"Length (μm)",double,g_half_length,.075)MUL(.5)LINK(CUBOID) PARAM(36,"Width (μm)",double,g_half_width,.075)MUL(.5)LINK(CUBOID) PARAM(37,"Height (μm)",double,g_half_height,.075)MUL(.5)LINK(CUBOID,CYLINDER,CONE,BICONE,PRISM,BIPYRAMID,PYRAMID) PARAM(38,"Base edges",int,g_edges,3)LINK(PRISM,BIPYRAMID,PYRAMID)TIP("Minimum of 3")
PARAM(9,"Direction of h vector",vector,g_unit_h,{0,0,-1})CODE(array[2]=param;)
PARAM(8,"X0",complex,g_chi_0,-.74984e-5+.88508e-7*I)
PARAM(7,"Xh",complex,g_chi_h,-.17751e-5+.66803e-7*I)
PARAM(32,"x displacement function",char256,g_xstr,{0})
PARAM(33,"y displacement function",char256,g_ystr,{0})
PARAM(34,"z displacement function",char256,g_zstr,{0})

GROUP("Beam")
PARAM(19,"Energy (keV)",double,g_energy,11.4)CODE(array[1]=param;)    MENU(5,"Beam type",g_beam_type,"Plane Wave","Gaussian","Fresnel Zone Plate")
PARAM(10,"Incident beam direction",vector,g_unit_s0,{0.598481,0,0.801137})CODE(direction->setParams((DoubleParam*)array[0],(DoubleParam*)array[1],(VectorParam*)array[2],(VectorParam*)param);group->line->addWidget(directionButton=new DialogButton("Calculate",direction,[=]{((VectorParam*)param)->set(direction->result());}));)
PARAM(11,"Beam center",vector,g_r0,{0,0,0})
PARAM(20,"Flux density (photons/mm²s)",double,g_flux_density,1e6)MUL(1e-6)LINK(PLANEWAVE,FZP)
PARAM(21,"Inner radius (μm)",double,g_inner_radius,20)LINK(FZP)       PARAM(24,"Outer radius (μm)",double,g_outer_radius,80)LINK(FZP)
PARAM(22,"Efficiency",double,g_efficiency,0.2)LINK(FZP)               PARAM(25,"Focal length (mm)",double,g_focal_length,147000)LINK(FZP)MUL(1000)
PARAM(23,"Width (μm)",double,g_w0,1)LINK(GAUSSIAN)                    PARAM(26,"Total flux (photons/s)",double,g_total_flux,1e10)LINK(GAUSSIAN)

#undef GROUP
#undef LINK
#undef MENU
#undef MUL
#undef PARAM
#undef TIP
#undef CODE
#undef FOREACH
#undef USE_ALL_MACROS
