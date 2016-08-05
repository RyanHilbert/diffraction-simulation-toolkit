#ifndef MAIN_H//main header file; contains various constants, structs, and macros; all configurable parameters should be defined at the end of this file
#define MAIN_H
#include"tinyexpr.h"//needed for displacement functions
#include"h5.h"

enum CrystalShape{CUBOID,SPHERE,CYLINDER,CONE,BICONE,PRISM,PYRAMID,BIPYRAMID};
enum BeamType{PLANEWAVE,GAUSSIAN,FZP};
typedef struct{double x,y,z;}vector;
typedef struct{double e11,e12,e13,e21,e22,e23,e31,e32,e33;}matrix;
typedef struct{int back,front,num;}boundary;
typedef struct{double scalar,x,y,z;}quaternion;
typedef struct{
	double pix_size; //pixel size in um
	double dist; //distance to origin
	int row; //number of pixels in horizontal direction
	int column; //number of pixels in horizontal direction
}detector;

static const double pi=3.1415926535897932384626433832795;
#undef complex

#ifdef __cplusplus
#include<complex>//treat complex in C++ as std::complex<double>
static const std::complex<double>I(0,1);
typedef std::complex<double>complex;
extern"C"{
#else//treat complex in C as double _Complex
typedef double _Complex complex;
#include<stdbool.h>
#endif

extern int g_max_cpu_cores;//global maximum number of CPU cores to use in calculation
extern int g_max_gpus;//global maximum number of GPUs to use in calculation
extern hid_t hdf5file;//global HDF5 file identifier

int get_num_procs_omp();

extern size_t getDeviceCountCuda();
extern size_t getMaxThreadsPerBlockCuda();
extern void setDeviceCuda(int);
extern void*mallocCuda(size_t);
extern void freeCuda(void*);
extern void memsetCuda(void*,int,size_t);
extern void kernelCuda(size_t,size_t,vector*,void*,double,double,vector,double,double);
extern void memcpyHostToDeviceCuda(void*,const void*,size_t);
extern void memcpyDeviceToHostCuda(void*,const void*,size_t);

vector vec_unit(vector);
vector vec_num(double, vector);
vector vec_plus(vector, vector);
vector vec_minus(vector, vector);
double vec_dist(vector, vector);
double vec_abs(vector);
double vec_dot(vector, vector);
vector vec_cross(vector, vector);
vector rotation(matrix, vector);
matrix rot_m(vector, double);
vector oblique_to_cartesian(double, vector, double, vector, double, vector);

vector rot_axis();
vector det_direction();
quaternion direction_to_quaternion(vector);

complex planewave(double,vector,vector,vector);
complex gaussian(double,vector,vector,vector);
complex fzp(double,vector,vector,vector);

bool cuboid(vector);
bool sphere(vector);
bool cylinder(vector);
bool cone(vector);
bool bicone(vector);
bool prism(vector);
bool pyramid(vector);
bool bipyramid(vector);

#ifdef __cplusplus
}
#endif
#endif//do all of the above only once per file

#ifdef USE_ALL_MACROS
#define CODE(...) __VA_ARGS__
#else//if USE_ALL_MACROS is not defined, process only macros used for declaring parameters
#define CODE(...)
#define LNK(...)
#define MUL(...)
#define GRP(...)
#define TIP(...)
#define MNU(name,var,...) PRM(name,int,var,0)
#endif

#ifdef FOREACH//if FOREACH is defined, apply it to each parameter
#define PRM FOREACH
#elif __cplusplus//otherwise, declare each of the parameters (standard header functionality, may be done repeatedly)
#define PRM(name,type,var,...) extern"C"{extern type var;}
#else
#define PRM(name,type,var,...) extern type var;
#endif

//declare all configurable global parameters below using the PRM(parameter) and MNU(menu) macros
//parameters declared on the same line will appear on the same line in the GUI
//the variable given to a MNU is set to the index of the currently selected item, beginning at zero
//all enums must be declared in the order they appear in their respective MNU
//follow a PRM with LNK(link) to link the parameter to the given values of the last added MNU, may be combined with MUL
//unlike MNU, LNK does not care about the order of its arguments
//follow a PRM with MUL(multiply) to multiply the user input by a conversion value before giving it to the computation (doubles only)
//use the GRP(group) macro to create a new group of parameters
//the TIP(tooltip) macro adds a helpful tooltip to a parameter
//the CODE(code) macro is used by the GUI implementation to insert irregular widgets into the parameter panels

//GRP syntax: GRP("display name of parameter group")
//MNU syntax: MNU("variable display name",variable type,variable name,"display names","of each","menu option"...)
//PRM syntax: PRM("variable display name",variable type,variable name,variable default value)
//LNK syntax: PRM(...)LNK(index0,index1,index2...)
//MUL syntax: PRM(...)MUL(conversion constant)
//TIP syntax: PRM(...)TIP("tooltip text")

GRP("Computation")
PRM("Save wave field slices?",bool,g_saving_wavefield,false)
PRM("Grid cell size (μm)",double,g_ds,0.005)MUL(1)                          PRM("Slice step size (μm)",double,g_dy,0.005)MUL(1)
PRM("Convergence threshold",double,g_tolerance,1e-6)
PRM("Rotation axis",vector,g_rot_axis,{0,0,0})TIP("An axis value will be calculated in place of all zeroes")
PRM("Start deviation angle (°)",double,g_dth_start,0)MUL(pi/180)
PRM("Angular step (°)",double,g_dth_step,0.001)MUL(pi/180)                  PRM("Angular scan points",int,g_num_angle,1)
PRM("Spiral scan radius (μm)",double,g_spiral_c,0.1)MUL(1)                  PRM("Spiral scan points",int,g_num_scan,1)

GRP("Detector")
PRM("Pixel size (μm)",double,g_pixel_size,55)MUL(1) PRM("Rotation angle (°)",double,g_alpha,0.0)MUL(pi/180) PRM("Rows",int,g_row,256)
PRM("Sample to detector distance (m)",double,g_det_dist,1.5)MUL(1e6)        PRM("Columns",int,g_column,256)

GRP("Crystal")
MNU("Shape",g_shape,"Cuboid","Sphere","Cylinder","Cone","Bicone","Prism","Pyramid","Bipyramid")
PRM("Radius (μm)",double,g_radius,.3)MUL(1)TIP("For shapes with base edges, this is the circumradius of the base")LNK(SPHERE,CYLINDER,CONE,BICONE,PRISM,BIPYRAMID,PYRAMID) PRM("Length (μm)",double,g_half_length,.15)MUL(.5)LNK(CUBOID) PRM("Width (μm)",double,g_half_width,.15)MUL(.5)LNK(CUBOID) PRM("Height (μm)",double,g_half_height,.15)MUL(.5)LNK(CUBOID,CYLINDER,CONE,BICONE,PRISM,BIPYRAMID,PYRAMID) PRM("Base edges",int,g_edges,3)LNK(PRISM,BIPYRAMID,PYRAMID)TIP("Minimum of 3")
PRM("d spacing (Å)",double,g_d_spacing,5.43/8.0)MUL(1e-4)
CODE(array[0]=param;group->line->addWidget(spacingButton=new DialogButton("Calculate",spacing,[=]{((DoubleParam*)param)->set(spacing->result());}));)
PRM("Direction of h vector",vector,g_unit_h,{0,0,-1})CODE(array[2]=param;)
PRM("X0",complex,g_chi_0,-0.74984e-5+0.88508e-7*I)
PRM("Xh",complex,g_chi_h,-0.17751e-5+0.66803e-7*I)
PRM("x displacement function",char*,g_xstr,0)
PRM("y displacement function",char*,g_ystr,0)
PRM("z displacement function",char*,g_zstr,0)

GRP("Beam")
PRM("Energy (keV)",double,g_energy,11.4)MUL(1)CODE(array[1]=param;)         MNU("Beam type",g_beam_type,"Plane Wave","Gaussian","Fresnel Zone Plate")
PRM("Incident beam direction",vector,g_unit_s0,{0.598481,0,0.801137})
CODE(direction->setParams((DoubleParam*)array[0],(DoubleParam*)array[1],(VectorParam*)array[2],(VectorParam*)param);group->line->addWidget(directionButton=new DialogButton("Calculate",direction,[=]{((VectorParam*)param)->set(direction->result());}));)
PRM("Beam center",vector,g_r0,{0,0,0})
PRM("Flux density (photons/mm²s)",double,g_flux_density,1e12)MUL(1e-6)LNK(PLANEWAVE,FZP)
PRM("Inner radius (μm)",double,g_inner_radius,20)MUL(1)LNK(FZP)             PRM("Outer radius (μm)",double,g_outer_radius,80)LNK(FZP)MUL(1)
PRM("Efficiency",double,g_efficiency,0.2)LNK(FZP)                           PRM("Focal length (mm)",double,g_focal_length,147)LNK(FZP)MUL(1000)
PRM("Width (μm)",double,g_w0,1)MUL(1)LNK(GAUSSIAN)                          PRM("Total flux (photons/s)",double,g_total_flux,1e10)LNK(GAUSSIAN)MUL(1)

#undef GRP
#undef LNK
#undef MNU
#undef MUL
#undef PRM
#undef TIP
#undef CODE
#undef FOREACH
#undef USE_ALL_MACROS