/*Numerical simulation code for x-ray dynamical diffraction from a particle,*/
/*based on the iterative solving algorithm published in PRB 89, 014104 (2014)*/
/*Written by Hanfei Yan, 2014*/

/// OpenMP, CUDA, & Qt functionality added by Ryan Hilbert
/// <--- Added comments preceded by triple slash

//V4.0

#undef __STRICT_ANSI__ ///allows non-standard j0 bessel function to be defined in MinGW's math.h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<complex.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include"main.h"

volatile bool halted=false;

#ifdef __cplusplus
extern"C"{
#endif///C functions for accessing Qt's signal-slot mechanism, defined in thread.cpp
extern void emit_output(char*);
extern void emit_maximum(int);
extern void emit_progress(int);
extern void emit_completed1st();
int g_max_cpu_cores=1;///maximum number of CPU cores to use, configurable by GUI
int g_max_gpus=0;///maximum number of GPUs to use, configurable by GUI
#ifdef __cplusplus
}
#endif

static const double spiral_phi = 2.4; //radian

///array of points located twice as far away from the polygon's center as each edge's midpoint
static complex*points=0;
///used for determining if a 2D point lies within a regular polygon for the PRISM and BIPYRAMID crystal() calculations

static inline double sq(double n){return n*n;}

#define printf(...) {char str[BUFSIZ];sprintf(str,__VA_ARGS__);emit_output(str);}
#ifdef __cplusplus
extern"C"{
#endif
void calc(){
	time_t now;
	time(&now);
	printf("\n________________________________\n");
	/*disabled due to large amount of compiler warnings generated
	///verify that these sizes are all different, otherwise the parameter printing will fail
	printf("sizeof(int) = %d",sizeof(int));
	printf("sizeof(double) = %d",sizeof(double));
	printf("sizeof(vector) = %d",sizeof(vector));
	///print the value of each parameter as the computation code sees it (disregard printed units; they do not take the conversion multiplier into account)
#define FOREACH(name,type,var,...) \
	if(sizeof(type)==sizeof(int)){printf("%s = %d",name,*(int*)&var);}\
	else if(sizeof(type)==sizeof(double)){printf("%s = %G",name,*(double*)&var);}\
	else if(sizeof(type)==sizeof(vector)){printf("%s = %Gx %+Gy %+Gz",name,((vector*)&var)->x,((vector*)&var)->y,((vector*)&var)->z);}
#include"main.h"
	*/
	complex(*beam)(double,vector,vector,vector);
	switch(g_beam_type){///choose beam function
	case PLANEWAVE:beam=planewave;break;
	case GAUSSIAN:beam=gaussian;break;
	case FZP:beam=fzp;break;
	}
	bool(*crystal)(vector);
	switch(g_shape){///choose crystal function
	case CUBOID:crystal=cuboid;break;
	case SPHERE:crystal=sphere;break;
	case CYLINDER:crystal=cylinder;break;
	case CONE:crystal=cone;break;
	case BICONE:crystal=bicone;break;
	case PRISM:crystal=prism;break;
	case PYRAMID:crystal=pyramid;break;
	case BIPYRAMID:crystal=bipyramid;break;
	}
	free(points);
	if(g_shape==PRISM||g_shape==PYRAMID||g_shape==BIPYRAMID){///initialize the list of comparison points for the calculation in crystal()
		complex vertices[g_edges];
		for(int i=0;i<g_edges;++i){
			vertices[i]=g_radius*sin(i*2*pi/g_edges)+g_radius*cos(i*2*pi/g_edges)*I;
		}
		points=(complex*)malloc(g_edges*sizeof*points);
		if(g_edges>0)points[0]=vertices[g_edges-1]+vertices[0];
		for(int i=1;i<g_edges;++i){
			points[i]=vertices[i-1]+vertices[i];
		}
	}
	const double lambda = 12.398/g_energy*1e-4;
	const double wave_number = 2.0*pi/lambda;
	const vector start_position = g_r0;

	const vector unit_s0 = vec_unit(g_unit_s0);
	const vector unit_s0_0 = unit_s0; // initial unit vector

	const vector unit_h = vec_unit(g_unit_h);
	const vector h = vec_num(2.0*pi / g_d_spacing, unit_h);

	//define intinal detector coordinate system and the rotation axis

	const vector det_z_0 = vec_unit(vec_plus(vec_num(wave_number, unit_s0_0), h)); //detector plane normal is along kh direction

	matrix det_R = rot_m(det_z_0, g_alpha);
	const vector det_x_0 = rotation(det_R,vec_unit(vec_cross(unit_h, unit_s0_0))); //detector column direction is perpendicular to the diffraction plane define by s0 and h
	//detector rotates in plane by an angle alpha

	const vector det_y_0 = vec_unit(vec_cross(det_z_0, det_x_0));

	///rot_axis is calculated only if zero is provided in gui
	const vector rot_axis = g_rot_axis.x||g_rot_axis.y||g_rot_axis.z?g_rot_axis:det_x_0; // rotation in the diffraction plane and a positive angle make the incidence angle to the diffracting lattice place larger.

	printf("det_x = (%G, %G, %G)\ndet_y = (%G, %G, %G)\ndet_z = (%G, %G, %G)\n", det_x_0.x,det_x_0.y,det_x_0.z,det_y_0.x,det_y_0.y,det_y_0.z,det_z_0.x,det_z_0.y,det_z_0.z);
	//create a log file
	/*char log_file_name[FILENAME_MAX];
	strcpy(log_file_name, dir_name);
	strcat(log_file_name, "/log.txt");
	FILE * log;
	log = fopen(log_file_name,"w+");*/
	const size_t thread_count =
		#ifdef _OPENMP
			g_max_cpu_cores<1?1:(g_max_cpu_cores>omp_get_num_procs()?omp_get_num_procs():g_max_cpu_cores);
	omp_set_dynamic(false);///disables automatic lowering of thread counts
	omp_set_nested(true);///allows nested parallelism
	printf("OpenMP Enabled: %d Threads\n", thread_count);
#else
			1;printf("OpenMP not enabled; limited to one GPU maximum!\n");
#endif
	const size_t device_count = g_max_gpus<0?0:(g_max_gpus>getDeviceCountCuda()?getDeviceCountCuda():g_max_gpus);
	const size_t cuda_block_size = getMaxThreadsPerBlockCuda();
	printf("%d CUDA Devices Enabled\n",device_count);
	////////////////////////////////////////////////////////
	//Note: xyz is crystal coordinates system
	const int num_iterations=g_num_scan*g_num_angle;
	emit_maximum(num_iterations);
	hid_t detector_group=h5group("Detector",hdf5file);///opens the HDF5 file's Detector group
	h5set_int(num_iterations,"Number of Images",detector_group);
	int iterations_completed=h5int("Iterations Completed",detector_group);
	if(iterations_completed<=0)iterations_completed=0;
	else emit_completed1st();
	emit_progress(iterations_completed);
	double minFrameX=INFINITY,maxFrameX=-INFINITY,minFrameZ=INFINITY,maxFrameZ=-INFINITY;
#pragma omp parallel for schedule(dynamic) num_threads(thread_count)
	for(int main_index = iterations_completed; main_index < num_iterations; ++main_index)if(!halted){///integer counter needed for OpenMP

		const double dth = main_index / g_num_scan * g_dth_step + g_dth_start;//deviation angle to the exact Bragg angle; scan position first
		const vector beam_center = {
			g_spiral_c*sqrt(main_index % g_num_scan)*cos(spiral_phi*(main_index % g_num_scan)) + start_position.x,
			g_spiral_c*sqrt(main_index % g_num_scan)*sin(spiral_phi*(main_index % g_num_scan)) + start_position.y,
			start_position.z
		};
		///fprintf(log, "%d\t%f\t%f\t%f\n", file_num, beam_center.x, beam_center.y, dth*180.0/pi);///Line causes error "Invalid parameter passed to C runtime function."

		//calculate the new k0 after a rotation
		const matrix R = rot_m(rot_axis, dth);
		const vector unit_s0 = rotation(R, unit_s0_0);
		printf("unit_s0 = (%G, %G, %G)\n", unit_s0.x, unit_s0.y, unit_s0.z);
		const vector k0 = vec_num(wave_number, unit_s0);

		//calculate the new kh after a roration
		const vector kh = vec_plus(k0, h);

		//calculate the corresponding unit vectors
		const vector unit_sh = vec_unit(kh);
		const vector unit_sy = vec_unit(vec_cross(unit_s0,unit_sh));

		//create a grid in oblique coordinates system

		////////////////////////////////////////////

		//initilize detector parameters
		const detector det_para = { g_pixel_size, g_det_dist, g_row, g_column };

		//rotate base vectors of detector
		const vector det_x = rotation(R, det_x_0);
		const vector det_y = rotation(R, det_y_0);
		const vector det_z = rotation(R, det_z_0);
		printf("After rotation \ndet_x = (%G, %G, %G)\ndet_y = (%G, %G, %G)\ndet_z = (%G, %G, %G)\n", det_x.x,det_x.y,det_x.z,det_y.x,det_y.y,det_y.z,det_z.x,det_z.y,det_z.z);
		complex *det = (complex*)calloc(det_para.row*det_para.column, sizeof(complex));
		vector *det_pix = (vector*)calloc(det_para.row*det_para.column, sizeof(vector)); //position of a pixel on detector
		for (int i = 0; i < det_para.row; i++){
			for (int j = 0; j < det_para.column; j++){
				det_pix[i*det_para.column+j] = vec_plus(vec_num((j - det_para.column/2)*det_para.pix_size, det_x), vec_num((i - det_para.row/2)*det_para.pix_size, det_y));
				det_pix[i*det_para.column+j] = vec_plus(det_pix[i*det_para.column+j], vec_num(det_para.dist, det_z));
				//                det_pix[i*det_para.column+j] = rotation(R, det_pix[i*det_para.column+j]);
			}
		}
		const int det_size = det_para.row*det_para.column;
		const int thread_id = 0
		#ifdef _OPENMP
				+omp_get_thread_num();
#endif
		const size_t cuda_grid_size = (det_size - 1)/cuda_block_size + 1;
		const int device_id = device_count>0?thread_id%device_count:0;
		setDeviceCuda(device_id);
		vector*cuda_det_pix=(vector*)mallocCuda(det_para.row*det_para.column*sizeof*det_pix);
		void*cuda_det=mallocCuda(det_para.row*det_para.column*sizeof*det);///actual type is cuDoubleComplex*
		memsetCuda(cuda_det, 0, det_para.row*det_para.column*sizeof*det);
		memcpyHostToDeviceCuda(cuda_det_pix, det_pix, det_para.row*det_para.column*sizeof*det_pix);

		///all printing in one statement eliminates concurrency issues
		printf("\nThread %d:\nBeam position = (%G, %G, %G)\tDeviation angle = %G\n", thread_id, beam_center.x, beam_center.y, beam_center.z, dth*180.0/pi);
		///////////////////////////////////////////
		//calculate the nmber of slide
		int flag_y_pos = 0, flag_y_neg = 0; //two status flags to show whether the positive or negative sy direction is searched
		int c = 0;
		double y, y_max = 0.0, y_min = 0.0;
		const double ds = g_ds;
		const double dy = g_dy;
		while (flag_y_pos == 0 || flag_y_neg == 0){	//search is not done
			if (flag_y_pos == 0){ //positive direction search is not done
				y = c*dy; //step toward positive sy axis
				c++;
			}
			else{
				y = -c*dy; //step toward negative sy direction
				c++;
			}
			//find the initial guess of the size of the grid in s0 direction
			int a = 0;
			int b = 0;

			vector oblique;
			oblique.x = b*ds;
			oblique.y = y;
			oblique.z = a*ds;
			vector cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			while (crystal(cartesian) == 1){
				b++;
				oblique.x = b*ds;
				cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			}
			int grid_size_s0 = 2 * b;
			//find the initial guess of the size of the grid in sh direction
			a = 0;
			b = 0;
			oblique.x = b*ds;
			oblique.y = y;
			oblique.z = a*ds;
			cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			while (crystal(cartesian) == 1){
				a++;
				oblique.z = a*ds;
				cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			}
			int grid_size_sh = 2*a;

			int status_s0 = 1;
			int status_sh = 1;
			//find the smallest grid that accomodate the particle
			//check four boundary lines of the grid, if there is an interception with the crystal, increase the size
			while (status_s0 == 1 || status_sh == 1){
				status_s0 = 0;
				oblique.x = (grid_size_s0 / 2)*ds;
				oblique.y = y;
				for (int i = -grid_size_sh / 2; i <= grid_size_sh / 2; i++){
					oblique.z = i*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_s0 = 1;
					}
				}
				oblique.x = -(grid_size_s0 / 2)*ds;
				oblique.y = y;
				for (int i = -grid_size_sh / 2; i <= grid_size_sh / 2; i++){
					oblique.z = i*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_s0 = 1;
					}
				}
				if (status_s0 == 1) {// there is an interception in sh direction
					grid_size_s0 += 2; // increase the size in s0 direction
				}
				status_sh = 0;
				oblique.z = (grid_size_sh / 2)*ds;
				oblique.y = y;
				for (int j = -grid_size_s0 / 2; j <= grid_size_s0 / 2; j++){
					oblique.x = j*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_sh = 1;
					}
				}
				oblique.z = -(grid_size_sh / 2)*ds;
				oblique.y = y;
				for (int j = -grid_size_s0 / 2; j <= grid_size_s0 / 2; j++){
					oblique.x = j*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_sh = 1;
					}
				}
				if (status_sh == 1){// there is an interception
					grid_size_sh += 2; // increase the size
				}
			}
			//take into account zero

			grid_size_sh++;
			grid_size_s0++;


			if (grid_size_s0 == 1 && grid_size_sh == 1){ // no interception point found
				if (flag_y_pos == 0){ //search on the positive direction
					flag_y_pos = 1; //postive y_max found
					y_max = y - dy; //y_max value
					c = 1; //reset the counter to step along the negative direction
				}
				else {
					flag_y_neg = 1; //negative y_min found
					y_min = y + dy;
				}
			}
		}
		char wavefield_group_name[32];
		hid_t wavefield_group;
		if(g_saving_wavefield){
			sprintf(wavefield_group_name,"det%dfield",main_index);
			wavefield_group=h5group(wavefield_group_name,detector_group);
		}
		const int num_slice = (y_max - y_min)/dy + 1;
		//printf("\nThread %d:\n%s\nBeam position = (%G, %G, %G)\tDeviation angle = %G\tNumber of slices: %d\n", thread_id, det_file_name,  beam_center.x, beam_center.y, beam_center.z, dth*180.0/pi, num_slice);
		//start calculation
		///printf("Iteration Number: %d\nThread Count: %d\n",main_index,(thread_count-1)/num_iterations+1);
#pragma omp parallel for schedule(dynamic) num_threads((thread_count-1)/num_iterations+1)///thread_count/remaining_iterations rounded-up
		for (int c = 0; c < num_slice; c++){//reuse c as a counter
			if(device_count)setDeviceCuda(thread_id%device_count);
#define DISPLACEMENT (vector){te_eval(xexpr),te_eval(yexpr),te_eval(zexpr)}
			vector v;///holds the variables used by displacement expressions, must be set before using DISPLACEMENT macro
			te_variable vars[]={{"x",&v.x},{"y",&v.y},{"z",&v.z}};
			te_expr*xexpr=te_compile(g_xstr?g_xstr:"0",vars,3,0);
			te_expr*yexpr=te_compile(g_ystr?g_ystr:"0",vars,3,0);
			te_expr*zexpr=te_compile(g_zstr?g_zstr:"0",vars,3,0);
			//find the initial guess of the size of the grid in s0 direction
			const double y = y_min + c*dy;
			int a = 0;
			int b = 0;

			vector oblique;
			oblique.x = b*ds;
			oblique.y = y;
			oblique.z = a*ds;
			vector cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			while (crystal(cartesian) == 1){
				b++;
				oblique.x = b*ds;
				cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			}
			int grid_size_s0 = 2 * b;
			//find the initial guess of the size of the grid in sh direction
			a = 0;
			b = 0;
			oblique.x = b*ds;
			oblique.y = y;
			oblique.z = a*ds;
			cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			while (crystal(cartesian) == 1){
				a++;
				oblique.z = a*ds;
				cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
			}
			int grid_size_sh = 2 * a;

			int status_s0 = 1;
			int status_sh = 1;
			//find the smallest grid that accomodates the particle
			//check four boundary lines of the grid, if there is an interception with the crystal, increase the size
			while (status_s0 == 1 || status_sh == 1){
				status_s0 = 0;
				oblique.x = (grid_size_s0 / 2)*ds;
				oblique.y = y;
				for (int i = -grid_size_sh / 2; i <= grid_size_sh / 2; i++){
					oblique.z = i*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_s0 = 1;
					}
				}
				oblique.x = -(grid_size_s0 / 2)*ds;
				oblique.y = y;
				for (int i = -grid_size_sh / 2; i <= grid_size_sh / 2; i++){
					oblique.z = i*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_s0 = 1;
					}
				}
				if (status_s0 == 1) {// there is an interception in sh direction
					grid_size_s0 += 2; // increase the size in s0 direction
				}
				status_sh = 0;
				oblique.z = (grid_size_sh / 2)*ds;
				oblique.y = y;
				for (int j = -grid_size_s0 / 2; j <= grid_size_s0 / 2; j++){
					oblique.x = j*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_sh = 1;
					}
				}
				oblique.z = -(grid_size_sh / 2)*ds;
				oblique.y = y;
				for (int j = -grid_size_s0 / 2; j <= grid_size_s0 / 2; j++){
					oblique.x = j*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){ //too small
						status_sh = 1;
					}
				}
				if (status_sh == 1){// there is an interception
					grid_size_sh += 2; // increase the size
				}
			}
			//take into account zero

			grid_size_sh++;
			grid_size_s0++;

			//allocate memory for complex matrix d0 and dh, the wavefields for incident and diffracted waves
			complex *dh = (complex*)malloc(sizeof(complex)*grid_size_s0*grid_size_sh);
			complex *d0 = (complex*)malloc(sizeof(complex)*grid_size_s0*grid_size_sh);
			double *phase = (double*)malloc(sizeof(double)*grid_size_s0*grid_size_sh);
			vector *poynting = (vector*)calloc(grid_size_s0*grid_size_sh, sizeof(vector));
			boundary *bound_d0 = (boundary*)malloc(sizeof(boundary)*grid_size_sh);
			boundary *bound_dh = (boundary*)malloc(sizeof(boundary)*grid_size_s0);

			//initilize the wavefield and find the boundary
			//found boundary conditions
			for (int i = 0; i < grid_size_sh; i++){
				status_s0 = 0;
				oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
				oblique.y = y;
				for (int j = 0; j < grid_size_s0; j++){
					oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){
						if (status_s0 == 0){
							bound_d0[i].back = j;
							status_s0++;
						}
						else if (status_s0 == 1){
							bound_d0[i].front = j;
							status_s0++;
						}
						else {
							bound_d0[i].front = j;
						}
					}
				}
				bound_d0[i].num = status_s0;
				if (status_s0 == 1){
					bound_d0[i].front = bound_d0[i].back;
				}
				//      printf("back = %d    front = %d\n", bound_d0[i].back, bound_d0[i].front);
			}
			for (int j = 0; j < grid_size_s0; j++){
				status_sh = 0;
				oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
				oblique.y = y;
				for (int i = 0; i < grid_size_sh; i++){
					oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (crystal(cartesian) == 1){
						if (status_sh == 0){
							bound_dh[j].back = i;
							status_sh++;
						}
						else if (status_sh == 1){
							bound_dh[j].front = i;
							status_sh++;
						}
						else {
							bound_dh[j].front = i;
						}
					}
				}
				bound_dh[j].num = status_sh;
				if (status_sh == 1){
					bound_dh[j].front = bound_dh[j].back;
				}
			}
			//  printf("boundary is defined\n");
			//initilize dh field and the phase due to a displacement
			for (int i = 0; i < grid_size_sh; i++){
				for (int j = 0; j < grid_size_s0; j++){
					dh[i*grid_size_s0 + j] = 0.0;
					oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
					oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
					oblique.y = y;
					v = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);///must use v for DISPLACEMENT macro
					phase[i*grid_size_s0 + j] = vec_dot(h,DISPLACEMENT);
				}
			}
			//  printf("dh and phase initilization is done\n");
			//initilize d0 field, calculated from the beam function
			for (int i = 0; i < grid_size_sh; i++){
				if (bound_d0[i].num == 0){
					for (int j = 0; j < grid_size_s0; j++){
						d0[i*grid_size_s0 + j] = 0.0; //outside the particle
					}
				}
				else {
					for (int j = 0; j < grid_size_s0; j++){
						if (j >= bound_d0[i].back && j <= bound_d0[i].front){
							oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
							oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
							oblique.y = y;
							cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
							d0[i*grid_size_s0 + j] = beam(wave_number, unit_s0, beam_center, cartesian);
						}
						else {
							d0[i*grid_size_s0 + j] = 0.0; //outside the particle
						}
					}
				}
			}

			//  printf("Initilization finished\n");
			//iteratively updating both fields

			complex *integral_s0 = (complex*)calloc(grid_size_s0, sizeof(complex));
			complex *integral_sh = (complex*)calloc(grid_size_sh, sizeof(complex));
			complex *dh_tmp = (complex*)calloc(grid_size_s0, sizeof(complex));

			const complex chi_h = g_chi_h;
			const complex chi_0 = g_chi_0;

			const complex ch = I*wave_number / 2.0*chi_h;
			const complex c0 = I*wave_number / 2.0*chi_0;
			const complex w0 = (1.0 - sq(vec_abs(kh) / wave_number)) + chi_0;
			const complex cw = I*wave_number*w0 / 2.0;
			complex cp, gp;
			double s0, sh, res = 1e10, avg = 0.0;
			int k = 0;

			//  printf("start iteration\n");
			while (res > g_tolerance){
				k++;
				res = 0.0;
				avg = 0.0;
				// updating d0
				for (int i = 0; i < grid_size_sh; ++i){
					if (bound_d0[i].num > 1){
						integral_s0[bound_d0[i].back] = 0.0;
						for (int j = bound_d0[i].back + 1; j <= bound_d0[i].front; ++j){
							s0 = -(grid_size_s0 - 1)*ds*0.5 + (j - 1)*ds;
							integral_s0[j] = integral_s0[j - 1] + ch*(dh[i*grid_size_s0 + j] + dh[i*grid_size_s0 + j - 1])*0.5*cexp(-c0*s0)*(1.0 - cexp(-c0*ds)) / c0;
							s0 += ds;
							d0[i*grid_size_s0 + j] = d0[i*grid_size_s0 + bound_d0[i].back] * cexp((j - bound_d0[i].back)*ds*c0) + cexp(c0*s0)*integral_s0[j];
						}
					}
				}
				//updating dh
				for (int j = 0; j < grid_size_s0; ++j){
					if (bound_dh[j].num > 1){
						integral_sh[bound_dh[j].back] = 0;
						for (int i = bound_dh[j].back + 1; i <= bound_dh[j].front; ++i){
							sh = -(grid_size_sh - 1)*ds*0.5 + (i - 1)*ds;
							cp = I*phase[(i - 1)*grid_size_s0 + j];
							gp = I*(phase[i*grid_size_s0 + j] - phase[(i - 1)*grid_size_s0 + j]) / ds;
							integral_sh[i] = integral_sh[i - 1] + (d0[i*grid_size_s0 + j] + d0[(i - 1)*grid_size_s0 + j])*0.5*cexp(-cp)*(1.0 - cexp(-(cw + gp)*ds))*cexp(-cw*sh) / (cw + gp);
							sh += ds;
							dh[i*grid_size_s0 + j] = ch*cexp(I*phase[i*grid_size_s0 + j] + cw*sh)*integral_sh[i];
						}
						res += cabs(cpow(dh_tmp[j] - dh[bound_dh[j].front*grid_size_s0 + j], 2));
						avg += cabs(cpow(dh[bound_dh[j].front*grid_size_s0 + j], 2));
						dh_tmp[j] = dh[bound_dh[j].front*grid_size_s0 + j];
					}
				}
				if (avg == 0.0){
					res = 0.0;
				}
				else {
					res = sqrt(res/avg);
				}
			}
			if(g_saving_wavefield){
				for (int i = 0; i < grid_size_sh; i++){
					if (bound_d0[i].num > 0){
						for (int j = bound_d0[i].back; j <= bound_d0[i].front; j++){
							poynting[i*grid_size_s0 + j] = vec_plus(vec_num(sq(cabs(dh[i*grid_size_s0 + j])), kh), vec_num(sq(cabs(d0[i*grid_size_s0 + j])), k0));
						}
					}
				}
				double minSliceX=INFINITY,maxSliceX=-INFINITY,minSliceZ=INFINITY,maxSliceZ=-INFINITY;
				float*array=(float*)poynting;///re-use poynting array to store wavefield data; assumes vector is 6x as large as float

				for (int i = 0; i < grid_size_sh; i++){
					if (bound_d0[i].num > 0){
						for (int j = bound_d0[i].back; j <= bound_d0[i].front; j++){
							oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
							oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
							oblique.y = y;
							cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);

							const int index = i*grid_size_s0+j;///
							*array++ = poynting[index].x;///first element MUST be poynting.x so we don't corrupt the first poynting element before we use it
							*array++ = poynting[index].z;
							*array++ = cartesian.x;
							*array++ = cartesian.z;
							*array++ = sq(cabs(dh[index]));
							*array++ = sq(cabs(d0[index]));

							if(cartesian.x<minSliceX)minSliceX=cartesian.x;
							else if(cartesian.x>maxSliceX)maxSliceX=cartesian.x;
							if(cartesian.z<minSliceZ)minSliceZ=cartesian.z;
							else if(cartesian.z>maxSliceZ)maxSliceZ=cartesian.z;

							///fprintf(output_dh, "%f      %f      %f      %f\n", cartesian.x, cartesian.y, cartesian.z, sq(cabs(dh[i*grid_size_s0 + j])));
							///fprintf(output_d0, "%f      %f      %f      %f\n", cartesian.x, cartesian.y, cartesian.z, sq(cabs(d0[i*grid_size_s0 + j])));
							///fprintf(output_poynting, "%f      %f      %f      %f      %f\n", cartesian.x, cartesian.y, cartesian.z, poynting[i*grid_size_s0 + j].x, poynting[i*grid_size_s0 + j].z);
						}
					}
				}
				if(minSliceX<minFrameX)minFrameX=minSliceX;
				if(maxSliceX>maxFrameX)maxFrameX=maxSliceX;
				if(minSliceZ<minFrameZ)minFrameZ=minSliceZ;
				if(maxSliceZ>maxFrameZ)maxFrameZ=maxSliceZ;
				char name[32];
				sprintf(name,"slice%d",c);
				h5write_dataset((float*)poynting,(vector*)array-poynting,6,name,wavefield_group);
				hid_t dataset_id=h5(name,wavefield_group);
				h5set_double(grid_size_s0,"x pixels",dataset_id);
				h5set_double(grid_size_sh,"z pixels",dataset_id);
				h5set_double(minSliceX,"x min",dataset_id);
				h5set_double(maxSliceX,"x max",dataset_id);
				h5set_double(minSliceZ,"z min",dataset_id);
				h5set_double(maxSliceZ,"z max",dataset_id);
				h5set_double_(y,"y",dataset_id);
			}
			/*
			for (int i = 0; i < grid_size_sh; i++){
				if (bound_d0[i].num > 0){
					oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + bound_d0[i].front*ds;
					oblique.z = -(grid_size_sh - 1)*ds / 2.0 + i*ds;
					oblique.y = y;
						cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					//fprintf(output_boundary_d0, "%f      %f      %f      %e\n", cartesian.x, cartesian.y, cartesian.z, pow(cabs(d0[i*grid_size_s0 + bound_d0[i].front]), 2.0));
				}
			}
			for (int j = 0; j < grid_size_s0; j++){
				if (bound_dh[j].num > 0){
					oblique.x = -(grid_size_s0 - 1)*ds / 2.0 + j*ds;
					oblique.z = -(grid_size_sh - 1)*ds / 2.0 + bound_dh[j].front*ds;
					oblique.y = y;
					sh = (bound_dh[j].front - bound_dh[j].back)*ds;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					//fprintf(output_boundary_dh, "%f      %f      %f      %e\n", cartesian.x, cartesian.y, cartesian.z, pow(cabs(dh[bound_dh[j].front*grid_size_s0 + j]), 2.0));
				}
			}
			*/
			//propagate to farfield detector

			int flag = 0;
			complex val_temp;
			vector pos_temp;
			for (int k = 0; k < grid_size_s0; k++){
				if (bound_dh[k].num > 0){
					oblique.x = -(grid_size_s0 - 1)*ds*0.5 + k*ds;
					oblique.z = -(grid_size_sh - 1)*ds*0.5 + bound_dh[k].front*ds;
					oblique.y = y;
					cartesian = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					if (flag){ //boundary points more than two
						v = vec_num(0.5, vec_plus(pos_temp, cartesian));///v = pos_avg for DISPLACEMENT macro
						const vector pos_dif = vec_minus(cartesian, pos_temp);
						const complex val_avg = (val_temp + dh[bound_dh[k].front*grid_size_s0 + k])*0.5;
						complex factor = g_dy*sqrt(sq(vec_abs(pos_dif)) - sq(vec_dot(pos_dif, kh)) / sq(vec_abs(kh)))*cexp(I*vec_dot(kh, v) - I*vec_dot(h, DISPLACEMENT))*val_avg;
						///if device_count > 0 offload the calculation to a GPU, otherwise perform it on the CPU
						if(device_count)kernelCuda(cuda_grid_size, cuda_block_size, cuda_det_pix, cuda_det, creal(factor), cimag(factor), v, wave_number, lambda);
						else for(int i = 0; i < det_size; ++i){
							const double distance = vec_dist(det_pix[i], v);
							double exponential[2];///manually perform the cexp() function on distance*wave_number*I for performance
							sincos(distance*wave_number,exponential+1,exponential);
							const complex temp = (I / (lambda*distance))*factor*(exponential[0]+exponential[1]*I);
							double*const type_punned_complex=(double*)(det+i);///treat the complex det[i] as an array of two doubles
#pragma omp atomic          ///atomic constructs don't work with complex numbers, so this type-punning is necessary
							type_punned_complex[0] -= creal(temp);///manually perform the real half of the complex subtraction atomically
#pragma omp atomic          ///these operations need to be atomic because multiple threads are iterating over this array simultaneously
							type_punned_complex[1] -= cimag(temp);///manually perform the imaginary half of the complex subtraction atomically
						}
						val_temp = dh[bound_dh[k].front*grid_size_s0 + k];
						pos_temp.x = cartesian.x;
						pos_temp.y = cartesian.y;
						pos_temp.z = cartesian.z;
					}
					else {
						val_temp = dh[bound_dh[k].front*grid_size_s0 + k];
						pos_temp.x = cartesian.x;
						pos_temp.y = cartesian.y;
						pos_temp.z = cartesian.z;
						flag = 1;
					}
				}
			}
			free(dh);
			free(d0);
			free(phase);
			free(poynting);
			free(bound_d0);
			free(bound_dh);
			free(integral_s0);
			free(integral_sh);
			te_free(xexpr);
			te_free(yexpr);
			te_free(zexpr);
		}
#ifdef CUDA///copy results from GPU and free used GPU memory
		memcpyDeviceToHostCuda(det, cuda_det, det_para.row*det_para.column*sizeof*det);
		freeCuda(cuda_det);
		freeCuda(cuda_det_pix);
#endif
		//write to detector file
		float*data = (float*)det;///re-use the det memory block to store results
		for(int i=0; i<det_size; ++i){
			data[i] = det_para.pix_size*det_para.pix_size*sq(cabs(det[i]));
		}
		char name[24];
		sprintf(name,"det%d",main_index);
		h5write_dataset(data,det_para.row,det_para.column,name,detector_group);
		hid_t id=h5(name,detector_group);
#pragma omp critical(OMP_CRIT_HDF5)///native HDF5 call must be synchronized
		{H5LTset_attribute_double(id,".","Beam center",&beam_center.x,3);}
		h5set_double_(dth,"Deviation angle",id);

		if(g_saving_wavefield){
			h5set_double(minFrameX,"x min",wavefield_group);
			h5set_double(maxFrameX,"x max",wavefield_group);
			h5set_double(minFrameZ,"z min",wavefield_group);
			h5set_double(maxFrameZ,"z max",wavefield_group);
			h5set_double_(num_slice,"slices",wavefield_group);
		}
#pragma omp critical(ITERATIONS_COMPLETED_UPDATE)
		{
			if(iterations_completed<=main_index){///if this is the highest-numbered iteration completed...
				iterations_completed=main_index+1;
				h5set_int(iterations_completed,"Iterations Completed",detector_group);
				emit_progress(iterations_completed);
			}
		}
		if(thread_id==0||main_index==0)h5save(detector_group);///the master thread saves progress
		if(main_index==0)emit_completed1st();///shows the first image when it is finished
		/*for (int i = 0; i < det_para.row; i++){
			for (int j = 0; j < det_para.column; j++){
				fprintf(output_det, "%e\t", pow(det_para.pix_size,2)*pow(cabs(det[i*det_para.column + j]), 2));
			}
			fprintf(output_det, "\n");
		}*/
		free(det);
		free(det_pix);
		/*
		fclose(output_dh);
		fclose(output_d0);
		fclose(output_poynting);
		fclose(output_boundary_d0);
		fclose(output_boundary_dh);
		*/
		///fclose(output_det);

		///file_increment++;

		//      printf("file_increment = %d\n", file_increment);
	}
	h5save_(detector_group);//saves and closes the HDF5 file's Detector group
	//fclose(log);
	time_t end;
	time(&end);
	const int seconds = (int)difftime(end, now);
	const int hrs = seconds / 3600;
	const int min = seconds % 3600 / 60;
	const int sec = seconds % 3600 % 60;
	printf("Calculation Time: %d hours %d minutes %d seconds\n", hrs, min, sec);
	halted=false;
}
#ifdef __cplusplus
}
#endif
int get_num_procs_omp(){
#ifdef _OPENMP
	return omp_get_num_procs();
#else
	return 1;
#endif
}
vector vec_unit(vector v){
	const double inverse_magnitude = 1/sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	v.x *= inverse_magnitude;
	v.y *= inverse_magnitude;
	v.z *= inverse_magnitude;
	return v;
}
vector vec_num(double num, vector vec){
	vector temp;
	temp.x = num*vec.x;
	temp.y = num*vec.y;
	temp.z = num*vec.z;
	return temp;
}
vector vec_plus(vector vec1, vector vec2){
	vector temp;
	temp.x = vec1.x + vec2.x;
	temp.y = vec1.y + vec2.y;
	temp.z = vec1.z + vec2.z;
	return temp;
}
vector vec_minus(vector vec1, vector vec2){
	vector temp;
	temp.x = vec1.x - vec2.x;
	temp.y = vec1.y - vec2.y;
	temp.z = vec1.z - vec2.z;
	return temp;
}
vector vec_cross(vector vec1, vector vec2){
	vector temp;
	temp.x = vec1.y*vec2.z - vec1.z*vec2.y;
	temp.y = vec1.z*vec2.x - vec1.x*vec2.z;
	temp.z = vec1.x*vec2.y - vec1.y*vec2.x;
	return temp;
}
double vec_dist(vector vec1, vector vec2){
	return sqrt(sq(vec1.x - vec2.x) + sq(vec1.y - vec2.y) + sq(vec1.z - vec2.z));
}
double vec_abs(vector vec){
	return sqrt(sq(vec.x) + sq(vec.y) + sq(vec.z));
}
double vec_dot(vector vec1, vector vec2){
	return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}
matrix rot_m(vector vec, double dth){
	double c, s, t;
	matrix R;
	c = cos(dth);
	s = sin(dth);
	t = 1.0 - cos(dth);
	R.e11 = t*vec.x*vec.x + c;
	R.e12 = t*vec.x*vec.y - s*vec.z;
	R.e13 = t*vec.x*vec.z + s*vec.y;
	R.e21 = t*vec.x*vec.y + s*vec.z;
	R.e22 = t*vec.y*vec.y + c;
	R.e23 = t*vec.y*vec.z - s*vec.x;
	R.e31 = t*vec.x*vec.z - s*vec.y;
	R.e32 = t*vec.y*vec.z + s*vec.x;
	R.e33 = t*vec.z*vec.z + c;
	return (R);
}
vector rotation(matrix R, vector vec){
	vector u;
	u.x = R.e11*vec.x + R.e12*vec.y + R.e13*vec.z;
	u.y = R.e21*vec.x + R.e22*vec.y + R.e23*vec.z;
	u.z = R.e31*vec.x + R.e32*vec.y + R.e33*vec.z;
	return u;
}
vector oblique_to_cartesian(double l1, vector s1, double l2, vector s2, double l3, vector s3){
	vector cartesian;
	cartesian = vec_plus(vec_num(l1, s1), vec_num(l2, s2));
	cartesian = vec_plus(cartesian, vec_num(l3, s3));
	return cartesian;
}
///converts vector pointing in v direction into a rotation on an arrow pointing in the y direction
quaternion direction_to_quaternion(vector v){
	double vmag=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
	v.x/=vmag;
	v.y/=vmag;
	v.z/=vmag;
	quaternion q={1,v.x,v.y,v.z};
	double qmag=sqrt(q.scalar*q.scalar+q.x*q.x+q.y*q.y+q.z*q.z);
	q.scalar/=qmag;
	q.x/=qmag;
	q.y/=qmag;
	q.z/=qmag;
	return q;
}
///calculates the rotation axis (returns user input if not all zeroes)
vector rot_axis(){
	if(g_rot_axis.x||g_rot_axis.y||g_rot_axis.z)return g_rot_axis;
	vector result;
	result.x=1;
	result.y=1;
	result.z=1;
	return result;
}
///calculates direction of detector's area vector
vector det_direction(){
	vector result;
	result.x=1;
	result.y=1;
	result.z=1;
	return result;
}
complex planewave(double k,vector z_axis,vector r0,vector pos){
	return csqrt(g_flux_density);
}
complex gaussian(double k,vector z_axis,vector r0,vector pos){
	const vector diff = vec_minus(pos, r0);
	const double zp = vec_dot(diff, z_axis);
	const double temp = sq(vec_abs(diff))-zp*zp;
	const double rp = temp<=0?0:sqrt(temp);
	const double ws = g_w0*g_w0*(1+sq(2*zp/(k*g_w0*g_w0)));
	const complex f = -sq(rp)/ws*(1-2*zp*I/(k*g_w0*g_w0));
	const complex g = 1/(1+2*zp*I/(k*g_w0*g_w0));
	return csqrt(g_total_flux/(.5*pi*g_w0*g_w0))*g*cexp(f);
}
complex fzp(double k,vector z_axis,vector r0,vector pos){
	const int n = 1024;
	const double dr = (g_outer_radius - g_inner_radius)/(n-1);
	const vector diff = vec_minus(pos, r0);
	const double zp = g_focal_length + vec_dot(diff, z_axis);
	const double temp = sq(vec_abs(diff)) - sq(vec_dot(diff, z_axis));
	const double rp = temp<=0?0:sqrt(temp);
	complex sum=0;
	for(int i=0;i<n;++i){
		const double r = g_inner_radius+i*dr;
		sum += cexp(.5*k*I*(1/zp-1/g_focal_length)*r*r)*j0(k*rp*r/zp)*r*dr;
	}
	return-csqrt(g_flux_density*g_efficiency)*sum*k*I*cexp(.5*k*rp*rp*I/zp)/zp;
}
bool cuboid(vector pos){
	return fabs(pos.z) <= g_half_height && fabs(pos.y) <= g_half_width && fabs(pos.x) <= g_half_length;
}
bool sphere(vector pos){
	return vec_abs(pos) <= g_radius;
}
bool cylinder(vector pos){
	return fabs(pos.z) <= g_half_height && pos.x*pos.x + pos.y*pos.y <= g_radius*g_radius;
}
bool cone(vector pos){
	if(fabs(pos.z)>g_half_height)return false;
	return pos.x*pos.x + pos.y*pos.y <= (pos.z+g_half_height)/(g_half_height*2)*g_radius*g_radius;
}
bool bicone(vector pos){
	const double z=fabs(pos.z);
	if(z>g_half_height)return false;
	return pos.x*pos.x + pos.y*pos.y <= (1-z/g_half_height)*g_radius*g_radius;
}
bool prism(vector pos){
	if(fabs(pos.z)>g_half_height)return false;
	const double squared_2d_center_distance = pos.x*pos.x + pos.y*pos.y;
	for(int i=0;i<g_edges;++i){
		if(sq(creal(points[i])-pos.x)+sq(cimag(points[i])-pos.y)<squared_2d_center_distance)return false;
	}
	return true;
}
bool pyramid(vector pos){
	if(fabs(pos.z)>g_half_height)return false;
	const double scalar = (pos.z+g_half_height)/(g_half_height*2);
	const double squared_2d_center_distance = pos.x*pos.x + pos.y*pos.y;
	for(int i=0;i<g_edges;++i){
		const complex point = scalar*points[i];
		if(sq(creal(point)-pos.x)+sq(cimag(point)-pos.y)<squared_2d_center_distance)return false;
	}
	return true;
}
bool bipyramid(vector pos){
	const double z=fabs(pos.z);
	if(z>g_half_height)return false;
	const double scalar = 1-z/g_half_height;
	const double squared_2d_center_distance = pos.x*pos.x + pos.y*pos.y;
	for(int i=0;i<g_edges;++i){
		const complex point = scalar*points[i];
		if(sq(creal(point)-pos.x)+sq(cimag(point)-pos.y)<squared_2d_center_distance)return false;
	}
	return true;
}
#ifndef CUDA///if not using CUDA, define all these functions as no-ops
size_t getDeviceCountCuda(){return 0;}
size_t getMaxThreadsPerBlockCuda(){return 1;}
void setDeviceCuda(int i){}
void*mallocCuda(size_t t){return 0;}
void freeCuda(void*_){}
void memsetCuda(void*_,int i,size_t t){}
void kernelCuda(size_t t,size_t _t,vector*_v,void*_,double d1,double d2,vector v,double d3,double d4){}
void memcpyHostToDeviceCuda(void*_,const void*_c,size_t t){}
void memcpyDeviceToHostCuda(void*_,const void*_c,size_t t){}
#endif