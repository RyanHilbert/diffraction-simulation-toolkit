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

#ifdef _MSC_VER
void sincos(double d,double*s,double*c){
	*s=sin(d);
	*c=cos(d);
}
#endif

#include"tif.h"

#ifdef __cplusplus
extern"C"{
#endif///C functions for accessing Qt's signal-slot mechanism, defined in thread.cpp
extern void emit_output(const char*);
extern void emit_maximum(int);
extern void emit_progress(int);
extern void emit_completed1st();
volatile bool g_halted=false;
int g_max_cpu_cores=1;///maximum number of CPU cores to use, configurable by GUI
int g_max_gpus=0;///maximum number of GPUs to use, configurable by GUI
#ifdef __cplusplus
}
#endif

complex planewave(double,vector,vector,vector);
complex gaussian(double,vector,vector,vector);
complex fzp(double,vector,vector,vector);

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
		case HEMISPHERE:crystal=hemisphere;break;
		case CYLINDER:crystal=cylinder;break;
		case CONE:crystal=cone;break;
		case BICONE:crystal=bicone;break;
		case PRISM:crystal=prism;break;
		case PYRAMID:crystal=pyramid;break;
		case BIPYRAMID:crystal=bipyramid;break;
	}
	free(points);
	points=0;
	if(g_edges>0&&(g_shape==PRISM||g_shape==PYRAMID||g_shape==BIPYRAMID)){///initialize the list of comparison points for the calculation in crystal()
		points=(complex*)malloc(g_edges*sizeof*points);
		for(int i=0;i<g_edges;++i){
			points[i]=g_radius*sin(i*2*PI/g_edges)+g_radius*cos(i*2*PI/g_edges)*I;
		}
		complex prev=points[g_edges-1];
		for(int i=0;i<g_edges;++i){
			const complex temp=points[i];
			points[i]=prev+temp;
			prev=temp;
		}
	}
	const double lambda = 12.398/g_energy*1e-4;
	const double wave_number = 2.0*PI/lambda;
	const vector start_position = g_r0;

	const vector unit_s0 = vec_unit(g_unit_s0);
	const vector unit_s0_0 = unit_s0; // initial unit vector

	const vector unit_h = vec_unit(g_unit_h);
	const vector h = vec_num(2.0*PI / g_d_spacing, unit_h);

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
	g_max_cpu_cores<1?1:(g_max_cpu_cores<=omp_get_num_procs()?g_max_cpu_cores:omp_get_num_procs());
	omp_set_dynamic(false);///disables automatic lowering of thread counts
	omp_set_nested(true);///allows nested parallelism
	printf("OpenMP Enabled: %d Threads\n", thread_count);
#else
	1;printf("OpenMP not enabled; limited to one GPU maximum!\n");
#endif
	const size_t device_count = g_max_gpus<=0?0:(g_max_gpus<=getDeviceCountCuda()?g_max_gpus:getDeviceCountCuda());
	const size_t cuda_block_size = device_count?getMaxThreadsPerBlockCuda():1024;
	printf("%d CUDA Devices Enabled\n",device_count);
	////////////////////////////////////////////////////////
	//Note: xyz is crystal coordinates system
	const size_t num_iterations = get_num_iterations();
	emit_maximum(num_iterations);
	size_t iterations_completed = tif_images();
	//hid_t detector_group=h5group("Detector",hdf5file);///opens the HDF5 file's Detector group
	//h5set_int(num_iterations,"Number of Images",detector_group);
	//int iterations_completed=h5int("Iterations Completed",detector_group);
	if(iterations_completed)emit_completed1st();
	emit_progress(iterations_completed);
#pragma omp parallel for schedule(dynamic) num_threads(thread_count)
	for(long long main_index = iterations_completed; main_index < num_iterations; ++main_index)if(!g_halted){///integer counter needed for OpenMP

		const double dth = main_index / (g_scan_type? g_scan_nx*g_scan_ny: g_num_scan) * g_dth_step + g_dth_start;//deviation angle to the exact Bragg angle; scan position first
		const vector beam_center = {
		g_scan_type?
		g_scan_x1 + (g_scan_x2 - g_scan_x1)/g_scan_nx*((main_index % (g_scan_nx*g_scan_ny))% g_scan_nx) + start_position.x:
		g_spiral_c*sqrt(main_index % g_num_scan)*cos(spiral_phi*(main_index % g_num_scan)) + start_position.x,
		g_scan_type?
		g_scan_y1 + (g_scan_y2 - g_scan_y1)/g_scan_ny*((main_index % (g_scan_nx*g_scan_ny))/g_scan_nx) + start_position.y:
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
		const detector det_para = {g_column, g_row, g_pixel_size, g_det_dist};

		//rotate base vectors of detector
		const vector det_x = rotation(R, det_x_0);
		const vector det_y = rotation(R, det_y_0);
		const vector det_z = rotation(R, det_z_0);
		printf("After rotation \ndet_x = (%G, %G, %G)\ndet_y = (%G, %G, %G)\ndet_z = (%G, %G, %G)\n", det_x.x,det_x.y,det_x.z,det_y.x,det_y.y,det_y.z,det_z.x,det_z.y,det_z.z);
		complex *det = (complex*)calloc_(det_para.row*det_para.col, sizeof(complex));
		vector *det_pix = (vector*)calloc_(det_para.row*det_para.col, sizeof(vector)); //position of a pixel on detector
		for (int i = 0; i < det_para.row; i++){
			for (int j = 0; j < det_para.col; j++){
				det_pix[i*det_para.col+j] = vec_plus(vec_num((j - det_para.col/2)*det_para.pix_size, det_x), vec_num((i - det_para.row/2)*det_para.pix_size, det_y));
				det_pix[i*det_para.col+j] = vec_plus(det_pix[i*det_para.col+j], vec_num(det_para.dist, det_z));
				//                det_pix[i*det_para.column+j] = rotation(R, det_pix[i*det_para.column+j]);
			}
		}
		const int det_size = det_para.row*det_para.col;
		const int thread_id = 0
#ifdef _OPENMP
		|omp_get_thread_num()
#endif
		;const size_t cuda_grid_size = (det_size - 1)/cuda_block_size + 1;
		const int device_id = device_count?thread_id%device_count:0;
		vector*cuda_det_pix=0;
		void*cuda_det=0;
		if(device_count){
			setDeviceCuda(device_id);
			cuda_det_pix=(vector*)mallocCuda(det_para.row*det_para.col*sizeof*det_pix);
			cuda_det=mallocCuda(det_para.row*det_para.col*sizeof*det);///actual type is cuDoubleComplex*
			memsetCuda(cuda_det, 0, det_para.row*det_para.col*sizeof*det);
			memcpyHostToDeviceCuda(cuda_det_pix, det_pix, det_para.row*det_para.col*sizeof*det_pix);
		}
		///all printing in one statement eliminates concurrency issues
		printf("\nThread %d:\nBeam position = (%G, %G, %G)\tDeviation angle = %G\n", thread_id, beam_center.x, beam_center.y, beam_center.z, dth);
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
		}/*
		char wavefield_group_name[32];
		hid_t wavefield_group;
		if(g_saving_wavefield){
			sprintf(wavefield_group_name,"det%dfield",main_index);
			wavefield_group=h5group(wavefield_group_name,detector_group);
		}*/
		const size_t num_slice = (y_max - y_min)/dy + 1;
		wavefield*const slices = g_saving_wavefield?(wavefield*)malloc_((num_slice+1)*sizeof*slices):0;
		if(slices){//null-terminate the array
			const wavefield wave = {0};
			slices[num_slice] = wave;
		}
		//printf("\nThread %d:\n%s\nBeam position = (%G, %G, %G)\tDeviation angle = %G\tNumber of slices: %d\n", thread_id, det_file_name,  beam_center.x, beam_center.y, beam_center.z, dth*180.0/pi, num_slice);
		//start calculation
		///printf("Iteration Number: %d\nThread Count: %d\n",main_index,(thread_count-1)/num_iterations+1);
#pragma omp parallel for schedule(dynamic) num_threads((thread_count-1)/num_iterations+1)///thread_count/remaining_iterations rounded-up
		for (long long slice = 0; slice < num_slice; ++slice){
			if(device_count)setDeviceCuda(thread_id%device_count);
			vector v;
			te_variable vars[]={{"x",&v.x,0,0},{"y",&v.y,0,0},{"z",&v.z,0,0}};
			te_expr*xexpr=te_compile(*g_xstr?g_xstr:"0",vars,3,0);
			te_expr*yexpr=te_compile(*g_ystr?g_ystr:"0",vars,3,0);
			te_expr*zexpr=te_compile(*g_zstr?g_zstr:"0",vars,3,0);
			//find the initial guess of the size of the grid in s0 direction
			const double y = y_min + slice*dy;
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
			++grid_size_sh;
			++grid_size_s0;

			//allocate memory for complex matrix d0 and dh, the wavefields for incident and diffracted waves
			complex *dh = (complex*)malloc_(sizeof(complex)*grid_size_s0*grid_size_sh);
			complex *d0 = (complex*)malloc_(sizeof(complex)*grid_size_s0*grid_size_sh);
			double *phase = (double*)malloc_(sizeof(double)*grid_size_s0*grid_size_sh);
			vector *poynting = (vector*)calloc_(grid_size_s0*grid_size_sh,sizeof(vector));
			boundary *bound_d0 = (boundary*)malloc_(sizeof(boundary)*grid_size_sh);
			boundary *bound_dh = (boundary*)malloc_(sizeof(boundary)*grid_size_s0);

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
					v = oblique_to_cartesian(oblique.x, unit_s0, oblique.y, unit_sy, oblique.z, unit_sh);
					const vector displacement = {te_eval(xexpr),te_eval(yexpr),te_eval(zexpr)};
					phase[i*grid_size_s0 + j] = vec_dot(h,displacement);
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

			complex *integral_s0 = (complex*)calloc_(grid_size_s0, sizeof(complex));
			complex *integral_sh = (complex*)calloc_(grid_size_sh, sizeof(complex));
			complex *dh_tmp = (complex*)calloc_(grid_size_s0, sizeof(complex));

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
				res = (avg?sqrt(res/avg):0);
			}
			if(g_saving_wavefield){
				for (int i = 0; i < grid_size_sh; i++){
					if (bound_d0[i].num > 0){
						for (int j = bound_d0[i].back; j <= bound_d0[i].front; j++){
							poynting[i*grid_size_s0 + j] = vec_plus(vec_num(sq(cabs(dh[i*grid_size_s0 + j])), kh), vec_num(sq(cabs(d0[i*grid_size_s0 + j])), k0));
						}
					}
				}
			}
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
						const vector displacement = {te_eval(xexpr),te_eval(yexpr),te_eval(zexpr)};
						complex factor = g_dy*sqrt(sq(vec_abs(pos_dif)) - sq(vec_dot(pos_dif, kh)) / sq(vec_abs(kh)))*val_avg*cexp(I*vec_dot(kh, v) - I*vec_dot(h, displacement));
						///if device_count > 0 offload the calculation to a GPU, otherwise perform it on the CPU
						if(device_count)kernelCuda(cuda_grid_size, cuda_block_size, cuda_det_pix, cuda_det, creal(factor), cimag(factor), v, wave_number, lambda);
						else for(int i = 0; i < det_size; ++i){
							const double distance = vec_dist(det_pix[i], v);
							double exponential[2];///manually perform the cexp() function on distance*wave_number*I for performance
							sincos(distance*wave_number,exponential+1,exponential);
							const complex temp = (I / (lambda*distance))*factor*(exponential[0]+exponential[1]*I);
							double*const type_punned_complex=(double*)(det+i);///treat the complex det[i] as an array of two doubles
#pragma omp atomic///atomic constructs don't work with complex numbers, so this type-punning is necessary
							type_punned_complex[0] -= creal(temp);///manually perform the real half of the complex subtraction atomically
#pragma omp atomic///these operations need to be atomic because multiple threads are iterating over this array simultaneously
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
			if(g_saving_wavefield){
				slices[slice] = rec_grid_interpolation(d0, dh, poynting, grid_size_s0, grid_size_sh, unit_s0, unit_sy, unit_sh, bound_d0, bound_dh, 2);
				if(BIG_ENDIAN)for(size_t i=0,n=slices[slice].col*slices[slice].row;i<n;++i){
					REVERSE_BYTES(slices[slice].d0[i]);
					REVERSE_BYTES(slices[slice].dh[i]);
					REVERSE_BYTES(slices[slice].poynting[i].x);
					REVERSE_BYTES(slices[slice].poynting[i].z);
				}
			}
			free(d0);
			free(dh);
			free(poynting);
			free(phase);
			free(bound_d0);
			free(bound_dh);
			free(integral_s0);
			free(integral_sh);
			te_free(xexpr);
			te_free(yexpr);
			te_free(zexpr);
		}
		///copy results from GPU and free used GPU memory
		if(device_count){
			memcpyDeviceToHostCuda(det,cuda_det,det_para.row*det_para.col*sizeof*det);
			freeCuda(cuda_det);
			freeCuda(cuda_det_pix);
		}
		//write to detector file
		float*data = (float*)det;///re-use the det memory block to store results
		for(int i=0; i<det_size; ++i) data[i] = sq(det_para.pix_size)*sq(cabs(det[i]));

		if(g_saving_wavefield){
			tif_write_wavefield(main_index,slices);
			for(size_t i=0;i<num_slice;++i){
				free(slices[i].d0);
				free(slices[i].dh);
				free(slices[i].poynting);
			}
		}
		tif_write_ifd(main_index,data);
#pragma omp critical(ITERATIONS_COMPLETED_UPDATE)
		{
			if(iterations_completed <= main_index){///if this is the highest-numbered iteration completed...
				emit_progress(iterations_completed = main_index + 1);
			}
		}
		if(main_index==0)emit_completed1st();///shows the first image when it is finished
		free(det);
		free(det_pix);
		free(slices);
	}
	time_t end;
	time(&end);
	const int seconds = (int)difftime(end, now);
	const int hrs = seconds / 3600;
	const int min = seconds % 3600 / 60;
	const int sec = seconds % 3600 % 60;
	printf("Calculation Time: %d hours %d minutes %d seconds\n", hrs, min, sec);
	g_halted=false;
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
size_t get_num_iterations(){
	return g_scan_type? g_scan_nx*(size_t)g_scan_ny*g_num_angle: g_num_scan*(size_t)g_num_angle;
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
	const double c = cos(dth);
	const double s = sin(dth);
	const double t = 1.0 - cos(dth);
	matrix R;
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
	k=k;z_axis=z_axis;r0=r0;pos=pos;//quiets compiler warnings
	return csqrt(g_flux_density);
}
complex gaussian(double k,vector z_axis,vector r0,vector pos){
	const vector diff = vec_minus(pos, r0);
	const double zp = vec_dot(diff, z_axis);
	const double temp = sq(vec_abs(diff))-zp*zp;
	const double rp = temp<=0?0:sqrt(temp);
	const double ws = g_w0*g_w0*(1+sq(2*zp/(k*g_w0*g_w0)));
	const complex f = -sq(rp)/ws*(1.-2*zp*I/(k*g_w0*g_w0));
	const complex g = 1./(1.+2*zp*I/(k*g_w0*g_w0));
	return csqrt(g_total_flux/(.5*PI*g_w0*g_w0))*g*cexp(f);
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
bool hemisphere(vector pos){
	return pos.z <= 0 && vec_abs(pos) <= g_radius;
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
/*if not using CUDA, define all these functions as no-ops
int getDeviceCountCuda(){return 0;}
int getMaxThreadsPerBlockCuda(){return 1;}
void setDeviceCuda(int i){i=i;}
void*mallocCuda(size_t t){t=t;return 0;}
void freeCuda(void*_){_=_;}
void memsetCuda(void*_,int i,size_t t){_=_;i=i;t=t;}
void kernelCuda(size_t t,size_t _t,vector*_v,void*_,double d1,double d2,vector v,double d3,double d4){t=t;_t=_t;_v=_v;_=_;d1=d1;d2=d2;v=v;d3=d3;d4=d4;}
void memcpyHostToDeviceCuda(void*_,const void*_c,size_t t){_=_;_c=_c;t=t;}
void memcpyDeviceToHostCuda(void*_,const void*_c,size_t t){_=_;_c=_c;t=t;}
*/
vector matrix_vector(matrix mat, vector vec){
	vector u;
	u.x = mat.e11*vec.x + mat.e12*vec.y + mat.e13*vec.z;
	u.y = mat.e21*vec.x + mat.e22*vec.y + mat.e23*vec.z;
	u.z = mat.e31*vec.x + mat.e32*vec.y + mat.e33*vec.z;
	return (u);
}
matrix assemble_matrix(vector v1, vector v2, vector v3, int mode){
	matrix mat;
	if(mode){
		mat.e11 = v1.x;
		mat.e12 = v1.y;
		mat.e13 = v1.z;
		mat.e21 = v2.x;
		mat.e22 = v2.y;
		mat.e23 = v2.z;
		mat.e31 = v3.x;
		mat.e32 = v3.y;
		mat.e33 = v3.z;
	}else{
		mat.e11 = v1.x;
		mat.e21 = v1.y;
		mat.e31 = v1.z;
		mat.e12 = v2.x;
		mat.e22 = v2.y;
		mat.e32 = v2.z;
		mat.e13 = v3.x;
		mat.e23 = v3.y;
		mat.e33 = v3.z;
	}
	return mat;
}
matrix inverse_matrix(matrix mat){
	const double a11 =  (mat.e22*mat.e33 - mat.e23*mat.e32);
	const double a21 = -(mat.e21*mat.e33 - mat.e23*mat.e31);
	const double a31 =  (mat.e21*mat.e32 - mat.e22*mat.e31);
	const double a12 = -(mat.e12*mat.e33 - mat.e13*mat.e32);
	const double a22 =  (mat.e11*mat.e33 - mat.e13*mat.e31);
	const double a32 = -(mat.e11*mat.e32 - mat.e12*mat.e31);
	const double a13 =  (mat.e12*mat.e23 - mat.e13*mat.e22);
	const double a23 = -(mat.e11*mat.e23 - mat.e13*mat.e21);
	const double a33 =  (mat.e11*mat.e22 - mat.e12*mat.e21);
	const double det = mat.e11*a11 + mat.e12*a21 + mat.e13*a31;
	matrix inv_mat;
	inv_mat.e11 = a11/det;
	inv_mat.e12 = a12/det;
	inv_mat.e13 = a13/det;
	inv_mat.e21 = a21/det;
	inv_mat.e22 = a22/det;
	inv_mat.e23 = a23/det;
	inv_mat.e31 = a31/det;
	inv_mat.e32 = a32/det;
	inv_mat.e33 = a33/det;
	return inv_mat;
}
pvector vec_p(vector v){
	const pvector p = {v.x,v.z};
	return p;
}
wavefield rec_grid_interpolation(complex*in_d0, complex*in_dh, vector*in_poynting, int size_x, int size_z, vector s0, vector sy, vector sh, boundary *boundary_d0, boundary *boundary_dh, int interp_order){
	boundary_dh = boundary_dh;//quiets compiler warning
	const vector zeroes = {0,0,0};
	wavefield rec_grid = {0};
	complex d0[4], dh[4];
	vector poynting[4]={{0}},corner[4]={{0}};
	printf("size_x = %d\tsize_z = %d\n", size_x, size_z);
	corner[0] = vec_num(0, s0);
	corner[1] = vec_num(1.0, s0);
	corner[2] = vec_num(1.0, sh);
	corner[3] = vec_plus(corner[1], corner[2]);
	double min_x = 0, max_x = 0, min_z = 0, max_z = 0;
	for (int i = 1; i < 4; i++){
		if (corner[i].x < min_x) min_x = corner[i].x;
		if (corner[i].x > max_x) max_x = corner[i].x;
		if (corner[i].z < min_z) min_z = corner[i].z;
		if (corner[i].z > max_z) max_z = corner[i].z;
	}
	const double dx = (max_x - min_x)/interp_order;
	const double dz = (max_z - min_z)/interp_order;
	printf("dx = %f\tdz = %f\n", dx, dz);
	//find the xz range in cartesian coordinates that accomodate the particle
	int flag = 0, num_x, num_z;
	for (int i = 0; i < size_z; i++){
		if (boundary_d0[i].num > 0){
			for (int j = boundary_d0[i].back; j <= boundary_d0[i].front; j++){
				const vector tmp = vec_plus(vec_num(j, s0), vec_num(i, sh));
				if (flag == 0){
					min_x = tmp.x;
					max_x = tmp.x;
					min_z = tmp.z;
					max_z = tmp.z;
					flag = 1;
				}else{
					if (tmp.x < min_x) min_x = tmp.x;
					if (tmp.x > max_x) max_x = tmp.x;
					if (tmp.z < min_z) min_z = tmp.z;
					if (tmp.z > max_z) max_z = tmp.z;
				}
			}
		}
	}
	if(!flag) return rec_grid;
	else{
		num_x = (max_x - min_x)/dx + 1;
		num_z = (max_z - min_z)/dz + 1;
	}
	printf("num_x = %d\tnum_z = %d\n", num_x, num_z);
	float *rec_grid_d0 = (float*)malloc_(num_x*num_z*sizeof*rec_grid_d0);
	float *rec_grid_dh = (float*)malloc_(num_x*num_z*sizeof*rec_grid_dh);
	pvector *rec_grid_poynting = (pvector*)malloc_(num_x*num_z*sizeof*rec_grid_poynting);
	if (rec_grid_d0 == NULL || rec_grid_dh == NULL || rec_grid_poynting == NULL){
		return rec_grid;
	}
	matrix mat = assemble_matrix(s0, sy, sh, 0);
	matrix inv_mat = inverse_matrix(mat);
	//printf("matrix inverse done!\n");
	vector oblique_vec, cartesian_vec;
	for (int i = 0; i < num_z; i++){
		for (int j = 0; j < num_x; j++){
			cartesian_vec.x = min_x + j*dx;
			cartesian_vec.y = 0.0;
			cartesian_vec.z = min_z + i*dz;
			oblique_vec = matrix_vector(inv_mat, cartesian_vec);
			const int index_x = (int) floor(oblique_vec.x);
			const int index_z = (int) floor(oblique_vec.z);
			const double rx = oblique_vec.x - index_x;
			const double rz = oblique_vec.z - index_z;
			if (index_x < size_x && index_z < size_z && index_x >= 0 && index_z >= 0){
				d0[0] = in_d0[index_z*size_x + index_x];
				dh[0] = in_dh[index_z*size_x + index_x];
				poynting[0] = in_poynting[index_z*size_x + index_x];
			}
			else{
				d0[0] = 0.0;
				dh[0] = 0.0;
				poynting[0] = zeroes;
			}
			if (index_x + 1 < size_x && index_z < size_z && index_x + 1 >= 0 && index_z >=0){
				d0[1] = in_d0[index_z*size_x + index_x + 1];
				dh[1] = in_dh[index_z*size_x + index_x + 1];
				poynting[1] = in_poynting[index_z*size_x + index_x + 1];
			}
			else{
				d0[1] = 0.0;
				dh[1] = 0.0;
				poynting[0] = zeroes;
			}
			if (index_z + 1 < size_z && index_x < size_x && index_z + 1 >= 0 && index_x >= 0){
				d0[2] = in_d0[(index_z + 1)*size_x + index_x];
				dh[2] = in_dh[(index_z + 1)*size_x + index_x];
				poynting[2] = in_poynting[(index_z + 1)*size_x + index_x];
			}
			else{
				d0[2] = 0.0;
				dh[2] = 0.0;
				poynting[0] = zeroes;
			}
			if (index_x + 1 < size_x && index_z + 1 < size_z && index_x + 1 >=0 && index_z + 1 >= 0){
				d0[3] = in_d0[(index_z + 1)*size_x + index_x + 1];
				dh[3] = in_dh[(index_z + 1)*size_x + index_x + 1];
				poynting[3] = in_poynting[(index_z + 1)*size_x + index_x + 1];
			}
			else{
				d0[3] = 0.0;
				dh[3] = 0.0;
				poynting[0] = zeroes;
			}
			rec_grid_d0[i*num_x + j] = sq(cabs((1-rx)*(1-rz)*d0[0] + rx*(1-rz)*d0[1] + (1-rx)*rz*d0[2] + rx*rz*d0[3]));
			rec_grid_dh[i*num_x + j] = sq(cabs((1-rx)*(1-rz)*dh[0] + rx*(1-rz)*dh[1] + (1-rx)*rz*dh[2] + rx*rz*dh[3]));
			rec_grid_poynting[i*num_x + j] = vec_p(vec_plus(vec_plus(vec_num((1-rx)*(1-rz), poynting[0]), vec_num(rx*(1-rz), poynting[1])), vec_plus(vec_num((1-rx)*rz, poynting[2]), vec_num(rx*rz, poynting[3]))));
		}
	}
	//printf("Calculation done!\n");
	rec_grid.d0 = rec_grid_d0;
	rec_grid.dh = rec_grid_dh;
	rec_grid.poynting = rec_grid_poynting;
	rec_grid.col = num_x;
	rec_grid.row = num_z;
	return rec_grid;
}
