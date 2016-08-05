#ifndef H5_H
#define H5_H
#include<hdf5.h>
#include<hdf5_hl.h>
#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<stdlib.h>
#include<stdbool.h>
#ifdef __cplusplus
#include<string>
extern"C"{
#endif
//main function for navigating HDF5 object hierarchies; accepts a group or file id
//returns the id of the HDF5 group or dataset at the given location with the given name
//will return a negative value and trigger HDF5 error-handling if the object does not exist
//ids returned by this function must eventually be closed to prevent resource leaks
//valid closing functions include H5Oclose and functions defined in this header ending in an underscore_
static inline hid_t h5(const char*name,hid_t loc){
	hid_t result=-1;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5Oopen(loc,name,H5P_DEFAULT);}
	return result;
}
//performs the above function and closes the passed group id (must not be a file id)
static inline hid_t h5_(const char*name,hid_t loc){
	hid_t result=h5(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}

//opens the HDF5 file with the given name or creates it if it does not exist
//returns either an id to the opened/created file or a negative value if a failure occurs
//valid ids returned by this function must eventually be closed with the H5Fclose(hid_t) function
//ids returned by this function may NOT be closed by functions in this header ending in an underscore_
static inline hid_t h5file(const char*name){
	FILE*f=fopen(name,"r+");
	hid_t result=-1;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=(f?fclose(f):1)?H5Fcreate(name,H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT):H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);}
	return result;
}

//returns an id to a new file with the given name, overwriting an identically named file if it exists
//the same closing restrictions as the above function apply
static inline hid_t h5overwrite(const char*name){
	hid_t result=-1;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);}
	return result;
}

//saves the file associated with the given id, which may identify any object in the file
static inline bool h5save(hid_t id){
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5Fflush(id,H5F_SCOPE_LOCAL)>=0;}
	return result;
}
//performs the above function and then closes the given id (must not be a file id)
static inline bool h5save_(hid_t id){
	bool result=h5save(id);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(id);}
	return result;
}

//checks if an object with the given name exists at the given location id
static inline bool h5exists(const char*name,hid_t loc){
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTpath_valid(loc,name,true);}
	return result;
}
//performs the above function and then closes the given location id (must be a group or dataset id)
static inline bool h5exists_(const char*name,hid_t loc){
	bool result=h5exists(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}

//returns an id for the group with the given name at the given location, creating it if it does not exist
static inline hid_t h5group(const char*name,hid_t loc){
	hid_t result=-1;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTpath_valid(loc,name,true)?H5Oopen(loc,name,H5P_DEFAULT):H5Gcreate1(loc,name,0);}
	return result;
}
//performs the above function and then closes the given group id (must not be a file id)
static inline hid_t h5group_(const char*name,hid_t loc){
	hid_t result=h5group(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Gclose(loc);}
	return result;
}


//toggles HDF5 error-handling on and off; returns true if handling is enabled after calling the function
static inline bool h5toggle_errors(){
	static H5E_auto1_t func;
	static void*data;
	static bool enabled=true;
#pragma omp critical(OMP_CRIT_HDF5)
	{
		if(enabled){
			H5Eget_auto1(&func,&data);
			H5Eset_auto1(0,0);
		}else{
			H5Eset_auto1(func,data);
		}
	}
	return enabled=!enabled;
}


//struct representing a 2D dataset of floats
//must free data pointer with free(dataset.data) if not using C++
struct h5dataset{
	float*const data;
	const size_t rows,columns;
#ifdef __cplusplus
	h5dataset():data(0),rows(0),columns(0){}
	h5dataset(float*const data,size_t rows,size_t columns):data(data),rows(rows),columns(columns){}
	float*operator[](size_t row){return row*columns+data;}//allows 2D array syntax in C++: dataset[row][col]
	~h5dataset(){free(data);}//allows automatic resource clean-up on destruction in C++; no need to free()
#endif
};
//retrieves the 2D dataset of floats with the given name at the given location (group or file id)
static inline struct h5dataset h5dataset(const char*name,hid_t loc){
	struct h5dataset zeroes={0,0,0};
	bool failure=false;
#pragma omp critical(OMP_CRIT_HDF5)
	{if(!H5LTfind_dataset(loc,name))failure=true;}
	if(failure)return zeroes;
	hsize_t dims[2]={0,0};
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTget_dataset_info(loc,name,dims,0,0)<0)failure=true;}
	if(failure)return zeroes;
	float*data=(float*)malloc(dims[0]*dims[1]*sizeof*data);
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTread_dataset_float(loc,name,data)<0)failure=true;}
	if(failure){free(data);return zeroes;}
	return
		#ifndef __cplusplus
			(struct h5dataset)
		#endif
	{data,(size_t)dims[0],(size_t)dims[1]};
}
//performs the above function and then closes the given group id (must not be a file id)
static inline struct h5dataset h5dataset_(const char*name,hid_t loc){
	struct h5dataset result=h5dataset(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Gclose(loc);}
	return result;
}
//creates a 2D dataset of floats with the given parameters at the given location (group or file id)
//returns true if successful or false if an error occured (the dataset cannot already exist)
static inline bool h5write_dataset(const float*data,size_t rows,size_t cols,const char*name,hid_t loc){
	hsize_t dims[2]={rows,cols};
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTmake_dataset_float(loc,name,2,dims,data)>=0;}
	return result;
}
//performs the above function and then closes the given group id (must not be a file id)
static inline bool h5write_dataset_(const float*data,size_t rows,size_t cols,const char*name,hid_t loc){
	bool result=h5write_dataset(data,rows,cols,name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}


//returns the int attribute with the given name attached to the given HDF5 file, group, or dataset
static inline int h5int(const char*name,hid_t loc){
	int result=INT_MIN;
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTfind_attribute(loc,name))H5LTget_attribute_int(loc,".",name,&result);}
	return result;
}
//performs the above function and then closes the given group or dataset id (must not be a file id)
static inline int h5int_(const char*name,hid_t loc){
	int result=h5int(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}
//sets the int attribute with the given name at the given location to the given value
static inline bool h5set_int(int value,const char*name,hid_t loc){
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTset_attribute_int(loc,".",name,&value,1)>=0;}
	return result;
}
//performs the above function and then closes the given group of dataset id (must not be a file id)
static inline bool h5set_int_(int value,const char*name,hid_t loc){
	bool result=h5set_int(value,name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}

//returns the double attribute with the given name attached to the given HDF5 file, group, or dataset
static inline double h5double(const char*name,hid_t loc){
	double result=NAN;
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTfind_attribute(loc,name))H5LTget_attribute_double(loc,".",name,&result);}
	return result;
}
//performs the above function and then closes the given group or dataset id (must not be a file id)
static inline double h5double_(const char*name,hid_t loc){
	double result=h5double(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}
//sets the double attribute with the given name at the given location to the given value
static inline bool h5set_double(double value,const char*name,hid_t loc){
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTset_attribute_double(loc,".",name,&value,1)>=0;}
	return result;
}
//performs the above function and then closes the given group or dataset id (must not be a file id)
static inline bool h5set_double_(double value,const char*name,hid_t loc){
	bool result=h5set_double(value,name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}

//sets the string attribute with the given name at the given location to the given value
static inline bool h5set_string(const char*value,const char*name,hid_t loc){
	bool result=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTset_attribute_string(loc,".",name,value);}
	return result;
}
//performs the above function and then closes the given group or dataset id (must not be a file id)
static inline bool h5set_string_(const char*value,const char*name,hid_t loc){
	bool result=h5set_string(value,name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}
#ifdef __cplusplus
}//end of extern "C" linkage; the following functions have no C implementation

//returns the string attribute with the given name attached to the given HDF5 file, group, or dataset
static inline std::string h5string(const char*name,hid_t loc){
	char result[BUFSIZ];
	result[0]=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTfind_attribute(loc,name))H5LTget_attribute_string(loc,".",name,result);}
	return result;
}
//performs the above function and then closes the given group or dataset id (must not be a file id)
static inline std::string h5string_(const char*name,hid_t loc){
	std::string result=h5string(name,loc);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(loc);}
	return result;
}

//returns the name of the file associated with the given id, which may be represent any object in the file
static inline std::string h5filename(hid_t object){
	char result[BUFSIZ];
	result[0]=0;
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Fget_name(object,result,BUFSIZ);}
	return result;
}
//performs the above function and then closes the given id (only valid for group and dataset ids)
static inline std::string h5filename_(hid_t object){
	std::string result=h5filename(object);
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(object);}
	return result;
}
#endif
#endif