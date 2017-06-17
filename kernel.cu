#include<cuComplex.h>

#ifdef __cplusplus
extern"C"{
#endif

typedef struct vector{double x,y,z;}vector;

void kernelCuda(size_t, size_t, vector*, void*, double, double, vector, double, double);
int getDeviceCountCuda();
void setDeviceCuda(int);
void freeCuda(void*);
int getMaxThreadsPerBlockCuda();
void*mallocCuda(size_t);
void memsetCuda(void*, int, size_t);
void memcpyHostToDeviceCuda(void*, const void*, size_t);
void memcpyDeviceToHostCuda(void*, const void*, size_t);

#ifdef __cplusplus
}
#endif

__device__  const cuDoubleComplex i = { 0, 1 };

__global__ void kernel(vector* input, cuDoubleComplex* output, cuDoubleComplex factor, vector avg, double wave, double lambda){
	const int id = blockDim.x*blockIdx.x + threadIdx.x;
	vector in = input[id];
	const double distance = sqrt((in.x - avg.x)*(in.x - avg.x) + (in.y - avg.y)*(in.y - avg.y) + (in.z - avg.z)*(in.z - avg.z));
	cuDoubleComplex exponential;
	sincos(wave*distance,&exponential.y,&exponential.x);
	output[id] = cuCsub(output[id], cuCmul(cuCdiv(i, make_cuDoubleComplex(lambda*distance, 0)), cuCmul(factor, exponential)));
}
void kernelCuda(size_t grid_size,size_t block_size,vector*input,void*output,double real,double imag,vector avg,double wave,double lambda){
	kernel<<<grid_size,block_size>>>(input,(cuDoubleComplex*)output,make_cuDoubleComplex(real,imag),avg,wave,lambda);
}
int getDeviceCountCuda(){
	int i;
	cudaGetDeviceCount(&i);
	return i;
}
void setDeviceCuda(int i) {
	cudaSetDevice(i);
}
void freeCuda(void*v) {
	cudaFree(v);
}
int getMaxThreadsPerBlockCuda() {
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,0);
	return prop.maxThreadsPerBlock;
}
void*mallocCuda(size_t size) {
	void*v;
	cudaMalloc(&v,size);
	return v;
}
void memsetCuda(void*ptr, int value, size_t size) {
	cudaMemset(ptr,value,size);
}
void memcpyHostToDeviceCuda(void*dst,const void*src,size_t count) {
	cudaMemcpy(dst,src,count,cudaMemcpyHostToDevice);
}
void memcpyDeviceToHostCuda(void*dst,const void*src,size_t count) {
	cudaMemcpy(dst,src,count,cudaMemcpyDeviceToHost);
}