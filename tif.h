#ifndef TIF_H
#define TIF_H

#include<stdio.h>
#include<stdint.h>
#include<limits.h>
#include"main.h"
#include"err.h"

ASSERT(CHAR_BIT==8)
ASSERT(sizeof(int)>=4)
ASSERT(sizeof(size_t)>=4)
ASSERT(sizeof(float)==4)
ASSERT(sizeof(double)==8)
ASSERT(sizeof(complex)==16)
ASSERT(sizeof(vector)==24)

#define IFD_SIZE (2+11*12+4)

#define TIF_READ(type,var) tif_read_##type(&var)
#define TIF_WRITE(file,id,type,var) tif_write_##type(file,id,var)
#define TIF_PRINT(type,str) tif_print_##type(str)

//used to separate values into bytes for use in char[] array initializations
#define C2(...) (unsigned char)(__VA_ARGS__),(unsigned char)((__VA_ARGS__)>>8),
#define C4(...) (unsigned char)(__VA_ARGS__),(unsigned char)((__VA_ARGS__)>>8),(unsigned char)((__VA_ARGS__)>>16),(unsigned char)((__VA_ARGS__)>>24),

//macros for determining machine endianness (not compile-time constants)
static const unsigned short one=1;
#define LITTLE_ENDIAN (*(unsigned char*)&one)
#define BIG_ENDIAN (!LITTLE_ENDIAN)

//macro to reverse the bytes of a given object (structs and arrays included)
#define REVERSE_BYTES(...) do for(size_t REVERSE_BYTES=0; REVERSE_BYTES < sizeof(__VA_ARGS__)>>1; ++REVERSE_BYTES)\
	((unsigned char*)&(__VA_ARGS__))[REVERSE_BYTES] ^= ((unsigned char*)&(__VA_ARGS__))[sizeof(__VA_ARGS__)-1-REVERSE_BYTES],\
	((unsigned char*)&(__VA_ARGS__))[sizeof(__VA_ARGS__)-1-REVERSE_BYTES] ^= ((unsigned char*)&(__VA_ARGS__))[REVERSE_BYTES],\
	((unsigned char*)&(__VA_ARGS__))[REVERSE_BYTES] ^= ((unsigned char*)&(__VA_ARGS__))[sizeof(__VA_ARGS__)-1-REVERSE_BYTES];\
while(0)

//TIFF data types - cannot use enum because Windows headers use some of these names
#define BYTE 1
#define ASCII 2
#define SHORT 3
#define LONG 4
#define RATIONAL 5
#define SBYTE 6
#define UNDEFINED 7
#define SSHORT 8
#define SLONG 9
#define SRATIONAL 10
#define FLOAT 11
#define DOUBLE 12

enum Tag{
	ImageWidth=256,   ImageLength, BitsPerSample, Compression,
	PhotometricInterpretation=262,
	StripOffsets=273,
	RowsPerStrip=278, StripByteCounts,
	XResolution=282,  YResolution,
	SampleFormat=339
};
#ifdef __cplusplus
extern"C"{
#define isnan std::isnan
#endif

extern FILE* tif;
extern long tif_offset;//offset of the first IFD
extern const char*const tmp_tif_name;

struct rational{
	unsigned long numerator;
	unsigned long denominator;
};
typedef struct{
	int32_t* data;
	int32_t* data_end;
	size_t total_wavefields;
	size_t current_index;
	wavefield current_wavefield;
}wavefield_iterator;

static struct rational recursionalize(float integral,unsigned int depth){
	struct rational result = {(unsigned long)integral,1};
	float fractional = modff(integral,&integral);
	if(!depth||fractional<=0)return result;
	result = recursionalize(1/fractional,depth-1);
	result.numerator   ^= result.denominator;
	result.denominator ^= result.numerator;
	result.numerator   ^= result.denominator;
	result.numerator   += integral*result.denominator;
	result.numerator   &= 0xFFFFFFFFUL;
	return result;
}
#define RETURN(n,d) return result.numerator=n,result.denominator=d,result
static struct rational rationalize(float x){
	struct rational result = {(unsigned long)x,1};
	if(isnan(x)) RETURN(0,0);
	if(x==0) RETURN(0,1);
	x = fabsf(x);
	if(x==INFINITY) RETURN(1,0);
	if(x>=0xFFFFFFFFUL) RETURN(0xFFFFFFFFUL,1);
	if(x<=1.f/0xFFFFFFFFUL) RETURN(1,0xFFFFFFFFUL);
	unsigned depth = 0;
	float difference = x-(unsigned long)x;
	while(true){
		struct rational test = recursionalize(x,++depth);
		float testdif = fabsf(x-test.numerator/(float)test.denominator);
		if(testdif >= difference)return result;
		difference = testdif;
		result = test;
	}
}
#undef RETURN
//functions for reading different data types from the file
//each accepts the location in memory to read into
//each returns zero on success - nonzero on a read error
static bool tif_read_char256(char256*var){
	int len = getc(tif);
	size_t read = fread(*var,1,len==EOF?0:(unsigned char)len,tif);
	(*var)[read] = 0;
	return len != (int)read;
}
static bool tif_read_char(char*var){
	int size = getc(tif);
	int c = getc(tif);
	*var = c;
	return size != 1 || c == EOF;
}
static bool tif_read_bool(bool*var){
	int size = getc(tif);
	int c = getc(tif);
	*var = c;
	return size != 1 || c == EOF;
}
static bool tif_read__Bool(bool*var){
	return tif_read_bool(var);
}
static bool tif_read_int(int*var){
	int size = getc(tif);
	size_t read = fread(var,1,4,tif);
	if(BIG_ENDIAN)REVERSE_BYTES(*var);
	return size != 4 || read != 4;
}
static bool tif_read_double(double*var){
	int size = getc(tif);
	size_t read = fread(var,1,8,tif);
	if(BIG_ENDIAN)REVERSE_BYTES(*var);
	return size != 8 || read != 8;
}
static bool tif_read_complex(complex*var){
	int size = getc(tif);
	size_t read = fread(var,1,16,tif);
	if(BIG_ENDIAN){
		REVERSE_BYTES(*(double*)var);
		REVERSE_BYTES(*(1+(double*)var));
	}
	return size != 16 || read != 16;
}
static bool tif_read_vector(vector*var){
	int size = getc(tif);
	size_t read = fread(var,1,24,tif);
	if(BIG_ENDIAN){
		REVERSE_BYTES(*(double*)var);
		REVERSE_BYTES(*(1+(double*)var));
		REVERSE_BYTES(*(2+(double*)var));
	}
	return size != 24 || read != 24;
}
static void tif_write_char256(FILE*file,unsigned char id,char256 var){
	const unsigned char len = strlen(var);
	fputc(id,file);
	fputc(len,file);
	fwrite(var,1,len,file);
}
static void tif_write_char(FILE*file,unsigned char id,char var){
	fputc(id,file);
	fputc(1,file);
	fputc(var,file);
}
static void tif_write_bool(FILE*file,unsigned char id,bool var){
	fputc(id,file);
	fputc(1,file);
	fputc(var,file);
}
static void tif_write__Bool(FILE*file,unsigned char id,bool var){
	tif_write_bool(file,id,var);
}
static void tif_write_int(FILE*file,unsigned char id,long var){
	if(BIG_ENDIAN)REVERSE_BYTES(var);
	fputc(id,file);
	fputc(4,file);
	fwrite(&var,1,4,file);
}
static void tif_write_double(FILE*file,unsigned char id,double var){
	if(BIG_ENDIAN)REVERSE_BYTES(var);
	fputc(id,file);
	fputc(sizeof var,file);
	fwrite(&var,1,sizeof var,file);
}
static void tif_write_complex(FILE*file,unsigned char id,complex var){
	if(BIG_ENDIAN){
		REVERSE_BYTES(*(double*)&var);
		REVERSE_BYTES(*(1+(double*)&var));
	}
	fputc(id,file);
	fputc(sizeof var,file);
	fwrite(&var,1,sizeof var,file);
}
static void tif_write_vector(FILE*file,unsigned char id,vector var){
	if(BIG_ENDIAN){
		REVERSE_BYTES(var.x);
		REVERSE_BYTES(var.y);
		REVERSE_BYTES(var.z);
	}
	fputc(id,file);
	fputc(sizeof var,file);
	fwrite(&var,1,sizeof var,file);
}
//renames the current temporary file, saving it from future overwriting
static bool tif_save_as(const char* filename){
	if(!tif)return-1;
#ifndef _WIN32
	return!rename(tmp_tif_name,filename);
#endif//things are more complicated on Windows since we can't rename a file while it's open or rename over an existing file
	remove(filename);
	bool result;
#pragma omp critical(TIF)
	{result = !rename(tmp_tif_name,filename) || !fclose(tif) && !rename(tmp_tif_name,filename) && (tif = fopen(filename,"r+b"));}
	return!result;
}
//attempts to open file with the given name and read its data
static FILE* tif_open(const char* filename){
	FILE* file = fopen(filename,"r+b");//attempt to open the file for update
	if(!file)return 0;//return NULL if it could not be opened
	unsigned char buf[256] = {0};//buffer large enough to hold any parameter
	if(fread(buf,1,16,file)!=16||(buf[0]!='I')|(buf[1]!='I')|(buf[2]!=42)|buf[3]){//read header, offset, & resolution
		fclose(file);//abort file opening if an unexpected header is found
		emit_error("Error reading file: ",filename,"\nInvalid header: ",buf);
		return 0;
	}
	if(tif)fclose(tif);//close the current file before replacing it
	tif = file;
	long offset = buf[4]|buf[5]<<8|buf[6]<<16|buf[7]<<24;//read in the 4-byte offset to the first IFD
	int tag = 0;
	while(tag != EOF && (tag = getc(tif)) && tag!=EOF){//read in all parameters
		switch(tag){//generate the following code for each parameter with a non-zero id
			#define FOREACH(id,name,type,var,...) IF(id)(\
			case id:\
				if(TIF_READ(type,var)){\
					tag = EOF;\
					emit_error("Error reading file: ",filename,"\nFailed to load parameter: " name);\
				}\
				continue;\
			)//^begin next loop iteration^
			#include"main.h"//compilation error should occur here on duplicate non-zero parameter IDs
		}//if the tag is unknown, simply read over and ignore it
		const int size = getc(tif);//all parameters begin with a byte indicating their remaining size
		if(size == EOF || fread(buf,1,size,tif) != (size_t)size){
			emit_error("Error reading file: ",filename,"\nFailed to load unknown parameter");
			break;//abort the loop
		}
	}
	if(offset) tif_offset = offset;//if the offset was non-zero, use it
	else{//otherwise calculate it from the current file position
		offset = 4+ftell(tif);//leave some space for our oblivious IFD-writing function
		if(offset&1) ++offset;//make sure the offset is even, as mandated by the TIFF spec
		tif_offset = offset;//finally, set the global offset variable
	}
	fseek(tif,tif_offset+get_num_iterations()*(IFD_SIZE+sizeof(float)*g_row*g_column),SEEK_SET);
	g_saving_wavefield = getc(tif) != EOF;//check if the file contains wavefield data
	return tif;//return the newly-read file
}
//truncate default TIF file and rewrite header, resolution, & parameters
static FILE* tif_update(){
	if(!freopen(tmp_tif_name,"w+b",tif))return 0;
	const char header[]={'I','I',42,0,0,0,0,0};
	fwrite(header,1,8,tif);
	struct rational resolution = rationalize(25400/g_pixel_size);
	if(BIG_ENDIAN){
		REVERSE_BYTES(resolution.numerator);
		REVERSE_BYTES(resolution.denominator);
	}
	fwrite(&resolution.numerator,1,4,tif);
	fwrite(&resolution.denominator,1,4,tif);
	for(size_t i=1;i<=0xFF;++i)switch(i){//write all parameters in order of id
		#define FOREACH(id,name,type,var,...) IF(id)(case id:TIF_WRITE(tif,id,type,var);continue;)
		#include"main.h"//compilation error should occur here on duplicate non-zero parameter IDs
	}
	putc(0,tif);//id of zero signals end of parameter list
	long offset = 4+ftell(tif);//add 4-bytes for junk offset before 1st IFD
	if(offset&1) ++offset;//make even if it would be odd
	tif_offset = offset;
	fflush(tif);
	return tif;
}
static FILE* tif_initialize(){//convenience initialization function to open the default file
	if(tif)return tif;//if a file is already open, do nothing
	if(tif_open(tmp_tif_name))return tif;//otherwise attempt to open the default file
	tif=fopen(tmp_tif_name,"w+b");//if that fails, create the default file
	return tif_update();//call the update function to populate the new file
}
static size_t tif_images(){//returns the number of IFDs (images in the file)
	size_t count = 0;
	unsigned char buf[255*255];//max possible size of parameter section
	#pragma omp critical(TIF)
	do{
		rewind(tif);
		if(fread(buf,1,8,tif)!=8)break;
		if(!buf[4]&!buf[5]&!buf[6]&!buf[7])break;//if the first offset is zero, stop here
		const size_t offset = buf[4]|buf[5]<<8|buf[6]<<16|buf[7]<<24;
		if(fread(buf,1,offset-8,tif)!=offset-8)break;
		count = 1;
		//while there are no read errors AND the next offset is non-zero...
		while(fread(buf,1,IFD_SIZE,tif)==IFD_SIZE && buf[IFD_SIZE-4]|buf[IFD_SIZE-3]|buf[IFD_SIZE-2]|buf[IFD_SIZE-1])++count;
	}while(0);
	return count;
}
static float* tif_read_ifd(size_t index,float*data){
	const size_t data_size = sizeof(float)*g_column*g_row;
	#pragma omp critical(TIF)
	{
		fseek(tif,tif_offset+get_num_iterations()*IFD_SIZE+index*data_size,SEEK_SET);
		if(!data) data = (float*)calloc_(data_size,1);//zero-initialized in-case the read fails
		fread(data,1,data_size,tif);
	}
	return data;
}
static void tif_write_ifd(size_t index,float*data){
	const size_t row = g_row;
	const size_t column = g_column;
	const size_t pixels = row*column;
	const size_t pixbytes = pixels*sizeof(float);
	const size_t ifds = get_num_iterations();
	const unsigned char ifd[IFD_SIZE]={C4(tif_offset+index*sizeof ifd)C2(sizeof ifd/12u)//previous offset, followed by count
		// Tag                Type       Count     Value/Offset
		C2(ImageWidth)     C2(LONG)      C4(1)     C4(column)
		C2(ImageLength)    C2(LONG)      C4(1)     C4(row)
		C2(BitsPerSample)  C2(SHORT)     C4(1)     C4(32)
		C2(Compression)    C2(SHORT)     C4(1)     C4(1)
		C2(PhotometricInterpretation)C2(SHORT)C4(1)C4(1)
		C2(StripOffsets)   C2(LONG)      C4(1)     C4(tif_offset+index*pixbytes+ifds*sizeof ifd)
		C2(RowsPerStrip)   C2(LONG)      C4(1)     C4(row)
		C2(StripByteCounts)C2(LONG)      C4(1)     C4(pixbytes)
		C2(XResolution)    C2(RATIONAL)  C4(1)     C4(8)
		C2(YResolution)    C2(RATIONAL)  C4(1)     C4(8)
		C2(SampleFormat)   C2(SHORT)     C4(1)     C4(3)
	};
	if(BIG_ENDIAN)for(size_t i=0;i<pixels;++i)REVERSE_BYTES(data[i]);
	#pragma omp critical(TIF)
	{
		fseek(tif,tif_offset+index*pixbytes+ifds*sizeof ifd,SEEK_SET);
		fwrite(data,1,pixbytes,tif);
		fseek(tif,tif_offset+index*sizeof ifd-4,SEEK_SET);
		fwrite(ifd,1,sizeof ifd,tif);
		if(!index){
			fseek(tif,4,SEEK_SET);
			fwrite(ifd,1,4,tif);
		}
		fflush(tif);
	}
}
static int tif_seek(unsigned long long pos){//needed for files >2GB on Windows
	if(pos <= LONG_MAX) return fseek(tif,pos,SEEK_SET);
	rewind(tif);
	do fseek(tif,LONG_MAX,SEEK_CUR); while((pos -= LONG_MAX) > LONG_MAX);
	return fseek(tif,pos,SEEK_CUR);
}
static wavefield_iterator get_wavefield(size_t index,wavefield_iterator in){
	if(index && index >= in.total_wavefields) index = in.total_wavefields - 1;
	if(index == in.current_index) return in;
	if(index > in.current_index){
		size_t current_index = in.current_index;
		int32_t*ptr = (int32_t*)in.current_wavefield.d0 - 2;
		while(current_index++ < index) ptr += 4 + 4ULL*ptr[0]*ptr[1];
		const size_t cols = ptr[0], rows = ptr[1];
		if(BIG_ENDIAN){
			REVERSE_BYTES(cols);
			REVERSE_BYTES(rows);
		}
		const wavefield current = {cols, rows, (float*)ptr+2, (float*)ptr+2+cols*rows, (pvector*)(ptr+2+2*cols*rows)};
		const wavefield_iterator result = {in.data, in.data_end, in.total_wavefields, index, current};
		return result;
	}else{
		size_t current_index = in.current_index;
		int32_t*ptr = (int32_t*)in.current_wavefield.d0 - 2;
		while(current_index-- > index) ptr -= 4 + 4ULL*ptr[-2]*ptr[-1];
		const size_t cols = ptr[0], rows = ptr[1];
		if(BIG_ENDIAN){
			REVERSE_BYTES(cols);
			REVERSE_BYTES(rows);
		}
		const wavefield current = {cols, rows, (float*)ptr+2, (float*)ptr+2+cols*rows, (pvector*)(ptr+2+2*cols*rows)};
		const wavefield_iterator result = {in.data, in.data_end, in.total_wavefields, index, current};
		return result;
	}
}
static wavefield_iterator tif_read_wavefield(size_t index){
	int32_t* data = 0;
	size_t count = 0;
	unsigned long long size = 0;
	#pragma omp critical(TIF)
	{
		unsigned long long pos=0;
		fseek(tif,tif_offset+20*index+get_num_iterations()*(IFD_SIZE+sizeof(float)*g_row*g_column),SEEK_SET);
		fread(&pos,1,8,tif);
		fread(&size,1,8,tif);
		fread(&count,1,4,tif);
		if(BIG_ENDIAN){
			REVERSE_BYTES(pos);
			REVERSE_BYTES(size);
			REVERSE_BYTES(count);
		}
		data = (int32_t*)malloc_(size);
		tif_seek(pos);
		fread(data,1,size,tif);
	}
	if(!size){
		const wavefield_iterator result = {0,0,0,0};
		return result;
	}
	int32_t*const end = (int32_t*)((char*)data+size);
	const size_t cols = data[0], rows = data[1];
	if(BIG_ENDIAN){
		REVERSE_BYTES(cols);
		REVERSE_BYTES(rows);
	}
	const wavefield first = {cols, rows, (float*)data+2, (float*)data+2+cols*rows, (pvector*)(data+2+2*cols*rows)};
	const wavefield_iterator result = {data, end, count, 0, first};
	return result;
}
static void tif_write_wavefield(const size_t index,const wavefield*waves){
#pragma omp critical(TIF)
{
	const size_t images = get_num_iterations();
	const long offset = tif_offset + images*(20+IFD_SIZE+sizeof(float)*g_row*g_column);
	unsigned long long pos;
	fseek(tif,offset,SEEK_SET);
	if(fread(&pos,1,8,tif)!=8) pos = offset + 8;
	else if(BIG_ENDIAN)REVERSE_BYTES(pos);
	tif_seek(pos);
	unsigned long long size = 0;
	size_t cols, rows;
	const wavefield*ptr = waves;
	while((cols=ptr->col)|(rows=ptr->row)){
		const size_t pixels = rows*cols;
		if(BIG_ENDIAN)REVERSE_BYTES(cols);
		if(BIG_ENDIAN)REVERSE_BYTES(rows);
		fwrite(&cols,1,4,tif);
		fwrite(&rows,1,4,tif);
		fwrite(ptr->d0,1,pixels*sizeof(float),tif);
		fwrite(ptr->dh,1,pixels*sizeof(float),tif);
		fwrite(ptr->poynting,1,pixels*2*sizeof(float),tif);
		fwrite(&cols,1,4,tif);
		fwrite(&rows,1,4,tif);
		size += 16 + 4*pixels*sizeof(float);
		++ptr;
	}
	unsigned long long newpos = pos + size;
	size_t count = ptr - waves;
	if(BIG_ENDIAN){
		REVERSE_BYTES(pos);
		REVERSE_BYTES(size);
		REVERSE_BYTES(newpos);
		REVERSE_BYTES(count);
	}
	fseek(tif,offset,SEEK_SET);
	fwrite(&newpos,1,8,tif);
	fseek(tif,offset-20*(images-index),SEEK_SET);
	fwrite(&pos,1,8,tif);
	fwrite(&size,1,8,tif);
	fwrite(&count,1,4,tif);
}}
static bool export_detector(const char*filename){
	FILE*file = fopen(filename,"wb");
	if(!file) return 1;
	unsigned char buf[BUFSIZ];
#pragma omp critical(TIF)
	{
	rewind(tif);
	for(size_t i=tif_offset+get_num_iterations()*(IFD_SIZE+sizeof(float)*g_row*g_column);i;){
		size_t bytes = fread(buf,1,i<sizeof buf?i:sizeof buf,tif);
		fwrite(buf,1,bytes,file);
		if(bytes != sizeof buf) break;
		i -= bytes;
	}}
	return fclose(file);
}
static bool export_wavefield(wavefield_iterator it,const char*d0filename,const char*dhfilename,const char*pfilename){
	FILE *d0file=0, *dhfile=0, *pfile=0;
	if(d0filename && !(d0file=fopen(d0filename,"wb"))) return 1;
	if(dhfilename && !(dhfile=fopen(dhfilename,"wb"))) return fclose(d0file), 1;
	if( pfilename && !( pfile=fopen( pfilename,"wb"))) return fclose(d0file), fclose(dhfile), 1;
	char header[]={'I','I',42,0,C4(tif_offset)};
	struct rational resolution = rationalize(25400/g_pixel_size);
	if(BIG_ENDIAN){
		REVERSE_BYTES(resolution.numerator);
		REVERSE_BYTES(resolution.denominator);
	}
	if(d0file)fwrite(header,1,8,d0file),fwrite(&resolution.numerator,1,4,d0file),fwrite(&resolution.denominator,1,4,d0file);
	if(dhfile)fwrite(header,1,8,dhfile),fwrite(&resolution.numerator,1,4,dhfile),fwrite(&resolution.denominator,1,4,dhfile);
	if( pfile)fwrite(header,1,8, pfile),fwrite(&resolution.numerator,1,4, pfile),fwrite(&resolution.denominator,1,4, pfile);
	for(size_t i=1;i<=0xFF;++i)switch(i){//write all parameters in order of id
		#define FOREACH(id,name,type,var,...) IF(id)(case id:if(d0file)TIF_WRITE(d0file,id,type,var);if(dhfile)TIF_WRITE(dhfile,id,type,var);if(pfile)TIF_WRITE(pfile,id,type,var);continue;)
		#include"main.h"//compilation error should occur here on duplicate non-zero parameter IDs
	}
	unsigned long long offset = 0;
	if(d0file)fwrite(&offset,1,tif_offset-ftell(d0file),d0file);
	if(dhfile)fwrite(&offset,1,tif_offset-ftell(dhfile),dhfile);
	if( pfile)fwrite(&offset,1,tif_offset-ftell( pfile), pfile);
	offset = tif_offset + it.total_wavefields*IFD_SIZE;
	for(it=get_wavefield(0,it); 1; it=get_wavefield(it.current_index+1,it)){
		const bool last = it.current_index == it.total_wavefields - 1;
		const wavefield wf = it.current_wavefield;
		const unsigned char ifd[IFD_SIZE]={C2(sizeof ifd/12u)
			// Tag                Type       Count     Value/Offset
			C2(ImageWidth)     C2(LONG)      C4(1)     C4(wf.col)
			C2(ImageLength)    C2(LONG)      C4(1)     C4(wf.row)
			C2(BitsPerSample)  C2(SHORT)     C4(1)     C4(32)
			C2(Compression)    C2(SHORT)     C4(1)     C4(1)
			C2(PhotometricInterpretation)C2(SHORT)C4(1)C4(1)
			C2(StripOffsets)   C2(LONG)      C4(1)     C4(offset)
			C2(RowsPerStrip)   C2(LONG)      C4(1)     C4(wf.row)
			C2(StripByteCounts)C2(LONG)      C4(1)     C4(sizeof(float)*wf.col*wf.row)
			C2(XResolution)    C2(RATIONAL)  C4(1)     C4(8)
			C2(YResolution)    C2(RATIONAL)  C4(1)     C4(8)
			C2(SampleFormat)   C2(SHORT)     C4(1)     C4(3)
			C4(last? 0: tif_offset+(it.current_index+1)*sizeof ifd)
		};
		offset += sizeof(float)*wf.col*wf.row;
		if(d0file)fwrite(ifd,1,sizeof ifd,d0file);
		if(dhfile)fwrite(ifd,1,sizeof ifd,dhfile);
		if( pfile)fwrite(ifd,1,sizeof ifd, pfile);
		if(last)break;
	}
	for(it=get_wavefield(0,it); 1; it=get_wavefield(it.current_index+1,it)){
		const bool last = it.current_index == it.total_wavefields - 1;
		const wavefield wf = it.current_wavefield;
		const size_t size = sizeof(float)*wf.col*wf.row;
		if(d0file)fwrite(wf.d0,1,size,d0file);
		if(dhfile)fwrite(wf.dh,1,size,dhfile);
		if( pfile)for(size_t i=0,n=wf.col*wf.row;i<n;++i){
			pvector pv = wf.poynting[i];
			if(BIG_ENDIAN){
				REVERSE_BYTES(pv.x);
				REVERSE_BYTES(pv.z);
			}
			float magnitude = sqrtf(pv.x*pv.x + pv.z*pv.z);
			if(BIG_ENDIAN)REVERSE_BYTES(magnitude);
			fwrite(&magnitude,1,sizeof(float),pfile);
		}
		if(last)return (d0file?fclose(d0file):0) | (dhfile?fclose(dhfile):0) | (pfile?fclose(pfile):0);
	}
}
/*
static size_t tif_print_char256(char*str){
	int len = getc(tif);
	size_t read = fread(str,1,len==EOF?0:(unsigned char)len,tif);
	return read;
}
static size_t tif_print_char(char*str){
	int size = getc(tif);
	int c = getc(tif);
	*str = '0'+c;
	return 1;
}
static size_t tif_print_bool(char*str){
	int size = getc(tif);
	int c = getc(tif);
	*str = c?'1':'0';
	return 1;
}
static size_t tif_print__Bool(char*str){
	return tif_print_bool(str);
}
static size_t tif_print_int(char*str){
	int size = getc(tif);
	int result = 0;
	size_t read = fread(&result,1,4,tif);
	if(BIG_ENDIAN)REVERSE_BYTES(read);
	return sprintf(str,"%d",result);
}
static size_t tif_print_double(char*str){
	int size = getc(tif);
	double result = 0;
	size_t read = fread(&result,1,8,tif);
	if(BIG_ENDIAN)REVERSE_BYTES(result);
	return sprintf(str,"%f",result);
}
static size_t tif_print_complex(char*str){
	int size = getc(tif);
	complex result = 0;
	size_t read = fread(&result,1,16,tif);
	if(BIG_ENDIAN){
		REVERSE_BYTES(*(double*)&result);
		REVERSE_BYTES(*(1+(double*)&result));
	}
	return sprintf(str,"%f%+fi",*(double*)&result,*(1+(double*)&result));
}
static size_t tif_print_vector(char*str){
	int size = getc(tif);
	vector result = {0,0,0};
	size_t read = fread(&result,1,24,tif);
	if(BIG_ENDIAN){
		REVERSE_BYTES(*(double*)&result);
		REVERSE_BYTES(*(1+(double*)&result));
		REVERSE_BYTES(*(2+(double*)&result));
	}
	return sprintf(str,"%fx %fy %fz",result.x,result.y,result.z);
}
static void tif_print(){
#pragma omp critical(TIF)
{
	char str[65536] = {0};
	size_t offset = 0;
	int tag = 0;
	fseek(tif,16,SEEK_SET);
	while(tag != EOF && (tag = getc(tif)) && tag != EOF){
		switch(tag){
			#define FOREACH(id,name,type,var,...) IF(id)(\
			case id:\
				offset += sprintf(str+offset,"<br>\n" #id " " #type " (" #var "): ");\
				offset += TIF_PRINT(type,str+offset);\
				continue;\
			)
			#include"main.h"
		}
		const int size = getc(tif);//all parameters begin with a byte indicating their remaining size
		char buf[256] = {0};
		if(size == EOF || fread(buf,1,size,tif) != (size_t)size){
			emit_error("Failed to read unknown parameter");
			break;
		}
		offset += sprintf(str+offset,"\n%d Unknown Parameter: %s",tag,buf);
	}
	emit_error(str);
}}
*/
#ifdef __cplusplus
#undef isnan
}
#endif
#endif