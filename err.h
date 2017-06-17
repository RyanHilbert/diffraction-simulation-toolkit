#ifndef ERR_H
#define ERR_H

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#ifdef __cplusplus
extern"C"{
#endif

extern void*emit_error(const char*str,...);

#define STR_(...) #__VA_ARGS__
#define STR(...) STR_(__VA_ARGS__)
#define LINE STR(__LINE__)

#define emit_error(...) emit_error(__VA_ARGS__,0)

static void*malloc_(size_t size,const char*file,const char*func,const char*line,const char*args){
	if(!size) return 0;
	void*result = malloc(size);
	if(result) return result;
	char buf[24] = {0};
	sprintf(buf,"%lu",(long)size);
	emit_error("Function ",func," failed to allocate ",buf," bytes of memory with malloc(",args,") on line ",line," of file ",file);
	return 0;
}
#define malloc_(...) malloc_(__VA_ARGS__,__FILE__,__func__,LINE,#__VA_ARGS__)

static void*calloc_(size_t num,size_t size,const char*file,const char*func,const char*line,const char*args){
	if(!num|!size) return 0;
	void*result = calloc(num,size);
	if(result) return result;
	char buf[24] = {0};
	sprintf(buf,"%lu",(long)(num*size));
	emit_error("Function ",func," failed to allocate ",buf," bytes of memory with calloc(",args,") on line ",line," of file ",file);
	return 0;
}
#define calloc_(...) calloc_(__VA_ARGS__,__FILE__,__func__,LINE,#__VA_ARGS__)

static void*realloc_(void*ptr,size_t size,const char*file,const char*func,const char*line,const char*args){
	if(!size){free(ptr);return 0;}
	void*result = realloc(ptr,size);
	if(result) return result;
	char buf[24] = {0};
	sprintf(buf,"%lu",(long)size);
	emit_error("Function ",func," failed to allocate ",buf," bytes of memory with realloc(",args,") on line ",line," of file ",file);
	return 0;
}
#define calloc_(...) calloc_(__VA_ARGS__,__FILE__,__func__,LINE,#__VA_ARGS__)

#ifdef __cplusplus
}
#endif
#endif