#include"thread.h"
#include"stdarg.h"
extern"C"{
extern volatile bool g_halted;
void calc();
void emit_warning(const char*str){emit calculation.warning(str);}
void emit_output(const char*str){emit calculation.output(str);}
void emit_maximum(int max){emit calculation.maximum(max);}
void emit_progress(int prg){emit calculation.progress(prg);}
void emit_completed1st(){emit calculation.completed1st();}
void*emit_error(const char*str,...){//accepts a null-terminated list of strings to output
	QString string = "";
	va_list v;
	va_start(v,str);
	while(str){
		string += str;
		str = va_arg(v,const char*);
	}
	va_end(v);
	bool calc = QThread::currentThread() == calculation.currentThread();
	if(calc) string += "<br>Halting Calculation";
	emit calculation.error(string);//connected to handling slot in window.cpp
	if(calc) calculation.sleep(~0UL);
	return 0;
}}
void Thread::begin(){
	if(!isRunning()){
		emit starting();
		start();
	}
}
void Thread::run(){calc();}
void Thread::stop(){g_halted=true;}
Thread calculation;