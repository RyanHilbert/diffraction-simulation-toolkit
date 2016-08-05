#include"thread.h"
extern"C"{
extern volatile bool halted;
void calc();
void emit_output(char*str){emit calculation.output(str);}
void emit_maximum(int max){emit calculation.maximum(max);}
void emit_progress(int prg){emit calculation.progress(prg);}
void emit_completed1st(){emit calculation.completed1st();}
}
void Thread::begin(){
	if(!isRunning()){
		emit starting();
		start(HighestPriority);
	}
}
void Thread::run(){
	calc();
}
void Thread::stop(){
	halted=true;
}
Thread calculation;