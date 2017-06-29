#include "dialogbutton.h"

DialogButton::DialogButton(QString str,QDialog*dlog,const std::function<void()>&func,QWidget*parent):QPushButton(str,parent),std::function<void()>(func),dialog(dlog){
	connect(this,&QPushButton::clicked,this,&DialogButton::call);
}
void DialogButton::call(){
	if(dialog->exec())(*this)();
}