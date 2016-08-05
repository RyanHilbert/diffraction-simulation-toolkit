#include<QApplication>
#include<QtPlugin>
#include"window.h"
#include"window3d.h"
#define FOREACH(name,type,var,...) type var = __VA_ARGS__;
#include"main.h"
hid_t hdf5file=0;
//Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin)
int main(int argc,char*argv[]){
	QApplication app(argc,argv);
	Window3D*window3D=new Window3D();
	Window window(window3D);
	window.show();
	QObject::connect(&window,&Window::edited,window3D,&Window3D::update);
	emit window.edited();
	return app.exec();
}