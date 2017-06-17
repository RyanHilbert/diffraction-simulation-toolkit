#include<QApplication>
#include<QtPlugin>
#include"window.h"
#include"window3d.h"
#define FOREACH(tag,name,type,var,...) type var = __VA_ARGS__;
#include"main.h"
static const char*_(){
	QString string=QStandardPaths::writableLocation(QStandardPaths::DataLocation)+"/.tif";
	QByteArray bytes=string.toUtf8();
	const char*data=bytes.constData();
	void*result=malloc(bytes.length()+1);
	memcpy(result,data,bytes.length()+1);
	return(const char*)result;
}
const char*const tmp_tif_name=_();
long tif_offset=0;
FILE*tif=0;
//Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin)
int main(int argc,char*argv[]){
	QApplication app(argc,argv);
	Q_INIT_RESOURCE(shaders);
	Window3D*window3D=new Window3D();
	Window window(window3D);
	window.show();
	QObject::connect(&window,&Window::edited,window3D,&Window3D::update);
	emit window.edited();
	return app.exec();
}