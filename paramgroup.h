#ifndef PARAMGROUP_H
#define PARAMGROUP_H
#include<QtWidgets>
#include"main.h"
#include"param.h"
class ParamGroup:public QGroupBox,private std::vector<Param*>{
	Q_OBJECT
public:
	QVBoxLayout*const layout=new QVBoxLayout();
	QHBoxLayout*line=new QHBoxLayout();
	ParamGroup(QString name="",QWidget*parent=0);
	using std::vector<Param*>::begin;
	using std::vector<Param*>::end;
	using std::vector<Param*>::push_back;
public slots:
	void enable();
	void disable();
	void add(Param*,bool inlined=false);
};
#endif