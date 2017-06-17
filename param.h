#ifndef PARAM_H
#define PARAM_H
#include<QtWidgets>
#include<initializer_list>
#include<complex>
#include"main.h"

//These classes are responsible for displaying a configurable parameter to the user
//Param itself cannot be instantiated and must be subclassed for each type of parameter

class Param:public QWidget{
	Q_OBJECT
	unsigned int links=0;//set of flags for linking Params to ComboParams
protected:
	void*data;//The piece of data to be updated by the interface in the upd8 slot
	QHBoxLayout*layout;//The layout that all components should be added to
	QLabel*label;
public:
	const signed char id;//should be between 1 and 127 inclusive
	const char*const name;
	Param(const signed char,const char*,void*,QWidget*parent=0);
signals:
	void edited();
public slots:
	void link(int);
	virtual void tip(QString);
	virtual void write()=0;
	virtual void read()=0;
	virtual void clear()=0;
};

class BoolParam:public Param{
	QCheckBox*edit=new QCheckBox(this);
public:
	BoolParam(const signed char,const char*,bool*,QWidget*parent=0);
	void read();
	void write();
	void clear();
	void tip(QString);
};

class IntParam:public Param{
	QSpinBox*edit=new QSpinBox(this);
public:
	IntParam(const signed char,const char*,int*,QWidget*parent=0);
	void read();
	void write();
	void clear();
};

class DoubleParam:public Param{
	QLineEdit*edit;
	double multiplier=1;//number to multiply by to convert input to desired units
public:
	DoubleParam(const signed char,const char*,double*,QWidget*parent=0);
	void set(QString);
	double setMultiplier(double);
	QString get();
	void read();
	void write();
	void clear();
};
class ComplexParam:public Param{
	QLineEdit*edit,*editI;
public:
	ComplexParam(const signed char,const char*,complex*,QWidget*parent=0);
	void read();
	void write();
	void clear();
};
class VectorParam:public Param{
	QLineEdit*editX,*editY,*editZ;
public:
	VectorParam(const signed char,const char*,vector*,QWidget*parent=0);
	void set(vector);
	vector get();
	void read();
	void write();
	void clear();
};

class FunctionParam:public Param{
	QLineEdit*edit;
public:
	FunctionParam(const signed char,const char*,char256*,QWidget*parent=0);
	void read();
	void write();
	void clear();
};

class ComboParam:public Param{
	Q_OBJECT
	QComboBox*box;
public:
	ComboParam(const signed char,const char*,char*,std::initializer_list<QString>,QWidget*parent=0);
	void read();
	void write();
	void clear();
signals:
	void currentIndexChanged(int);
};
//Overloaded convenience methods for generating Param widgets
BoolParam*newParam(const signed char,const char*,bool*,QWidget*parent=0);
IntParam*newParam(const signed char,const char*,int*,QWidget*parent=0);
DoubleParam*newParam(const signed char,const char*,double*,QWidget*parent=0);
ComplexParam*newParam(const signed char,const char*,complex*,QWidget*parent=0);
VectorParam*newParam(const signed char,const char*,vector*,QWidget*parent=0);
FunctionParam*newParam(const signed char,const char*,char256*,QWidget*parent=0);
#endif