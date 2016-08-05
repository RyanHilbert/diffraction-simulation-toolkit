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
	int links=0;//set of flags for linking Params to ComboParams
protected:
	void*data;//The piece of data to be updated by the interface in the upd8 slot
	QHBoxLayout*layout;//The layout that all components should be added to
	QLabel*label;
public:
	const char*const name;
	Param(const char*,void*,QWidget*parent=0);
signals:
	void edited();
public slots:
	void link(int);
	void tip(QString);
	virtual void upd8()=0;
	virtual void clear()=0;
	virtual bool hdf5read(hid_t)=0;
	virtual bool hdf5write(hid_t)=0;
};

class BoolParam:public Param{
	QCheckBox*edit=new QCheckBox(this);
public:
	BoolParam(const char*,bool*,QWidget*parent=0);
	void upd8();
	void clear();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};

class IntParam:public Param{
	QSpinBox*edit=new QSpinBox(this);
public:
	IntParam(const char*,int*,QWidget*parent=0);
	void upd8();
	void clear();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};

class DoubleParam:public Param{
	QLineEdit*edit;
	double multiplier=1;//number to multiply by to convert input to desired units
public:
	DoubleParam(const char*,double*,QWidget*parent=0);
	void upd8();
	void clear();
	void set(QString);
	double setMultiplier(double);
	QString get();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};
class ComplexParam:public Param{
	QLineEdit*edit,*editI;
public:
	ComplexParam(const char*,complex*,QWidget*parent=0);
	void upd8();
	void clear();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};
class VectorParam:public Param{
	QLineEdit*editX,*editY,*editZ;
public:
	VectorParam(const char*,vector*,QWidget*parent=0);
	void upd8();
	void clear();
	void set(vector);
	vector get();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};

class FunctionParam:public Param{
	QLineEdit*edit;
public:
	FunctionParam(const char*,char**,QWidget*parent=0);
	void upd8();
	void clear();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
};

class ComboParam:public Param{
	Q_OBJECT
	QComboBox*box;
public:
	ComboParam(const char*,int*,std::initializer_list<QString>,QWidget*parent=0);
	void upd8();
	void clear();
	bool hdf5read(hid_t);
	bool hdf5write(hid_t);
signals:
	void currentIndexChanged(int);
};

//Overloaded convenience methods for generating Param widgets
BoolParam*newParam(const char*,bool*,QWidget*parent=0);
IntParam*newParam(const char*,int*,QWidget*parent=0);
DoubleParam*newParam(const char*,double*,QWidget*parent=0);
ComplexParam*newParam(const char*,complex*,QWidget*parent=0);
VectorParam*newParam(const char*,vector*,QWidget*parent=0);
FunctionParam*newParam(const char*,char**,QWidget*parent=0);
#endif