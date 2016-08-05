#ifndef DIRECTIONCALCULATOR_H
#define DIRECTIONCALCULATOR_H
#include"dialog.h"
#include"param.h"
#include"main.h"
class DirectionCalculator:public Dialog{
	DoubleParam*spacingParam,*energyParam;
	VectorParam*hParam,*beamParam;
	QVBoxLayout*vertical=new QVBoxLayout(frame);
	QLineEdit
	*spacingBox=new QLineEdit(this),
	*energyBox=new QLineEdit(this),
	*braggBox=new QLineEdit(this),
	*hx=new QLineEdit(this),
	*hy=new QLineEdit(this),
	*hz=new QLineEdit(this),
	*userx=new QLineEdit(this),
	*usery=new QLineEdit(this),
	*userz=new QLineEdit(this),
	*resx=new QLineEdit(this),
	*resy=new QLineEdit(this),
	*resz=new QLineEdit(this);
public:
	DirectionCalculator(QWidget*parent=0);
	vector result();
	void setParams(DoubleParam*,DoubleParam*,VectorParam*,VectorParam*);
	void setVisible(bool);
private slots:
	void calculate();
	void onAccept();
};
#endif