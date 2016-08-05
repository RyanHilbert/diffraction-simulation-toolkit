#ifndef SPACINGCALCULATOR_H
#define SPACINGCALCULATOR_H
#include"dialog.h"
class SpacingCalculator:public Dialog{
#define _ new QLineEdit(this)
	QLineEdit*spacing=_,*edita=_,*editb=_,*editc=_,*editA=_,*editB=_,*editC=_;
#undef _
	QSpinBox*edith=new QSpinBox(this),*editk=new QSpinBox(this),*editl=new QSpinBox(this);
public:
	SpacingCalculator(QWidget*parent=0);
	QString result();
	void hdf5read();
	void hdf5write();
private slots:
	void calculate();
};
#endif