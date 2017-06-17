#ifndef SPACINGCALCULATOR_H
#define SPACINGCALCULATOR_H
#include"dialog.h"
class SpacingCalculator:public Dialog{
	QLineEdit*spacing=new QLineEdit(this),*edita=new QLineEdit(this),*editb=new QLineEdit(this),*editc=new QLineEdit(this),*editA=new QLineEdit(this),*editB=new QLineEdit(this),*editC=new QLineEdit(this);
	QSpinBox*edith=new QSpinBox(this),*editk=new QSpinBox(this),*editl=new QSpinBox(this);
public:
	SpacingCalculator(QWidget*parent=0);
	QString result();
public slots:
	void read();
	void write();
private slots:
	void calculate();
};
#endif