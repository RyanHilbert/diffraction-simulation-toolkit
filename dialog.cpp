#include "dialog.h"

Dialog::Dialog(QWidget*parent,Qt::WindowFlags flags):QDialog(parent,flags){
	QVBoxLayout*base=new QVBoxLayout(this);
	QDialogButtonBox*buttons=new QDialogButtonBox(QDialogButtonBox::Apply|QDialogButtonBox::Cancel,this);
	base->addWidget(frame);
	base->addWidget(buttons);
	connect(buttons->button(QDialogButtonBox::Apply),&QPushButton::clicked,this,&QDialog::accept);
	connect(buttons,&QDialogButtonBox::accepted,this,&QDialog::accept);
	connect(buttons,&QDialogButtonBox::rejected,this,&QDialog::reject);
}