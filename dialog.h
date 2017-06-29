#ifndef DIALOG_H
#define DIALOG_H
#include<QtWidgets>
class Dialog:public QDialog{
protected:
	QFrame*frame=new QFrame(this);
public:
	Dialog(QWidget*parent=0,Qt::WindowFlags flags=0);
};
#endif