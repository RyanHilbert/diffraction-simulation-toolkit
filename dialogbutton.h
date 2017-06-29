#ifndef DIALOGBUTTON_H
#define DIALOGBUTTON_H
#include<functional>
#include<QtWidgets>
class DialogButton:public QPushButton,public std::function<void()>{
	Q_OBJECT
	QDialog*dialog;
public:
	DialogButton(QString,QDialog*,const function&,QWidget*parent=0);
private slots:
	void call();
};
#endif