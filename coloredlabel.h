#ifndef COLOREDLABEL_H
#define COLOREDLABEL_H

#include<QtWidgets>

class ColoredLabel:public QLabel{

public:
	ColoredLabel(const QColor&color,const QString&text,QLayout*parent=0);
};
#endif