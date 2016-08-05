#include"coloredlabel.h"

ColoredLabel::ColoredLabel(const QColor&color,const QString&text,QLayout*parent):QLabel(text){
	setAlignment(Qt::AlignCenter);
	QPalette pal=palette();
	pal.setColor(QPalette::WindowText,color);
	setPalette(pal);
	QFont fon=this->font();
	fon.setBold(true);
	setFont(fon);
	parent->addWidget(this);
}