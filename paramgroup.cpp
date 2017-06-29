#include "paramgroup.h"
ParamGroup::ParamGroup(QString name,QWidget*parent):QGroupBox(name,parent){
	QLayout*outer=new QVBoxLayout(this);
	outer->setContentsMargins(0,0,0,0);
	QWidget*frame=new QWidget(this);
	outer->addWidget(frame);
	frame->setLayout(layout);
	layout->addLayout(line);
	setStyleSheet("font-weight:bold");
	frame->setStyleSheet("font-weight:normal");
	setCheckable(true);
	connect(this,&QGroupBox::toggled,frame,&QWidget::setVisible);
}
void ParamGroup::enable(){
	for(auto param:*this)param->setEnabled(true);
}
void ParamGroup::disable(){
	for(auto param:*this)param->setEnabled(false);
}
void ParamGroup::add(Param*param,bool inlined){
	if(!inlined)layout->addLayout(line=new QHBoxLayout());
	line->addWidget(param);
	push_back(param);
}