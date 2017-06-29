#include"param.h"
#include"tif.h"

Param::Param(const signed char id,const char*name,void*dfault,QWidget*parent):QWidget(parent),data(dfault),layout(new QHBoxLayout(this)),label(new QLabel(name)),id(id),name(name){
	layout->setMargin(0);
	layout->addWidget(label);
	connect(this,&Param::edited,this,&Param::write);
}
void Param::link(int flags){
	if(links)setVisible(links&(1<<flags));
	else links=flags;
}
void Param::tip(QString string){
	label->setToolTip(string);
	QFont font=label->font();
	font.setUnderline(true);
	label->setFont(font);
}
BoolParam::BoolParam(const signed char id,const char*name,bool*dfault,QWidget*parent):Param(id,name,dfault,parent),edit(new QCheckBox(name,this)){
	layout->removeWidget(label);
	delete label;
	layout->addWidget(edit);
	edit->setChecked(*dfault);
	connect(edit,&QCheckBox::toggled,this,&Param::edited);
}
void BoolParam::read(){
	const bool blocked = blockSignals(true);
	edit->setChecked(*(bool*)data);
	blockSignals(blocked);
}
void BoolParam::write(){
	*(bool*)data=edit->isChecked();
}
void BoolParam::clear(){
	edit->setChecked(false);
}
void BoolParam::tip(QString string){
	edit->setToolTip(string);
	QFont font=edit->font();
	font.setUnderline(true);
	edit->setFont(font);
}

IntParam::IntParam(const signed char id,const char*name,int*dfault,QWidget*parent):Param(id,name,dfault,parent){
	edit->setMinimum(1);
	edit->setMaximum(INT_MAX);
	edit->setValue(*dfault);
	layout->addWidget(edit);
	connect(edit,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,this,&Param::edited);
}
void IntParam::read(){
	const bool blocked = blockSignals(true);
	edit->setValue(*(int*)data);
	blockSignals(blocked);
}
void IntParam::write(){
	*(int*)data=edit->text().toInt();
}
void IntParam::clear(){
	edit->setValue(edit->minimum());
}

DoubleParam::DoubleParam(const signed char id,const char*name,double*dfault,QWidget*parent):Param(id,name,dfault,parent),edit(new QLineEdit(QString::number(*dfault),this)){
	edit->setValidator(new QDoubleValidator());
	layout->addWidget(edit);
	connect(edit,&QLineEdit::textChanged,this,&Param::edited);
}
void DoubleParam::read(){
	const bool blocked = blockSignals(true);
	edit->setText(QString::number(*(double*)data/multiplier));
	blockSignals(blocked);
}
void DoubleParam::write(){
	*(double*)data=edit->text().toDouble()*multiplier;
}
void DoubleParam::clear(){
	edit->setText("");
}
void DoubleParam::set(QString str){
	edit->setText(str);
}
double DoubleParam::setMultiplier(double mul){
	double temp = multiplier;
	multiplier = mul;
	return temp;
}
QString DoubleParam::get(){
	return edit->text();
}

ComplexParam::ComplexParam(const signed char id,const char*name,complex*dfault,QWidget*parent):Param(id,name,dfault,parent),edit(new QLineEdit(QString::number(dfault->real()),this)),editI(new QLineEdit(QString::number(dfault->imag()),this)){
	edit->setValidator(new QDoubleValidator());
	editI->setValidator(new QDoubleValidator());
	layout->addWidget(edit);
	layout->addWidget(new QLabel("+"));
	layout->addWidget(editI);
	layout->addWidget(new QLabel("i"));
	connect(edit,&QLineEdit::textChanged,this,&Param::edited);
	connect(editI,&QLineEdit::textChanged,this,&Param::edited);
}
void ComplexParam::read(){
	const complex c = *(complex*)data;
	const bool blocked = blockSignals(true);
	edit->setText(QString::number(c.real()));
	editI->setText(QString::number(c.imag()));
	blockSignals(blocked);
}
void ComplexParam::write(){
	*(complex*)data=edit->text().toDouble()+I*editI->text().toDouble();
}
void ComplexParam::clear(){
	edit->setText("");
	editI->setText("");
}

VectorParam::VectorParam(const signed char id,const char*name,vector*dfault,QWidget*parent):Param(id,name,dfault,parent),editX(new QLineEdit(QString::number(dfault->x),this)),editY(new QLineEdit(QString::number(dfault->y),this)),editZ(new QLineEdit(QString::number(dfault->z),this)){
	editX->setValidator(new QDoubleValidator(editX));
	editY->setValidator(new QDoubleValidator(editY));
	editZ->setValidator(new QDoubleValidator(editZ));
	layout->addWidget(editX);
	layout->addWidget(new QLabel("x"));
	layout->addWidget(editY);
	layout->addWidget(new QLabel("y"));
	layout->addWidget(editZ);
	layout->addWidget(new QLabel("z"));
	connect(editX,&QLineEdit::textChanged,this,&Param::edited);
	connect(editY,&QLineEdit::textChanged,this,&Param::edited);
	connect(editZ,&QLineEdit::textChanged,this,&Param::edited);
}
void VectorParam::read(){
	vector v = *(vector*)data;
	const bool blocked = blockSignals(true);
	editX->setText(QString::number(v.x));
	editY->setText(QString::number(v.y));
	editZ->setText(QString::number(v.z));
	blockSignals(blocked);
}
void VectorParam::write(){
	*(vector*)data={editX->text().toDouble(),editY->text().toDouble(),editZ->text().toDouble()};
}
void VectorParam::clear(){
	editX->setText("");
	editY->setText("");
	editZ->setText("");
}
void VectorParam::set(vector v){
	editX->setText(QString::number(v.x));
	editY->setText(QString::number(v.y));
	editZ->setText(QString::number(v.z));
}
vector VectorParam::get(){
	return{editX->text().toDouble(),editY->text().toDouble(),editZ->text().toDouble()};
}

FunctionParam::FunctionParam(const signed char id,const char*name,char256*dfault,QWidget*parent):Param(id,name,dfault,parent),edit(new QLineEdit("0",this)){
	edit->setMaxLength(255);
	layout->addWidget(edit);
	connect(edit,&QLineEdit::editingFinished,[=]{//set our own validator to ensure the function is valid and only uses variables "x" "y" "z"
		double x,y,z;
		te_variable vars[]={{"x",&x,0,0},{"y",&y,0,0},{"z",&z,0,0}};
		te_expr*expr=te_compile(edit->text().toUtf8().data(),vars,3,0);
		if(expr)te_free(expr);
		else edit->setText("0");
		if(edit->isModified())emit edited();
		edit->setModified(false);
	});
}
void FunctionParam::read(){
	const char*str = *(char256*)data;
	double x,y,z;
	te_variable vars[]={{"x",&x,0,0},{"y",&y,0,0},{"z",&z,0,0}};
	te_expr*expr=te_compile(str,vars,3,0);
	const bool blocked = blockSignals(true);
	if(expr){
		te_free(expr);
		edit->setText(str);
	}else edit->setText("0");
	blockSignals(blocked);
}
void FunctionParam::write(){
	double x,y,z;
	te_variable vars[]={{"x",&x,0,0},{"y",&y,0,0},{"z",&z,0,0}};
	QByteArray array=edit->text().toUtf8();
	char*str=array.data();
	te_expr*expr=te_compile(str,vars,3,0);
	if(expr){
		te_free(expr);
		memcpy(*(char256*)data,str,array.length()+1);
	}else{
		edit->setText("0");
		**(char256*)data=0;
	}
}
void FunctionParam::clear(){
	edit->setText("0");
}

ComboParam::ComboParam(const signed char id,const char*name,char*dfault,std::initializer_list<QString>list,QWidget*parent):Param(id,name,dfault,parent),box(new QComboBox(this)){
	QList<QString>qlist = QList<QString>();
	for (QString string : list)qlist.append(string);
	box->insertItems(0,qlist);
	layout->addWidget(box);
	connect(box,(void(QComboBox::*)(int))&QComboBox::currentIndexChanged,this,&ComboParam::currentIndexChanged);
	connect(box,(void(QComboBox::*)(int))&QComboBox::currentIndexChanged,this,&ComboParam::edited);
}
void ComboParam::read(){
	const bool blocked = blockSignals(true);
	box->setCurrentIndex(*(char*)data);
	blockSignals(blocked);
	emit currentIndexChanged(box->currentIndex());
}
void ComboParam::write(){
	*(char*)data=box->currentIndex();
}
void ComboParam::clear(){
	box->setCurrentIndex(-1);
}
BoolParam*newParam(const signed char id,const char*name,bool*dfault,QWidget*parent){return new BoolParam(id,name,dfault,parent);}
IntParam*newParam(const signed char id,const char*name,int*dfault,QWidget*parent){return new IntParam(id,name,dfault,parent);}
DoubleParam*newParam(const signed char id,const char*name,double*dfault,QWidget*parent){return new DoubleParam(id,name,dfault,parent);}
ComplexParam*newParam(const signed char id,const char*name,complex*dfault,QWidget*parent){return new ComplexParam(id,name,dfault,parent);}
VectorParam*newParam(const signed char id,const char*name,vector*dfault,QWidget*parent){return new VectorParam(id,name,dfault,parent);}
FunctionParam*newParam(const signed char id,const char*name,char256*dfault,QWidget*parent){return new FunctionParam(id,name,dfault,parent);}