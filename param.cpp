#include"param.h"

Param::Param(const char*name,void*dfault,QWidget*parent):QWidget(parent),data(dfault),layout(new QHBoxLayout(this)),label(new QLabel(name)),name(name){
	layout->setMargin(0);
	layout->addWidget(label);
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
BoolParam::BoolParam(const char*name,bool*dfault,QWidget*parent):Param(name,dfault,parent),edit(new QCheckBox(name,this)){
	layout->removeWidget(label);
	delete label;
	layout->addWidget(edit);
	edit->setChecked(*dfault);
	connect(edit,&QCheckBox::toggled,this,&Param::edited);
}
void BoolParam::upd8(){
	*(bool*)data=edit->isChecked();
}
void BoolParam::clear(){
	edit->setChecked(false);
}
bool BoolParam::hdf5read(hid_t loc){
	int attribute=h5int(name,loc);
	bool blocked=blockSignals(true);
	edit->setChecked(attribute);
	blockSignals(blocked);
	return attribute!=INT_MIN;
}
bool BoolParam::hdf5write(hid_t loc){
	return h5set_int(edit->isChecked(),name,loc);
}

IntParam::IntParam(const char*name,int*dfault,QWidget*parent):Param(name,dfault,parent){
	edit->setMinimum(1);
	edit->setMaximum(INT_MAX);
	edit->setValue(*dfault);
	layout->addWidget(edit);
	connect(edit,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,this,&Param::edited);
}
void IntParam::upd8(){
	*(int*)data=edit->text().toInt();
}
void IntParam::clear(){
	edit->setValue(edit->minimum());
}
bool IntParam::hdf5read(hid_t loc){
	int attribute=h5int(name,loc);
	bool blocked=blockSignals(true);
	edit->setValue(attribute);
	blockSignals(blocked);
	return attribute!=INT_MIN;
}
bool IntParam::hdf5write(hid_t loc){
	return h5set_int(edit->value(),name,loc);
}

DoubleParam::DoubleParam(const char*name,double*dfault,QWidget*parent):Param(name,dfault,parent),edit(new QLineEdit(QString::number(*dfault),this)){
	edit->setValidator(new QDoubleValidator());
	layout->addWidget(edit);
	connect(edit,&QLineEdit::textChanged,this,&Param::edited);
}
void DoubleParam::upd8(){
	*(double*)data=edit->text().toDouble()*multiplier;
}
void DoubleParam::clear(){
	edit->setText("");
}
void DoubleParam::set(QString str){
	edit->setText(str);
}
double DoubleParam::setMultiplier(double mul){
	double temp=multiplier;
	multiplier=mul;
	return temp;
}
QString DoubleParam::get(){
	return edit->text();
}
bool DoubleParam::hdf5read(hid_t loc){
	double attribute=h5double(name,loc);
	bool blocked=blockSignals(true);
	edit->setText(QString::number(attribute));
	blockSignals(blocked);
	return attribute==attribute;//will return false if attribute is NaN
}
bool DoubleParam::hdf5write(hid_t loc){
	return h5set_double(edit->text().toDouble(),name,loc);
}
ComplexParam::ComplexParam(const char*name,complex*dfault,QWidget*parent):Param(name,dfault,parent),edit(new QLineEdit(QString::number(dfault->real()),this)),editI(new QLineEdit(QString::number(dfault->imag()),this)){
	edit->setValidator(new QDoubleValidator());
	editI->setValidator(new QDoubleValidator());
	layout->addWidget(edit);
	layout->addWidget(new QLabel("+"));
	layout->addWidget(editI);
	layout->addWidget(new QLabel("i"));
	connect(edit,&QLineEdit::textChanged,this,&Param::edited);
	connect(editI,&QLineEdit::textChanged,this,&Param::edited);
}
void ComplexParam::upd8(){
	*(complex*)data=edit->text().toDouble()+I*editI->text().toDouble();
}
void ComplexParam::clear(){
	edit->setText("");
	editI->setText("");
}
bool ComplexParam::hdf5read(hid_t loc){
	double buffer[2]={NAN,NAN};
	bool result=false,blocked=blockSignals(true);
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTfind_attribute(loc,name))result=H5LTget_attribute_double(loc,".",name,buffer)>=0;}
	edit->setText(QString::number(buffer[0]));
	editI->setText(QString::number(buffer[1]));
	blockSignals(blocked);
	return result;
}
bool ComplexParam::hdf5write(hid_t loc){
	bool result=false;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTset_attribute_double(loc,".",name,(double*)data,2);}
	return result;
}
VectorParam::VectorParam(const char*name,vector*dfault,QWidget*parent):Param(name,dfault,parent),editX(new QLineEdit(QString::number(dfault->x),this)),editY(new QLineEdit(QString::number(dfault->y),this)),editZ(new QLineEdit(QString::number(dfault->z),this)){
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
void VectorParam::upd8(){
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
bool VectorParam::hdf5read(hid_t loc){
	double buffer[3]={NAN,NAN,NAN};
	bool result=false,blocked=blockSignals(true);
#pragma omp critical(OMP_CRIT_HDF5)
	{if(H5LTfind_attribute(loc,name))result=H5LTget_attribute_double(loc,".",name,buffer)>=0;}
	editX->setText(QString::number(buffer[0]));
	editY->setText(QString::number(buffer[1]));
	editZ->setText(QString::number(buffer[2]));
	blockSignals(blocked);
	return result;
}
bool VectorParam::hdf5write(hid_t loc){
	bool result=false;
#pragma omp critical(OMP_CRIT_HDF5)
	{result=H5LTset_attribute_double(loc,".",name,(double*)data,3)>=0;}
	return result;
}

FunctionParam::FunctionParam(const char*name,char**dfault,QWidget*parent):Param(name,dfault,parent),edit(new QLineEdit("0",this)){
	layout->addWidget(edit);
	connect(edit,&QLineEdit::textChanged,this,&Param::edited);
	connect(edit,&QLineEdit::editingFinished,[=]{//set our own validator to ensure the function is valid and only uses variables "x" "y" "z"
		double x,y,z;
		te_variable vars[]={{"x",&x},{"y",&y},{"z",&z}};
		te_expr*expr=te_compile(edit->text().toUtf8().data(),vars,3,0);
		if(expr)te_free(expr);
		else edit->setText("0");
	});
}
void FunctionParam::upd8(){
	double x,y,z;
	te_variable vars[]={{"x",&x},{"y",&y},{"z",&z}};
	QByteArray array=edit->text().toUtf8();
	char*str=array.data();
	te_expr*expr=te_compile(str,vars,3,0);
	if(expr){
		te_free(expr);
		int length=array.length()+1;
		free(*(char**)data);
		*(char**)data=(char*)malloc(length);
		memcpy(*(char**)data,str,length);
	}
	else{
		edit->setText("0");
		free(*(char**)data);
		*(char**)data=0;
	}
}
void FunctionParam::clear(){
	edit->setText("0");
}
bool FunctionParam::hdf5read(hid_t loc){
	std::string attribute=h5string(name,loc);
	bool blocked=blockSignals(true);
	edit->setText(attribute.c_str());
	blockSignals(blocked);
	return!attribute.empty();
}
bool FunctionParam::hdf5write(hid_t loc){
	QByteArray array=edit->text().toUtf8();
	return h5set_string(array.data(),name,loc);
}

ComboParam::ComboParam(const char*name,int*dfault,std::initializer_list<QString>list,QWidget*parent):Param(name,dfault,parent),box(new QComboBox(this)){
	QList<QString>qlist = QList<QString>();
	for (QString string : list)qlist.append(string);
	box->insertItems(0,qlist);
	layout->addWidget(box);
	connect(box,(void(QComboBox::*)(int))&QComboBox::currentIndexChanged,this,&ComboParam::currentIndexChanged);
	connect(box,(void(QComboBox::*)(int))&QComboBox::currentIndexChanged,this,&ComboParam::edited);
}
void ComboParam::upd8(){
	*(int*)data=box->currentIndex();
}
void ComboParam::clear(){
	box->setCurrentIndex(-1);
}
bool ComboParam::hdf5read(hid_t loc){
	std::string attribute=h5string(name,loc);
	bool blocked=blockSignals(true);
	box->setCurrentIndex(box->findText(QString::fromStdString(attribute)));
	blockSignals(blocked);
	emit currentIndexChanged(box->currentIndex());
	return!attribute.empty();
}
bool ComboParam::hdf5write(hid_t loc){
	QByteArray array=box->currentText().toUtf8();
	return h5set_string(array.data(),name,loc);
}

BoolParam*newParam(const char*name,bool*dfault,QWidget*parent){return new BoolParam(name,dfault,parent);}
IntParam*newParam(const char*name,int*dfault,QWidget*parent){return new IntParam(name,dfault,parent);}
DoubleParam*newParam(const char*name,double*dfault,QWidget*parent){return new DoubleParam(name,dfault,parent);}
ComplexParam*newParam(const char*name,complex*dfault,QWidget*parent){return new ComplexParam(name,dfault,parent);}
VectorParam*newParam(const char*name,vector*dfault,QWidget*parent){return new VectorParam(name,dfault,parent);}
FunctionParam*newParam(const char*name,char**dfault,QWidget*parent){return new FunctionParam(name,dfault,parent);}