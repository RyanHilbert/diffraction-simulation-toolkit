#include<QtWidgets>

#include"window.h"
#include"coloredlabel.h"
#include"paramgroup.h"
#include"thread.h"
#include"main.h"

static int flagset(std::initializer_list<unsigned int>list){//used for linking Params with the LNK macro
	unsigned int flags=0;
	for(unsigned int i:list)flags|=1<<i;
	return flags;
}
Window::Window(QWindow*scene3D,QWidget*parent):QMainWindow(parent){
	setWindowTitle(defaultTitle);
	QPlainTextEdit*console=new QPlainTextEdit("Console initialized!\n");

	cpus->setSuffix(" CPU cores");
	cpus->setMinimum(1);
	cpus->setMaximum(get_num_procs_omp());
	connect(cpus,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,[](int value){g_max_cpu_cores=value;});
	cpus->setValue(cpus->maximum());

	gpus->setSuffix(" GPUs");
	gpus->setMaximum(getDeviceCountCuda());
	connect(gpus,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,[](int value){g_max_gpus=value;});
	gpus->setValue(gpus->maximum());

	QMenu*menu=menuBar()->addMenu("File");
	menu->addAction("New",this,SLOT(clear()),QKeySequence::New);
	menu->addAction("Open",this,SLOT(load()),QKeySequence::Open);
	menu->addAction("Save",this,SLOT(save()),QKeySequence::Save);
	menu->addAction("Export",this,SLOT(exportiff()),QKeySequence::AddTab);

	QWidget*base=new QWidget(this);
	setCentralWidget(base);
	QHBoxLayout*baseLayout=new QHBoxLayout(base);
	{//leftLayout scope
		QVBoxLayout*leftLayout=new QVBoxLayout;
		baseLayout->addLayout(leftLayout);
		{//upperLeftLayout scope
			QHBoxLayout*upperLeftLayout=new QHBoxLayout();
			leftLayout->addLayout(upperLeftLayout);
			upperLeftLayout->addWidget(editable);
			upperLeftLayout->addWidget(start);
			upperLeftLayout->addWidget(pause);
			upperLeftLayout->addWidget(gpus);
			upperLeftLayout->addWidget(cpus);
			upperLeftLayout->addWidget(check3D);
			connect(editable,&QCheckBox::stateChanged,[=](int checked){
				if(checked){
					spacingButton->setEnabled(true);
					directionButton->setEnabled(true);
					for(ParamGroup*group:*this)group->enable();
				}
				else{
					spacingButton->setEnabled(false);
					directionButton->setEnabled(false);
					for(ParamGroup*group:*this)group->disable();
				}
			});
		}
		QScrollArea*scroll=new QScrollArea();
		scroll->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Preferred);
		leftLayout->addWidget(scroll);
		QWidget*frame=new QWidget(scroll);
		QVBoxLayout*groupLayout=new QVBoxLayout(frame);
		ParamGroup*group;
		Param*param;
		ComboParam*combo;
		Param*array[3]={0};
		int line=0;
#define USE_ALL_MACROS//Excessive use of macro sorcery below to make main.h readable. Edit main.h itself to alter parameter panels.
#define GRP(...) group=new ParamGroup(__VA_ARGS__);groupLayout->addWidget(group);push_back(group);
#define MNU(name,var,...) param=combo=new ComboParam(name,&var,{__VA_ARGS__},group);group->add(param,line==__LINE__);connect(param,&Param::edited,this,&Window::onEdit);line=__LINE__;
#define FOREACH(name,type,var,...) param=newParam(name,&var,group);group->add(param,line==__LINE__);connect(param,&Param::edited,this,&Window::onEdit);line=__LINE__;
#define LNK(...) param->link(flagset({__VA_ARGS__}));param->link(0);connect(combo,&ComboParam::currentIndexChanged,param,&Param::link);
#define MUL(...) ((DoubleParam*)param)->setMultiplier(__VA_ARGS__);
#define TIP(...) param->tip(__VA_ARGS__);
#include"main.h"
		editable->setChecked(true);
		editable->setChecked(false);
		groupLayout->addWidget(console);
		scroll->setWidgetResizable(true);
		scroll->setWidget(frame);
		scroll->setContentsMargins(0,0,0,0);
		leftLayout->addWidget(progress);
	}//end of leftLayout scope
	{//rightLayout scope
		QStackedLayout*stack=new QStackedLayout(baseLayout);
		QFrame*frame3D=new QFrame();
		frame3D->setAutoFillBackground(true);
		QPalette palette=frame3D->palette();
		palette.setColor(QPalette::Window,Qt::black);
		frame3D->setPalette(palette);
		stack->addWidget(manager);
		stack->addWidget(new QWidget(this));
		stack->addWidget(frame3D);
		QVBoxLayout*layout3D=new QVBoxLayout(frame3D);
		QWidget*widget3D=QWidget::createWindowContainer(scene3D,frame3D);
		widget3D->setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
		layout3D->addWidget(widget3D);
		QHBoxLayout*labels=new QHBoxLayout();
		layout3D->addLayout(labels);
		new ColoredLabel(Qt::red,"x axis",labels);
		new ColoredLabel(Qt::green,"y axis",labels);
		new ColoredLabel(Qt::blue,"z axis",labels);
		new ColoredLabel(Qt::yellow,"rotation axis",labels);
		new ColoredLabel(Qt::magenta,"h vector",labels);
		new ColoredLabel(Qt::cyan,"beam direction",labels);
		connect(&calculation,&Thread::completed1st,manager,&Spectromanager::load1st);
		connect(&calculation,&Thread::completed1st,check3D,&QCheckBox::setChecked);
		connect(&calculation,&Thread::progress,manager,&Spectromanager::setFrameCount);
		connect(check3D,&QCheckBox::stateChanged,stack,&QStackedLayout::setCurrentIndex);
		check3D->setChecked(true);
	}//end of rightLayout scope
	connect(start,&QPushButton::clicked,[=]{
		hdf5write();
		h5save(hdf5file);
		start->setEnabled(false);
		pause->setEnabled(true);
		cpus->setEnabled(false);
		gpus->setEnabled(false);
		menuBar()->setEnabled(false);
		editable->setEnabled(false);
		editable->setChecked(false);
		upd8();
		calculation.begin();
	});
	connect(&calculation,&QThread::finished,this,&Window::onFinish);
	connect(pause,&QPushButton::clicked,[=]{
		pause->setEnabled(false);
		calculation.stop();
	});
	pause->setEnabled(false);

	QPalette palette=console->palette();
	palette.setColor(QPalette::Base,Qt::black);
	palette.setColor(QPalette::Text,Qt::white);
	console->setPalette(palette);
	console->setFont(QFont("Courier"));
	console->setReadOnly(true);
	connect(&calculation,&Thread::output,console,&QPlainTextEdit::appendPlainText);

	progress->setTextVisible(true);
	progress->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Fixed);
	connect(&calculation,&Thread::maximum,progress,&QProgressBar::setMaximum);
	connect(&calculation,&Thread::progress,progress,&QProgressBar::setValue);

	FILE*f=fopen(defaultFilename,"r");
	bool exists=f?!fclose(f):false;
	hdf5file=h5file(defaultFilename);
	if(exists)hdf5read();
	else upd8();
	calculation.stop();
	calculation.begin();
	if(scene3D)scene3D->show();
}
Window::~Window(){
	H5Fclose(hdf5file);
	hdf5file=h5file(defaultFilename);
	hdf5write();
	H5Fclose(hdf5file);
}
void Window::clear(){
	if(QMessageBox::Ok==QMessageBox::warning(this,"Confirm","This will clear your current parameters. Are you sure?",QMessageBox::Ok|QMessageBox::Cancel,QMessageBox::Cancel)){
		H5Fclose(hdf5file);
		hdf5file=h5file(defaultFilename);
		for(ParamGroup*group:*this)for(Param*param:*group)param->clear();
		manager->clear();
	}
}
bool Window::load(){
	QString name=QFileDialog::getOpenFileName(this,"Open","","HDF5 (*.h5);;All Files (*)");
	if(name.isNull())return false;
	setWindowTitle(defaultTitle+(" - "+QFileInfo(name).fileName()));
	H5Fclose(hdf5file);
	QByteArray array=name.toUtf8();
	hdf5file=h5file(array.data());
	bool result=hdf5read();
	calculation.stop();
	calculation.begin();
	manager->setFrame();
	return result;
}
bool Window::save(){
	if(windowTitle()!=defaultTitle)return false;
	QString newName=QFileDialog::getSaveFileName(this,"Save As","","HDF5 (*.h5)");
	if(newName.isNull())return false;
	setWindowTitle(defaultTitle+(" - "+QFileInfo(newName).fileName()));
	QByteArray array=newName.toUtf8();
	char*utf8name=array.data();
	H5Fclose(hdf5file);
	remove(utf8name);
	rename(defaultFilename,utf8name);
	hdf5file=h5file(utf8name);
	return hdf5write();
}
bool Window::hdf5read(){
	bool result=true,blocked=blockSignals(true);
	for(ParamGroup*group:*this)result&=group->hdf5read(hdf5file);
	blockSignals(blocked);
	spacing->hdf5read();
	upd8();
	return result;
}
bool Window::hdf5write(){
	bool result=true;
	for(ParamGroup*group:*this)result&=group->hdf5write(hdf5file);
	spacing->hdf5write();
	h5save(hdf5file);
	return result;
}
void Window::lock(){
	for(auto group:*this)group->disable();
}
void Window::unlock(){
	for(auto group:*this)group->enable();
}
void Window::upd8(){
	for(auto group:*this)for(auto param:*group)param->upd8();
}
void Window::onEdit(){
	setWindowTitle(defaultTitle);
	H5Fclose(hdf5file);
	hdf5file=h5overwrite(defaultFilename);
	hdf5write();
	h5save(hdf5file);
	progress->setValue(progress->minimum());
	manager->clear();
	upd8();
	emit edited();
}
void Window::onFinish(){
	menuBar()->setEnabled(true);
	editable->setEnabled(true);
	pause->setEnabled(false);
	start->setEnabled(true);
	cpus->setEnabled(true);
	gpus->setEnabled(true);
}
//all code beyond this line is for generating TIFF files
enum Tag{
	ImageWidth=256,ImageLength,BitsPerSample,Compression,
	PhotometricInterpretation=262,
	StripOffsets=273,
	RowsPerStrip=278,StripByteCounts,
	XResolution=282,YResolution,
	ResolutionUnit=296,
	SampleFormat=339
};
bool Window::exportiff(){
	QString exportname=QFileDialog::getSaveFileName(this,"Export As","","TIFF (*.tif)");
	if(exportname.isNull())return false;
	FILE*f=fopen(exportname.toUtf8().data(),"wb");
	if(!f)return false;

	hid_t detector=h5group("Detector",hdf5file);
	const int rows=h5int("Rows",detector);
	const int columns=h5int("Columns",detector);
	const int images=h5int("Iterations Completed",detector);

	uint32_t u4=16;
	uint16_t u2=42;

	fwrite((*(char*)&u2)?"II":"MM",1,2,f);//endianness
	fwrite(&u2,2,1,f);//magic number 42
	fwrite(&u4,4,1,f);//offset to first ImageFileDirectory(IFD)
	fwrite(&(u4=1),4,1,f);//numerator of rational resolution fields
	fwrite(&u4,4,1,f);//denominator of rational resolution fields

#define FIELD_COUNT 12//number of fields per ImageFileDirectory(IFD)
#define IFD_SIZE (12*FIELD_COUNT+2+4)//size of each IFD in bytes
#define HEADER_SIZE (8+8)//size of header and resolution field info in bytes

#define SHORT_FIELD(tag,value) fwrite(&(u2=tag),2,1,f);fwrite(&(u2=3),2,1,f);fwrite(&(u4=1),4,1,f);fwrite(&(u2=value),2,1,f);fwrite(&(u2=0),2,1,f);
#define LONG_FIELD(tag,value) fwrite(&(u2=tag),2,1,f);fwrite(&(u2=4),2,1,f);fwrite(&(u4=1),4,1,f);fwrite(&(u4=value),4,1,f);
#define RATIONAL_FIELD(tag) fwrite(&(u2=tag),2,1,f);fwrite(&(u2=5),2,1,f);fwrite(&(u4=1),4,1,f);fwrite(&(u4=8),4,1,f);

	for(int i=0;i<images;++i){//write each IFD to TIFF file
		fwrite(&(u2=FIELD_COUNT),2,1,f);//field count, followed by each 12-byte field
		LONG_FIELD(ImageWidth,columns)
				LONG_FIELD(ImageLength,rows)
				SHORT_FIELD(BitsPerSample,32)
				SHORT_FIELD(Compression,1)
				SHORT_FIELD(PhotometricInterpretation,1)
				LONG_FIELD(RowsPerStrip,rows)
				LONG_FIELD(StripOffsets,HEADER_SIZE+IFD_SIZE*images+i*rows*columns*sizeof(float))
				LONG_FIELD(StripByteCounts,rows*columns*sizeof(float))
				RATIONAL_FIELD(XResolution)
				RATIONAL_FIELD(YResolution)
				SHORT_FIELD(ResolutionUnit,1)
				SHORT_FIELD(SampleFormat,3)
				fwrite(&(u4=(i+1==images)?0:(HEADER_SIZE+(i+1)*IFD_SIZE)),4,1,f);//offset of next IFD
	}
	for(int i=0;i<images;++i){
		struct h5dataset data=h5dataset(("det"+std::to_string(i)).c_str(),detector);
		fwrite(data.data,sizeof(float),rows*columns,f);
	}
#pragma omp critical(OMP_CRIT_HDF5)//synchronize raw HDF5 call
	{H5Oclose(detector);}
	fclose(f);
	return true;
}