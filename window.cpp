#include<QtWidgets>

#include"window.h"
#include"coloredlabel.h"
#include"paramgroup.h"
#include"thread.h"
#include"tif.h"
extern"C"{
	extern int omp_get_num_procs();
}
static int flagset(std::initializer_list<size_t>list){//used for linking Params with the LNK macro
	size_t flags=0;
	for(size_t i:list)flags|=1<<i;
	return flags;
}
Window::Window(QWindow*scene3D,QWidget*parent):QMainWindow(parent){
	setWindowTitle(defaultTitle);
	QPlainTextEdit*console=new QPlainTextEdit();
	QPalette palette=console->palette();
	palette.setColor(QPalette::Base,Qt::black);
	palette.setColor(QPalette::Text,Qt::white);
	console->setPalette(palette);
	console->setFont(QFont("Courier"));
	console->setReadOnly(true);
	connect(&calculation,&Thread::output,console,&QPlainTextEdit::appendPlainText);
	editable->setToolTip("Allows parameter editing; changing any parameters will clear the current results");

	cpus->setSuffix(" CPU cores");
	cpus->setMinimum(1);
	cpus->setMaximum(omp_get_num_procs());
	connect(cpus,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,[](int value){g_max_cpu_cores=value;});
	cpus->setValue(cpus->maximum());
	cpus->setToolTip("Each core will be allocated one OpenCL device, so make sure this is at least as high as your GPUs value");

	gpus->setSuffix(" GPUs");
	gpus->setMaximum(getDeviceCountOpenCL());
	gpus->setMinimum(-gpus->maximum());
	connect(gpus,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,[](int value){g_max_gpus=value;});
	gpus->setValue(gpus->maximum());
	gpus->setToolTip("A negative value indicates that the calculation should use higher-indexed OpenCL devices; useful for machines with multiple types of OpenCL devices");

	QMenu*menu=menuBar()->addMenu("File");
	menu->addAction("New",this,SLOT(clear()),QKeySequence::New);
	menu->addAction("Open",this,SLOT(load()),QKeySequence::Open);
	menu->addAction("Save",this,SLOT(save()),QKeySequence::Save);

	QWidget*base=new QWidget(this);
	setCentralWidget(base);
	QHBoxLayout*baseLayout=new QHBoxLayout(base);
	{//leftLayout scope
		QVBoxLayout*leftLayout=new QVBoxLayout;
		baseLayout->addLayout(leftLayout);
		{//upperLeftLayout scope
			QHBoxLayout*upperLeftLayout=new QHBoxLayout();
			leftLayout->addLayout(upperLeftLayout);
			//QPushButton*dbug = new QPushButton("D");
			//connect(dbug,&QPushButton::clicked,[=]{tif_print();});
			//upperLeftLayout->addWidget(dbug);
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
		ParamGroup*group=0;
		Param*param=0;
		ComboParam*combo=0;
		Param*array[3]={0};
		size_t line=0;
#define USE_ALL_MACROS//Excessive use of macro sorcery below to make main.h readable. Edit main.h itself to alter parameter panels.
#define GROUP(...) group=new ParamGroup(__VA_ARGS__);groupLayout->addWidget(group);push_back(group);
#define MENU(tag,name,var,...) param=combo=new ComboParam(tag,name,&var,{__VA_ARGS__},group);group->add(param,line==__LINE__);connect(param,&Param::edited,this,&Window::onEdit);line=__LINE__;
#define FOREACH(id,name,type,var,...) if(group){param=newParam(id,name,&var,group);group->add(param,line==__LINE__);connect(param,&Param::edited,this,&Window::onEdit);line=__LINE__;}
#define LINK(...) param->link(flagset({__VA_ARGS__}));param->link(0);connect(combo,&ComboParam::currentIndexChanged,param,&Param::link);
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
		scene3D->show();
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
		error->setModal(true);
		connect(&calculation,&Thread::error,error,(void(QErrorMessage::*)(const QString&))&QErrorMessage::showMessage);
		connect(&calculation,&Thread::completed1st,manager,&Spectromanager::load1st);
		connect(&calculation,&Thread::completed1st,check3D,&QCheckBox::setChecked);
		connect(&calculation,&Thread::progress,manager,&Spectromanager::setFrameCount);
		connect(check3D,&QCheckBox::stateChanged,stack,&QStackedLayout::setCurrentIndex);
		check3D->setChecked(true);
	}//end of rightLayout scope
	connect(start,&QPushButton::clicked,[=]{
		start->setEnabled(false);
		pause->setEnabled(true);
		cpus->setEnabled(false);
		gpus->setEnabled(false);
		menuBar()->setEnabled(false);
		editable->setEnabled(false);
		editable->setChecked(false);
		write();
		calculation.begin();
	});
	connect(&calculation,&QThread::finished,this,&Window::onFinish);
	connect(pause,&QPushButton::clicked,[=]{
		pause->setEnabled(false);
		calculation.stop();
	});
	pause->setEnabled(false);
	progress->setTextVisible(true);
	progress->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Fixed);
	connect(&calculation,&Thread::maximum,progress,&QProgressBar::setMaximum);
	connect(&calculation,&Thread::progress,progress,&QProgressBar::setValue);
#ifdef _WIN32
	button->progress()->setVisible(true);
	connect(&calculation,&Thread::maximum,button->progress(),&QWinTaskbarProgress::setMaximum);
	connect(&calculation,&Thread::progress,button->progress(),&QWinTaskbarProgress::setValue);
#endif
	if(!tif_initialize())calculation.error("Could not open file!");
	else read();
	calculation.stop();
	calculation.begin();
}
void Window::clear(){
	if(QMessageBox::Ok==QMessageBox::warning(this,"Confirm","This will clear your current parameters. Are you sure?",QMessageBox::Ok|QMessageBox::Cancel,QMessageBox::Cancel)){
		for(ParamGroup*group:*this)for(Param*param:*group)param->clear();
		manager->clear();
		tif_update();
	}
}
bool Window::load(){
	QString name=QFileDialog::getOpenFileName(this,"Open","","TIFF (*.tif);;All Files (*)");
	if(name.isNull())return false;
	setWindowTitle(defaultTitle+(" - "+QFileInfo(name).fileName()));
	if(!tif_open(name.toUtf8().data())){
		error->showMessage("Failed to open file: "+name);
		return false;
	}
	read();
	calculation.stop();
	calculation.begin();
	manager->setFrame();
	return true;
}
bool Window::save(){
	if(windowTitle()!=defaultTitle)return false;
	QString newName=QFileDialog::getSaveFileName(this,"Save As","","TIFF (*.tif)");
	if(newName.isNull())return false;
	QByteArray array=newName.toUtf8();
	char*utf8name=array.data();
	remove(utf8name);
	if(tif_save_as(utf8name)){
		emit_error("Error: Unable to save file! ",strerror(errno));
		return false;
	}
	setWindowTitle(defaultTitle+(" - "+QFileInfo(newName).fileName()));
	return true;
}
#ifdef _WIN32
void Window::showEvent(QShowEvent*e){
	e=e;
	button->setWindow(windowHandle());
}
#endif
void Window::lock(){
	for(auto group:*this)group->disable();
}
void Window::unlock(){
	for(auto group:*this)group->enable();
}
void Window::read(){
	spacing->read();
	for(auto group:*this)for(auto param:*group)param->read();
}
void Window::write(){
	spacing->write();
	for(auto group:*this)for(auto param:*group)param->write();
}
void Window::onEdit(){
	setWindowTitle(defaultTitle);
	write();
	tif_update();
	progress->setValue(progress->minimum());
	manager->clear();
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