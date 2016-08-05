#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include<QMainWindow>
#include"param.h"
#include"paramgroup.h"
#include"spectromanager.h"
#include"dialogbutton.h"
#include"directioncalculator.h"
#include"spacingcalculator.h"
class Window:public QMainWindow,private std::vector<ParamGroup*>{
	Q_OBJECT
	Spectromanager*manager=new Spectromanager(this);
	QProgressBar*progress=new QProgressBar(this);
	SpacingCalculator*spacing=new SpacingCalculator(this);
	DirectionCalculator*direction=new DirectionCalculator(this);
	DialogButton*spacingButton;
	DialogButton*directionButton;
	QCheckBox*editable=new QCheckBox("Allow Editing");
	QPushButton*start=new QPushButton("Start Calculation");
	QPushButton*pause=new QPushButton("Pause Calculation");
	QSpinBox*cpus=new QSpinBox(this);
	QSpinBox*gpus=new QSpinBox(this);
	QCheckBox*check3D=new QCheckBox("3D",this);
public:
	const char*defaultTitle="xDiffraction";
	const char*defaultFilename=".~";//simplest name for hidden temporary file
	Window(QWindow*scene3D=0,QWidget*parent=0);
	~Window();
	using std::vector<ParamGroup*>::begin;
	using std::vector<ParamGroup*>::end;
signals:
	void edited();
public slots:
	void clear();
	bool exportiff();
	bool load();
	bool save();
	bool hdf5read();
	bool hdf5write();
	void lock();
	void unlock();
	void upd8();
private slots:
	void onEdit();
	void onFinish();
};
#endif