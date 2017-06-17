#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#ifdef _WIN32
#include<QtWinExtras>
#endif
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
	QErrorMessage*error=QErrorMessage::qtHandler();
#ifdef _WIN32
	QWinTaskbarButton*button=new QWinTaskbarButton(this);
	void showEvent(QShowEvent*)override;
#endif
public:
	const char*defaultTitle="Diffraction Simulation Toolkit";
	Window(QWindow*scene3D=0,QWidget*parent=0);
	using std::vector<ParamGroup*>::begin;
	using std::vector<ParamGroup*>::end;
signals:
	void edited();
public slots:
	void clear();
	bool load();
	bool save();
	void lock();
	void unlock();
	void read();
	void write();
private slots:
	void onEdit();
	void onFinish();
};
#endif