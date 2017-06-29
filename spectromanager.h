#ifndef SPECTROMANAGER_H
#define SPECTROMANAGER_H
#include<QFrame>
#include"spectrogram.h"
#include"tif.h"
class Spectromanager:public QFrame{
	QGridLayout*const grid=new QGridLayout(this);
	QSlider*const frameSlider=new QSlider(Qt::Horizontal,this),*const sliceSlider=new QSlider(Qt::Horizontal,this);
	QLabel*const beamLabel=new QLabel("Beam Position: 0x, 0y, 0z"),*const angleLabel=new QLabel("Deviation Angle: 0");
	Spectrogram*det,*poynting,*dh,*d0;
	wavefield_iterator waves={0,0,0,0};
	bool focused=false;
	bool setSliceCount(size_t count=0);
public:
	Spectromanager(QWidget*parent=0);
public slots:
	int getFrameMax();
	int getSliceMax();
	int getFrame();
	int getSlice();
	bool load1st();
	bool setFrame(int i=0);
	bool setSlice(int i=0);
	bool setFrameCount(size_t count=0);
	bool clear();
private slots:
	bool focus(QWidget*widget=0);
};
#endif