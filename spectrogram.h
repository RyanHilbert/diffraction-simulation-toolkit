#ifndef SPECTROGRAM_H
#define SPECTROGRAM_H
#include"qcustomplot.h"
class Spectrogram:public QCustomPlot{
	Q_OBJECT
	QDialog*dialog=new QDialog(this);
protected:
	QCheckBox*autoscale=new QCheckBox("Auto",dialog);
	QCPColorMap*map=new QCPColorMap(xAxis,yAxis);
	QCPColorMapData*data=map->data();
	QCPColorScale*scale=new QCPColorScale(this);
public:
	Spectrogram(QString name="",QWidget*parent=0);
	virtual void load(size_t,size_t,float*)=0;
public slots:
	virtual bool tiffExport()=0;
	void adjustScale();
	void clear();
	bool saveImage();
	void clipboard();
	void resetScale();
	void resetView();
	void setSize(int x,int y);
	void setGradient(int gradient);
};
#endif