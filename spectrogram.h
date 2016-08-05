#ifndef SPECTROGRAM_H
#define SPECTROGRAM_H
#include"qcustomplot.h"
#include"h5.h"
class Spectrogram:public QCustomPlot{
	Q_OBJECT
	QDialog*dialog=new QDialog(this);
protected:
	QCPColorMap*map=new QCPColorMap(xAxis,yAxis);
	QCPColorMapData*data=map->data();
public:
	Spectrogram(QString name="",QWidget*parent=0);
	virtual void load(struct h5dataset&data)=0;
signals:
	void loaded(int);
public slots:
	void clear();
	bool saveImage();
	void clipboard();
	void resetScale();
	void resetView();
	double setRangeLower(double);
	double setRangeUpper(double);
	QCPRange setRange(QCPRange);
	void setRanges(double xmin,double xmax,double ymin,double ymax);
	void setSize(int x,int y);
	void setGradient(int gradient);
};
#endif