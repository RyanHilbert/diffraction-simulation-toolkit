#include "spectrogram.h"
Spectrogram::Spectrogram(QString name,QWidget*parent):QCustomPlot(parent){
	setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
	QCPColorScale*scale=new QCPColorScale(this);
	map->setColorScale(scale);
	map->setInterpolate(false);
	map->setTightBoundary(true);
	addPlottable(map);
	plotLayout()->addElement(0,1,scale);
	plotLayout()->insertRow(0);
	plotLayout()->addElement(0,0,new QCPPlotTitle(this,name));

	dialog->setWindowTitle(name);
	QLineEdit*upper=new QLineEdit,*lower=new QLineEdit;
	QComboBox*gradients=new QComboBox(dialog);

	QCPAxis*axis=scale->axis();
	connect(axis,(void(QCPAxis::*)(const QCPRange&))&QCPAxis::rangeChanged,[=](const QCPRange&range){
		upper->setText(QString::number(range.upper));
		lower->setText(QString::number(range.lower));
	});
	connect(upper,&QLineEdit::editingFinished,[=]{
		axis->setRangeUpper(upper->text().toDouble());
		replot();
	});
	connect(lower,&QLineEdit::editingFinished,[=]{
		axis->setRangeLower(lower->text().toDouble());
		replot();
	});
	QCheckBox*logBox=new QCheckBox("Logarithmic",dialog);
	connect(logBox,&QCheckBox::stateChanged,[=](int state){
		axis->setScaleType(state?QCPAxis::stLogarithmic:QCPAxis::stLinear);
		replot();
	});
	QFormLayout*layout=new QFormLayout(dialog);
	layout->addRow("Scale Max:",upper);
	layout->addRow("Scale Min:",lower);
	QHBoxLayout*bottom=new QHBoxLayout();
	layout->addRow(bottom);
	bottom->addWidget(logBox);
	bottom->addWidget(gradients);

	connect(gradients,(void(QComboBox::*)(int))(&QComboBox::currentIndexChanged),this,&Spectrogram::setGradient);
	gradients->addItem("Grayscale");
	gradients->addItem("Hot");
	gradients->addItem("Cold");
	gradients->addItem("Night");
	gradients->addItem("Candy");
	gradients->addItem("Geography");
	gradients->addItem("Ion");
	gradients->addItem("Thermal");
	gradients->addItem("Polar");
	gradients->addItem("Spectrum");
	gradients->addItem("Jet");
	gradients->addItem("Hues");
	gradients->setMaxVisibleItems(gradients->count());
	gradients->setCurrentIndex(8);

	QMenu*menu=new QMenu(this);
	menu->addAction("Show Controls",dialog,SLOT(show()));
	menu->addAction("Reset View",this,SLOT(resetView()));
	menu->addAction("Reset Scale",this,SLOT(resetScale()));
	menu->addAction("Copy to Clipboard",this,SLOT(clipboard()));
	menu->addAction("Save Image As...",this,SLOT(saveImage()));
	setContextMenuPolicy(Qt::CustomContextMenu);
	connect(this,&QWidget::customContextMenuRequested,[=](const QPoint&pos){menu->popup(mapToGlobal(pos));});

	QLabel*label=new QLabel("0",this);
	label->setMinimumWidth(127);
	connect(this,&QCustomPlot::mouseMove,[=](QMouseEvent*event){
		label->setText(QString::number(data->data(xAxis->pixelToCoord(event->x()),yAxis->pixelToCoord(event->y()))));
	});
}
void Spectrogram::clear(){
	map->clearData();
	replot();
}
bool Spectrogram::saveImage(){
	QString formats="";
	for(QByteArray format:QImageWriter::supportedImageFormats())formats+=format.toUpper()+" (*."+format+");;";
	formats.chop(2);
	QString format="PNG (*.png)";
	const QString name=QFileDialog::getSaveFileName(this,"Save Image","",formats,&format);
	if(name.isNull())return false;
	saveRastered(name,0,0,1,QFileInfo(name).suffix().toUtf8().data());
	return true;
}
void Spectrogram::clipboard(){
	QApplication::clipboard()->setPixmap(toPixmap());
}
void Spectrogram::resetScale(){
	map->rescaleDataRange(true);
	replot();
}
void Spectrogram::resetView(){
	rescaleAxes();
	replot();
}
double Spectrogram::setRangeLower(const double range){
	const double oldRange=yAxis2->range().lower;
	yAxis2->setRangeLower(range);
	replot();
	return oldRange;
}
double Spectrogram::setRangeUpper(const double range){
	const double oldRange=yAxis2->range().upper;
	yAxis2->setRangeLower(range);
	replot();
	return oldRange;
}
QCPRange Spectrogram::setRange(const QCPRange range){
	return QCPRange(setRangeLower(range.lower),setRangeUpper(range.upper));
}
void Spectrogram::setRanges(double xmin, double xmax, double ymin, double ymax){
	data->setKeyRange(QCPRange(xmin,xmax));
	data->setValueRange(QCPRange(ymin,ymax));
}
void Spectrogram::setSize(int x,int y){
	data->setSize(x,y);
}
void Spectrogram::setGradient(int gradient){
	map->setGradient((QCPColorGradient::GradientPreset)gradient);
	replot();
}