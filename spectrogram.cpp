#include "spectrogram.h"
Spectrogram::Spectrogram(QString name,QWidget*parent):QCustomPlot(parent){
	setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);
	map->setColorScale(scale);
	map->setInterpolate(false);
	map->setTightBoundary(true);
	plotLayout()->addElement(0,1,scale);
	plotLayout()->insertRow(0);
	data->setRange(QCPRange(0,1),QCPRange(0,1));
	QFont titleFont;
	titleFont.setWeight(99);
	titleFont.setPointSize(titleFont.pointSize()*1.5);
	plotLayout()->addElement(0,0,new QCPTextElement(this,name,titleFont));

	dialog->setWindowTitle(name);
	QLineEdit*upper=new QLineEdit,*lower=new QLineEdit;
	QComboBox*gradients=new QComboBox(dialog);

	QCPAxis*axis=scale->axis();
	connect(axis,(void(QCPAxis::*)(const QCPRange&))&QCPAxis::rangeChanged,[=](const QCPRange&range){
		upper->setText(QString::number(range.upper));
		lower->setText(QString::number(range.lower));
	});
	connect(upper,&QLineEdit::editingFinished,[=]{
		autoscale->setChecked(false);
		axis->setRangeUpper(upper->text().toDouble());
		replot();
	});
	connect(lower,&QLineEdit::editingFinished,[=]{
		autoscale->setChecked(false);
		axis->setRangeLower(lower->text().toDouble());
		replot();
	});
	QCheckBox*logBox=new QCheckBox("Log",dialog);
	connect(logBox,&QCheckBox::stateChanged,[=](int state){
		axis->setScaleType(state?QCPAxis::stLogarithmic:QCPAxis::stLinear);
		replot();
	});
	logBox->setToolTip("Enables logarithmic scaling");
	autoscale->setToolTip("Attempts to automatically detect reasonable max/min scale values");
	autoscale->setChecked(true);
	QFormLayout*layout=new QFormLayout(dialog);
	layout->addRow("Scale Max:",upper);
	layout->addRow("Scale Min:",lower);
	QHBoxLayout*bottom=new QHBoxLayout();
	layout->addRow(bottom);
	bottom->addWidget(autoscale);
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
	menu->addAction("Export to TIFF...",this,SLOT(tiffExport()));
	menu->addAction("Save Image As...",this,SLOT(saveImage()));
	setContextMenuPolicy(Qt::CustomContextMenu);
	connect(this,&QWidget::customContextMenuRequested,[=](const QPoint&pos){menu->popup(mapToGlobal(pos));});

	QLabel*label=new QLabel("0",this);
	label->setMinimumWidth(127);
	connect(this,&QCustomPlot::mouseMove,[=](QMouseEvent*event){
		label->setText(QString::number(data->data(xAxis->pixelToCoord(event->x()),yAxis->pixelToCoord(event->y()))));
	});
}
static int comp(const void*v1,const void*v2){
	if(*(float*)v1 < *(float*)v2) return -1;
	else if(*(float*)v1 > *(float*)v2) return 1;
	else return 0;
}
void Spectrogram::adjustScale(){
	const int keySize = data->keySize(), valueSize = data->valueSize();
	float*values = (float*)malloc(keySize*valueSize*sizeof*values);
	for(int key=0; key<keySize; ++key){
		for(int value=0; value<valueSize; ++value){
			values[key*valueSize+value] = data->cell(key,value);
		}
	}
	qsort(values, keySize*valueSize, sizeof*values, &comp);
	int n = 0;//filter out zeroes
	while(n < keySize*valueSize && values[n] <= 0) ++n;
	scale->axis()->setRange((values+n)[(keySize*valueSize-n)>>3],values[keySize*valueSize-1]);
	free(values);
	replot();
}
void Spectrogram::clear(){
	data->clear();
	replot();
}
bool Spectrogram::saveImage(){
	QString formats="";
	for(QByteArray format:QImageWriter::supportedImageFormats()){
		formats += format.toUpper() + " (*."+format+");;";
	}
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
	if(autoscale->isChecked())adjustScale();
	else{
		map->rescaleDataRange(true);
		replot();
	}
}
void Spectrogram::resetView(){
	rescaleAxes();
	replot();
}
void Spectrogram::setSize(int x,int y){
	data->setSize(x,y);
}
void Spectrogram::setGradient(int gradient){
	map->setGradient((QCPColorGradient::GradientPreset)gradient);
	replot();
}