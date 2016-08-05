#include"spectromanager.h"
#include"h5.h"
#include"main.h"
#include<QtWidgets>
extern"C"{extern hid_t hdf5file;}

static void fill(QCPColorMapData*data){//corrects gaps in wave field graphs
	const int xsize=data->keySize(),ysize=data->valueSize();
	float floats[xsize][ysize];
	for(int x=0;x<xsize;++x)for(int y=0;y<ysize;++y)
		if(data->cell(x,y)==0)floats[x][y]=fmax(fmax(data->cell(x-1,y),data->cell(x+1,y)),fmax(data->cell(x,y-1),data->cell(x,y+1)));
		else floats[x][y]=0;
	for(int x=0;x<xsize;++x)for(int y=0;y<ysize;++y)if(floats[x][y]!=0){
		data->setCell(x,y,floats[x][y]);
	}
}
Spectromanager::Spectromanager(QWidget*parent):QFrame(parent){
	enum{xpoynting_index, ypoynting_index, x_index, y_index, dh_index, d0_index};//index of each piece of data

	class D:public Spectrogram{
		using Spectrogram::Spectrogram;
		void load(struct h5dataset&dataset){
			//yAxis->setRangeReversed(true);
			data->setSize(0,0);
			data->setSize(dataset.columns,dataset.rows);
			data->setRange(QCPRange(0,dataset.columns),QCPRange(0,dataset.rows));
			for(size_t row=0;row<dataset.rows;++row)for(size_t column=0;column<dataset.columns;++column)data->setCell(column,row,dataset[row][column]);
			replot();
		}
	};
	class P:public Spectrogram{
		using Spectrogram::Spectrogram;
		void load(struct h5dataset&dataset){
			yAxis->setRangeReversed(true);
			for(size_t row=0;row<dataset.rows;++row){
				const double x=dataset[row][x_index],y=dataset[row][y_index],xvalue=dataset[row][xpoynting_index],yvalue=dataset[row][ypoynting_index],value=sqrt(xvalue*xvalue+yvalue*yvalue);
				int cellx,celly;
				data->coordToCell(x,y,&cellx,&celly);
				if(data->cell(cellx,celly)<value)data->setCell(cellx,celly,value);
			}
			fill(data);
			replot();
		}
	};
	class H:public Spectrogram{
		using Spectrogram::Spectrogram;
		void load(struct h5dataset&dataset){
			yAxis->setRangeReversed(true);
			for(size_t row=0;row<dataset.rows;++row){
				const double x=dataset[row][x_index],y=dataset[row][y_index],value=dataset[row][dh_index];
				int cellx,celly;
				data->coordToCell(x,y,&cellx,&celly);
				if(data->cell(cellx,celly)<value)data->setCell(cellx,celly,value);
			}
			fill(data);
			replot();
		}
	};
	class O:public Spectrogram{
		using Spectrogram::Spectrogram;
		void load(struct h5dataset&dataset){
			yAxis->setRangeReversed(true);
			for(size_t row=0;row<dataset.rows;++row){
				const double x=dataset[row][x_index],y=dataset[row][y_index],value=dataset[row][d0_index];
				int cellx,celly;
				data->coordToCell(x,y,&cellx,&celly);
				if(data->cell(cellx,celly)<value)data->setCell(cellx,celly,value);
			}
			fill(data);
			replot();
		}
	};
	det=new D("Detector");
	poynting=new P("Poynting");
	dh=new H("Diffraction");
	d0=new O("Transmission");
	connect(det,&Spectrogram::mouseDoubleClick,[=]{focus(det);});
	connect(poynting,&Spectrogram::mouseDoubleClick,[=]{focus(poynting);});
	connect(dh,&Spectrogram::mouseDoubleClick,[=]{focus(dh);});
	connect(d0,&Spectrogram::mouseDoubleClick,[=]{focus(d0);});

	QSpinBox*frameSpinner=new QSpinBox(this);
	QSpinBox*sliceSpinner=new QSpinBox(this);
	frameSpinner->setMinimumWidth(50);
	sliceSpinner->setMinimumWidth(50);
	connect(frameSlider,&QSlider::rangeChanged,frameSpinner,&QSpinBox::setRange);
	connect(sliceSlider,&QSlider::rangeChanged,sliceSpinner,&QSpinBox::setRange);
	connect(frameSlider,&QSlider::valueChanged,frameSpinner,&QSpinBox::setValue);
	connect(sliceSlider,&QSlider::valueChanged,sliceSpinner,&QSpinBox::setValue);
	connect(frameSpinner,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,frameSlider,&QSlider::setValue);
	connect(sliceSpinner,(void(QSpinBox::*)(int))&QSpinBox::valueChanged,sliceSlider,&QSlider::setValue);

	frameSlider->setMaximum(0);
	sliceSlider->setMaximum(0);
	frameSlider->setTickInterval(1);
	sliceSlider->setTickInterval(1);
	frameSlider->setTickPosition(QSlider::TicksBothSides);
	sliceSlider->setTickPosition(QSlider::TicksBothSides);
	connect(frameSlider,&QSlider::valueChanged,this,&Spectromanager::setFrame);
	connect(sliceSlider,&QSlider::valueChanged,this,&Spectromanager::setSlice);

	QHBoxLayout*frameLayout=new QHBoxLayout;
	QHBoxLayout*sliceLayout=new QHBoxLayout;
	frameLayout->addWidget(frameSpinner);
	sliceLayout->addWidget(sliceSpinner);
	frameLayout->addWidget(new QLabel("Frame\t"));
	sliceLayout->addWidget(new QLabel("Slice\t"));
	frameLayout->addWidget(frameSlider);
	sliceLayout->addWidget(sliceSlider);
	setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
	beamLabel->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Maximum);
	angleLabel->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Maximum);
	grid->addWidget(beamLabel,2,0,1,1);
	grid->addWidget(angleLabel,2,1,1,1);
	grid->addLayout(frameLayout,3,0,1,2);
	grid->addLayout(sliceLayout,4,0,1,2);
	setFrameShape(QFrame::Box);
	focus(det);
}
int Spectromanager::getFrame(){return frameSlider->value();}
int Spectromanager::getSlice(){return sliceSlider->value();}
int Spectromanager::getFramemax(){return frameSlider->maximum();}
int Spectromanager::getSliceMax(){return sliceSlider->maximum();}

bool Spectromanager::load1st(){
	setFrame(0);
	const bool result=setSlice(0);
	det->resetScale();
	det->resetView();
	poynting->resetScale();
	poynting->resetView();
	dh->resetScale();
	dh->resetView();
	d0->resetScale();
	d0->resetView();
	if(g_saving_wavefield){
		if(focused)focus();
	}
	else{
		if(!focused)focus(det);
	}
	return result;
}
bool Spectromanager::setFrame(int frame){
	if(frame<0)return setFrame(0);
	int max=getFramemax();
	if(frame>max)return setFrame(max);
	frameSlider->setValue(frame);
	char dataset_name[32];
	char slice_group_name[64];
	sprintf(dataset_name,"det%i",frame);
	sprintf(slice_group_name,"det%ifield",frame);
	hid_t det_group=h5group("Detector",hdf5file);
	if(!h5exists(dataset_name,det_group)){
#pragma omp critical(OMP_CRIT_HDF5)
		{H5Oclose(det_group);}
		return clear();
	}
	hid_t data_id=h5(dataset_name,det_group);
	double position[3];
#pragma omp critical(OMP_CRIT_HDF5)
	{H5LTget_attribute_double(data_id,".","Beam center",position);}
	beamLabel->setText("Beam position: "+QString::number(position[0])+" x, "+QString::number(position[1])+" y, "+QString::number(position[2])+" z");
	angleLabel->setText("Deviation angle: "+QString::number(h5double_("Deviation angle",data_id)));
	struct h5dataset data=h5dataset(dataset_name,det_group);
	det->load(data);
	if(h5exists(slice_group_name,det_group)){
		hid_t slice_group=h5(slice_group_name,det_group);
		double xmin=h5double("x min",slice_group);
		double xmax=h5double("x max",slice_group);
		double ymin=h5double("z min",slice_group);
		double ymax=h5double("z max",slice_group);
		int slices=h5int_("slices",slice_group);
		setSliceCount(slices);
		poynting->setRanges(xmin,xmax,ymin,ymax);
		dh->setRanges(xmin,xmax,ymin,ymax);
		d0->setRanges(xmin,xmax,ymin,ymax);
		setSlice(getSlice());
	}
#pragma omp critical(OMP_CRIT_HDF5)
	{H5Oclose(det_group);}
	return true;
}
bool Spectromanager::setSlice(int slice){
	int max=getSliceMax();
	if(slice>max)return setSlice(max);
	const int frame=getFrame();
	char name[64];
	sprintf(name,"Detector/det%ifield/slice%i",frame,slice);
	if(slice<0||!h5exists(name,hdf5file)){
		poynting->clear();
		dh->clear();
		d0->clear();
		return false;
	}
	char group_name[32];
	char slice_name[32];
	sprintf(group_name,"det%ifield",frame);
	sprintf(slice_name,"slice%i",slice);
	hid_t slice_group=h5_(group_name,h5("Detector",hdf5file));
	hid_t slice_id=h5(slice_name,slice_group);
	const int xpixels=h5int("x pixels",slice_id);
	const int zpixels=h5int_("z pixels",slice_id);
	struct h5dataset data=h5dataset_(slice_name,slice_group);
	poynting->setSize(0,0);
	dh->setSize(0,0);
	d0->setSize(0,0);
	poynting->setSize(xpixels,zpixels);
	dh->setSize(xpixels,zpixels);
	d0->setSize(xpixels,zpixels);
	poynting->load(data);
	dh->load(data);
	d0->load(data);
	sliceSlider->setValue(slice);
	return true;
}
bool Spectromanager::setFrameCount(size_t count){
	frameSlider->setMaximum(count<1?0:count-1);
	return true;
}
bool Spectromanager::setSliceCount(size_t count){
	sliceSlider->setMaximum(count<1?0:count-1);
	return true;
}
bool Spectromanager::focus(QWidget*widget){
	if(focused){
		grid->addWidget(det,0,0);
		grid->addWidget(poynting,0,1);
		grid->addWidget(dh,1,0);
		grid->addWidget(d0,1,1);
	}else{
		det->setParent(0);
		poynting->setParent(0);
		dh->setParent(0);
		d0->setParent(0);
		grid->addWidget(widget,0,0,2,2);
	}
	return focused=!focused;
}
bool Spectromanager::clear(){
	det->clear();
	poynting->clear();
	dh->clear();
	d0->clear();
	setFrameCount(0);
	setSliceCount(0);
	return true;
}