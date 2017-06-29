#include"spectromanager.h"
#include<QtWidgets>

Spectromanager::Spectromanager(QWidget*parent):QFrame(parent){
    enum{xpoynting_index, ypoynting_index, x_index, y_index, dh_index, d0_index};//index of each piece of data

    class D:public Spectrogram{
        using Spectrogram::Spectrogram;
        void load(size_t columns,size_t rows,float*f){
            //yAxis->setRangeReversed(true);
            data->setSize(0,0);
            data->setSize(columns,rows);
            data->setRange(QCPRange(0,columns),QCPRange(0,rows));
            for(int row=0;row<rows;++row){
                const size_t index=row*columns;
                for(size_t column=0;column<columns;++column){
                    float cell = f[index+column];
                    if(BIG_ENDIAN)REVERSE_BYTES(cell);
                    data->setCell(column,row,cell);
                }
            }
            if(autoscale->isChecked()) adjustScale();
            else replot();
        }
        bool tiffExport(){
            const QString name=QFileDialog::getSaveFileName(this,"Save Image","","TIFF (*.tif)");
            if(name.isNull())return false;
            return export_detector(name.toUtf8().data());
        }
    };
#define ARROW_DENSITY 8
    class P:public Spectrogram{
        QCPItemLine*arrows[ARROW_DENSITY+1][ARROW_DENSITY+1];
        void load(size_t cols,size_t rows,float*f){
            yAxis->setRangeReversed(true);
            size_t count = 0;
            float min = INFINITY;
            float max = -INFINITY;
            for(size_t row=0;row<rows;++row){
                const float thisRow = row*ARROW_DENSITY/(float)(rows-1);
                const float prevRow = row? (row-1)*ARROW_DENSITY/(float)(rows-1): thisRow;
                const float nextRow = row<rows-1? (row+1)*ARROW_DENSITY/(float)(rows-1): thisRow;
                const size_t closestRow = roundf(thisRow);
                const float rowDif = fabsf(thisRow-closestRow);
                const bool isClosestRow = rowDif<=fabsf(prevRow-closestRow) && rowDif<=fabsf(nextRow-closestRow);
                const size_t index=row*cols;
                for(size_t col=0;col<cols;++col){
                    pvector pv = ((pvector*)f)[index+col];
                    if(BIG_ENDIAN){
                        REVERSE_BYTES(pv.x);
                        REVERSE_BYTES(pv.z);
                    }
                    const float val = sqrtf(pv.x*pv.x + pv.z*pv.z);
                    data->setCell(col,row,val);
                    if(isClosestRow){
                        const float thisCol = col*ARROW_DENSITY/(float)(cols-1);
                        const float prevCol = col? (col-1)*ARROW_DENSITY/(float)(cols-1): thisCol;
                        const float nextCol = col<cols-1? (col+1)*ARROW_DENSITY/(float)(cols-1): thisRow;
                        const size_t closestCol = roundf(thisCol);
                        const float colDif = fabsf(thisCol - closestCol);
                        if(colDif<=fabsf(prevCol-closestCol) && colDif<=fabsf(nextCol-closestCol)){
                            const float mag = sqrtf(pv.x*pv.x+pv.z*pv.z);
                            const float x = col/(float)(cols-1);
                            const float y = row/(float)(rows-1);
                            QCPItemLine*arrow = arrows[closestRow][closestCol];
                            arrow->start->setCoords(x-pv.x/mag*.03,y-pv.z/mag*.03);
                            arrow->end->setCoords(x+pv.x/mag*.03,y+pv.z/mag*.03);
                        }
                    }
                }
            }
            if(autoscale->isChecked()) adjustScale();
            else replot();
        }
        bool tiffExport(){
            const QString name=QFileDialog::getSaveFileName(this,"Save Image","","TIFF (*.tif)");
            if(name.isNull())return false;
            return export_wavefield(((Spectromanager*)this->parentWidget())->waves,0,0,name.toUtf8().data());
        }
        public:P(QString name="",QWidget*parent=0):Spectrogram(name,parent){
            for(size_t x=0;x<=ARROW_DENSITY;++x){
                for(size_t y=0;y<=ARROW_DENSITY;++y){
                    arrows[x][y] = new QCPItemLine(this);
                    arrows[x][y]->setHead(QCPLineEnding::esSpikeArrow);
                    arrows[x][y]->start->setType(QCPItemPosition::ptPlotCoords);
                    arrows[x][y]->end->setType(QCPItemPosition::ptPlotCoords);
                }
            }
        }
    };
#undef ARROW_DENSITY
    class H:public Spectrogram{
        using Spectrogram::Spectrogram;
        void load(size_t columns,size_t rows,float*f){
            yAxis->setRangeReversed(true);
            for(size_t row=0;row<rows;++row){
                const size_t index=row*columns;
                for(size_t column=0;column<columns;++column){
                    float cell = f[index+column];
                    if(BIG_ENDIAN)REVERSE_BYTES(cell);
                    data->setCell(column,row,cell);
                }
            }
            if(autoscale->isChecked()) adjustScale();
            else replot();
        }
        bool tiffExport(){
            const QString name=QFileDialog::getSaveFileName(this,"Save Image","","TIFF (*.tif)");
            if(name.isNull())return false;
            return export_wavefield(((Spectromanager*)this->parentWidget())->waves,0,name.toUtf8().data(),0);
        }
    };
    class O:public Spectrogram{
        using Spectrogram::Spectrogram;
        void load(size_t columns,size_t rows,float*f){
            yAxis->setRangeReversed(true);
            for(size_t row=0;row<rows;++row){
                const size_t index=row*columns;
                for(size_t column=0;column<columns;++column){
                    float cell = f[index+column];
                    if(BIG_ENDIAN)REVERSE_BYTES(cell);
                    data->setCell(column,row,cell);
                }
            }
            if(autoscale->isChecked()) adjustScale();
            else replot();
        }
        bool tiffExport(){
            const QString name=QFileDialog::getSaveFileName(this,"Save Image","","TIFF (*.tif)");
            if(name.isNull())return false;
            return export_wavefield(((Spectromanager*)this->parentWidget())->waves,name.toUtf8().data(),0,0);
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
int Spectromanager::getFrameMax(){return frameSlider->maximum();}
int Spectromanager::getSliceMax(){return sliceSlider->maximum();}

bool Spectromanager::load1st(){
    setFrame(0);
    setSlice(0);
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
    return true;
}
bool Spectromanager::setFrame(int frame){
    int max = getFrameMax();
    if(frame < 0) frame = 0;
    if(frame > max) frame = max;
    frameSlider->setValue(frame);
    const vector beam_center = {
        g_scan_type?
            g_scan_x1 + (g_scan_x2 - g_scan_x1)/g_scan_nx*((frame % (g_scan_nx*g_scan_ny))% g_scan_nx) + g_r0.x:
            g_spiral_c*sqrt(frame % g_num_scan)*cos(2.4*(frame % g_num_scan)) + g_r0.x,
        g_scan_type?
            g_scan_y1 + (g_scan_y2 - g_scan_y1)/g_scan_ny*floor((frame % (g_scan_nx*g_scan_ny))/g_scan_nx) + g_r0.y:
            g_spiral_c*sqrt(frame % g_num_scan)*sin(2.4*(frame % g_num_scan)) + g_r0.y,
        g_r0.z
    };
    beamLabel->setText("Beam position: "+QString::number((float)beam_center.x)+" x, "+QString::number((float)beam_center.y)+" y, "+QString::number((float)beam_center.z)+" z");
    angleLabel->setText("Deviation angle: "+QString::number(frame/g_num_scan*g_dth_step+g_dth_start));
    float*data = tif_read_ifd(frame,0);
    det->load(g_column,g_row,data);
    free(data);
    if(g_saving_wavefield){
        free(waves.data);
        waves = tif_read_wavefield(frame);
        setSliceCount(waves.total_wavefields);
        setSlice(getSlice());
    }
    return true;
}
bool Spectromanager::setSlice(int slice){
    if(!g_saving_wavefield) return false;
    int max = getSliceMax();
    if(slice < 0) slice = 0;
    if(slice > max) slice = max;
    waves = get_wavefield(slice,waves);
    const wavefield wave = waves.current_wavefield;
    const size_t cols = wave.col;
    const size_t rows = wave.row;
    poynting->setSize(0,0);
    dh->setSize(0,0);
    d0->setSize(0,0);
    poynting->setSize(cols,rows);
    dh->setSize(cols,rows);
    d0->setSize(cols,rows);
    poynting->load(cols,rows,(float*)wave.poynting);
    dh->load(cols,rows,wave.dh);
    d0->load(cols,rows,wave.d0);
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