QT += widgets printsupport 3dcore 3drender 3dextras

QMAKE_CFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp

#library source files
SOURCES += tinyexpr.c qcustomplot.cpp
HEADERS += tinyexpr.h qcustomplot.h

#source files
SOURCES += calc.c thread.cpp param.cpp window.cpp coloredlabel.cpp dialog.cpp dialogbutton.cpp paramgroup.cpp directioncalculator.cpp spacingcalculator.cpp spectrogram.cpp spectromanager.cpp main.cpp window3d.cpp vector3d.cpp
HEADERS += h5.h main.h thread.h param.h window.h coloredlabel.h dialog.h dialogbutton.h paramgroup.h directioncalculator.h spacingcalculator.h spectrogram.h spectromanager.h window3d.h vector3d.h

#source files to compile separately with NVCC
OTHER_FILES += kernel.cu
unix:OBJECTS += kernel.o

LIBS += -lhdf5 -lhdf5_hl -fopenmp

#specify the location of 32-bit HDF5 library files on Windows
win32:INCLUDEPATH += "C:\Program Files (x86)\HDF_Group\HDF5\1.8.17\include"
win32:LIBS += -L"C:\Program Files (x86)\HDF_Group\HDF5\1.8.17\lib"

unix:INCLUDEPATH += /usr/include/hdf5/serial
unix:LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial

#unix compilation assumes CUDA support
unix:LIBS += -lcuda -lcudart
unix:DEFINES += CUDA

DISTFILES += GPL.txt HDF5.txt
