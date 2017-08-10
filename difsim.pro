CONFIG += c++11
QT += widgets printsupport
win32:QT += winextras 3dcore 3drender 3dextras

win32:QMAKE_CXXFLAGS +=  /openmp
unix: QMAKE_CXXFLAGS += -fopenmp
unix: QMAKE_CFLAGS   += -fopenmp -std=c99

RESOURCES += shaders.qrc
win32:HEADERS += vector3d.h window3d.h
HEADERS += coloredlabel.h dialog.h dialogbutton.h directioncalculator.h err.h main.h param.h paramgroup.h qcustomplot.h spacingcalculator.h spectrogram.h spectromanager.h thread.h tif.h tinyexpr.h window.h
SOURCES += coloredlabel.cpp dialog.cpp dialogbutton.cpp directioncalculator.cpp main.cpp param.cpp paramgroup.cpp qcustomplot.cpp spacingcalculator.cpp spectrogram.cpp spectromanager.cpp thread.cpp tinyexpr.c window.cpp
win32:SOURCES += calc.cpp vector3d.cpp window3d.cpp
unix:SOURCES += shaders/calc.c

DISTFILES += configure README.md LICENSE.txt debian/source/format debian/*
win32:INCLUDEPATH += $$(CUDA_PATH)/include $$(AMDAPPSDKROOT)/include $$(INTELOCLSDKROOT)/include
win32:LIBS += -L$$(CUDA_PATH)/lib/x64 -L$$(AMDAPPSDKROOT)/lib/x86_64 -L$$(INTELOCLSDKROOT)/lib64
unix:LIBS += -fopenmp -lrt -ldl
LIBS += -lOpenCL
