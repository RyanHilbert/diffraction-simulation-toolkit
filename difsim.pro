CONFIG += c++11
QT += widgets printsupport 3dcore 3drender 3dextras
win32:QT += winextras

win32:QMAKE_CXXFLAGS +=  /openmp
unix: QMAKE_CXXFLAGS += -fopenmp
unix: QMAKE_CFLAGS   += -fopenmp

RESOURCES += shaders.qrc
HEADERS += *.h
win32:SOURCES += calc.cpp
unix:SOURCES += shaders/calc.c
SOURCES += coloredlabel.cpp dialog.cpp dialogbutton.cpp directioncalculator.cpp main.cpp param.cpp paramgroup.cpp qcustomplot.cpp spacingcalculator.cpp spectrogram.cpp spectromanager.cpp thread.cpp tinyexpr.c vector3d.cpp window.cpp window3d.cpp

DISTFILES += configure README.md LICENSE.txt debian/*
win32:LIBS += -L$$(CUDA_PATH)/lib/x64 -L$$(AMDAPPSDKROOT)/lib/x86_64 -L$$(INTELOCLSDKROOT)/lib64 -lOpenCL
unix:LIBS += -fopenmp -lrt -ldl

win32:INCLUDEPATH += $$(CUDA_PATH)/include $$(AMDAPPSDKROOT)/include $$(INTELOCLSDKROOT)/include