CONFIG += c++11
QT += widgets printsupport 3dcore 3drender 3dextras
win32:QT += winextras

win32:QMAKE_CXXFLAGS +=  /openmp
unix: QMAKE_CXXFLAGS += -fopenmp
unix: QMAKE_CFLAGS   += -fopenmp -std=c99

RESOURCES += shaders.qrc
HEADERS += *.h
unix:SOURCES += shaders/calc.c
win32:SOURCES += calc.cpp
SOURCES += coloredlabel.cpp dialog.cpp dialogbutton.cpp directioncalculator.cpp main.cpp param.cpp paramgroup.cpp qcustomplot.cpp spacingcalculator.cpp spectrogram.cpp spectromanager.cpp thread.cpp tinyexpr.c vector3d.cpp window.cpp window3d.cpp

#magic needed to make CUDA work
DISTFILES += configure kernel.cu README.md LICENSE.txt debian/*
LIBS += $$OUT_PWD/kernel$$QMAKE_EXT_OBJ -lcuda -lcudart_static
win32:LIBS += -L$$(CUDA_PATH)/lib/x64
unix:LIBS += -fopenmp -lrt -ldl

kernel.target += $$OUT_PWD/kernel$$QMAKE_EXT_OBJ
kernel.commands += nvcc -c -o $$OUT_PWD/kernel$$QMAKE_EXT_OBJ $$PWD/kernel.cu
QMAKE_EXTRA_TARGETS += kernel
PRE_TARGETDEPS += $$OUT_PWD/kernel$$QMAKE_EXT_OBJ
