CONFIG += c++11
QT += widgets printsupport 3dcore 3drender 3dextras
win32:QT += winextras

win32:QMAKE_CXXFLAGS +=  /openmp
unix: QMAKE_CXXFLAGS += -fopenmp
unix: QMAKE_CFLAGS   += -fopenmp -std=c99

RESOURCES += shaders.qrc
SOURCES += *.c *.cpp
HEADERS += *.h

#magic needed to make CUDA work
DISTFILES += configure kernel.cu README.md LICENSE.txt debian/*
LIBS += -fopenmp -lcuda -lcudart_static $$OUT_PWD/kernel$$QMAKE_EXT_OBJ
win32:LIBS += -L$$(CUDA_PATH)/lib/x64

kernel.target += $$OUT_PWD/kernel$$QMAKE_EXT_OBJ
kernel.commands += nvcc -c -o $$OUT_PWD/kernel$$QMAKE_EXT_OBJ $$PWD/kernel.cu
#kernel.commands += -Xcompiler /MD # Remove this line for most builds
QMAKE_EXTRA_TARGETS += kernel
PRE_TARGETDEPS += $$OUT_PWD/kernel$$QMAKE_EXT_OBJ