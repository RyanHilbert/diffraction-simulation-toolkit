#ifndef VECTOR3D_H
#define VECTOR3D_H

#include<Qt3DCore>
#include<Qt3DRender>
#include<Qt3DExtras>

using namespace Qt3DCore;
using namespace Qt3DRender;
using namespace Qt3DExtras;

enum{AXIS=1,FLIP=2};

class Vector3D:public QEntity{
	Q_OBJECT
	static const char tailLength=8;
	static const char headLength=2;
	const bool flipped;
	Qt3DCore::QTransform
	*transform=new Qt3DCore::QTransform(this),//handles rotations
	*tailTransform=new Qt3DCore::QTransform(),//handles head translation
	*headTransform=new Qt3DCore::QTransform();//handles tail translation
public:
	Vector3D(QColor color,unsigned char type=0,QNode*parent=0);
public slots:
	void setDirection(float x=0,float y=0,float z=0);
	void setDistance(float distance=0);
	void setLength(float length=1);
};
#endif