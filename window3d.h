#ifndef WINDOW3D_H
#define WINDOW3D_H

#include"vector3d.h"

#include<Qt3DCore>
#include<Qt3DRender>
#include<Qt3DExtras>

using namespace Qt3DCore;
using namespace Qt3DRender;
using namespace Qt3DExtras;

class Window3D:public Qt3DWindow{
	Q_OBJECT

	QPoint lastMousePosition=QPoint(0,0);

	QEntity*root=new QEntity();
	QEntity*crystals[9];
	Vector3D*vectors[6];
	QPhongAlphaMaterial*material=new QPhongAlphaMaterial();
	Qt3DCore::QTransform*transform=new Qt3DCore::QTransform();

	Qt3DCore::QTransform*detransform=new Qt3DCore::QTransform();
	Qt3DCore::QTransform*botransform=new Qt3DCore::QTransform();
	Qt3DCore::QTransform*toptransform=new Qt3DCore::QTransform();
	QEntity*det=plane3D();

	QEntity*plane3D();
	QEntity*cuboid3D();
	QEntity*sphere3D();
	QEntity*hemisphere3D();
	QEntity*cylinder3D();
	QEntity*cone3D();
	QEntity*bicone3D();
	QEntity*prism3D();
	QEntity*pyramid3D();
	QEntity*bipyramid3D();

	Vector3D*xaxis=new Vector3D(Qt::red,AXIS,root),*yaxis=new Vector3D(Qt::green,AXIS,root),*zaxis=new Vector3D(Qt::blue,AXIS,root);
	Vector3D*rotationAxis=new Vector3D(Qt::yellow,0,root),*hVector=new Vector3D(Qt::magenta,0,root),*incidentBeam=new Vector3D(Qt::cyan,FLIP,root);
protected:
	void mouseMoveEvent(QMouseEvent*)override;
	void mousePressEvent(QMouseEvent*)override;
	void mouseReleaseEvent(QMouseEvent*)override;
	void wheelEvent(QWheelEvent*)override;
public:
	Window3D();
signals:
	void updated(int);
public slots:
	void update();
};
#endif