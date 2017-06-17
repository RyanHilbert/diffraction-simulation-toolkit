#include"vector3d.h"
#include"main.h"

Vector3D::Vector3D(QColor color,unsigned char type,QNode*parent):QEntity(parent),flipped(type&FLIP){
	addComponent(transform);
	QPhongMaterial*material=new QPhongMaterial(this);
	material->setShareable(true);
	material->setAmbient(color);
	QEntity*tail=new QEntity(this);
	tail->addComponent(material);
	QCylinderMesh*tailMesh=new QCylinderMesh(tail);
	tailMesh->setLength(tailLength);
	tailMesh->setRadius(.01f);
	tail->addComponent(tailMesh);
	tail->addComponent(tailTransform);
	if(type&AXIS){
		tailMesh->setRadius(.009f);
		tailMesh->setLength(9999);
		tailTransform->setTranslation(QVector3D(0,tailMesh->length()/2,0));
	}else{
		QEntity*head=new QEntity(this);
		head->addComponent(material);
		QConeMesh*headMesh=new QConeMesh(head);
		headMesh->setLength(headLength);
		headMesh->setBottomRadius(tailMesh->radius()*8);
		head->addComponent(headMesh);
		head->addComponent(headTransform);
		setDistance();
	}
}
void Vector3D::setDirection(float x,float y,float z){
	quaternion q=direction_to_quaternion(vector{x,y,z});
	transform->setRotation(QQuaternion(q.scalar,q.x,q.y,q.z));
}
void Vector3D::setDistance(float distance){
	distance*=tailLength+headLength;
	if(flipped){
		tailTransform->setTranslation(QVector3D(0,tailLength/-2.-headLength-distance,0));
		headTransform->setTranslation(QVector3D(0,headLength/-2.-distance,0));
	}else{
		tailTransform->setTranslation(QVector3D(0,tailLength/2.+distance,0));
		headTransform->setTranslation(QVector3D(0,tailLength+headLength/2.+distance,0));
	}
}
void Vector3D::setLength(float length){
	transform->setScale3D(QVector3D(length,length/(tailLength+headLength),length));
}