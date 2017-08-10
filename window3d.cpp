#include"window3d.h"
#include"main.h"
QEntity*Window3D::plane3D(){
	QEntity*plane=new QEntity(root);
	QEntity*top=new QEntity(plane);
	QEntity*bottom=new QEntity(plane);
	QPlaneMesh*mesh=new QPlaneMesh();
	QPhongMaterial*material=new QPhongMaterial();
	material->setShareable(true);
	material->setAmbient(Qt::lightGray);
	mesh->setShareable(true);
	botransform->setRotationX(180);
	plane->addComponent(detransform);
	top->addComponent(mesh);
	top->addComponent(material);
	top->addComponent(toptransform);
	bottom->addComponent(mesh);
	bottom->addComponent(material);
	bottom->addComponent(botransform);
	return plane;
}
QEntity*Window3D::cuboid3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QCuboidMesh*mesh=new QCuboidMesh(crystal);
	crystal->addComponent(mesh);
	return crystal;
}
QEntity*Window3D::sphere3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QSphereMesh*mesh=new QSphereMesh(crystal);
	crystal->addComponent(mesh);
	return crystal;
}
QEntity*Window3D::hemisphere3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QSphereMesh*mesh=new QSphereMesh(crystal);
	crystal->addComponent(mesh);
	return crystal;
}
QEntity*Window3D::cylinder3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QCylinderMesh*mesh=new QCylinderMesh(crystal);
	crystal->addComponent(mesh);
	return crystal;
}
QEntity*Window3D::cone3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QConeMesh*mesh=new QConeMesh(crystal);
	crystal->addComponent(mesh);
	mesh->setLength(2);
	return crystal;
}
QEntity*Window3D::bicone3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(transform);
	QEntity*top=new QEntity(crystal);
	QEntity*bottom=new QEntity(crystal);
	top->addComponent(material);
	bottom->addComponent(material);
	QConeMesh*mesh=new QConeMesh();
	mesh->setShareable(true);
	top->addComponent(mesh);
	bottom->addComponent(mesh);
	mesh->setLength(2);
	auto*topTransform=new Qt3DCore::QTransform(top);
	auto*bottomTransform=new Qt3DCore::QTransform(bottom);
	top->addComponent(topTransform);
	bottom->addComponent(bottomTransform);
	bottomTransform->setRotationX(180);
	topTransform->setTranslation(QVector3D(0,1,0));
	bottomTransform->setTranslation(QVector3D(0,-1,0));
	connect(crystal,&QEntity::enabledChanged,top,&QEntity::setEnabled);
	connect(crystal,&QEntity::enabledChanged,bottom,&QEntity::setEnabled);
	return crystal;
}
QEntity*Window3D::prism3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QCylinderMesh*mesh=new QCylinderMesh(crystal);
	crystal->addComponent(mesh);
	connect(this,&Window3D::updated,mesh,&QCylinderMesh::setSlices);
	return crystal;
}
QEntity*Window3D::pyramid3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(material);
	crystal->addComponent(transform);
	QConeMesh*mesh=new QConeMesh(crystal);
	crystal->addComponent(mesh);
	mesh->setLength(2);
	connect(this,&Window3D::updated,mesh,&QConeMesh::setSlices);
	return crystal;
}
QEntity*Window3D::bipyramid3D(){
	QEntity*crystal=new QEntity(root);
	crystal->addComponent(transform);
	QEntity*top=new QEntity(crystal);
	QEntity*bottom=new QEntity(crystal);
	top->addComponent(material);
	bottom->addComponent(material);
	QConeMesh*mesh=new QConeMesh();
	mesh->setShareable(true);
	top->addComponent(mesh);
	bottom->addComponent(mesh);
	mesh->setLength(2);
	auto*topTransform=new Qt3DCore::QTransform(top);
	auto*bottomTransform=new Qt3DCore::QTransform(bottom);
	top->addComponent(topTransform);
	bottom->addComponent(bottomTransform);
	bottomTransform->setRotationX(180);
	topTransform->setTranslation(QVector3D(0,1,0));
	bottomTransform->setTranslation(QVector3D(0,-1,0));
	connect(crystal,&QEntity::enabledChanged,top,&QEntity::setEnabled);
	connect(crystal,&QEntity::enabledChanged,bottom,&QEntity::setEnabled);
	connect(this,&Window3D::updated,mesh,&QConeMesh::setSlices);
	return crystal;
}
Window3D::Window3D(){
	defaultFrameGraph()->setClearColor(Qt::black);
	setTitle("3D");
	setRootEntity(root);

	QRenderSettings*render=new QRenderSettings(root);
	render->setRenderPolicy(QRenderSettings::OnDemand);
	render->setActiveFrameGraph(activeFrameGraph());
	root->addComponent(render);

	detransform->setShareable(true);
	transform->setShareable(true);
	transform->setRotationX(90);
	material->setShareable(true);
	material->setAlpha(.9f);
	material->setAmbient(Qt::darkGray);

	QCamera*cam=camera();
	cam->setPosition(QVector3D(4,4,4));
	cam->setUpVector(QVector3D(0,0,1));
	cam->setViewCenter(QVector3D(0,0,0));

	crystals[CUBOID]=cuboid3D();
	crystals[SPHERE]=sphere3D();
	crystals[HEMISPHERE]=hemisphere3D();
	crystals[CYLINDER]=cylinder3D();
	crystals[CONE]=cone3D();
	crystals[BICONE]=bicone3D();
	crystals[PRISM]=prism3D();
	crystals[PYRAMID]=pyramid3D();
	crystals[BIPYRAMID]=bipyramid3D();

	vectors[0]=xaxis;
	vectors[1]=yaxis;
	vectors[2]=zaxis;
	vectors[3]=rotationAxis;
	vectors[4]=hVector;
	vectors[5]=incidentBeam;

	incidentBeam->setDistance(.4f);
	xaxis->setDirection(1,0,0);
	yaxis->setDirection(0,1,0);
	zaxis->setDirection(0,0,1);

	connect(this,&QWindow::activeChanged,this,&Window3D::update);
	connect(this,&QWindow::windowStateChanged,this,&Window3D::update);
}
void Window3D::update(){
	for(size_t i=0;i<sizeof(crystals)/sizeof(*crystals);++i)crystals[i]->setEnabled(false);
	crystals[(size_t)g_shape]->setEnabled(true);

	vector raxis=rot_axis();
	rotationAxis->setDirection(raxis.x,raxis.y,raxis.z);
	hVector->setDirection(g_unit_h.x,g_unit_h.y,g_unit_h.z);
	incidentBeam->setDirection(g_unit_s0.x,g_unit_s0.y,g_unit_s0.z);
	quaternion q=direction_to_quaternion(det_direction());
	detransform->setRotation(QQuaternion(q.scalar,q.x,q.y,q.z));

	float max_dim=0;
	switch(g_shape){
	case CUBOID:
		max_dim=2*std::max(std::max(g_half_width,g_half_height),g_half_length);
		transform->setScale3D(QVector3D(2*g_half_length,2*g_half_height,2*g_half_width));
		break;
	case SPHERE:case HEMISPHERE:
		max_dim=2*g_radius;
		transform->setScale3D(QVector3D(g_radius,g_radius,g_radius));
		break;
	case CONE:case PYRAMID:
		max_dim=2*std::max(g_radius,g_half_height);
		transform->setScale3D(QVector3D(g_radius,g_half_height,g_radius));
		break;
	case BICONE:case BIPYRAMID:
		max_dim=2*std::max(g_radius,g_half_height);
		transform->setScale3D(QVector3D(g_radius,g_half_height/2,g_radius));
		break;
	default:
		max_dim=2*std::max(g_radius,g_half_height);
		transform->setScale3D(QVector3D(g_radius,2*g_half_height,g_radius));
	}
	for(size_t i=0;i<sizeof(vectors)/sizeof(*vectors);++i)vectors[i]->setLength(max_dim*2);
	detransform->setScale3D(QVector3D(max_dim*2,1,max_dim*2));
	botransform->setTranslation(QVector3D(0,max_dim*2.5,0));
	toptransform->setTranslation(QVector3D(0,max_dim*2.5,0));

	const int edges=g_edges;
	emit updated(edges>16?16:edges);
}
void Window3D::mouseMoveEvent(QMouseEvent*e){
	if(!lastMousePosition.isNull()){
		QPoint pos=e->pos();
		QCamera*cam=camera();
		cam->panAboutViewCenter(lastMousePosition.x()-pos.x(),QVector3D(0,0,1));
		cam->tiltAboutViewCenter(pos.y()-lastMousePosition.y());
		lastMousePosition=pos;
	}
}
void Window3D::mousePressEvent(QMouseEvent*e){
	lastMousePosition=e->pos();
}
void Window3D::mouseReleaseEvent(QMouseEvent*e){
	e=e;
	lastMousePosition=QPoint(0,0);
}
void Window3D::wheelEvent(QWheelEvent*e){
	QCamera*cam=camera();
	int angle=e->angleDelta().y();
	if(angle<0||cam->position().length()>1.1)cam->translate(QVector3D(0,0,angle?(angle<0?-1:1):0),QCamera::DontTranslateViewCenter);
}