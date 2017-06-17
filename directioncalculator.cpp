#include"directioncalculator.h"

DirectionCalculator::DirectionCalculator(QWidget*parent):Dialog(parent,Qt::MSWindowsFixedSizeDialogHint){
	vertical->setDirection(QBoxLayout::BottomToTop);
	QPushButton*calc=new QPushButton("Calculate");
	connect(calc,&QPushButton::clicked,this,&DirectionCalculator::calculate);
	connect(this,&QDialog::accepted,this,&DirectionCalculator::onAccept);
	userx->setValidator(new QDoubleValidator(userx));
	usery->setValidator(new QDoubleValidator(usery));
	userz->setValidator(new QDoubleValidator(userz));
	resx->setReadOnly(true);
	resy->setReadOnly(true);
	resz->setReadOnly(true);
	braggBox->setReadOnly(true);
	QLabel*label=new QLabel("Define the diffraction plane by specifying a vector non-parallel to h");
	QHBoxLayout*upper=new QHBoxLayout,*middle=new QHBoxLayout,*lower=new QHBoxLayout;
	vertical->addLayout(lower);
	vertical->addLayout(middle);
	vertical->addSpacing(spacingBox->height());
	vertical->addLayout(upper);
	vertical->addWidget(label);
	upper->addWidget(userx);
	upper->addWidget(new QLabel("x"));
	upper->addWidget(usery);
	upper->addWidget(new QLabel("y"));
	upper->addWidget(userz);
	upper->addWidget(new QLabel("z"));
	middle->addWidget(new QLabel("Bragg angle(°) ="));
	middle->addWidget(braggBox);
	middle->addWidget(calc);
	lower->addWidget(new QLabel("Incident beam direction ="));
	lower->addWidget(resx);
	lower->addWidget(new QLabel("x"));
	lower->addWidget(resy);
	lower->addWidget(new QLabel("y"));
	lower->addWidget(resz);
	lower->addWidget(new QLabel("z"));

	QHBoxLayout*higher=new QHBoxLayout,*highest=new QHBoxLayout;
	highest->addWidget(new QLabel("d spacing(Å)"));
	highest->addWidget(spacingBox);
	highest->addWidget(new QLabel("Energy(keV)"));
	highest->addWidget(energyBox);
	higher->addWidget(new QLabel("h vector"));
	higher->addWidget(hx);
	higher->addWidget(new QLabel("x"));
	higher->addWidget(hy);
	higher->addWidget(new QLabel("y"));
	higher->addWidget(hz);
	higher->addWidget(new QLabel("z"));
	vertical->addLayout(higher);
	vertical->addLayout(highest);

	setMaximumSize(minimumSize());
}
void DirectionCalculator::calculate(){
	vector v={userx->text().toDouble(),usery->text().toDouble(),userz->text().toDouble()};
	double spacing=spacingBox->text().toDouble();
	double energy=energyBox->text().toDouble();
	vector h={hx->text().toDouble(),hy->text().toDouble(),hz->text().toDouble()};
	//user input captured in vector v, defined above

	//calculation code goes here
	double wave_number=2.0*PI/(12.398/energy), x1, x2, c, bragg_angle;
	vector k0;

	//calculate the reciprocal lattice vector
	h = vec_num(2.0*PI/(vec_abs(h)*spacing),h);

	//calculate Bragg angle
	if ((12.398/energy)/(2.0*spacing)>1.0){
		//no solution
		resx->setText(QString::number(NAN));
		resy->setText(QString::number(NAN));
		resz->setText(QString::number(NAN));
		braggBox->setText(QString::number(NAN));
		return;
	}
	else {
		//calculated Bragg angle in degree
		bragg_angle=asin((12.398/energy)/(2.0*spacing))*180.0/PI;
	}

	//calculate vector k0 at Bragg condition
	if (vec_dot(v, h) == 0.0){
		//v and h are perpendicular to each other
		x2=-0.5;
		if (pow(wave_number,2.0)-vec_dot(h, h)/4.0<0.0){
			//no solution
			resx->setText(QString::number(NAN));
			resy->setText(QString::number(NAN));
			resz->setText(QString::number(NAN));
			braggBox->setText(QString::number(NAN));
			return;
		}
		else {
			x1=sqrt((pow(wave_number,2.0)-vec_dot(h, h)/4.0)/vec_dot(v,v));
			k0=vec_plus(vec_num(x1,v),vec_num(x2,h));
		}
	}
	else {
		//v and h are not perpendicular
		c=(pow(vec_dot(h,h),2.0)*vec_dot(v,v)/(4.0*pow(vec_dot(v,h),2.0))-pow(wave_number,2.0))/(pow(vec_dot(h,h),2.0)*vec_dot(v,v)/pow(vec_dot(v,h),2.0) - vec_dot(h,h));
		if (1.0-4.0*c<0.0){
			//no solution
			resx->setText(QString::number(NAN));
			resy->setText(QString::number(NAN));
			resz->setText(QString::number(NAN));
			braggBox->setText(QString::number(NAN));
			return;
		}
		else {
			//force incident beam direction from left to right (non-negative x component of k0)
			if (-(1.0+sqrt(1.0-4.0*c))/2.0>0.0)
				x2=-(1.0+sqrt(1.0-4.0*c))/2.0;
			else
				x2=-(1.0-sqrt(1.0-4.0*c))/2.0;
			x1=-(0.5+x2)*vec_dot(h,h)/vec_dot(v,h);
			k0=vec_plus(vec_num(x1,v),vec_num(x2,h));
		}
	}
	//convert it to a unit vector
	k0=vec_num(1.0/vec_abs(k0),k0);

	//replace each NAN with the corresponding calculated vector component
	resx->setText(QString::number(k0.x));
	resy->setText(QString::number(k0.y));
	resz->setText(QString::number(k0.z));
	braggBox->setText(QString::number(bragg_angle));
}
vector DirectionCalculator::result(){
	return{resx->text().toDouble(),resy->text().toDouble(),resz->text().toDouble()};
}
void DirectionCalculator::setParams(DoubleParam*spacing,DoubleParam*energy,VectorParam*h,VectorParam*beam){
	spacingParam=spacing;
	energyParam=energy;
	hParam=h;
	beamParam=beam;
}
void DirectionCalculator::setVisible(bool visible){
	if(visible){
		spacingBox->setText(spacingParam->get());
		energyBox->setText(energyParam->get());
		vector h=hParam->get();
		hx->setText(QString::number(h.x));
		hy->setText(QString::number(h.y));
		hz->setText(QString::number(h.z));
	}
	Dialog::setVisible(visible);
}
void DirectionCalculator::onAccept(){
	spacingParam->set(spacingBox->text());
	energyParam->set(energyBox->text());
	hParam->set({hx->text().toDouble(),hy->text().toDouble(),hz->text().toDouble()});
	beamParam->set({resx->text().toDouble(),resy->text().toDouble(),resz->text().toDouble()});
}