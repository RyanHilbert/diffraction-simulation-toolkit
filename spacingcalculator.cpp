#include"spacingcalculator.h"
#include"main.h"

SpacingCalculator::SpacingCalculator(QWidget*parent):Dialog(parent,Qt::MSWindowsFixedSizeDialogHint){
	QPushButton*calc=new QPushButton("Calculate");
	connect(calc,&QPushButton::clicked,this,&SpacingCalculator::calculate);
	spacing->setReadOnly(true);
	edita->setValidator(new QDoubleValidator(edita));
	editb->setValidator(new QDoubleValidator(editb));
	editc->setValidator(new QDoubleValidator(editc));
	editA->setValidator(new QDoubleValidator(editA));
	editB->setValidator(new QDoubleValidator(editB));
	editC->setValidator(new QDoubleValidator(editC));
	QVBoxLayout*vertical=new QVBoxLayout(frame);
	QHBoxLayout*upper=new QHBoxLayout;
	QHBoxLayout*lower=new QHBoxLayout;
	vertical->addLayout(upper);
	vertical->addLayout(lower);
	lower->addWidget(calc);
	lower->addWidget(new QLabel("d (Å) ="));
	lower->addWidget(spacing);
	QFormLayout*form1=new QFormLayout,*form2=new QFormLayout,*form3=new QFormLayout;
	upper->addLayout(form1);
	upper->addLayout(form2);
	upper->addLayout(form3);
	form1->addRow("a (Å)",edita);
	form2->addRow("b (Å)",editb);
	form3->addRow("c (Å)",editc);
	form1->addRow("α (°)",editA);
	form2->addRow("β (°)",editB);
	form3->addRow("γ (°)",editC);
	form1->addRow("h",edith);
	form2->addRow("k",editk);
	form3->addRow("l",editl);
	setMaximumSize(minimumSize());
}
void SpacingCalculator::calculate(){
	double  a = edita->text().toDouble(),
			b = editb->text().toDouble(),
			c = editc->text().toDouble(),
			alpha = editA->text().toDouble(),
			beta = editB->text().toDouble(),
			gamma = editC->text().toDouble();
	int     h = edith->value(),
			k = editk->value(),
			l = editl->value();
	//above defines ints h, k, l and doubles a, b, c, alpha, beta, gamma

	vector a1, a2, a3, b1, b2, b3, diff_h;
	double c1, c2, tx, ty, vol, result;

	//convert degree to radian
	alpha = alpha*pi/180.0;
	beta = beta*pi/180.0;
	gamma = gamma*pi/180.0;

	//assume vector a1 along x axis, and a2 in xy plane
	a1 = {a, 0, 0};
	a2 = {b*cos(gamma), b*sin(gamma), 0.0};
	c1 = cos(alpha)-cos(beta)*cos(gamma);
	c2 = pow(sin(gamma)*sin(beta),2.0)-pow(c1,2.0);

	if (c2 > 0.0){
		//calculate a vector that makes angles alpha and beta with a2 and a1, respectively
		ty = c1/sqrt(c2);
		tx = cos(beta)*sqrt(1.0+ty*ty)/sin(beta);
	}
	else {
		//no solution
		spacing->setText(QString::number(NAN));
		return;//no result; set spacing box to NAN and return
	}
	//assume as has a nonzero z component, which is always true
	a3 = {tx, ty, 1.0};

	//make the length of vector a3 equal to c
	a3 = vec_num(c/vec_abs(a3), a3);

	//calculate volume of the unit cell
	vol = vec_dot(vec_cross(a2, a3), a1);

	//calculate the base reciprocal lattice vectors b1, b2 and b3
	b1 = vec_num(1.0/vol, vec_cross(a2, a3));
	b2 = vec_num(1.0/vol, vec_cross(a3, a1));
	b3 = vec_num(1.0/vol, vec_cross(a1, a2));

	//calculate the reciprocal lattice vector corresponding to the hkl reflection
	diff_h = vec_plus(vec_num((double) h, b1), vec_num((double) k, b2));
	diff_h = vec_plus(diff_h, vec_num((double) l, b3));

	//calculate the corresponding d-spacing
	result = 1.0/vec_abs(diff_h);

	spacing->setText(QString::number(result));
}
QString SpacingCalculator::result(){
	return spacing->text();
}
void SpacingCalculator::hdf5read(){
	hid_t group=h5group_("Calculator",h5group("Crystal",hdf5file));
	editA->setText(QString::number(h5double("α (°)",group)));
	editB->setText(QString::number(h5double("β (°)",group)));
	editC->setText(QString::number(h5double("γ (°)",group)));
	edita->setText(QString::number(h5double("a (Å)",group)));
	editb->setText(QString::number(h5double("b (Å)",group)));
	editc->setText(QString::number(h5double("c (Å)",group)));
	edith->setValue(h5int("h",group));
	editk->setValue(h5int("k",group));
	editl->setValue(h5int_("l",group));
}
void SpacingCalculator::hdf5write(){
	hid_t group=h5group_("Calculator",h5group("Crystal",hdf5file));
	h5set_double(editA->text().toDouble(),"α (°)",group);
	h5set_double(editB->text().toDouble(),"β (°)",group);
	h5set_double(editC->text().toDouble(),"γ (°)",group);
	h5set_double(edita->text().toDouble(),"a (Å)",group);
	h5set_double(editb->text().toDouble(),"b (Å)",group);
	h5set_double(editc->text().toDouble(),"c (Å)",group);
	h5set_int(edith->value(),"h",group);
	h5set_int(editk->value(),"k",group);
	h5set_int_(editl->value(),"l",group);
}