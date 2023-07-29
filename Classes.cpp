#include<iostream>
#include<cmath>

#include "Classes.h"
//#include "Point.h"

double toln = 0.01 / 50., toli = 0.1*toln; //toln should be~0.01*(h1,h2)

//Curvebase
Curvebase::Curvebase(const double pmin, const double pmax) : pmin_(pmin), pmax_(pmax)
{
	if (pmin > pmax) {
		rev = true;
	}
	else {
		rev = false;
	}
}

const double Curvebase::x(const double s) {
	if (s < 0 || s > 1) exit(EXIT_FAILURE);
	double p = newton((1 - s)*pmin_ + s * pmax_, toln, s); // linear initial guess
	return xp(p);
}

const double Curvebase::y(const double s) {
	if (s < 0 || s > 1) exit(EXIT_FAILURE);
	double p = newton((1 - s)*pmin_ + s * pmax_, toln, s);
	return yp(p);
}

// find p corresponding s
const double Curvebase::newton(const double p0_, const double toln, const double s) {
	double p0, p1;
	length = ASI(&Curvebase::func, pmin_, pmax_, toli); // arc length
	p1 = p0_; // to be substituted to p0
	p0 = p1 + toln + 1.; // to go in loop
	while (abs(p1 - p0) > toln) {
		p0 = p1;
		//p(i+1)=p(i)-F(p(i))/F'(p(i)); F(p)=l(p)-s*length
		p1 = p0 - (ASI(&Curvebase::func, pmin_, p0, toli) - s * length) / sqrt(dxp(p0)*dxp(p0) + dyp(p0)*dyp(p0));
	};
	if (p1 < pmin_ || p1 > pmax_) exit(EXIT_FAILURE);
	return p1;
}

const double Curvebase::ASI(const FunctionPointer f, const double& a, const double& b, const double& tol) const{
	double I1, I2, errest;
	if (f == nullptr) //null check
	{
		std::cout << "error" << std::endl;
		exit(EXIT_FAILURE);
	}
	I1 = I(f, a, b);
	I2 = I(f, a, (a + b) / 2.) + I(f, (a + b) / 2., b); //I(alpha,gamma)+I(gamma,beta)
	errest = abs(I1 - I2);
	if (errest < 15.*tol) return I2;
	return ASI(f, a, (a + b) / 2., tol / 2.) + ASI(f, (a + b) / 2., b, tol / 2.); //Recursion
}

const double Curvebase::func(const double& x) const {
	return sqrt(dxp(x)*dxp(x) + dyp(x)*dyp(x)); //
}

const double Curvebase::I(const FunctionPointer f, const double& a, const double& b) const{
	if (f == nullptr) //null check
	{
		std::cout << "error" << std::endl;
		exit(EXIT_FAILURE);
	}
	return (b - a) / 6.*((this->*f)(a) + 4.*(this->*f)((a + b) / 2.) + (this->*f)(b));
}

//Lines
Lines::Lines(const double pmin, const double pmax, const double val, const bool dir) : Curvebase(pmin, pmax), val_(val), dir_(dir) {
	length = abs(pmin - pmax);
}

const double Lines::x(const double s){ //overwrite base
	if (s < 0 || s > 1) exit(EXIT_FAILURE);
	if (!dir_) {
		if (!rev) {//normal orientatinon
			return xp(pmin_ + s * length);
		}
		else {//reverse orientatinon
			return xp(pmin_ - (1 - s) * length);
		}
	}
	else {
		return val_;
	}
}

const double Lines::y(const double s) { //overwrite base
	if (s < 0 || s > 1) exit(EXIT_FAILURE);
	if (dir_) {
		if (!rev) {//normal orientatinon
			return yp(pmin_ + s * length);
		}
		else {//reverse orientatinon
			return yp(pmin_ - (1 - s) * length);
		}
	}
	else {
		return val_;
	}
}

const double Lines::xp(const double p)const {
	if (!dir_) {
		return p;
	}
	else {
		return val_;
	}
}

const double Lines::yp(const double p)const {
	if (!dir_) {
		return val_;
	}
	else {
		return p;
	}
}

const double Lines::dxp(const double p)const {
	if (!dir_) {
		return 1.;
	}
	else{
		return 0.;
	}
}

const double Lines::dyp(const double p)const {
	if (!dir_) {
		return 0.;
	}
	else{
		return 1.;
	}
}

//Bcurve
const double Bcurve::xp(const double p) const{
	return p;
}

const double Bcurve::yp(const double p) const {
	if (p < -3) {
		return 0.5 / (1. + std::exp(-3.*(p + 6.)));
	}
	else {
		return 0.5 / (1. + std::exp(3.*p));
	}
}

const double Bcurve::dxp(const double p) const {
	return 1.;
}

const double Bcurve::dyp(const double p) const {
	if (p < -3) {
		return 1.5*std::exp(-3.*(p + 6.)) / (1. + std::exp(-3.*(p + 6.))) / (1. + std::exp(-3.*(p + 6.)));
	}
	else {
		return -1.5 *std::exp(3.*p) / (1. + std::exp(3.*p)) / (1. + std::exp(3.*p));
	}
}

//Domain
Domain::Domain(Curvebase& s0, Curvebase& s1, Curvebase& s2, Curvebase& s3) {
	sides[0] = &s0; //dynamic type changed
	sides[1] = &s1;
	sides[2] = &s2;
	sides[3] = &s3;
	if (!check_consistency())
		sides[0] = sides[1] = sides[2] = sides[3] = nullptr;
	m_ = n_ = 0; //no grid
	x_ = y_ = nullptr;
}

void Domain::generate_grid(const int m, const int n) {
	if (m < 1 || n < 1) exit(EXIT_FAILURE);
	if (m_ > 0) { // already exist grid
		delete[] x_;
		delete[] y_;
	}
	m_ = m; n_ = n;
	x_ = new double[(m_ + 1)*(n_ + 1)];
	y_ = new double[(m_ + 1)*(n_ + 1)];
	double h1 = 1. / m, h2 = 1. / n;
	for (int i = 0; i <= m_; i++) {
		for (int j = 0; j <= n_; j++) {
			int ind = i + j * (m_ + 1); //1D array
			x_[ind] = phi1(i*h1)*sides[3]->x(j*h2) + phi2(i*h1)*sides[1]->x(j*h2)
					+ phi1(j*h2)*sides[0]->x(i*h1) + phi2(j*h2)*sides[2]->x(i*h1)
					- phi1(i*h1)*phi1(j*h2)*sides[0]->x(0)
					- phi2(i*h1)*phi1(j*h2)*sides[1]->x(0)
					- phi1(i*h1)*phi2(j*h2)*sides[3]->x(1)
					- phi2(i*h1)*phi2(j*h2)*sides[2]->x(1);
			y_[ind] = phi1(i*h1)*sides[3]->y(j*h2) + phi2(i*h1)*sides[1]->y(j*h2)
					+ phi1(j*h2)*sides[0]->y(i*h1) + phi2(j*h2)*sides[2]->y(i*h1)
					- phi1(i*h1)*phi1(j*h2)*sides[0]->y(0)
					- phi2(i*h1)*phi1(j*h2)*sides[1]->y(0)
					- phi1(i*h1)*phi2(j*h2)*sides[3]->y(1)
					- phi2(i*h1)*phi2(j*h2)*sides[2]->y(1);
		}
	}
}

void Domain::write_grid()const {
	FILE *fp;
	fopen_s(&fp, "grid.bin", "wb");
	fwrite(x_, sizeof(double), (m_ + 1)*(n_ + 1), fp);
	fwrite(y_, sizeof(double), (m_ + 1)*(n_ + 1), fp);
	fclose(fp);
}

Domain::~Domain() {
	if (m_ > 0) {
		delete[] x_;
		delete[] y_;
	}
}

bool Domain::check_consistency() const {
	double tol = 1.e-5; //permit small non-consistensy
	if (abs(sides[0]->x(1) - sides[1]->x(0)) < tol && abs(sides[0]->y(1) - sides[1]->y(0)) < tol &&
		abs(sides[1]->x(1) - sides[2]->x(1)) < tol && abs(sides[1]->y(1) - sides[2]->y(1)) < tol &&
		abs(sides[2]->x(0) - sides[3]->x(1)) < tol && abs(sides[2]->y(0) - sides[3]->y(1)) < tol &&
		abs(sides[3]->x(0) - sides[0]->x(0)) < tol && abs(sides[3]->y(0) - sides[0]->y(0)) < tol)
	{
		return true;
	}
	else {
		std::cout << "check_consistency false" << std::endl;
		return false;
	}
}

//const Point Domain::operator()(const int i, const int j) const{
//	if (i < 0 || i > m_ || j < 0 || j > n_) {
//		exit(EXIT_FAILURE);
//	}
//	int ind = i + j * (m_ + 1);
//	return Point(x_[ind],y_[ind]);
//}