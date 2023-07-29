#ifndef CLASSES_H
#define CLASSES_H

#include<stdlib.h>
#include<cmath>

//#include "Point.h"

class Curvebase {
public:
	Curvebase(const double pmin = 0.0, const double pmax = 1.0); //constructor
	virtual const double x(const double s); //arc length parametrization 0~1
	virtual const double y(const double s); //arc length parametrization 0~1
	//virtual ~Curvebase();
protected:
	// parameters
	const double pmin_, pmax_; // minimam & maximam value for p, a & b, user provided
	bool rev; // orientation of the curve, reverse
	double length; //whole arc length
	// functions
	virtual const double xp(const double p) const = 0; //pure virtual function
	virtual const double yp(const double p) const = 0;
	virtual const double dxp(const double p) const = 0;
	virtual const double dyp(const double p) const = 0;
	const double newton(const double a, const double toln, const double s);
	// integration functions from Project 1
	typedef const double(Curvebase::*FunctionPointer)(const double&)const; //define type of "FunctionPointer" (const double)(Curvebase::*)(const double) const
	const double ASI(const FunctionPointer f, const double& pmin, const double& p, const double& toli)const; //Adaptive Simpson Integration
	const double func(const double& x)const; //Integrated function
	const double I(const FunctionPointer f, const double& a, const double& b)const; //Simpson rule
};

class Lines : public Curvebase {
public:
	//constructor
	Lines(const double pmin, const double pmax, const double val, const bool dir);
	const double x(const double s); //overwrite base class
	const double y(const double s); //overwrite base class
private:
	//parameters
	const bool dir_; // 0 for y=val_, 1 for x=val_
	const double val_;
	//functions(overwrite)
	const double xp(const double p)const;
	const double yp(const double p)const;
	const double dxp(const double p)const;
	const double dyp(const double p)const;
};

//Bottom curve
class Bcurve : public Curvebase {
public:
	//constructor
	Bcurve(const double pmin, const double pmax) : Curvebase(pmin, pmax){ 
		if (rev) std::cout << "input error at construction of Bcurve" << std::endl;
	}
private:
	//functions(overwrite)
	const double xp(const double p)const;
	const double yp(const double p)const;
	const double dxp(const double p)const;
	const double dyp(const double p)const;
};


class Domain {
public:
	//constructor
	Domain(Curvebase& s0, Curvebase& s1, Curvebase& s2, Curvebase& s3);
	void generate_grid(const int m, const int n);
	void write_grid()const;
	~Domain();
	//const Point operator()(const int i, const int j) const;
protected:
	Curvebase *sides[4]; //pointers for only Curvebase (static bind)
	int m_ = 0, n_ = 0;
	double *x_, *y_; //grid data
	const double phi1(const double x) const { return 1-x; }
	const double phi2(const double x) const { return x; }
	bool check_consistency()const; //check corner
};

#endif // !CLASSES_H