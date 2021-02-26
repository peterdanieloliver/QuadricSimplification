#pragma once

class icVector2;
class icVector3;
class icVector4;

extern "C"
{
#include <math.h>
#include <stdlib.h>
}

// Start class icVector2
class icVector2{
public:
  inline icVector2();
  inline icVector2(double d);
  inline icVector2(double d0,double d1);

  inline icVector2(const icVector2& a);
  inline icVector2(const double*    a);

  inline icVector2& set(double d);
  inline icVector2& set(double d0, double d1);

  inline icVector2& set(const icVector2& a);
  inline icVector2& set(const double*    a);

  inline icVector2& operator=(double d);
  inline icVector2& operator=(const icVector2& a);
  inline icVector2& operator=(const double*    a);

  inline int operator==(const icVector2& a) const;
  inline int operator!=(const icVector2& a) const;

  inline int operator==(double d) const;
  inline int operator!=(double d) const;

  inline icVector2& operator+=(double d);
  inline icVector2& operator-=(double d);
  inline icVector2& operator*=(double d);
  inline icVector2& operator/=(double d);

  inline icVector2& operator+=(const icVector2& a);
  inline icVector2& operator-=(const icVector2& a);
  inline icVector2& operator*=(const icVector2& a);
  inline icVector2& operator/=(const icVector2& a);
  inline double length(const icVector2& a);
  inline void   normalize(icVector2& a);

  inline double    dot(const icVector2& a,const icVector2& b);
  inline icVector2 cross(const icVector2& a);

public:
	union {
		double entry[2]; /// all the vector entries
		double array[2]; /// all the vector entries (an additional alias added later)
		struct {
			double x;
			double y;
		};
	};
};

inline icVector2 operator-(const icVector2& a);

inline icVector2 operator+(const icVector2& a, const icVector2& b);
inline icVector2 operator-(const icVector2& a, const icVector2& b);

inline icVector2 operator+(const icVector2& a, double b);
inline icVector2 operator-(const icVector2& a, double b);
inline icVector2 operator*(const icVector2& a, double b);

inline icVector2 operator+(double a, const icVector2& b);
inline icVector2 operator-(double a, const icVector2& b);
inline icVector2 operator*(double a, const icVector2& b);

inline double length(const icVector2& a);
inline void   normalize(icVector2& a);


inline double dot(const icVector2& a,const icVector2& b);
inline icVector2 cross(const icVector2& a);

inline icVector2::icVector2() {
  entry[0] = entry[1] = 0.0;
}
inline icVector2::icVector2(double d) {
  entry[0] = entry[1] = d;
}

inline icVector2::icVector2(double d0,double d1) {
  entry[0] = d0;
  entry[1] = d1;
}

inline icVector2::icVector2(const icVector2& a) {
  entry[0] = a.entry[0];
  entry[1] = a.entry[1];
}

inline icVector2::icVector2(const double* a) {
  entry[0] = a[0];
  entry[1] = a[1];
}

inline icVector2& icVector2::set(double d) {
  entry[0] = d;
  entry[1] = d;
  return (*this);
}

inline icVector2& icVector2::set(double d0, double d1) {
  entry[0] = d0;
  entry[1] = d1;
  return (*this);
}

inline icVector2& icVector2::set(const icVector2& a) {
  entry[0] = a.entry[0];
  entry[1] = a.entry[1];
  return (*this);
}

inline icVector2& icVector2::set(const double* a) {
  entry[0] = a[0];
  entry[1] = a[1];
  return (*this);
}

inline icVector2 operator-(const icVector2& a) {
  return icVector2(-a.entry[0],-a.entry[1]);
}

inline icVector2& icVector2::operator=(double d) {
  return set(d);
}

inline icVector2& icVector2::operator=(const icVector2& a) {
  return set(a);
}

inline icVector2& icVector2::operator=(const double* a) {
  return set(a);
}

//-------------------------------------------------------------------

inline int icVector2::operator==(const icVector2& a) const {
  return ((entry[0] == a.entry[0]) && (entry[1] == a.entry[1]));
}

inline int icVector2::operator!=(const icVector2& a) const {
  return ((entry[0] != a.entry[0]) || (entry[1] != a.entry[1]));
}

inline int icVector2::operator==(double d) const {
  return ((entry[0] == d) && (entry[1] == d));
}

inline int icVector2::operator!=(double d) const {
  return ((entry[0] != d) || (entry[1] != d));
}

//-------------------------------------------------------------------

inline icVector2& icVector2::operator+=(double d) {
  entry[0] += d;
  entry[1] += d;
  return (*this);
}

inline icVector2& icVector2::operator-=(double d) {
  entry[0] -= d;
  entry[1] -= d;
  return (*this);
}

inline icVector2& icVector2::operator*=(double d) {
  entry[0] *= d;
  entry[1] *= d;
  return (*this);
}

inline icVector2& icVector2::operator+=(const icVector2& a) {
  entry[0] += a.entry[0];
  entry[1] += a.entry[1];
  return (*this);
}

inline icVector2& icVector2::operator-=(const icVector2& a) {
  entry[0] -= a.entry[0];
  entry[1] -= a.entry[1];
  return (*this);
}

inline icVector2& icVector2::operator*=(const icVector2& a) {
  entry[0] *= a.entry[0];
  entry[1] *= a.entry[1];
  return (*this);
}

//-------------------------------------------------------------------

inline icVector2 operator+(const icVector2& a,const icVector2& b) {
  return icVector2(a.entry[0] + b.entry[0], a.entry[1] + b.entry[1]);
}

inline icVector2 operator-(const icVector2& a,const icVector2& b) {
  return icVector2(a.entry[0] - b.entry[0], a.entry[1] - b.entry[1]);
}

inline icVector2 operator+(const icVector2& a,double b){
  return icVector2(a.entry[0] + b, a.entry[1] + b);
}

inline icVector2 operator-(const icVector2& a,double b){
  return icVector2(a.entry[0] - b, a.entry[1] - b);
}

inline icVector2 operator*(const icVector2& a,double b){
  return icVector2(a.entry[0] * b, a.entry[1] * b);
}

inline icVector2 operator+(double a,const icVector2& b){
  return icVector2(a + b.entry[0], a + b.entry[1]);
}

inline icVector2 operator-(double a,const icVector2& b){
  return icVector2(a - b.entry[0], a - b.entry[1]);
}

inline icVector2 operator*(double a,const icVector2& b){
  return icVector2(a * b.entry[0], a * b.entry[1]);
}

inline double length(const icVector2& a) {
  return sqrt(a.entry[0] * a.entry[0] + a.entry[1] * a.entry[1]);
}

inline void normalize(icVector2& a) {
  register double m = length(a);
  if (m != 0) a *= (1/m);
}

inline double dot(const icVector2& a, const icVector2& b) {
  return (a.entry[0] * b.entry[0] + a.entry[1] * b.entry[1]);
}

inline icVector2 cross(const icVector2& a) {
  return icVector2(-a.entry[1], a.entry[0]);
}

//-------------------------------------------------------------------

// Start class icVector3
class icVector3{
public:
	inline icVector3();
	inline icVector3(double d);
	inline icVector3(double d0,double d1,double d2);

	inline icVector3(const icVector3& a);
	inline icVector3(const double*    a);

	inline icVector3& set(double d);
	inline icVector3& set(double d0, double d1,double d2);

	inline icVector3& set(const icVector3& a);
	inline icVector3& set(const double*    a);

	inline icVector3& operator=(double d);
	inline icVector3& operator=(const icVector3& a);
	inline icVector3& operator=(const double*    a);

	inline int operator==(const icVector3& a) const;
	inline int operator!=(const icVector3& a) const;

	inline int operator==(double d) const;
	inline int operator!=(double d) const;

	inline icVector3& operator+=(double d);
	inline icVector3& operator-=(double d);
	inline icVector3& operator*=(double d);
	inline icVector3& operator/=(double d);

	inline icVector3& operator+=(const icVector3& a);
	inline icVector3& operator-=(const icVector3& a);
	inline icVector3& operator*=(const icVector3& a);
	inline icVector3& operator/=(const icVector3& a);

public:
	union {
		double entry[3]; /// all the vector entries
		double array[3]; /// all the vector entries (an additional alias added later)
		struct {
			double x;
			double y;
			double z;
		};
	};
};

inline icVector3 operator-(const icVector3& a);

inline icVector3 operator+(const icVector3& a, const icVector3& b);
inline icVector3 operator-(const icVector3& a, const icVector3& b);

inline icVector3 operator+(const icVector3& a, double b);
inline icVector3 operator-(const icVector3& a, double b);
inline icVector3 operator*(const icVector3& a, double b);

inline icVector3 operator+(double a, const icVector3& b);
inline icVector3 operator-(double a, const icVector3& b);
inline icVector3 operator*(double a, const icVector3& b);

inline double length(const icVector3& a);
inline void normalize(icVector3& a);


inline double dot(const icVector3& a,const icVector3& b);
inline icVector3 cross(const icVector3& a, const icVector3& b);

inline icVector3::icVector3() {
  entry[0] = entry[1] = entry[2] = 0.0;
}
inline icVector3::icVector3(double d) {
  entry[0] = entry[1] = entry[2] = d;
}

inline icVector3::icVector3(double d0,double d1,double d2) {
  entry[0] = d0;
  entry[1] = d1;
  entry[2] = d2;
}

inline icVector3::icVector3(const icVector3& a) {
  entry[0] = a.entry[0];
  entry[1] = a.entry[1];
  entry[2] = a.entry[2];
}

inline icVector3::icVector3(const double* a) {
  entry[0] = a[0];
  entry[1] = a[1];
  entry[2] = a[2];
}

//-------------------------------------------------------------------

inline icVector3& icVector3::set(double d) {
  entry[0] = d;
  entry[1] = d;
  entry[2] = d;
  return (*this);
}

inline icVector3& icVector3::set(double d0, double d1, double d2) {
  entry[0] = d0;
  entry[1] = d1;
  entry[2] = d2;
  return (*this);
}

inline icVector3& icVector3::set(const icVector3& a) {
  entry[0] = a.entry[0];
  entry[1] = a.entry[1];
  entry[2] = a.entry[2];
  return (*this);
}

inline icVector3& icVector3::set(const double* a) {
  entry[0] = a[0];
  entry[1] = a[1];
  entry[2] = a[2];
  return (*this);
}

inline icVector3 operator-(const icVector3& a) {
  return icVector3(-a.entry[0],-a.entry[1],-a.entry[2]);
}

inline icVector3& icVector3::operator=(double d) {
  return set(d);
}

inline icVector3& icVector3::operator=(const icVector3& a) {
  return set(a);
}

inline icVector3& icVector3::operator=(const double* a) {
  return set(a);
}

//-------------------------------------------------------------------

inline int icVector3::operator==(const icVector3& a) const {
  return (
	  (entry[0] == a.entry[0]) &&
	  (entry[1] == a.entry[1]) &&
	  (entry[2] == a.entry[2]));
}

inline int icVector3::operator!=(const icVector3& a) const {
  return (
	  (entry[0] != a.entry[0]) ||
	  (entry[1] != a.entry[1]) ||
	  (entry[2] != a.entry[2]));
}

inline int icVector3::operator==(double d) const {
  return (
	  (entry[0] == d) &&
	  (entry[1] == d) &&
	  (entry[2] == d));
}

inline int icVector3::operator!=(double d) const {
  return (
	  (entry[0] != d) ||
	  (entry[1] != d) ||
	  (entry[2] != d));
}

//-------------------------------------------------------------------

inline icVector3& icVector3::operator+=(double d) {
  entry[0] += d;
  entry[1] += d;
  entry[2] += d;
  return (*this);
}

inline icVector3& icVector3::operator-=(double d) {
  entry[0] -= d;
  entry[1] -= d;
  entry[2] -= d;
  return (*this);
}

inline icVector3& icVector3::operator*=(double d) {
  entry[0] *= d;
  entry[1] *= d;
  entry[2] *= d;
  return (*this);
}

inline icVector3& icVector3::operator/=(double d) {
  entry[0] /= d;
  entry[1] /= d;
  entry[2] /= d;
  return (*this);
}

inline icVector3& icVector3::operator+=(const icVector3& a) {
  entry[0] += a.entry[0];
  entry[1] += a.entry[1];
  entry[2] += a.entry[2];
  return (*this);
}

inline icVector3& icVector3::operator-=(const icVector3& a) {
  entry[0] -= a.entry[0];
  entry[1] -= a.entry[1];
  entry[2] -= a.entry[2];
  return (*this);
}

inline icVector3& icVector3::operator*=(const icVector3& a) {
  entry[0] *= a.entry[0];
  entry[1] *= a.entry[1];
  entry[2] *= a.entry[2];
  return (*this);
}

inline icVector3& icVector3::operator/=(const icVector3& a) {
  entry[0] /= a.entry[0];
  entry[1] /= a.entry[1];
  entry[2] /= a.entry[2];
  return (*this);
}

//-------------------------------------------------------------------

inline icVector3 operator+(const icVector3& a,const icVector3& b) {
  return icVector3(a.entry[0] + b.entry[0], a.entry[1] + b.entry[1], 
a.entry[2] + b.entry[2]);
}

inline icVector3 operator-(const icVector3& a,const icVector3& b) {
  return icVector3(a.entry[0] - b.entry[0], a.entry[1] - b.entry[1], 
a.entry[2] - b.entry[2]);
}

inline icVector3 operator+(const icVector3& a,double b){
  return icVector3(a.entry[0] + b, a.entry[1] + b, a.entry[2] + b);
}

inline icVector3 operator-(const icVector3& a,double b){
  return icVector3(a.entry[0] - b, a.entry[1] - b, a.entry[2] - b);
}

inline icVector3 operator*(const icVector3& a,double b){
  return icVector3(a.entry[0] * b, a.entry[1] * b, a.entry[2] * b);
}

inline icVector3 operator+(double a,const icVector3& b){
  return icVector3(a + b.entry[0], a + b.entry[1], a + b.entry[2]);
}

inline icVector3 operator-(double a,const icVector3& b){
  return icVector3(a - b.entry[0], a - b.entry[1], a - b.entry[2]);
}

inline icVector3 operator*(double a,const icVector3& b){
  return icVector3(a * b.entry[0], a * b.entry[1], a * b.entry[2]);
}

inline double length(const icVector3& a) {
  return sqrt(a.entry[0] * a.entry[0] + a.entry[1] * a.entry[1] + a.entry[2] 
* a.entry[2]);
}

inline void normalize(icVector3& a) {
  register double m = length(a);
  if (m != 0) a *= (1/m);
}

inline double dot(const icVector3& a,const icVector3& b) {
  return (a.entry[0] * b.entry[0] + a.entry[1] * b.entry[1] + a.entry[2] * 
b.entry[2]);
}

inline icVector3 cross(const icVector3& a, const icVector3& b) {
  return icVector3(
	a.entry[1] * b.entry[2] - a.entry[2] * b.entry[1],
	a.entry[2] * b.entry[0] - a.entry[0] * b.entry[2],
	a.entry[0] * b.entry[1] - a.entry[1] * b.entry[0]);
}

//-------------------------------------------------------------------

// Start class icVector4
class icVector4
{
public:
	inline icVector4();
	inline icVector4(double d);
	inline icVector4(double d0, double d1, double d2, double d3);

	inline icVector4(const icVector4& a);
	inline icVector4(const double* a);

	inline icVector4& set(double d);
	inline icVector4& set(double d0, double d1, double d2, double d3);

	inline icVector4& set(const icVector4& a);
	inline icVector4& set(const double* a);

	inline icVector4& operator=(double d);
	inline icVector4& operator=(const icVector4& a);
	inline icVector4& operator=(const double* a);

	inline int operator==(const icVector4& a) const;
	inline int operator!=(const icVector4& a) const;

	inline int operator==(double d) const;
	inline int operator!=(double d) const;

	inline icVector4& operator+=(double d);
	inline icVector4& operator-=(double d);
	inline icVector4& operator*=(double d);
	inline icVector4& operator/=(double d);

	inline icVector4& operator+=(const icVector4& a);
	inline icVector4& operator-=(const icVector4& a);
	inline icVector4& operator*=(const icVector4& a);
	inline icVector4& operator/=(const icVector4& a);

public:
	union
	{
		double entry[4];
		double array[4];
		struct
		{
			double x;
			double y;
			double z;
		};
	};
};

inline icVector4 operator-(const icVector4& a);

inline icVector4 operator+(const icVector4& a, const icVector4& b);
inline icVector4 operator-(const icVector4& a, const icVector4& b);

inline icVector4 operator+(const icVector4& a, double b);
inline icVector4 operator-(const icVector4& a, double b);
inline icVector4 operator*(const icVector4& a, double b);

inline icVector4 operator+(double a, const icVector4& b);
inline icVector4 operator-(double a, const icVector4& b);
inline icVector4 operator*(double a, const icVector4& b);

inline double length(const icVector4& a);
inline void normalize(icVector4& a);
inline double dot(const icVector4& a, const icVector4& b);

//-----------------------------------------------------------------

inline icVector4::icVector4()
{
	entry[0] = entry[1] = entry[2] = entry[3] = 0.0;
}

inline icVector4::icVector4(double d)
{
	entry[0] = entry[1] = entry[2] = entry[3] = d;
}

inline icVector4::icVector4(const icVector4& a)
{
	entry[0] = a.entry[0];
	entry[1] = a.entry[1];
	entry[2] = a.entry[2];
	entry[3] = a.entry[3];
}

inline icVector4::icVector4(const double* a)
{
	entry[0] = a[0];
	entry[1] = a[1];
	entry[2] = a[2];
	entry[2] = a[3];
}

//-------------------------------------------------------------------

inline icVector4& icVector4::set(double d)
{
	entry[0] = d;
	entry[1] = d;
	entry[2] = d;
	entry[3] = d;
	return (*this);
}

inline icVector4& icVector4::set(double d0, double d1, double d2, double d3)
{
	entry[0] = d0;
	entry[1] = d1;
	entry[2] = d2;
	entry[3] = d3;
	return (*this);
}

inline icVector4& icVector4::set(const icVector4& a)
{
	entry[0] = a.entry[0];
	entry[1] = a.entry[1];
	entry[2] = a.entry[2];
	entry[3] = a.entry[3];
	return (*this);
}

inline icVector4& icVector4::set(const double* a)
{
	entry[0] = a[0];
	entry[1] = a[1];
	entry[2] = a[2];
	entry[3] = a[3];
	return (*this);
}

inline icVector4 operator-(const icVector4& a) 
{
	return icVector4(-a.entry[0], -a.entry[1], -a.entry[2], -a.entry[3]);
}

inline icVector4& icVector4::operator=(double d)
{
	return set(d);
}

inline icVector4& icVector4::operator=(const icVector4& a) 
{
	return set(a);
}

inline icVector4& icVector4::operator=(const double* a) 
{
	return set(a);
}

//-------------------------------------------------------------------

inline int icVector4::operator==(const icVector4& a) const 
{
	return 
	(
		(entry[0] == a.entry[0]) &&
		(entry[1] == a.entry[1]) &&
		(entry[2] == a.entry[2]) &&
		(entry[3] == a.entry[3])
	);
}

inline int icVector4::operator!=(const icVector4& a) const 
{
	return
	(
		(entry[0] != a.entry[0]) ||
		(entry[1] != a.entry[1]) ||
		(entry[2] != a.entry[2]) ||
		(entry[3] != a.entry[3])
	);
}

inline int icVector4::operator==(double d) const
{
	return
	(
		(entry[0] == d) &&
		(entry[1] == d) &&
		(entry[2] == d) &&
		(entry[3] == d)
	);
}

inline int icVector4::operator!=(double d) const
{
	return
	(
		(entry[0] != d) ||
		(entry[1] != d) ||
		(entry[2] != d) ||
		(entry[3] != d)
	);
}

//-------------------------------------------------------------------

inline icVector4& icVector4::operator+=(double d)
{
	entry[0] += d;
	entry[1] += d;
	entry[2] += d;
	entry[3] += d;
	return (*this);
}

inline icVector4& icVector4::operator-=(double d)
{
	entry[0] -= d;
	entry[1] -= d;
	entry[2] -= d;
	entry[3] -= d;
	return (*this);
}

inline icVector4& icVector4::operator*=(double d)
{
	entry[0] *= d;
	entry[1] *= d;
	entry[2] *= d;
	entry[3] *= d;
	return (*this);
}

inline icVector4& icVector4::operator/=(double d) 
{
	entry[0] /= d;
	entry[1] /= d;
	entry[2] /= d;
	entry[3] /= d;
	return (*this);
}

inline icVector4& icVector4::operator+=(const icVector4& a)
{
	entry[0] += a.entry[0];
	entry[1] += a.entry[1];
	entry[2] += a.entry[2];
	entry[3] += a.entry[3];
	return (*this);
}

inline icVector4& icVector4::operator-=(const icVector4& a)
{
	entry[0] -= a.entry[0];
	entry[1] -= a.entry[1];
	entry[2] -= a.entry[2];
	entry[3] -= a.entry[3];
	return (*this);
}

inline icVector4& icVector4::operator*=(const icVector4& a)
{
	entry[0] *= a.entry[0];
	entry[1] *= a.entry[1];
	entry[2] *= a.entry[2];
	entry[3] *= a.entry[3];
	return (*this);
}

inline icVector4& icVector4::operator/=(const icVector4& a)
{
	entry[0] /= a.entry[0];
	entry[1] /= a.entry[1];
	entry[2] /= a.entry[2];
	entry[3] /= a.entry[3];
	return (*this);
}

//-------------------------------------------------------------------

inline icVector4 operator+(const icVector4& a, const icVector4& b) 
{
	return icVector4(a.entry[0] + b.entry[0],
					 a.entry[1] + b.entry[1],
					 a.entry[2] + b.entry[2],
					 a.entry[3] + b.entry[3]);
}

inline icVector4 operator-(const icVector4& a, const icVector4& b)
{
	return icVector4(a.entry[0] - b.entry[0], 
					 a.entry[1] - b.entry[1],
					 a.entry[2] - b.entry[2],
					 a.entry[3]	- b.entry[3]);
}

inline icVector4 operator+(const icVector4& a, double b)
{
	return icVector4(a.entry[0] + b, a.entry[1] + b, a.entry[2] + b, a.entry[3] + b);
}

inline icVector4 operator-(const icVector4& a, double b)
{
	return icVector4(a.entry[0] - b, a.entry[1] - b, a.entry[2] - b, a.entry[3] + b);
}

inline icVector4 operator*(const icVector4& a, double b) 
{
	return icVector4(a.entry[0] * b, a.entry[1] * b, a.entry[2] * b, a.entry[3] * b);
}

inline icVector4 operator+(double a, const icVector4& b) {
	return icVector4(a + b.entry[0], a + b.entry[1], a + b.entry[2], a + b.entry[3]);
}

inline icVector4 operator-(double a, const icVector4& b)
{
	return icVector4(a - b.entry[0], a - b.entry[1], a - b.entry[2], a - b.entry[3]);
}

inline icVector4 operator*(double a, const icVector4& b)
{
	return icVector4(a * b.entry[0], a * b.entry[1], a * b.entry[2], a * b.entry[3]);
}

inline double length(const icVector4& a)
{
	return sqrt(a.entry[0] * a.entry[0] + 
				a.entry[1] * a.entry[1] + 
				a.entry[2] * a.entry[2] +
				a.entry[3] * a.entry[3]);
}

inline void normalize(icVector4& a)
{
	register double m = length(a);
	if (m != 0) { a *= (1 / m); }
}

inline double dot(const icVector4& a, const icVector4& b)
{
	return (a.entry[0] * b.entry[0] +
			a.entry[1] * b.entry[1] +
			a.entry[2] * b.entry[2] +
			a.entry[3] * b.entry[3]);
}