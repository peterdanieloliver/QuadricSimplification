#pragma once
//#include <QtOpenGL>
#include <QVector2D>
#include <QVector3D>
#include <QQuaternion>
#include <QtMath>

class Trackball
{
private:
	
	float TRACKBALLSIZE = 0.8;

	// project x,y pair onto a sphere of radius r OR a hyperbolic sheet
	// if far away from the center of the sphere. Returs x,y,z as a vector
	float projectToSphere(float radius, QVector2D pos)
	{
		float d, t, z;

		d = pos.length();
		if (d < radius * sqrt(0.5))
		{
			// inside sphere
			z = sqrt(radius * radius - d * d);
		}
		else
		{
			// on hyperbola
			t = radius / sqrt(2.0);
			z = t * t / d;
		}

		return z;
	}

public:
	
	Trackball(float size)
	{
		TRACKBALLSIZE = size;
	}

	QQuaternion getRotationQuat(QVector2D pold, QVector2D pnew)
	{
		QVector3D axis;							// axis of rotation
		float phi;								// rotation amount
		QVector3D proj_new, proj_old, diff;		// projections
		float t;								// something
		QQuaternion rotquat = QQuaternion();	// rotation quaternion

		// check for repeated values
		if (pnew.x() == pold.x() && pnew.y() == pold.y())
		{
			return rotquat;
		}

		// find z-coordinate projections for pnew and pold
		proj_new = QVector3D(pnew.x(), pnew.y(), projectToSphere(TRACKBALLSIZE, pnew));
		proj_old = QVector3D(pold.x(), pold.y(), projectToSphere(TRACKBALLSIZE, pold));

		axis = QVector3D::crossProduct(proj_new, proj_old);

		// determine rotation amount around axis
		diff = proj_new - proj_old;
		t = diff.length() / (2.0 * TRACKBALLSIZE);
		// limit rotation amount
		if (t > 1.0) { t = 1.0; }
		else if (t < -1.0) { t = -1.0; }
		// compute rotation
		phi = 2.0 * asin(t);
		phi = qRadiansToDegrees(phi);

		// create quaternion from angle and rotation
		rotquat = QQuaternion::fromAxisAndAngle(axis, phi);
		return rotquat;
	}
};