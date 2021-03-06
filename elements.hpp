#pragma once
#include "icVector.hpp"
#include "icMatrix.hpp"
#include <vector>

// forward declarations
class Face;
class Edge;
class Corner;
class PairContraction;
class PairCompare;

//////////////////////////////////////////////////////////
// Contraction Classes
//////////////////////////////////////////////////////////

class Vertex
{
	// fields
public:

	int index;
	double x, y, z;
	double texcoord[2];

	int nfaces;
	int pcount;
	Face** faces;

	int nedges;
	int ecount;
	Edge** edges;
	double* eweights;

	// number of corners is the same as number of faces
	int ccount;
	Corner** corners;
	double total_angle;

	icVector3 normal;
	int forward;

	icMatrix4x4 error_quad;
	std::vector<PairContraction*> pairs;

	void* other_props;

	// methods
public:

	Vertex(double xx, double yy, double zz)
	{
		x = xx;
		y = yy;
		z = zz;
	}

	icVector3 pos()
	{
		return icVector3(x, y, z);
	}

	void set_pos(double xx, double yy, double zz)
	{
		x = xx;
		y = yy;
		z = zz;
	}

	icVector4 pos4()
	{
		return icVector4(x, y, z, 1);
	}

	// computes geometric error of the current vertex position
	double error()
	{
		return dot((this->pos4() * this->error_quad),this->pos4());
	}
};

class Edge
{
public:

	int index;
	Vertex* verts[2];

	int nfaces;
	Face **faces;

	// number of corners is the same as number of faces
	int ccount;
	Corner **corners;

	double length;

	Vertex* getOtherVert(Vertex* vert)
	{
		if (vert == verts[0])
		{
			return verts[1];
		}
		else if (vert == verts[1])
		{
			return verts[0];
		}
		else
		{
			return nullptr;
		}
	}

	Face* getOtherFace(Face* face)
	{
		if (nfaces == 2)
		{
			if (face == faces[0])
			{
				return faces[1];
			}
			else if (face == faces[1])
			{
				return faces[0];
			}
			else
			{
				return nullptr;
			}
		}
		else
		{
			return nullptr;
		}
	}
};

class Corner
{
public:

	int index;
	double angle;
	double angle_cot;

	Vertex* vert;
	Face* face;
	Edge* edge;
	Corner* next;
	Corner* prev;
	Corner* oppo;

	Corner(Face* face_in, Edge* edge_in, Vertex* vert_in)
	{
		face = face_in;
		edge = edge_in;
		vert = vert_in;
	}
};

class Face
{
public:

	int index;

	Vertex **verts;
	int nverts;

	Edge **edges;
	int nedges;

	Corner* corners[3];

	icVector3 normal;
	icVector4 plane;
	
	icVector3 center;
	double area;
	void* other_props;

	Face(int i)
	{
		index = i;
	}
};

//////////////////////////////////////////////////////////
// Contraction Classes
//////////////////////////////////////////////////////////

class PairContraction
{
public:

	Vertex* v1;
	Vertex* v2;
	Edge* edge = nullptr;

	icMatrix4x4 error_quad;
	icVector4 target;
	double error;
	bool allowed = true;

public:

	PairContraction(Vertex* vert1, Vertex* vert2, Edge* edge_in)
	{
		// set up field variables
		v1 = vert1;
		v2 = vert2;
		edge = edge_in;
		error_quad = v1->error_quad + v2->error_quad;

		// compute optimum contraction target
		icVector4 tvect(0.0, 0.0, 0.0, 1.0);
		icMatrix4x4 temp(error_quad);
		temp.entry[3][0] = 0;
		temp.entry[3][1] = 0;
		temp.entry[3][2] = 0;
		temp.entry[3][3] = 1;

		if (determinant(temp) != 0.0)
		{
			target = inverse(temp) * tvect;
		}
		else
		{
			// if the matrix is not invertible, choose between midpoint and endpoints
			icVector4 end1(v1->x, v1->y, v1->z, 1);
			icVector4 end2(v2->x, v2->y, v2->z, 1);
			icVector4 mid((v1->x + v2->x) / 2, (v1->y + v2->y) / 2, (v1->z + v2->z) / 2, 1);

			double err1 = dot((end1 * error_quad), end1);
			double err2 = dot((end2 * error_quad), end2);
			double emid = dot((mid * error_quad), mid);

			if (err1 < err2 && err1 < emid)
			{
				target = end1;
			}
			else if (err2 < err1 && err2 < emid)
			{
				target = end2;
			}
			else
			{
				target = mid;
			}
		}

		error = dot((target * error_quad), target);
	}

	void updateError()
	{
		// compute optimum contraction target
		icVector4 tvect(0.0, 0.0, 0.0, 1.0);
		icMatrix4x4 temp(error_quad);
		temp.entry[3][0] = 0;
		temp.entry[3][1] = 0;
		temp.entry[3][2] = 0;
		temp.entry[3][3] = 1;

		if (determinant(temp) != 0.0)
		{
			target = inverse(temp) * tvect;
		}
		else
		{
			// if the matrix is not invertible, choose between midpoint and endpoints
			icVector4 end1(v1->x, v1->y, v1->z, 1);
			icVector4 end2(v2->x, v2->y, v2->z, 1);
			icVector4 mid((v1->x + v2->x) / 2, (v1->y + v2->y) / 2, (v1->z + v2->z) / 2, 1);

			double err1 = dot((end1 * error_quad), end1);
			double err2 = dot((end2 * error_quad), end2);
			double emid = dot((mid * error_quad), mid);

			if (err1 < err2 && err1 < emid)
			{
				target = end1;
			}
			else if (err2 < err1 && err2 < emid)
			{
				target = end2;
			}
			else
			{
				target = mid;
			}
		}

		error = dot((target * error_quad), target);
	}

	Vertex* getOtherVert(Vertex* vert)
	{
		if (vert == v1) { return v2; }
		else if (vert == v2) { return v1; }
		else { return nullptr; }
	}
};

class PairCompare
{
public:
	int operator() (const PairContraction* p1, const PairContraction* p2)
	{
		return (p1->error > p2->error);
	}
};

//////////////////////////////////////////////////////////
// Drawing utility classes
//////////////////////////////////////////////////////////

class LineSegment
{
public:
	LineSegment(icVector3 s, icVector3 e) {
		start = s;
		end = e;
		len = length(end - start);
	};
	LineSegment(double sx, double sy, double sz, double ex, double ey, double ez) {
		start = icVector3(sx, sy, sz);
		end = icVector3(ex, ey, ez);
		len = length(end - start);
	};

	icVector3 start, end;  // all in local coordinate systems
	double len;
};

typedef std::vector<LineSegment> PolyLine;