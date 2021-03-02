#pragma once
#include "icVector.hpp"
#include "icMatrix.hpp"
#include <vector>

// forward declarations
class Face;
class Edge;
class Corner;

// classes
class Vertex
{
// fields
public:

	int index;
	double x, y, z;
	double texcoord[2];

	int nfaces;
	int pcount;
	Face **faces;

	int nedges;
	int ecount;
	Edge **edges;
	double * eweights;

	// number of corners is the same as number of faces
	int ccount;
	Corner **corners;
	double total_angle;

	icVector3 normal;
	int forward;

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
	
	icVector3 center;
	double area;
	void* other_props;

	Face(int i)
	{
		index = i;
	}
};

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
