#pragma once

#include <math.h>
#include <fstream>
#include <float.h>
#include <qmath.h>
#include "ply.h"
#include "ply_io.h"
#include "elements.hpp"
#include "icVector.hpp"
#include <iostream>
#include <vector>
#include <queue>

#define WEIGHT_UNIF			1
#define WEIGHT_CORD			2
#define WEIGHT_FLOW			3
#define WEIGHT_MEAN			4

class Polyhedron
{
	// private fields
private:

	PlyFile* in_ply;
	bool TRI_MESH = true;

	// fields
public:

	// faces
	Face** flist;
	int nfaces;
	int selected;

	// edges
	Edge** elist;
	int nedges;

	// vertices
	Vertex** vlist;
	int nverts;
	
	// corners
	Corner** clist;
	int ncorners;

	icVector3 center;
	double radius;
	double mean_elength;
	double area;
	unsigned char orientation; // 0 for ccw, 1 for cw

	PlyOtherProp* vert_other;
	PlyOtherProp* face_other;

	// contraction containers
	std::priority_queue<PairContraction*, std::vector<PairContraction*>, PairCompare> pair_queue;

// private helper functions
private:

	//////////////////////////////////////////////////////////
	// Pointer and Edge creation functions
	//////////////////////////////////////////////////////////

	// create pointers between verts, edges, and faces
	void create_pointers()
	{
		// index verts and quads
		for (int i = 0; i < nverts; i++)
		{
			vlist[i]->index = i;
		}
		for (int i = 0; i < nfaces; i++)
		{
			flist[i]->index = i;
		}

		// create pointers from verts to faces
		vertex_to_face_ptrs();

		// make edges
		create_edges();

		// create pointers from verts to edges
		vertex_to_edge_ptrs();

		// index edges
		for (int i = 0; i < nedges; i++)
		{
			elist[i]->index = i;
		}

		if (TRI_MESH)
		{
			// make corners
			create_corners();
			corner_to_corner_ptrs();
			calculate_corner_angles();
			calculate_vertex_angles();
		}
	}

	// create vert to face pointers
	void vertex_to_face_ptrs()
	{
		Face* face;
		Vertex* vert;

		// zero out vertex Face count
		for (int i = 0; i < nverts; i++)
		{
			vlist[i]->pcount = 0;
		}

		// count face pointers needed for each vert
		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				face->verts[j]->pcount++;
			}
		}

		// allocate memory for face pointers in each vert
		for (int i = 0; i < nverts; i++)
		{
			vlist[i]->faces = new Face * [vlist[i]->pcount];
			vlist[i]->nfaces = 0;
		}

		// create face pointers in each vertex
		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				vert = face->verts[j];
				vert->faces[vert->nfaces] = face;
				vert->nfaces++;
			}
		}
	}

	void create_edges()
	{
		Face* face;
		Vertex* v1;
		Vertex* v2;
		double edge_estimate = 0;

		// get rough estimate for number of edges
		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				v1 = face->verts[j];
				v2 = face->verts[(j + 1) % (face->nverts)];
				double adjacent_faces = adjacent_face_count(face, v1, v2);
				edge_estimate += 1 / (adjacent_faces + 1);
			}
		}

		// allocate memory for edges with a little extra
		elist = new Edge * [(int)(edge_estimate * 1.1)];

		// zero out all pointers from faces to edges
		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			face->nedges = 0;
			for (int j = 0; j < face->nverts; j++)
			{
				face->edges[j] = nullptr;
			}
		}

		// create edges for each Face

		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				v1 = face->verts[j];
				v2 = face->verts[(j + 1) % (face->nverts)];
				// check for duplicate edges
				bool duplicate = false;
				Edge* edge = nullptr;
				for (int k = 0; k < face->nedges; k++)
				{
					edge = face->edges[k];
					if ((v1 == edge->verts[0] && v2 == edge->verts[1]) ||
						(v2 == edge->verts[0] && v1 == edge->verts[1]))
					{
						duplicate = true;
					}
				}
				if (!duplicate)
				{
					create_edge(v1, v2);
				}
			}
		}

		std::cout << std::to_string(nedges) << " edges created" << std::endl;
	}

	void vertex_to_edge_ptrs()
	{
		Edge* edge;
		Vertex* vert;

		// zero out the number of edges for each vertex
		for (int i = 0; i < nverts; i++)
		{
			vlist[i]->ecount = 0;
		}

		// count edge pointers needed for each vertex
		for (int i = 0; i < nedges; i++)
		{
			edge = elist[i];
			edge->verts[0]->ecount++;
			edge->verts[1]->ecount++;
		}

		// allocate memory for edge pointers in vertices
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			vert->edges = new Edge * [vert->ecount];
			vert->nedges = 0;
		}

		// create edge pointers in vertices
		for (int i = 0; i < nedges; i++)
		{
			edge = elist[i];
			for (int j = 0; j < 2; j++)
			{
				vert = edge->verts[j];
				vert->edges[vert->nedges] = edge;
				vert->nedges++;
			}
		}
	}

	// counts the number of Faces sharing an edge with a given face
	int adjacent_face_count(Face* face, Vertex* v1, Vertex* v2)
	{
		Face* adjacent = nullptr;
		int count = 0;

		// loop through faces of first vertex
		for (int i = 0; i < v1->nfaces; i++)
		{
			adjacent = v1->faces[i];
			if (adjacent == face) { continue; }
			// look for a match to the second vertex
			for (int j = 0; j < adjacent->nverts; j++)
			{
				if (adjacent->verts[j] == v2) { count++; }
			}
		}

		return count;
	}

	// creates an edge between the given vertices
	void create_edge(Vertex* v1, Vertex* v2)
	{
		Face* face;

		// create the edge
		elist[nedges] = new Edge;
		Edge* edge = elist[nedges];
		edge->index = nedges;
		edge->verts[0] = v1;
		edge->verts[1] = v2;
		nedges++;

		// zero out edge face number
		edge->nfaces = 0;

		// count all quads associated with this edge
		for (int i = 0; i < v1->nfaces; i++)
		{
			face = v1->faces[i];
			for (int j = 0; j < face->nverts; j++)
			{
				if (face->verts[j] == v2)
				{
					edge->nfaces++;
					break;
				}
			}
		}

		// allocate memory for Face pointers
		edge->faces = new Face * [edge->nfaces];

		// create pointers from edges to faces and vice versa
		edge->nfaces = 0;
		for (int i = 0; i < v1->nfaces; i++)
		{
			face = v1->faces[i];

			for (int j = 0; j < face->nverts; j++)
			{
				if (face->verts[j] == v2)
				{
					edge->faces[edge->nfaces] = face;
					edge->nfaces++;

					face->edges[face->nedges] = edge;
					face->nedges++;

					break;
				}
			}
		}
	}

	//////////////////////////////////////////////////////////
	// Corner creation and pointer functions 
	// (only for tri meshes)
	//////////////////////////////////////////////////////////

	void create_corners()
	{
		Face* face;
		Edge* edge;
		Vertex* vert;
		Corner* corn;

		// allocate memory in corner list, edges, verts
		clist = new Corner * [(nfaces * 3)];
		ncorners = 0;
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			vert->corners = new Corner * [vert->nfaces];
			vert->ccount = 0;
		}
		for (int i = 0; i < nedges; i++)
		{
			edge = elist[i];
			edge->corners = new Corner * [edge->nfaces];
			edge->ccount = 0;
		}

		// create corners for each face
		for (int i = 0; i < nfaces; i++)
		{
			face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				vert = face->verts[j];
				for (int k = 0; k < face->nedges; k++)
				{
					// find opposite edge
					if (face->edges[k]->verts[0] != vert && 
						face->edges[k]->verts[1] != vert)
					{
						edge = face->edges[k];
					}
				}

				// create corner
				corn = new Corner(face, edge, vert);
				corn->index = ncorners;
				clist[ncorners] = corn;
				ncorners++;

				// create pointers in face, edge, vert
				face->corners[j] = corn;
				edge->corners[edge->ccount] = corn;
				edge->ccount++;
				vert->corners[vert->ccount] = corn;
				vert->ccount++;
			}
		}

		std::cout << std::to_string(ncorners) << " corners created" << std::endl;
	}

	void corner_to_corner_ptrs()
	{
		Corner* corner;
		for (int i = 0; i < nfaces; i++)
		{
			for (int j = 0; j < flist[i]->nverts; j++)
			{
				corner = flist[i]->corners[j];
				corner->next = flist[i]->corners[(j + 1) % 3];
				corner->prev = flist[i]->corners[(j + 2) % 3];
				if (corner->edge->ccount != 2)
				{
					corner->oppo = nullptr;
				}
				else
				{
					if (corner->edge->corners[0] != corner)
					{
						corner->oppo = corner->edge->corners[0];
					}
					else if (corner->edge->corners[1] != corner)
					{
						corner->oppo = corner->edge->corners[1];
					}
					else
					{
						corner->oppo = nullptr;
					}
				}
			}
		}
	}

	void calculate_corner_angles()
	{
		Vertex* v0;
		Vertex* v1;
		Vertex* v2;
		icVector3 vect1;
		icVector3 vect2;

		for (int i = 0; i < ncorners; i++)
		{
			v0 = clist[i]->vert;
			v1 = clist[i]->next->vert;
			v2 = clist[i]->prev->vert;

			vect1 = icVector3(v1->x - v0->x, v1->y - v0->y, v1->z - v0->z);
			vect2 = icVector3(v2->x - v0->x, v2->y - v0->y, v2->z - v0->z);

			clist[i]->angle = acos(dot(vect1, vect2) / (length(vect1) * length(vect2)));
			clist[i]->angle_cot = 1.0 / tan(clist[i]->angle);
		}
	}

	void calculate_vertex_angles()
	{
		Vertex* vert;
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			vert->total_angle = 0.0;

			// sum angles of all associated corners
			for (int j = 0; j < vert->nfaces; j++)
			{
				vert->total_angle += vert->corners[j]->angle;
			}
		}
	}

	//////////////////////////////////////////////////////////
	// Dimension and normal creation functions
	//////////////////////////////////////////////////////////

	void calculate_dimensions()
	{
		// calculate edge lengths
		Edge* edge;
		double sum = 0.0;
		for (int i = 0; i < nedges; i++)
		{
			edge = elist[i];
			//calculate the length of the edge
			edge->length = sqrt((edge->verts[0]->x - edge->verts[1]->x) * (edge->verts[0]->x - edge->verts[1]->x)
				+ (edge->verts[0]->y - edge->verts[1]->y) * (edge->verts[0]->y - edge->verts[1]->y)
				+ (edge->verts[0]->z - edge->verts[1]->z) * (edge->verts[0]->z - edge->verts[1]->z));

			sum += edge->length;
		}
		mean_elength = (sum / (double)nedges);

		icVector3 min, max;
		Vertex* vert;
		// find min and max values to find overall dimensions
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			if (i == 0)
			{
				min.set(vert->x, vert->y, vert->z);
				max.set(vert->x, vert->y, vert->z);
			}
			else
			{
				if (vert->x < min.x) { min.x = vert->x; }
				if (vert->x > max.x) { max.x = vert->x; }
				if (vert->y < min.y) { min.y = vert->y; }
				if (vert->y > max.y) { max.y = vert->y; }
				if (vert->z < min.z) { min.z = vert->z; }
				if (vert->z > max.z) { max.z = vert->z; }
			}
		}
		center = (min + max) * 0.5;
		radius = length(center - min);
	}

	void create_normals()
	{
		icVector3 v0, v1, v2;
		icVector3 poly_center = center;
		Face* face;
		Vertex* vert;
		double signed_volume = 0.0;

		area = 0.0; // zero out mesh area

		// get normal vectors for faces and calculate face areas
		for (int i = 0; i < nfaces; i++)
		{
			// get face normal vectors
			face = flist[i];
			v0.set(face->verts[0]->x, face->verts[0]->y, face->verts[0]->z);
			v1.set(face->verts[1]->x, face->verts[1]->y, face->verts[1]->z);
			v2.set(face->verts[2]->x, face->verts[2]->y, face->verts[2]->z);
			face->normal = cross(v0 - v1, v2 - v1);
			normalize(face->normal);

			// calculate area of each face
			calc_face_area(face);
			area += face->area;

			signed_volume += dot(poly_center - v0, face->normal) * face->area;
		}

		// determine face orientations
		signed_volume /= area;
		if (signed_volume < 0)
		{
			orientation = 0;
		}
		else
		{
			orientation = 1;
			for (int i = 0; i < nfaces; i++)
			{
				flist[i]->normal *= -1;
			}
		}

		// get vertex normals by averaging face normals
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			vert->normal.set(0.0);
			for (int j = 0; j < vert->nfaces; j++)
			{
				face = vert->faces[j];
				vert->normal += face->normal;
			}
			normalize(vert->normal);
		}
	}

	void calc_face_area(Face* face)
	{
		icVector3 v1, v2, sum, center;
		
		for (int i = 0; i < face->nverts; i++)
		{
			v1.set(face->verts[i]->x, face->verts[i]->y, face->verts[i]->z);
			v2.set(face->verts[(i + 1) % face->nverts]->x, 
				   face->verts[(i + 1) % face->nverts]->y,
				   face->verts[(i + 1) % face->nverts]->z);

			sum = sum + cross(v1, v2);

			center.x += face->verts[i]->x;
			center.y += face->verts[i]->y;
			center.z += face->verts[i]->z;
		}
		
		face->area = 0.5 * length(sum);
		center /= (double)face->nverts;
		face->center = center;
	}

	//////////////////////////////////////////////////////////
	// Functions for quadric error simplification
	//////////////////////////////////////////////////////////

	// calculate the error quadric for each vertex. Only call 
	// this once during mesh initialization
	void calc_error_quadrics()
	{
		Face* face;
		Vertex* vert;
		icVector4 plane;

		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			vert->error_quad = 0.0;
			for (int j = 0; j < vert->nfaces; j++)
			{
				face = vert->faces[j];
				plane.set
				(
					face->normal.x,
					face->normal.y,
					face->normal.z,
					-(	(face->normal.x * face->verts[0]->x) +
						(face->normal.y * face->verts[0]->y) + 
						(face->normal.z * face->verts[0]->z) )
				);

				vert->error_quad += square(plane);
			}
		}
	}

	// find all eligable pair contractions and put them in a min heap
	void find_valid_pairs(double threshold)
	{
		Vertex* vert;
		Vertex* oppo;
		Edge* edge;
		icVector3 temp;
		Edge* pair_edge = nullptr;
		PairContraction* contraction;

		// search through all vertex pairs
		for (int i = 0; i < nverts; i++)
		{
			vert = vlist[i];
			

			for (int j = (i+1); j < nverts; j++)
			{
				oppo = vlist[j];
				pair_edge = nullptr;

				// check to see if it is an edge
				for (int k = 0; k < vert->nedges; k++)
				{
					edge = vert->edges[k];
					if (edge->getOtherVert(vert) == oppo)
					{
						pair_edge = edge;
						break;
					}
				}
				
				// check the distance between points
				temp = (oppo->pos() - vert->pos());
				if ((length(temp) < threshold) || (pair_edge != nullptr))
				{
					// create contraction and insert into priority queue
					contraction = new PairContraction(vert, oppo, pair_edge);
					pair_queue.push(contraction);
					vert->pairs.push_back(contraction);
					oppo->pairs.push_back(contraction);
				}
			}
		}
	}

	// contract a vertex pair and update mesh
	void contract_pair(PairContraction* pair)
	{
		// move both vertices to the new target position
		pair->v1->set_pos(pair->target.x, pair->target.y, pair->target.z);
		pair->v2->set_pos(pair->target.x, pair->target.y, pair->target.z);

		// if pair is an edge, delete degenerate faces and edge and update pointers
		if (pair->edge != nullptr)
		{
			Face* face;
			Face* face_update;
			Edge* edge;
			Edge* edge_update;
			Corner* corn;
			Corner* oppv1;
			Corner* oppv2;
			Vertex* vert;

			if (pair->edge->index == 2641)
			{
				int test = 0;
				face = pair->edge->faces[0];
				face = pair->edge->faces[1];
				corn = pair->edge->corners[0];
				corn = pair->edge->corners[1];
			}
			// remove degenerate faces
			for (int i = 0; i < pair->edge->nfaces; i++)
			{
				face = pair->edge->faces[i];
				// remove degenerate corners
				for (int j = 0; j < 3; j++)
				{
					corn = face->corners[j];
					if (corn->vert == pair->v1)
					{ 
						oppv1 = corn->oppo;
					}
					if (corn->vert == pair->v2)
					{
						oppv2 = corn->oppo;
					}
					// remove all references to corner and delete
					remove_corner(corn);
				}
				for (int j = 0; j < 3; j++)
				{
					vert = face->verts[j];
				}
				// remove degenerate edges
				for (int j = 0; j < 3; j++)
				{
					edge = face->edges[j];
					if (edge == pair->edge)
					{
						continue;
					}
					else if ((edge->verts[0] == pair->v1) || (edge->verts[1] == pair->v1))
					{
						edge_update = edge;
					}
					else if ((edge->verts[0] == pair->v2) || (edge->verts[1] == pair->v2))
					{
						face_update = edge->getOtherFace(face);
						remove_edge(edge);
					}
				}

				// remove all references to pair edge and delete
				if (i == (pair->edge->nfaces - 1))
				{
					remove_edge(pair->edge);
				}
				
				// remove all references to face and delete
				remove_face(face);

				// update pointers
				oppv1->edge = edge_update;
				oppv1->oppo = oppv2;
				oppv2->oppo = oppv1;
				for (int j = 0; j < edge_update->nfaces; j++)
				{
					if (edge_update->corners[j] == nullptr)
					{
						edge_update->corners[j] = oppv1;
					}
				}
				for (int j = 0; j < edge_update->nfaces; j++)
				{
					if (edge_update->faces[j] == nullptr)
					{
						edge_update->faces[j] = face_update;
					}
				}
				for (int j = 0; j < 3; j++)
				{
					if (face_update->edges[j] == nullptr)
					{
						face_update->edges[j] = edge_update;
					}
				}
			}
		}

		// disallow pair for future contraction and transfer all elements attached to v2 to v1 and delete v2
		pair->allowed = false;
		remove_vert(pair->v2, pair->v1);
		delete(pair);
	}

	//////////////////////////////////////////////////////////
	// Element deletion functions
	//////////////////////////////////////////////////////////

	void remove_corner(Corner* corn)
	{
		// make all references in other elements nullptr
		if (corn->next != nullptr) { corn->next->prev = nullptr; }
		if (corn->prev != nullptr) { corn->prev->next = nullptr; }
		if (corn->oppo != nullptr) { corn->oppo->oppo = nullptr; }
		if (corn->vert != nullptr)
		{
			for (int i = 0; i < (corn->vert->nfaces); i++)
			{
				if (corn->vert->corners[i] == corn)
				{
					corn->vert->corners[i] = corn->vert->corners[corn->vert->nfaces-1];
					break;
				}
			}
		}
		if (corn->edge != nullptr)
		{
			for (int i = 0; i < corn->edge->nfaces; i++)
			{
				if (corn->edge->corners[i] == corn)
				{
					
					corn->edge->corners[0] = nullptr;
					break;
				}
			}
		}
		if (corn->face != nullptr)
		{
			for (int i = 0; i < 3; i++)
			{
				if (corn->face->corners[i] == corn)
				{
					corn->face->corners[i] = nullptr;
					break;
				}
			}
		}
		
		// remove from corner list and delete object
		ncorners--;
		clist[corn->index] = clist[ncorners];
		clist[corn->index]->index = corn->index;
		delete(corn);
	}

	void remove_edge(Edge* edge)
	{
		// make all references in other elements nullptr
		for (int i = 0; i < 2; i++)
		{
			if (edge->verts[i] != nullptr)
			{
				for (int j = 0; j < edge->verts[i]->nedges; j++)
				{
					if (edge->verts[i]->edges[j] == edge)
					{
						edge->verts[i]->nedges--;
						edge->verts[i]->edges[j] = edge->verts[i]->edges[edge->verts[i]->nedges];
						break;
					}
				}
			}
		}
		for (int i = 0; i < edge->nfaces; i++)
		{
			if (edge->faces[i] != nullptr)
			{
				for (int j = 0; j < 3; j++)
				{
					if (edge->faces[i]->edges[j] == edge)
					{
						edge->faces[i]->edges[j] = nullptr;
						break;
					}
				}
			}
		}
		for (int i = 0; i < edge->nfaces; i++)
		{
			if (edge->corners[i] != nullptr)
			{
				edge->corners[i]->edge = nullptr;
			}
		}
		for (int i = 0; i < edge->verts[0]->pairs.size(); i++)
		{
			if (edge->verts[0]->pairs[i]->edge == edge)
			{
				edge->verts[0]->pairs[i]->allowed = false;
			}
		}
		for (int i = 0; i < edge->verts[1]->pairs.size(); i++)
		{
			if (edge->verts[1]->pairs[i]->edge == edge)
			{
				edge->verts[1]->pairs[i]->allowed = false;
			}
		}

		// remove from edge list and delete object
		nedges--;
		elist[edge->index] = elist[nedges];
		elist[edge->index]->index = edge->index;
		delete(edge);
	}

	void remove_face(Face* face)
	{
		// make all references in other elements nullptr
		for (int i = 0; i < 3; i++)
		{
			if (face->verts[i] != nullptr)
			{
				for (int j = 0; j < face->verts[i]->nfaces; j++)
				{
					if (face->verts[i]->faces[j] == face)
					{
						face->verts[i]->nfaces--;
						face->verts[i]->faces[j] = face->verts[i]->faces[face->verts[i]->nfaces];
						break;
					}
				}
			}
		}
		for (int i = 0; i < 3; i++)
		{
			if (face->edges[i] != nullptr)
			{
				for (int j = 0; j < face->edges[i]->nfaces; j++)
				{
					if (face->edges[i]->faces[j] == face)
					{
						face->edges[i]->faces[j] = nullptr;
						break;
					}
				}
			}
		}
		for (int i = 0; i < 3; i++)
		{
			if (face->corners[i] != nullptr)
			{
				face->corners[i]->face = nullptr;
			}
		}

		// remove from face list and delete object
		nfaces--;
		flist[face->index] = flist[nfaces];
		flist[face->index]->index = face->index;
		delete(face);
	}

	void remove_vert(Vertex* vrem, Vertex* vnew)
	{
		// initialize new arrays for vertex faces, edges, and corners
		int ccount, ecount, fcount;
		ccount = ecount = fcount = 0;
		Corner** ctemp = new Corner * [vrem->nfaces + vnew->nfaces];
		Edge** etemp = new Edge * [vrem->nedges + vnew->nedges];
		Face** ftemp = new Face * [vrem->nfaces + vnew->nfaces];
		Vertex* vtemp;

		// replace all references to vrem with references to vnew
		for (int i = 0; i < vrem->nfaces; i++)
		{
			if (vrem->corners[i] != nullptr)
			{
				vrem->corners[i]->vert = vnew;
				ctemp[ccount] = vrem->corners[i];
				ccount++;
			}
		}
		for (int i = 0; i < vrem->nedges; i++)
		{
			if (vrem->edges[i] != nullptr)
			{
				if (vrem->edges[i]->verts[0] == vrem) { vrem->edges[i]->verts[0] = vnew; }
				if (vrem->edges[i]->verts[1] == vrem) { vrem->edges[i]->verts[1] = vnew; }
				etemp[ecount] = vrem->edges[i];
				ecount++;
			}
		}
		for (int i = 0; i < vrem->nfaces; i++)
		{
			if (vrem->faces[i] != nullptr)
			{
				for (int j = 0; j < 3; j++)
				{
					if (vrem->faces[i]->verts[j] == vrem)
					{
						vrem->faces[i]->verts[j] = vnew;
						ftemp[fcount] = vrem->faces[i];
						fcount++;
					}
				}
			}
		}
		for (PairContraction* pair : vrem->pairs)
		{
			if (pair->v1 == vrem) { pair->v1 = vnew; }
			else if (pair->v2 == vrem) { pair->v2 = vnew; }
			if (pair->allowed)
			{
				vnew->pairs.push_back(pair);
			}
		}

		// carry over all references originally from vnew
		for (int i = 0; i < vnew->nedges; i++)
		{
			if (vnew->edges[i] != nullptr)
			{
				etemp[ecount] = vnew->edges[i];
				ecount++;
			}
		}
		for (int i = 0; i < vnew->nfaces; i++)
		{
			if (vnew->faces[i] != nullptr)
			{
				ftemp[fcount] = vnew->faces[i];
				fcount++;
			}
		}
		for (int i = 0; i < vnew->nfaces; i++)
		{
			if (vnew->corners[i] != nullptr)
			{
				ctemp[ccount] = vnew->corners[i];
				ccount++;
			}
		}

		// add error quadrics
		vnew->error_quad += vrem->error_quad;

		// delete old vertex lists and replace with temp lists
		delete[](vnew->corners);
		delete[](vnew->edges);
		delete[](vnew->faces);
		vnew->corners = ctemp;
		vnew->edges = etemp;
		vnew->faces = ftemp;
		vnew->nedges = ecount;
		vnew->nfaces = fcount;

		// remove disallowed contractions and update values for others
		for (int i = 0; i < vnew->pairs.size(); )
		{
			if (!vnew->pairs[i]->allowed)
			{
				vnew->pairs.erase(vnew->pairs.begin() + i);
			}
			else
			{
				// update pair costs
				vnew->pairs[i]->updateError();
				i++;
			}
		}

		// remove from vert list and delete object
		nverts--;
		vlist[vrem->index] = vlist[nverts];
		vlist[vrem->index]->index = vrem->index;
		delete(vrem);
	}

	//////////////////////////////////////////////////////////
	// For printing information about the mesh
	//////////////////////////////////////////////////////////

	double total_error()
	{
		double error_sum = 0;
		for (int i = 0; i < nverts; i++)
		{
			error_sum += vlist[i]->error();
		}
		return error_sum;
	}

	void print_vef()
	{
		std::cout << "Geometric Error = " << std::to_string(total_error()) << std::endl;
		int v_e_f = nverts - nedges + nfaces;
		std::cout << "V-E+F = " << std::to_string(v_e_f) << std::endl << std::endl;
	}

	void print_stats()
	{
		std::cout << "# Verts = " << std::to_string(nverts) << std::endl;
		std::cout << "# Edges =	" << std::to_string(nedges) << std::endl;
		std::cout << "# Faces = " << std::to_string(nfaces) << std::endl;
		std::cout << "# Corners = " << std::to_string(ncorners) << std::endl;
		std::cout << "Geometric Error = " << std::to_string(total_error()) << std::endl;
		print_vef();
	}

// methods
public:

	// full constructor
	Polyhedron(FILE* file)
	{
		int elem_count;
		char* elem_name;
		in_ply = read_ply(file);

		for (int i = 0; i < in_ply->num_elem_types; i++)
		{
			elem_name = setup_element_read_ply(in_ply, i, &elem_count);

			if (equal_strings("vertex", elem_name))
			{
				// vertex list
				nverts = elem_count;
				vlist = new Vertex * [nverts];

				// setup for extracting vertex elements
				setup_property_ply(in_ply, &vert_props[0]);
				setup_property_ply(in_ply, &vert_props[1]);
				setup_property_ply(in_ply, &vert_props[2]);

				vert_other = get_other_properties_ply(in_ply, offsetof(Vertex_io, other_props));

				// extract all vertex elements
				for (int j = 0; j < nverts; j++)
				{
					Vertex_io vert;
					get_element_ply(in_ply, (void*)&vert);

					// copy info from vert structure
					vlist[j] = new Vertex(vert.x, vert.y, vert.z);
					vlist[j]->other_props = vert.other_props;
				}

				std::cout << std::to_string(nverts) << " vertices created" << std::endl;
			}
			else if (equal_strings("face", elem_name))
			{
				// face list
				nfaces = elem_count;
				flist = new Face * [nfaces];

				// setup for extracting face elements
				setup_property_ply(in_ply, &face_props[0]);
				face_other = get_other_properties_ply(in_ply, offsetof(Face_io, other_props));

				// extract face elements
				for (int j = 0; j < nfaces; j++)
				{
					Face_io face;
					get_element_ply(in_ply, (void*)&face);

					flist[j] = new Face(j);
					flist[j]->nverts = face.nverts;
					flist[j]->nedges = 0;
					flist[j]->verts = new Vertex* [face.nverts];
					flist[j]->edges = new Edge * [face.nverts];

					int vert;
					for (int k = 0; k < face.nverts; k++)
					{
						flist[j]->verts[k] = (Vertex*)face.verts[k];
					}

					flist[j]->other_props = face.other_props;

					// check if mesh is a tri mesh
					if (face.nverts != 3) { TRI_MESH = false; }
				}

				std::cout << std::to_string(nfaces) << " faces created" << std::endl;
			}
			else
			{
				get_other_element_ply(in_ply);
			}
		}

		// close ply file
		close_ply(in_ply);

		// match up pointers from faces to verts
		for (int i = 0; i < nfaces; i++)
		{
			Face* face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				face->verts[j] = vlist[(int)flist[i]->verts[j]];
			}
		}

		// eliminate Faces that use the same vertex more than once make sure pointers match up
		for (int i = nfaces-1; i > 0; i--)
		{
			Face* face = flist[i];
			for (int j = 0; j < face->nverts; j++)
			{
				Vertex* vert = face->verts[j];
				for (int k = 0; k < face->nverts; k++)
				{
					if (j != k && vert == face->verts[k])
					{
						free(flist[i]);
						nfaces--;
						flist[i] = flist[nfaces];
						flist[i]->index = i;
					}
				}
			}
		}

		// initialize edge count
		nedges = 0;
	}

	// empty constructor
	Polyhedron()
	{
		nverts = nedges = nfaces = ncorners = 0;
		
		vlist = new Vertex * [64];
		flist = new Face * [64];
		elist = new Edge * [64];
		clist = new Corner * [64];
	}

	// set up pointers, create edges, calculate normals, create corners
	void initialize()
	{
		create_pointers();
		calculate_dimensions();
		create_normals();
		calc_error_quadrics();

		if(TRI_MESH) { print_vef(); }
	}

	void finalize()
	{
		for (int i = 0; i < nfaces; i++)
		{
			//free(flist[i]->other_props);
			delete[](flist[i]->edges);
			delete[](flist[i]->verts);
			delete(flist[i]);
		}
		for (int i = 0; i < nedges; i++) {
			delete[](elist[i]->faces);
			delete(elist[i]);
		}
		for (int i = 0; i < ncorners; i++)
		{
			delete(clist[i]);
		}
		for (int i = 0; i < nverts; i++) {
			delete[](vlist[i]->faces);
			delete[](vlist[i]->edges);
			if(TRI_MESH){ delete[](vlist[i]->corners); }
			//free(vlist[i]->other_props);
			delete(vlist[i]);
		}
		while (!pair_queue.empty())
		{
			delete(pair_queue.top());
			pair_queue.pop();
		}

		delete[](flist);
		delete[](elist);
		if (TRI_MESH) { delete[](clist); }
		delete[](vlist);
		if (!vert_other) { delete(vert_other); }
		if (!face_other) { delete(face_other); }
	}

	void smooth(int weight_scheme, double step, int iterations)
	{
		if (!TRI_MESH) { return; }

		// allocate vector array for holding holding old vert positions
		icVector3** vold = new icVector3 * [nverts];
		Vertex* vert;
		for (int i = 0; i < nverts; i++)
		{
			vold[i] = new icVector3();
		}

		// Vertex operations
		double weight_sum;
		double weight;
		icVector3 vmove = icVector3(0.0, 0.0, 0.0);
		Vertex* avert;
		Edge* edge;
		Corner* corn;
		for (int k = 0; k < iterations; k++)
		{
			// copy vert positions into old positions list
			for (int i = 0; i < nverts; i++)
			{
				vert = vlist[i];
				vold[i]->set(vert->x, vert->y, vert->z);
			}

			// update positions of all vertices
			for (int i = 0; i < nverts; i++)
			{
				weight_sum = 0.0;
				vmove.set(0.0, 0.0, 0.0);
				vert = vlist[i];
				switch (weight_scheme)
				{
				case WEIGHT_UNIF:
				{
					// calculate weight normalizing sum
					weight_sum = (double)vert->nedges;
					// calculate move vector sum
					for (int j = 0; j < vert->nedges; j++)
					{
						edge = vert->edges[j];
						avert = edge->getOtherVert(vert);
						weight = 1.0;
						vmove += ((weight / weight_sum) * ((*vold[avert->index]) - (*vold[vert->index])));
					}
				}
				break;

				case WEIGHT_CORD:
				{
					// calculate weight normalizing sum
					for (int j = 0; j < vert->nedges; j++)
					{
						edge = vert->edges[j];
						weight_sum += (1.0 / edge->length);
					}

					// calculate move vector sum
					for (int j = 0; j < vert->nedges; j++)
					{
						edge = vert->edges[j];
						avert = edge->getOtherVert(vert);
						weight = 1.0 / edge->length;
						vmove += ((weight / weight_sum) * ((*vold[avert->index]) - (*vold[vert->index])));
					}
				}
				break;

				case WEIGHT_FLOW:
				{
					// calculate weight normalizing sum
					for (int j = 0; j < vert->nedges; j++)
					{
						corn = vert->corners[j];
						weight_sum += (((1.0 / tan(corn->next->angle)) + (1.0 / tan(corn->prev->angle))) / 2.0);
					}

					// calculate move vector sum
					for (int j = 0; j < vert->nedges; j++)
					{
						edge = vert->edges[j];
						avert = edge->getOtherVert(vert);
						weight = (((1.0 / tan(edge->corners[0]->angle)) + (1.0 / tan(edge->corners[1]->angle))) / 2.0);
						vmove += ((weight / weight_sum) * ((*vold[avert->index]) - (*vold[vert->index])));
					}
				}
				break;

				case WEIGHT_MEAN:
				{
					// calculate weight normalizing sum
					for (int j = 0; j < vert->nedges; j++)
					{
						corn = vert->corners[j];
						weight_sum += tan(corn->angle / 2.0);
					}

					// calculate move vector sum
					for (int j = 0; j < vert->nedges; j++)
					{
						corn = vert->corners[j];
						edge = corn->prev->edge;
						avert = edge->getOtherVert(vert);
						weight = ((tan(corn->angle / 2.0) + tan(corn->prev->oppo->prev->angle / 2.0)) / 2.0);
						vmove += ((weight / weight_sum) * ((*vold[avert->index]) - (*vold[vert->index])));
					}
				}
				break;
				}

				// update vertex position
				vmove *= step;
				vert->x += vmove.x;
				vert->y += vmove.y;
				vert->z += vmove.z;
			}

			// update edge lengths, face areas, and normals
			calculate_dimensions();
			create_normals();
		}

		// delete old vertex locations
		for (int i = 0; i < nverts; i++)
		{
			delete(vold[i]);
		}
		delete[](vold);
	}

	void pairSimplify(int face_target, double error_tolerance, int max_contractions)
	{
		find_valid_pairs(0.0);

		// iteratively contract pairs
		int count = 0;
		while (!pair_queue.empty())
		{
			// disallow pairs with deleted edges (find a better fix for this)
			if (pair_queue.top()->edge != nullptr)
			{
				if (pair_queue.top()->edge->index < 0 || pair_queue.top()->edge->index > nedges)
				{
					pair_queue.top()->allowed = false;
				}
			}

			// check if it is allowed
			if (pair_queue.top()->allowed)
			{
				contract_pair(pair_queue.top());
				count++;
			}
			pair_queue.pop();

			// check if a stopping point has been reached
			if ((nfaces <= face_target) || (total_error() > error_tolerance) || (count >= max_contractions))
			{
				break;
			}
		}

		calculate_dimensions();
		create_normals();
		calculate_corner_angles();
		calculate_vertex_angles();
		print_stats();
	}

};