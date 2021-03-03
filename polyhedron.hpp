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

	// curvature variables
	double gmean, gstdev;
	double mmean, mstdev;

	PlyOtherProp* vert_other;
	PlyOtherProp* face_other;

	std::vector<LineSegment*> silhouette;
	std::vector<LineSegment*> major_streamlines;
	std::vector<LineSegment*> minor_streamlines;

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
		// and compute vertex error quadrics
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

	// YOU ARE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// find all eligable pair contractions and put them in a priority queue

	//////////////////////////////////////////////////////////
	// For printing information about the mesh
	//////////////////////////////////////////////////////////

	void print_stats()
	{
		int v_e_f = nverts - nedges + nfaces;
		std::cout << "V-E+F = " << std::to_string(v_e_f) << std::endl << std::endl;
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

		// eliminate Faces that use the same vertex more than once
		// make sure pointers match up
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

		if(TRI_MESH) { print_stats(); }
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

};