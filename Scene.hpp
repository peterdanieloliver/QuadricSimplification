#pragma once

#include "polyhedron.hpp"
#include <QtOpenGL>
#include <GL/GLU.h>
#include <qmath.h>
#include "stb_image.h"

// display mode definitions
#define DISPLAY_SOLID				1
#define DISPLAY_WIREFRAME			2
#define DISPLAY_CHECKER				3
#define DISPLAY_TEXTURE				4

class Scene
{
private:
	
	// for defining checkerboard texture pattern
	double L = .25;

	// for scaling coloring schemes
	double factor = 2.0;

	// texture
	GLuint texID;

	// flag for needing update
	bool needs_update;

public:

	Polyhedron* mesh;
	icVector3 scene_center;
	double scene_radius;
	unsigned char scene_orientation;
	float color[4] = { 0.0,0.0,0.0,1.0 };

	// constructor
	Scene(Polyhedron* mesh_in)
	{
		mesh = mesh_in;
		scene_center.set(mesh->center.x, mesh->center.y, mesh->center.z);
		scene_radius = mesh->radius;
		scene_orientation = mesh->orientation;
	}

	Scene()
	{
		mesh = new Polyhedron;
		scene_center.set(0.0, 0.0, 0.0);
		scene_radius = 0.0;
		scene_orientation = 0;
	}

	~Scene()
	{
		if (mesh != nullptr)
		{
			mesh->finalize();
		}
		delete(mesh);
	}

	void setNeedsUpdate()
	{
		needs_update = true;
	}

	void draw(int display_mode, GLenum mode)
	{
		glClearColor(0.2, 0.2, 0.2, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// check for mesh
		if (mesh == nullptr) { return; }

		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glShadeModel(GL_SMOOTH);

		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.67, 0.67, 0.67, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		glColorMaterial(GL_FRONT, GL_DIFFUSE);
		glEnable(GL_COLOR_MATERIAL);

		switch (display_mode)
		{
		case DISPLAY_SOLID:	// color control mode
		{
			Face* face;
			Vertex* vert;
			for (int i = 0; i < mesh->nfaces; i++)
			{
				if (mode == GL_SELECT) { glLoadName(i + 1); }

				face = mesh->flist[i];
				glBegin(GL_POLYGON);
				for (int j = 0; j < face->nverts; j++)
				{
					vert = face->verts[j];
					glNormal3d(vert->normal.x, vert->normal.y, vert->normal.z);
					if (i == mesh->selected) { glColor3f(1.0, 0.0, 0.0); }
					else { glColor3f(color[0], color[1], color[2]); }

					glVertex3d(vert->x, vert->y, vert->z);
				}
				glEnd();
			}
		}
		break;

		case DISPLAY_WIREFRAME: // wireframe mode
		{
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			Face* face;
			Vertex* vert;
			for (int i = 0; i < mesh->nfaces; i++)
			{
				face = mesh->flist[i];
				glBegin(GL_POLYGON);
				for (int j = 0; j < face->nverts; j++)
				{
					vert = face->verts[j];
					glNormal3d(vert->normal.x, vert->normal.y, vert->normal.z);
					glColor3f(color[0], color[1], color[2]);
					glVertex3d(vert->x, vert->y, vert->z);
				}
				glEnd();
			}
		}
		break;

		case DISPLAY_CHECKER:	// checkerboard color pattern
		{
			Face* face;
			Vertex* vert;
			for (int i = 0; i < mesh->nfaces; i++)
			{
				face = mesh->flist[i];
				glBegin(GL_POLYGON);
				for (int j = 0; j < face->nverts; j++)
				{
					vert = face->verts[j];
					glNormal3d(vert->normal.x, vert->normal.y, vert->normal.z);
					int r = !(qFloor(vert->x / L) % 2);
					int g = !(qFloor(vert->y / L) % 2);
					int b = !(qFloor(vert->z / L) % 2);
					glColor3f(r, g, b);
					glVertex3d(vert->x, vert->y, vert->z);
				}
				glEnd();
			}
		}
		break;

		case DISPLAY_TEXTURE:
		{
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, texID);

			Face* face;
			Vertex* vert;
			for (int i = 0; i < mesh->nfaces; i++)
			{
				face = mesh->flist[i];
				glBegin(GL_POLYGON);
				for (int j = 0; j < face->nverts; j++)
				{
					vert = face->verts[j];
					glNormal3d(vert->normal.x, vert->normal.y, vert->normal.z);
					glTexCoord2d(vert->texcoord[0], vert->texcoord[1]);
					glVertex3d(vert->x, vert->y, vert->z);
				}
				glEnd();
			}
		}
		break;

		// mark as updated
		needs_update = false;
		}
	}

	// initialize a texture to be applied to the model
	void initTexture()
	{
		// initialize char array for storing texture locally and other variables
		unsigned char* texBuffer;
		int texWidth, texHeight, texBpp;

		// load in image usig stb library
		stbi_set_flip_vertically_on_load(1);
		texBuffer = stbi_load("../12.jpg", &texWidth, &texHeight, &texBpp, 4);

		// generate texture
		glEnable(GL_TEXTURE_2D);
		glGenTextures(1, &texID);
		glBindTexture(GL_TEXTURE_2D, texID);

		// texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

		// attach image data to texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, texBuffer);

		if (texBuffer)
		{
			stbi_image_free(texBuffer);
		}
	}

	private:

		// draw a sphere
		// x, y, z are the coordiate of the dot
		// radius of the sphere
		// R: the red channel of the color, ranges [0, 1]
		// G: the green channel of the color, ranges [0, 1]
		// B: the blue channel of the color, ranges [0, 1]
		void drawDot(double x, double y, double z, double radius = 0.1)
		{
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0, 1.0);
			glEnable(GL_DEPTH_TEST);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glShadeModel(GL_SMOOTH);

			GLUquadric* quad = gluNewQuadric();

			glPushMatrix();
			glTranslatef(x, y, z);
			gluSphere(quad, radius, 50, 50);
			glPopMatrix();

			gluDeleteQuadric(quad);
			glDisable(GL_POLYGON_OFFSET_FILL);
		}

		// draw a line segment
		// width: the width of the line, should bigger than 0
		// R: the red channel of the color, ranges [0, 1]
		// G: the green channel of the color, ranges [0, 1]
		// B: the blue channel of the color, ranges [0, 1]
		void drawLineSegment(LineSegment ls, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {

			glDisable(GL_LIGHTING);
			glEnable(GL_LINE_SMOOTH);
			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glLineWidth(width);

			glBegin(GL_LINES);
			glColor3f(R, G, B);
			glVertex3f(ls.start.x, ls.start.y, ls.start.z);
			glVertex3f(ls.end.x, ls.end.y, ls.end.z);
			glEnd();

			glDisable(GL_BLEND);
		}

		// convert scalar to rainbow color value
		void rainbowColor(double val, double min, double max)
		{
			double scaled_val = (val - min) / (max - min);

			if (scaled_val < 0.25)
			{
				color[0] = 0.0;
				color[1] = scaled_val / 0.25;
				color[2] = 1.0;
			}
			else if (0.25 <= scaled_val && scaled_val < 0.5)
			{
				color[0] = 0.0;
				color[1] = 1.0;
				color[2] = (scaled_val - 0.5) / -0.25;
			}
			else if (0.5 <= scaled_val && scaled_val < 0.75)
			{
				color[0] = (scaled_val - 0.5) / 0.25;
				color[1] = 1.0;
				color[2] = 0.0;
			}
			else if (0.75 <= scaled_val)
			{
				color[0] = 1.0;
				color[1] = (scaled_val - 1.0) / -0.25;
				color[2] = 0.0;
			}
			else
			{
				color[0] = 0.0;
				color[1] = 0.0;
				color[2] = 0.0;
			}
		}
};