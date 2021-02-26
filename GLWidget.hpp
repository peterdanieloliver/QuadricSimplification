#pragma once

#include <qopenglwidget.h>
#include <qopenglfunctions.h>
#include <qsurfaceformat.h>
#include <QtOpenGL>
#include <GL/GLU.h>
#include <qcolor.h>
#include "Scene.hpp"
#include "trackball.hpp"

#define MOUSE_NONE              1
#define MOUSE_TRANSLATE         2
#define MOUSE_ROTATE            3
#define MOUSE_SELECT            4

#define SELECT_NONE             11
#define SELECT_INDIVIDUAL       12

#define VIEW_PERSPECTIVE        21
#define VIEW_ORTHOGRAPHIC       22

class GLWidget : public QOpenGLWidget, public QOpenGLFunctions
{

private:
	Scene* scene;
	Trackball* trackball;

	// mouse control variables
	QVector2D pnew, pold;
	int mouse_mode = MOUSE_NONE;
	QMatrix4x4 rotmat;
	QVector2D transvect;
	float aspectratio;
	float zoom = 1.0;
	float zoom_speed = 0.9;

    // display control variable
    int display_mode = DISPLAY_SOLID;

    // selection control variable
    int selection_mode = SELECT_INDIVIDUAL;

    // view control variable
    int view_mode = VIEW_PERSPECTIVE;

public: // public methods

	// constructor
	GLWidget(QWidget* parent = Q_NULLPTR)
		: QOpenGLWidget(parent)
	{
        setScene(new Scene);
	}

	void setScene(Scene* scene_in)
	{
		scene = scene_in;
		trackball = new Trackball(0.8);
        this->update();
	}

    Scene* getScene()
    {
        return scene;
    }
    
    void deleteScene()
    {
        delete(scene);
        delete(trackball);
    }

    void setDisplayMode(int mode)
    {
        display_mode = mode;
        this->update();
    }

    void setSelectMode(int mode)
    {
        selection_mode = mode;
    }

    void setViewMode(int mode)
    {
        view_mode = mode;
        this->update();
    }

	QSize sizeHint() const override
	{
		return QSize(600, 600);
	}

	QSize minimumSizeHint() const override
	{
		return QSize(400, 400);
	}

public slots:

    void setColor(QColor color)
    {
        qreal r, g, b, a;
        color.getRgbF(&r, &g, &b, &a);
        scene->color[0] = (double)r;
        scene->color[1] = (double)g;
        scene->color[2] = (double)b;
        scene->color[3] = (double)a;
        this->update();
    }

private: // private helper functions

    void setView(GLenum mode)
    {
        // do some lighting stuff
        glMatrixMode(GL_PROJECTION);
        if (mode == GL_RENDER)
        {
            glLoadIdentity();
        }

        // view setup
        if (view_mode == VIEW_PERSPECTIVE)
        {
            double fov = 2 * atan(zoom / 2);
            fov = qRadiansToDegrees(fov);
            gluPerspective(fov, aspectratio, 0.1, 10000);
        }
        if (view_mode == VIEW_ORTHOGRAPHIC)
        {
            glOrtho(-zoom, zoom, -zoom, zoom, 0.1, 10000);
        }

        // Initialize model view matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // more lighting stuff
    }

    void setTransforms(GLenum mode) // note to self: transforms are applied in reverse order
    {
        // input translations
        glTranslatef(transvect.x(), transvect.y(), -2 * scene->scene_radius);

        // input rotations
        int index = 0;
        GLfloat mat[16];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                mat[index++] = rotmat(i, j);
            }
        }
        glMultMatrixf(mat);

        // initial positioning and rotation
        glRotatef(45, 1.0, 0.0, 0.0);
        glRotatef(-45, 0.0, 1.0, 0.0);

        glTranslatef(-scene->scene_center.x,
            -scene->scene_center.y,
            -scene->scene_center.z);
    }

    // used for trackball inputs
    QVector2D scaleWinCoords(QVector2D point)
    {
        QVector2D scaled_point;
        scaled_point.setX((2.0 * point.x() / width()) - 1.0);
        scaled_point.setY((2.0 * (height() - point.y()) / height()) - 1.0);
        return scaled_point;
    }

    // used for selecting a vertex
    void selectFace(int x, int y)
    {
        if (selection_mode == SELECT_NONE)
        {
            return;
        }

        GLuint selectBuf[1024];
        glSelectBuffer(1024, selectBuf);
        GLint hits;
        GLint viewport[4];

        glGetIntegerv(GL_VIEWPORT, viewport);
        glRenderMode(GL_SELECT);

        glInitNames();
        glPushName(0);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();

        // create a 5x5 pixel selection region
        gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1, 1, viewport);

        setView(GL_SELECT);
        setTransforms(GL_SELECT);
        scene->draw(display_mode, GL_SELECT);

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glFlush();

        glMatrixMode(GL_MODELVIEW);

        hits = glRenderMode(GL_RENDER);

        int selected = processHits(hits, selectBuf);
        if (selected != -1)
        {
            switch (selection_mode)
            {
        
                case SELECT_INDIVIDUAL:
                {
                    scene->mesh->selected = selected;
                    std::cout << "Selected face #" << std::to_string(selected) << std::endl;
                }
                break;

            }
        }
        
        this->update();
    }

    // process the hits from picking matrix
    int processHits(GLint hits, GLuint buffer[])
    {
        GLuint names;
        GLuint* ptr;
        double shallowest = 1.0e+20;
        double current;
        int seed_id = -1;
        bool needs_update = false;

        ptr = (GLuint*)buffer;
        for (int i = 0; i < hits; i++)
        {
            needs_update = false;
            names = *ptr;
            ptr++;

            current = (double)*ptr / 0x7fffffff;
            if (current < shallowest)
            {
                shallowest = current;
                needs_update = true;
            }
            ptr++;
            current = (double)*ptr / 0x7fffffff;
            if (current < shallowest)
            {
                shallowest = current;
                needs_update = true;
            }
            ptr++;
            for (int j = 0; j < names; j++)
            {
                if (needs_update)
                {
                    seed_id = *ptr - 1;
                }
                ptr++;
            }
        }
        return seed_id;
    }

protected: // protected methods

	void initializeGL() override
	{
		initializeOpenGLFunctions();
		glClearColor(0.2, 0.2, 0.2, 1.0);
        glShadeModel(GL_FLAT);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_NORMALIZE);
        if (scene->scene_orientation == 0) { glFrontFace(GL_CW); }
        else { glFrontFace(GL_CCW); }

        resizeGL(this->width(), this->height());
	}

	void paintGL() override
	{   
        setView(GL_RENDER);
        setTransforms(GL_RENDER);

        scene->draw(display_mode,GL_RENDER);
        glFlush();
	}

	void resizeGL(int w, int h) override
	{
        glViewport((-w/2) , (-h/2), (w/2), (h/w));
        aspectratio = (float)w / (float)h;
        setView(GL_RENDER);
	}

    void mousePressEvent(QMouseEvent* event) override
    {
        resizeGL(this->width(), this->height());

        Qt::MouseButton button = event->button();
        pnew = scaleWinCoords(QVector2D(event->pos()));
        if (button == Qt::RightButton)
        {
            if (mouse_mode == MOUSE_NONE)
            {
                mouse_mode = MOUSE_TRANSLATE; // translation mode
            }
        }
        else if (button == Qt::MiddleButton)
        {
            if (mouse_mode == MOUSE_NONE)
            {
                mouse_mode = MOUSE_ROTATE; // rotation mode
            }
        }
        else if (button == Qt::LeftButton)
        {
            if (mouse_mode == MOUSE_NONE)
            {
                mouse_mode = MOUSE_SELECT; // selection mode
                selectFace(event->pos().x(),event->pos().y());
            }
        }
        else
        {
            mouse_mode = 0; // no action mode
        }
    }

    void mouseMoveEvent(QMouseEvent* event) override
    {
        resizeGL(this->width(), this->height());

        if (mouse_mode == MOUSE_TRANSLATE || mouse_mode == MOUSE_ROTATE)
        {
            pold = pnew;
            pnew = scaleWinCoords(QVector2D(event->pos()));
            if (mouse_mode == MOUSE_TRANSLATE)
            {
                transvect += QVector2D(2 * (pnew - pold));
                this->update();
            }
            else if (mouse_mode == MOUSE_ROTATE)  // rotate
            {
                QQuaternion rotquat = trackball->getRotationQuat(pold, pnew);
                rotmat.rotate(rotquat);

                this->update();
            }
        }
    }

    void mouseReleaseEvent(QMouseEvent* event) override
    {
        resizeGL(this->width(), this->height());

        Qt::MouseButton button = event->button();
        if (button == Qt::RightButton && mouse_mode == MOUSE_TRANSLATE)
        {
            mouse_mode = MOUSE_NONE;
        }
        else if (button == Qt::MiddleButton && mouse_mode == MOUSE_ROTATE)
        {
            mouse_mode = MOUSE_NONE;
        }
        else if (button == Qt::LeftButton && mouse_mode == MOUSE_SELECT)
        {
            mouse_mode = MOUSE_NONE;
        }
    }

    void wheelEvent(QWheelEvent* event) override
    {
        float direction = event->angleDelta().y();
        if (direction < 0)
        {
            zoom /= zoom_speed;
            this->update();
        }
        else if (direction > 0)
        {
            zoom *= zoom_speed;
            this->update();
        }
    }
};