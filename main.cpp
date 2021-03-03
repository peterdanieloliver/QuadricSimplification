#include "mainwindow.h"
#include "GLWidget.hpp"
#include "window.hpp"
#include "polyhedron.hpp"
#include "Scene.hpp"
#include "messagestream.hpp"
#include <QtWidgets/QApplication>

MainWindow* mainWindow;
Window* window;

int main(int argc, char *argv[])
{
    icMatrix4x4 test0(2.0, -3.0, 4.0, 1.0,
                     -2.0, -3.0, 0.0, 2.0,
                     -4.0, 3.0, 1.0, 2.0,
                      1.0, 0.0, 3.0, -3.0);

    icMatrix4x4 test1(1.0, -2.0, 3.0, 3.0,
                     -2.0, -4.0, 2.0, 0.0,
                     -1.0, 1.0, 1.0, 2.0,
                      3.0, 2.0, -3.0, 0.0);
    
    icMatrix4x4 test2;

    icVector4 vect0(1.0, -2.0, 3.0, -4.0);
    icVector4 vect2;

    test0.leftMultiply(test1);
    vect2 = test1 * vect0;
    vect2 = vect0 * test1;
    test2 = square(vect0);

    QApplication a(argc, argv);
    mainWindow = new MainWindow;
    window = new Window(mainWindow);
    MessageStream qout(std::cout, window->getMessageBox());
    mainWindow->setWindow(window);
    mainWindow->resize(1300, 800);
    mainWindow->show();
    return a.exec();
}
