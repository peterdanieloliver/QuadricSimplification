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
    QApplication a(argc, argv);
    mainWindow = new MainWindow;
    window = new Window(mainWindow);
    MessageStream qout(std::cout, window->getMessageBox());
    mainWindow->setWindow(window);
    mainWindow->resize(1300, 800);
    mainWindow->show();
    return a.exec();
}
