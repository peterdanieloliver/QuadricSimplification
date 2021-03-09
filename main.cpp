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
    //Testy** tlist = new Testy * [7];
    //for (int i = 0; i < 7; i++)
    //{
    //    tlist[i] = new Testy(i);
    //}

    //tlist[0] = nullptr;

    //Testy* ttemp;
    //for (int i = 0; i < 7; i++)
    //{
    //    ttemp = tlist[i];
    //}
    
    QApplication a(argc, argv);
    mainWindow = new MainWindow;
    window = new Window(mainWindow);
    MessageStream qout(std::cout, window->getMessageBox());
    mainWindow->setWindow(window);
    mainWindow->resize(1300, 800);
    mainWindow->show();
    return a.exec();
}
