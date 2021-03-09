#pragma once

#include <QtWidgets/QMainWindow>
#include <qmenu.h>
#include <qmenubar.h>
#include "window.hpp"
#include <qstring.h>
#include <iostream>
#include "dialogs.hpp"
#include "messagestream.hpp"

class MainWindow : public QMainWindow
{
	Window* window;
	FILE* open_file;

public:
	MainWindow(QWidget* parent = Q_NULLPTR)
		: QMainWindow(parent)
	{
		// setting up menu bar
		QMenuBar* menuBar = new QMenuBar;

		// File tab
		QMenu* fileMenu = menuBar->addMenu(tr("&File"));

		QAction* fileOpen = new QAction(fileMenu);
		fileOpen->setText(tr("Open"));
		fileMenu->addAction(fileOpen);
		connect(fileOpen, &QAction::triggered, this, &MainWindow::onFileOpen);

		QAction* fileClose = new QAction(fileMenu);
		fileClose->setText(tr("Close"));
		fileMenu->addAction(fileClose);
		connect(fileClose, &QAction::triggered, this, &MainWindow::onFileClose);

		QAction* fileExit = new QAction(fileMenu);
		fileExit->setText(tr("Exit"));
		fileMenu->addAction(fileExit);
		connect(fileExit, &QAction::triggered, this, &MainWindow::onFileExit);

		// view tab
		QMenu* viewMenu = menuBar->addMenu(tr("&View"));

		QAction* viewPerspective = new QAction(viewMenu);
		viewPerspective->setText("Perspective");
		viewMenu->addAction(viewPerspective);
		connect(viewPerspective, &QAction::triggered, this, &MainWindow::onViewPerspective);

		QAction* viewOrtho = new QAction(viewMenu);
		viewOrtho->setText("Orthographic");
		viewMenu->addAction(viewOrtho);
		connect(viewOrtho, &QAction::triggered, this, &MainWindow::onViewOrtho);

		// display tab
		QMenu* displayMenu = menuBar->addMenu(tr("&Display"));

		QAction* displaySolid = new QAction(displayMenu);
		displaySolid->setText(tr("Solid Color"));
		displayMenu->addAction(displaySolid);
		connect(displaySolid, &QAction::triggered, this, &MainWindow::onDisplaySolid);

		QAction* displayWireframe = new QAction(displayMenu);
		displayWireframe->setText(tr("Wireframe"));
		displayMenu->addAction(displayWireframe);
		connect(displayWireframe, &QAction::triggered, this, &MainWindow::onDisplayWireframe);

		QAction* displayChecker = new QAction(displayMenu);
		displayChecker->setText(tr("Checkerboard"));
		displayMenu->addAction(displayChecker);
		connect(displayChecker, &QAction::triggered, this, &MainWindow::onDisplayChecker);

		QAction* displayTexture = new QAction(displayMenu);
		displayTexture->setText(tr("Texture"));
		displayMenu->addAction(displayTexture);
		connect(displayTexture, &QAction::triggered, this, &MainWindow::onDisplayTexture);

		// mesh tab
		QMenu* meshMenu = menuBar->addMenu(tr("&Mesh"));

		QAction* meshSmooth = new QAction(meshMenu);
		meshSmooth->setText(tr("Smooth"));
		meshMenu->addAction(meshSmooth);
		connect(meshSmooth, &QAction::triggered, this, &MainWindow::onMeshSmooth);

		QAction* meshPairSimp = new QAction(meshMenu);
		meshPairSimp->setText(tr("Pair Simplify"));
		meshMenu->addAction(meshPairSimp);
		connect(meshPairSimp, &QAction::triggered, this, &MainWindow::onMeshPairSimp);

		setMenuBar(menuBar);
	}

	void setWindow(Window* w)
	{
		window = w;
		setCentralWidget(window);
	}

private slots:

	void onFileOpen()
	{
		QString fileName = QFileDialog::getOpenFileName(this);
		if (!fileName.isEmpty())
		{
			QByteArray strArray = fileName.toLocal8Bit();
			const char* fileStr = strArray.constData();
			open_file = fopen(fileStr, "r");

			std::cout << endl << "Opening " << fileStr << std::endl;

			Polyhedron* mesh = new Polyhedron(open_file);
			mesh->initialize();
			Scene* scene = new Scene(mesh);

			window->getGLWidget()->deleteScene();
			window->getGLWidget()->setScene(scene);
			window->resendColorData();
		}
	}

	void onFileClose()
	{
		window->getGLWidget()->deleteScene();
		window->getGLWidget()->setScene(new Scene);
	}

	void onFileExit()
	{
		window->getGLWidget()->deleteScene();
		this->close();
	}

	void onViewPerspective()
	{
		window->getGLWidget()->setViewMode(VIEW_PERSPECTIVE);
	}

	void onViewOrtho()
	{
		window->getGLWidget()->setViewMode(VIEW_ORTHOGRAPHIC);
	}

	void onDisplaySolid()
	{
		window->getGLWidget()->setSelectMode(SELECT_INDIVIDUAL);
		window->getGLWidget()->setDisplayMode(DISPLAY_SOLID);
	}

	void onDisplayWireframe()
	{
		window->getGLWidget()->setSelectMode(SELECT_NONE);
		window->getGLWidget()->setDisplayMode(DISPLAY_WIREFRAME);
	}

	void onDisplayChecker()
	{
		window->getGLWidget()->setSelectMode(SELECT_NONE);
		window->getGLWidget()->setDisplayMode(DISPLAY_CHECKER);
	}

	void onDisplayTexture()
	{
		window->getGLWidget()->setSelectMode(SELECT_NONE);
		window->getGLWidget()->setDisplayMode(DISPLAY_TEXTURE);
	}

	void onMeshSmooth()
	{
		SmoothDialog smoothDialog("Smooth Mesh", this);
		smoothDialog.setModal(true);
		int result = smoothDialog.exec();
		if (result == QDialog::Accepted)
		{
			// grab values from dialog
			int weightScheme = smoothDialog.weightSelector->currentIndex();
			QString weightText = smoothDialog.weightSelector->currentText();
			double stepSize = smoothDialog.stepSelector->value();
			int iterations = smoothDialog.iterSelector->value();
			
			// Smoothing initiated message
			QString message = QString("Initiated mesh smoothing with %1 weighting, %2 time step, %3 iterations:")
				.arg(weightText, QString::number(stepSize), QString::number(iterations));
			std::cout << std::endl << message.toStdString().c_str() << std::endl;

			// execute smoothing and time execution
			QTime timer;
			timer.start();
			window->getGLWidget()->getScene()->mesh->smooth((weightScheme + 1), stepSize, iterations);
			window->getGLWidget()->getScene()->setNeedsUpdate();
			int time_elapsed = timer.elapsed();

			// smoothing finished message
			message = QString("\tMesh smoothing complete, %1ms elapsed").arg(QString::number(time_elapsed));
			std::cout << message.toStdString().c_str() << std::endl;
		}
	}

	void onMeshPairSimp()
	{
		PairSimpDialog pairDialog(this);
		pairDialog.setModal(true);
		int result = pairDialog.exec();
		if (result == QDialog::Accepted)
		{
			// grab values from dialog
			int faceTarget = pairDialog.targetSelector->value();
			double errorTolerance = pairDialog.toleranceSelector->value();
			int maxContractions = pairDialog.maxSelector->value();

			// Simplification initiated message
			QString message = QString("Initiated pair contraction mesh simplification with %1 target face count, %2 geometric error tolerance, %3 max contractions:")
				.arg(QString::number(faceTarget), QString::number(errorTolerance), QString::number(maxContractions));
			printLine(message.toStdString(),true);

			// execute simplification and time execution
			QTime timer;
			timer.start();
			window->getGLWidget()->getScene()->mesh->pairSimplify(faceTarget, errorTolerance, maxContractions);
			window->getGLWidget()->getScene()->setNeedsUpdate();
			int time_elapsed = timer.elapsed();

			// Simplification finished message
			message = QString("\tPair contraction simplification complete, %1 elapsed").arg(QString::number(time_elapsed));
			printLine(message.toStdString());
		}
	}

};
