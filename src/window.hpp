#pragma once

#include <qwidget.h>
#include <qslider.h>
#include <qboxlayout.h>
#include <qdesktopwidget.h>
#include <qapplication.h>
#include <qtextedit.h>
#include "GLWidget.hpp"
#include <iostream>
#include <qtoolbox.h>
#include <qcolordialog.h>

class MainWindow;

class Window : public QWidget
{
private:

	MainWindow* mainWindow;
	GLWidget* glWidget;
	QTextEdit* messageBox;
	QColorDialog* colorDialog;
	

public:

	Window(MainWindow* mw)
		: mainWindow(mw)
	{
		// initialize glwidget
		glWidget = new GLWidget;

		// create color selection dialog
		colorDialog = new QColorDialog(this);
		colorDialog->setWindowFlags(Qt::Widget);
		colorDialog->setOptions(QColorDialog::DontUseNativeDialog | QColorDialog::NoButtons);
		connect(colorDialog, &QColorDialog::currentColorChanged, glWidget, &GLWidget::setColor);

		// create message box
		messageBox = new QTextEdit(this);

		// set up side panel toolbox
		QToolBox* sidePanel = new QToolBox;
		sidePanel->setFixedWidth(520);
		sidePanel->addItem(colorDialog, tr("Color"));
		sidePanel->addItem(messageBox, tr("Output"));
		sidePanel->setCurrentIndex(sidePanel->count()-1);
		
		// set up main grid layout
		QGridLayout* mainLayout = new QGridLayout;
		mainLayout->addWidget(glWidget,0,0);
		mainLayout->addWidget(sidePanel,0,1);
		setLayout(mainLayout);
	}

	GLWidget* getGLWidget()
	{
		return glWidget;
	}

	QTextEdit* getMessageBox()
	{
		return messageBox;
	}

	void resendColorData()
	{
		colorDialog->currentColorChanged(colorDialog->currentColor());
	}

};