#pragma once

#include <qdialog.h>

class SmoothDialog : public QDialog
{
public:

	QComboBox* weightSelector;
	QDoubleSpinBox* stepSelector;
	QSpinBox* iterSelector;

public:

	explicit SmoothDialog(QString title, QWidget* parent = nullptr)
	{
		QLabel* weightLabel = new QLabel("Weight Sceme:");
		QLabel* stepLabel = new QLabel("Step Size:");
		QLabel* iterLabel = new QLabel("Iterations:");

		weightSelector = new QComboBox;
		weightSelector->addItem("uniform");
		weightSelector->addItem("cord");
		weightSelector->addItem("mean curvature flow");
		weightSelector->addItem("mean values coordinates");
		weightSelector->setCurrentIndex(0);

		iterSelector = new QSpinBox;
		iterSelector->setRange(1, 100);
		iterSelector->setSingleStep(1);
		iterSelector->setValue(1);

		stepSelector = new QDoubleSpinBox;
		stepSelector->setRange(0.01, 10);
		stepSelector->setSingleStep(0.1);
		stepSelector->setValue(0.1);

		QPushButton* okButton = new QPushButton("&OK",this);
		connect(okButton, &QPushButton::clicked, this, &SmoothDialog::accept);
		QPushButton* cancelButton = new QPushButton("&Cancel",this);
		connect(cancelButton, &QPushButton::clicked, this, &SmoothDialog::reject);
		
		QGridLayout* layout = new QGridLayout;
		layout->addWidget(weightLabel, 0, 0);
		layout->addWidget(weightSelector, 0, 1);
		layout->addWidget(stepLabel, 1, 0);
		layout->addWidget(stepSelector, 1, 1);
		layout->addWidget(iterLabel, 2, 0);
		layout->addWidget(iterSelector, 2, 1);
		layout->addWidget(okButton, 3, 0);
		layout->addWidget(cancelButton, 3, 1);

		setLayout(layout);

		setWindowTitle(title);
		Qt::WindowFlags flags = windowFlags();
		Qt::WindowFlags helpFlag = Qt::WindowContextHelpButtonHint;
		flags = flags & (~helpFlag);
		setWindowFlags(flags);

	}
};