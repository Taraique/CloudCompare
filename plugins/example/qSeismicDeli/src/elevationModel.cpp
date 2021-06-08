#include "elevationmodel.h"

//Qt
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
#include <QMainWindow>
#include <QSettings>
#include <QCloseEvent>


//qCC_fbo
#include <ccGlFilter.h>
#include <ccGLWidget.h>
#include <ccPointCloud.h>
#include <ccScalarField.h>
#include <ScalarField.h>

//à remplacer

#include <ccScalarField.h>
#include <ScalarField.h>
#include "Jacobi.h"

//System
#include <assert.h>


#include "ccStdPluginInterface.h"
#include "ccMainAppInterface.h"

elevationModel::elevationModel(float hU, float lU, float hV, float lV, ccMainAppInterface* app) :
	QDialog(app ? app->getMainWindow() : nullptr, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint),
	Ui::elevationModel(),
	m_glWindow(nullptr),
	m_app(app)
{

	setupUi(this);

	cloud = static_cast<ccPointCloud*>(m_app->getSelectedEntities().front())->cloneThis();
	cloud->showColors(true);

	QWidget* glWidget = nullptr;
	
	m_app->createGLWindow(m_glWindow, glWidget);
	assert(m_glWindow && glWidget);
	
	frame->setLayout(new QHBoxLayout());
	frame->layout()->addWidget(glWidget);
	m_glWindow->setPerspectiveState(false, true);
	m_glWindow->displayOverlayEntities(true);
	m_glWindow->setInteractionMode(ccGLWindow::MODE_TRANSFORM_CAMERA);
	m_glWindow->setPickingMode(ccGLWindow::NO_PICKING);
	m_glWindow->setViewportParameters(cloud->getDisplay()->getViewportParameters());

	m_glWindow->addToOwnDB(cloud);
	m_glWindow->redraw();
	
	highestU = hU;
	lowestU = lU;
	highestV = hV;
	lowestV = lV;
	lengthU = hU - lU;
	lengthV = hV - lV;

	m_app->dispToConsole(QString::number(hU) + "   " + QString::number(lU) + "   " + QString::number(lengthU));


	//PENSER A FAIRE LE DESTRUCTEUR
}


elevationModel::~elevationModel()
{
	if (m_glWindow)
	{
		m_glWindow->getOwnDB()->removeAllChildren();
		if (m_app)
		{
			m_app->destroyGLWindow(m_glWindow);
			m_glWindow = nullptr;
		}
	}
}

void elevationModel::on_cmBOX_valueChanged(double arg1) {
	pxBOX->setValue(arg1 / 2);
}

void elevationModel::on_pxBOX_valueChanged(double arg1) {

	double step;

	step = arg1 / 10;

	generateGrid(step);
	//cmBOX->setValue(truc * 2);
}

void elevationModel::on_buttonBox_accepted()
{
	pixelValue = pxBOX->value();
	if (!FilePath->text().isEmpty())
		exampleValue = 1;

	close();
}

void elevationModel::on_buttonBox_rejected()
{
	
	close();
}

void elevationModel::getUVH() {

}

void elevationModel::generateGrid(double step) {

	delete[] grid;



	//FAIRE UNE FONCTION !!!!
	if (lowestU < 0 && highestU < 0) {
		float lowestUSubstitute = lowestU;
		lowestU = abs(highestU);
		highestU = abs(lowestUSubstitute);
	}
	if (lowestV < 0 && highestV < 0) {
		float lowestVSubstitute = lowestV;
		lowestV = abs(highestV);
		highestV = abs(lowestVSubstitute);
	}

	//faire un assert ici par rapport à u & v
	CCCoreLib::ScalarField* u = cloud->getScalarField(cloud->getScalarFieldIndexByName("u"));
	CCCoreLib::ScalarField* v = cloud->getScalarField(cloud->getScalarFieldIndexByName("v"));

	int nbOfPoints = cloud->size();

	float uStep = lengthU / step;
	float vStep = lengthV / step;

	/*if (pxBOX->value() == 0) {
		m_glWindow->displayNewMessage("Put a value in px or cm", ccGLWindow::UPPER_CENTER_MESSAGE, false, 3600, ccGLWindow::CUSTOM_MESSAGE);
		m_glWindow->redraw();
		return;
	}*/

	int nbOfColumns = lengthV / step + 1 ;
	int nbOfLines = lengthU / step + 1 ;
	int numberOfBoxes = nbOfLines * nbOfColumns;

	m_app->dispToConsole(QString::number(nbOfColumns) + "  " + QString::number(nbOfLines) + "   " + QString::number(numberOfBoxes) + "  " + QString::number(lengthU) + "  " + QString::number(lengthV));

	grid = new std::vector<const CCVector3*>[numberOfBoxes];

	float uCoord = lowestU;
	float vCoord = lowestV;

	int currentPointLine;
	int currentPointColumn;

	int pos;


	for (int i = 0; i < nbOfPoints; i++) {

		currentPointLine = (abs(u->getValue(i)) - lowestU)/step;
		currentPointColumn = ( abs(v->getValue(i)) - lowestV)/step;
		pos = currentPointLine * nbOfLines + currentPointColumn;


		

		
		if (currentPointColumn % 2 == 0 && currentPointLine % 2 == 0 || currentPointColumn % 2 != 0 && currentPointLine % 2 != 0)
			cloud->setPointColor(i, ccColor::orangeRGB);
		else
			cloud->setPointColor(i, ccColor::blueRGB);

		//grid[currentPointLine * nbOfLines + currentPointColumn].push_back(cloud->getPoint(i));
	}

	m_app->dispToConsole("BBBBBBBB " + QString::number(abs(u->getValue(30))) + "  " + QString::number(abs(u->getValue(30)) - lowestU) + "  " + QString::number(step));

	m_glWindow->redraw();
}


void elevationModel::on_browseToolButton_clicked() {
	generateGrid(0.1);
}