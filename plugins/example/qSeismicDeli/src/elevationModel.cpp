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

	pxBOX->setSingleStep(0.01);
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
	step = arg1;

	generateGrid();
	//cmBOX->setValue(truc * 2);
}

void elevationModel::on_Generate_clicked()
{
	Variogram();
}

void elevationModel::on_Close_clicked()
{
	close();
}

void elevationModel::generateGrid() {

	if (step == 0)
		return;

	delete[] grid;

	//faire un assert ici par rapport à u & v
	CCCoreLib::ScalarField* u = cloud->getScalarField(cloud->getScalarFieldIndexByName("u"));
	CCCoreLib::ScalarField* v = cloud->getScalarField(cloud->getScalarFieldIndexByName("v"));

	int nbOfPoints = cloud->size();

	float uStep = lengthU / step;
	float vStep = lengthV / step;

	nbOfColumns = lengthV / step + 1 ;
	nbOfLines = lengthU / step + 1 ;
	numberOfBoxes = nbOfLines * nbOfColumns;

	m_app->dispToConsole( "   " + QString::number(numberOfBoxes));

	grid = new std::vector<int>[numberOfBoxes];

	float uCoord = lowestU;
	float vCoord = lowestV;

	int currentPointLine;
	int currentPointColumn;

	int pos;

	int l = 0;
	int ll = 0;

	test = nbOfColumns + 80;

	for (int i = 0; i < nbOfPoints; i++) {

		currentPointLine =  trunc((u->getValue(i) - lowestU) / step);
		currentPointColumn =  trunc((v->getValue(i) - lowestV)  / step);

		pos = currentPointLine * nbOfColumns + currentPointColumn;

		if (currentPointColumn % 2 == 0 && currentPointLine % 2 == 0 || currentPointColumn % 2 != 0 && currentPointLine % 2 != 0)
			cloud->setPointColor(i, ccColor::orangeRGB);
		else
			cloud->setPointColor(i, ccColor::blueRGB);

		grid[pos].push_back(i);
	}

	m_glWindow->redraw();
}

void elevationModel::on_browseToolButton_clicked() {
	QString filePath =
		QFileDialog::getSaveFileName(this, "Where to save the results", "results",
			"CSV(*.csv);;Text files (*.txt)");
	FilePath->setText(filePath);
}

float elevationModel::findT(std::vector<int> box) {
	int nbOfPointsInBox = box.size();

	float T = 0;
	float currentMin;
	float tempDistance;
	const CCVector3* currentPoint;
	for (int i = 0; i < nbOfPointsInBox; i++) {
		currentMin = INFINITY;
		currentPoint = (cloud->getPoint(box.at(i)));
		for (int j = 0; j < nbOfPointsInBox; j++) {
			tempDistance = distance(currentPoint, cloud->getPoint(box.at(j)));
			if (currentMin > tempDistance && tempDistance != 0)
				currentMin = tempDistance;
		}
		T += currentMin;
	}

	T /= nbOfPointsInBox;
	return T;
}

double elevationModel::distance(const CCVector3* point1, const CCVector3* point2) {
	return sqrt(
		pow(point1->x - point2->x, 2) +
		pow(point1->y - point2->y, 2) +
		pow(point1->z - point2->z, 2)
	);
}

void elevationModel::Variogram() {
	int a;
	CCCoreLib::ScalarField* h = cloud->getScalarField(cloud->getScalarFieldIndexByName("h"));

	float T = findT(grid[test]);
	int nbOfT = round((step * 2) / (4 * T));


	m_app->dispToConsole(QString::number(nbOfT) + "  " + QString::number(T) + "  ");

	int nbOfCells = numberOfBoxes * nbOfT;
	
	auto couplesDirectory = new std::vector<coupleOfPoints>[nbOfCells];
	double* gammaResults = new double[nbOfCells];

	coupleOfPoints currentCouple;
	int TPosition;
	

	int nbOfPointsInCell;

	for (int currentCell = 0; currentCell < numberOfBoxes; currentCell++) {
		nbOfPointsInCell = grid[currentCell].size();

		for (int i = 0; i < nbOfPointsInCell; i++) {
			currentCouple.firstSpouse = grid[currentCell].at(i);
			for (int j = i; j < nbOfPointsInCell; j++) {
				currentCouple.secondSpouse = grid[currentCell].at(j);
				TPosition = round(distance(cloud->getPoint(currentCouple.firstSpouse), cloud->getPoint(currentCouple.secondSpouse)) / T);
				if (TPosition <= nbOfT && TPosition > 0) {																						//probablement inutile
					TPosition -=  1;
					couplesDirectory[currentCell * nbOfT + TPosition].push_back(currentCouple);
					a++;
				}
			}
		}
	}
	m_app->dispToConsole(QString::number(a));
	delete[] grid;

	float currentGammaResult;
	float machin;
	int nbOfCouplesInCell;
	int TMultiplier = 1;

	for (int currentCell = 0; currentCell < nbOfCells; currentCell++)
	{
		nbOfCouplesInCell = couplesDirectory[currentCell].size();
		currentGammaResult = 0;

		for (int i = 0; i < nbOfCouplesInCell; i++)
		{
			currentCouple = couplesDirectory[currentCell].at(i);
			currentGammaResult += pow(
				h->getValue(currentCouple.firstSpouse) - (h->getValue(currentCouple.secondSpouse) + T*TMultiplier)
				, 2);
		}
		if (nbOfCouplesInCell == 0)
			currentGammaResult = -1;
		else {

			currentGammaResult *= (1. / (2. * nbOfCouplesInCell));
		}

		if (currentCell == 3 || currentCell == 12 || currentCell == 3000 || currentCell == 800) {
			m_app->dispToConsole(QString::number(nbOfCouplesInCell) + "   "  + QString::number(1./(2.* nbOfCouplesInCell))
				+ "   " + QString::number(currentGammaResult));

		}

		TMultiplier++;
		gammaResults[currentCell] = currentGammaResult;
		if (TMultiplier > nbOfT)
			TMultiplier = 1;
	}

	//Writing the result

	QFile file(FilePath->text()); 
	if (file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		// We're going to streaming text to the file
		QTextStream stream(&file);
		stream << "T;  ";

		for (size_t i = 0; i < nbOfLines;i++)
		{
			for (size_t j = 0; j < nbOfColumns; j = j + 1)
			{
				stream << "C" << i<<"_"<< j << ";  ";
			}
		}
		stream << "\n";
		
		for (size_t currentT = 0; currentT < nbOfT; currentT++)
		{
			stream << T << currentT <<"= " << T * (currentT + 1) << ";  ";

			for (size_t i = currentT; i < nbOfCells; i = i + nbOfT)
			{
				stream << gammaResults[i] << ";  ";
			}
			stream << ";\n";
		}

		file.close();
		qDebug() << "Writing finished";
	}




	delete[] couplesDirectory;
	delete[] gammaResults;
} 