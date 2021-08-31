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

elevationModel::elevationModel(float hU, float lU, float hV, float lV, const CCVector3* center, ccMainAppInterface* app) :
	QDialog(app ? app->getMainWindow() : nullptr, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint),
	Ui::elevationModel(),
	m_glWindow(nullptr),
	m_app(app)
{

	setupUi(this);


	cloud = static_cast<ccPointCloud*> (m_app->getSelectedEntities().front())->cloneThis();

	cloud->showSF(false);
	cloud->showColors(true);

	superPXSlider->setMaximum(50);


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

	centroid->x = center->x;
	centroid->y = center->y;
	centroid->z = center->z;


	highestU = hU;
	lowestU = lU;
	highestV = hV;
	lowestV = lV;

	lengthU = hU - lU;
	lengthV = hV - lV;

	nbOfPoints = cloud->size();
	grayValues = new float[nbOfPoints];
	originalColors = new ccColor::Rgba[nbOfPoints];
	gridColors = new ccColor::Rgba[nbOfPoints];

	//generates the cloud point in shades of Gray for the LBP computation
	for (size_t i = 0; i < nbOfPoints; i++)
	{
		grayValues[i] = cloud->getPointColor(i).r * 0.299 + cloud->getPointColor(i).g * 0.587 + cloud->getPointColor(i).b * 0.114;
		originalColors[i] = cloud->getPointColor(i);
	}

	pxBOX->setSingleStep(0.01);
}


elevationModel::~elevationModel()
{
	delete[] grid;
	delete[] grayValues;
	delete[] originalColors;
	delete[] gridColors;
	delete[] localBinaryPattern;
	delete[] alpha;
	delete[] gammaResults;

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

//Separates the point cloud into many subsets to reduce the program complexity
void elevationModel::generateGrid() {

	//interface update and data clearing
	Generate->setText("Generate");
	Generate2->setText("Generate");

	VariogramGenerated = false;
	LBPGenerated = false;
	SPXGenerated = false;

	superPXSlider->setValue(0);

	{
		int alphaIndex = cloud->getScalarFieldIndexByName("alpha");
		if (alphaIndex != -1)
			cloud->deleteScalarField(alphaIndex);

		int superPXsIndex = cloud->getScalarFieldIndexByName("SuperPXsField");
		if (superPXsIndex != -1)
			cloud->deleteScalarField(superPXsIndex);
	}

	showColors->setEnabled(true);
	showAlpha->setEnabled(false);
	showGrid->setEnabled(false);

	cloud->showSF(false);

	/////////////

	if (stepU == 0 || stepV == 0)
		return;

	delete[] grid;
	delete[] localBinaryPattern;

	CCCoreLib::ScalarField* u = cloud->getScalarField(cloud->getScalarFieldIndexByName("u"));
	CCCoreLib::ScalarField* v = cloud->getScalarField(cloud->getScalarFieldIndexByName("v"));


	nbOfColumns = lengthV / stepV + 1;
	nbOfLines = lengthU / stepU + 1;
	numberOfBoxes = nbOfLines * nbOfColumns;

	grid = new std::vector<int>[numberOfBoxes];
	localBinaryPattern = new int[numberOfBoxes];

	for (int i = 0; i < nbOfPoints; i++) {

		int currentPointLine = trunc((u->getValue(i) - lowestU) / stepU);
		int currentPointColumn = trunc((v->getValue(i) - lowestV) / stepV);

		int pos = currentPointLine * nbOfColumns + currentPointColumn;

		if (currentPointLine < nbOfLines && currentPointColumn < nbOfColumns) {

			if (currentPointColumn % 2 == 0 && currentPointLine % 2 == 0 || currentPointColumn % 2 != 0 && currentPointLine % 2 != 0) {
				cloud->setPointColor(i, ccColor::orange);
				gridColors[i] = ccColor::orange;
			}
			else {
				cloud->setPointColor(i, ccColor::blue);
				gridColors[i] = ccColor::blue;
			}
			grid[pos].push_back(i);
		}
	}
	m_glWindow->redraw();
}

void elevationModel::on_browseToolButton_clicked() {
	QString filePath =
		QFileDialog::getSaveFileName(this, "Where to save the results", "resultsVARIOGRAM",
			"CSV(*.csv);;Text files (*.txt)");
	FilePath->setText(filePath);

	if (!VariogramGenerated)
		if (!filePath.isEmpty())
			Generate->setText("Generate and Save");
}

void elevationModel::on_browseToolButton2_clicked() {
	QString filePath =
		QFileDialog::getSaveFileName(this, "Where to save the results", "resultsLBP",
			"CSV(*.csv);;Text files (*.txt)");
	FilePath2->setText(filePath);

	if (!LBPGenerated)
		if (!filePath.isEmpty())
			Generate2->setText("Generate and Save");
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

//Euclidean computation method
double elevationModel::distance(const CCVector3* point1, const CCVector3* point2) {
	return sqrt(
		(point1->x - point2->x) * (point1->x - point2->x) +
		(point1->y - point2->y) * (point1->y - point2->y) +
		(point1->z - point2->z) * (point1->z - point2->z)
	);
}

//Uses a Variogram to compute the alpha values of each cell
void elevationModel::Variogram() {
	if (!VariogramGenerated) {
		{
			int alphaIndex = cloud->getScalarFieldIndexByName("alpha");
			if (alphaIndex != -1)
				cloud->deleteScalarField(alphaIndex);
		}

		CCCoreLib::ScalarField* h = cloud->getScalarField(cloud->getScalarFieldIndexByName("h"));

		t = findT(grid[numberOfBoxes / 2]);

		nbOfT = round((stepU + stepV) / (4 * t));
		int nbOfCellsWithDistances = numberOfBoxes * nbOfT;

		auto couplesDirectory = new std::vector<coupleOfPoints>[nbOfCellsWithDistances];

		delete[] gammaResults;
		gammaResults = new double[nbOfCellsWithDistances];

		coupleOfPoints currentCouple;

		//searching couple of points for each t
		for (int currentCell = 0; currentCell < numberOfBoxes; currentCell++) {
			int nbOfPointsInCell = grid[currentCell].size();
			int TPosition;
			for (int i = 0; i < nbOfPointsInCell; i++) {
				currentCouple.firstSpouse = grid[currentCell].at(i);
				for (int j = i; j < nbOfPointsInCell; j++) {
					currentCouple.secondSpouse = grid[currentCell].at(j);
					TPosition = round(distance(cloud->getPoint(currentCouple.firstSpouse), cloud->getPoint(currentCouple.secondSpouse)) / t);
					if (TPosition <= nbOfT && TPosition > 0) {																						//probablement inutile
						TPosition -= 1;
						couplesDirectory[currentCell * nbOfT + TPosition].push_back(currentCouple);
					}
				}
			}
		}

		float currentGammaResult;
		int nbOfCouplesInCell;
		int TMultiplier = 1;

		//Variogram computation
		for (int currentCell = 0; currentCell < nbOfCellsWithDistances; currentCell++)
		{
			nbOfCouplesInCell = couplesDirectory[currentCell].size();
			currentGammaResult = 0;

			for (int i = 0; i < nbOfCouplesInCell; i++)
			{
				currentCouple = couplesDirectory[currentCell].at(i);
				currentGammaResult +=
					(h->getValue(currentCouple.firstSpouse) - (h->getValue(currentCouple.secondSpouse) + t * TMultiplier)) *
					(h->getValue(currentCouple.firstSpouse) - (h->getValue(currentCouple.secondSpouse) + t * TMultiplier));
			}
			if (nbOfCouplesInCell == 0)
				currentGammaResult = 0;
			else {
				currentGammaResult *= (1. / (2. * nbOfCouplesInCell));
			}

			TMultiplier++;
			gammaResults[currentCell] = currentGammaResult;
			if (TMultiplier > nbOfT)
				TMultiplier = 1;
		}

		//linear regression after logarithm transformation
		double avgTln = 0;
		alpha = new double[numberOfBoxes];

		for (size_t i = 1; i < nbOfT + 1; i++)
		{
			avgTln += log(t * i);
		}

		avgTln /= nbOfT;

		double* avgGammaTln = new double[numberOfBoxes];
		double lnGammaT;
		int position = 0;

		for (size_t i = 0; i < nbOfCellsWithDistances; i = i + nbOfT)
		{
			lnGammaT = 0;
			for (size_t j = i; j < i + nbOfT; j++)
			{
				lnGammaT += log(gammaResults[i]);
			}
			avgGammaTln[position++] = lnGammaT / nbOfT;
		}


		//fractal dimension computation
		double currentBoxAlpha;
		double div;
		ccScalarField* alphaField = new ccScalarField("alpha");

		for (size_t currentBox = 0; currentBox < numberOfBoxes; currentBox++)
		{
			currentBoxAlpha = 0;
			div = 0;
			for (size_t i = 0; i < nbOfT; i++)
			{
				currentBoxAlpha += (log(t * (i + 1)) - avgTln) * (log(gammaResults[currentBox * nbOfT + i]) - avgGammaTln[currentBox]);
				div += pow(log(t * (i + 1)) - avgTln, 2);
			}
			currentBoxAlpha /= div;
			alpha[currentBox] = currentBoxAlpha;
		}

		alphaField->resize(nbOfPoints);


		for (size_t currentBox = 0; currentBox < numberOfBoxes; currentBox++)
		{
			//When the cases are too small the alpha sometimes contains -inf as a value, we don't want that to happen as it messes with the ScalarField
			if (alpha[currentBox] == -INFINITY)
				alpha[currentBox] = NAN;

			for (size_t i = 0; i < grid[currentBox].size(); i++)
				alphaField->at(grid[currentBox].at(i)) = alpha[currentBox];

		}

		alphaField->computeMinAndMax();

		alphaField->resizeSafe(nbOfPoints, true, CCCoreLib::NAN_VALUE);


		cloud->addScalarField(alphaField);
		cloud->setCurrentDisplayedScalarField(cloud->getScalarFieldIndexByName("alpha"));
		cloud->showSF(true);
		showGrid->setEnabled(true);

		Generate->setText(QString("Save"));
		m_glWindow->redraw();
		delete[] couplesDirectory;
	}

	//Writing the result

	if (!FilePath->text().isEmpty()) {

		QFile file(FilePath->text());
		if (file.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			QTextStream stream(&file);

			stream << cloud->getName() << "\n" << "case dimensions : " << floor(stepU * 100.) << "cm²" << "\n" << "centroid : x = " << centroid->x << " y = " << centroid->y << " z = " << centroid->z << "\n Variogram\n";

			for (size_t currentLine = 0; currentLine < nbOfLines; currentLine++)
			{
				stream << "T;  ";
				for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
				{
					stream << "C " << currentLine << "_" << currentColumn << ";";
				}
				stream << "\n";

				for (size_t currentT = 0; currentT < nbOfT; currentT++) {
					stream << t * (currentT + 1) << ";  ";

					for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
					{
						stream << gammaResults[currentLine * nbOfColumns * nbOfT + nbOfT * currentColumn + currentT] << "; ";
					}
					stream << "\n";
				}
			}
			file.close();
		}
	

		//writing the average fractal dimension of each line
		if (AvgCSV->isChecked()) {
			QFile fileAVG(FilePath->text().insert(FilePath->text().size() - 4, "AVG"));
			QTextStream stream(&fileAVG);

			stream << cloud->getName() << "\n" << "case dimensions : " << floor(stepU * 100.) << "cm²" << "\n" << "centroid : x = " << centroid->x << " y = " << centroid->y << " z = " << centroid->z << "\n Variogram AVG\n";

			if (fileAVG.open(QIODevice::WriteOnly | QIODevice::Text)) {
				for (size_t currentLine = 0; currentLine < nbOfLines; currentLine++)
				{
					double avgAlpha = 0;

					for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
					{
						if (!isnan(alpha[currentLine * nbOfColumns + currentColumn]))
							avgAlpha += alpha[currentLine * nbOfColumns + currentColumn];
					}

					avgAlpha /= nbOfColumns;
					stream << "ligne " << currentLine << "; " << avgAlpha << " ;" << "\n";
				}
				fileAVG.close();
			}
		}
	}
	
	VariogramGenerated = true;
}

//converts binary numbers to decimal
int elevationModel::BinaryToDecimal(int n) {
	int decimalNumber = 0;
	int base = 1;
	int temp = n;
	for (temp; temp != 0; temp /= 10) {
		int lastDigit = temp % 10;
		decimalNumber += lastDigit * base;
		base = base * 2;
	}
	return decimalNumber;
}

//Regroups the cases into subsets using their distances with eachother and their alpha value
void elevationModel::generateSuperPixelsWithSLIC() {

	int numberOfSuperPx = superPXSlider->value();

	if (numberOfSuperPx == 0)
		return;
	{
		int superPXsIndex = cloud->getScalarFieldIndexByName("SuperPXsField");
		if (superPXsIndex != -1)
			cloud->deleteScalarField(superPXsIndex);
	}

	float S = numberOfBoxes / numberOfSuperPx;
	superPXs.clear();
	superPXs.resize(numberOfSuperPx);
	double* superPXsAVG = new double[numberOfSuperPx];
	cells = new cellAttributes[numberOfBoxes];
	int* superPXsNewPosition = new int[numberOfSuperPx];


	for (size_t i = 0; i < numberOfSuperPx; i++)
	{
		int pos = i * S;
		int superPXLine = pos / nbOfColumns;
		int superPXColumn = pos - superPXLine * nbOfColumns;
		superPXs[i].column = superPXColumn;
		superPXs[i].line = superPXLine;
	}

	for (size_t i = 0; i < numberOfBoxes; i++)
	{
		//Here 500000 is used as an equivalent to INFINITY, INFINITY has not been used because a finite number was needed
		cells[i].distance = 500000;
		cells[i].superPX = 0;

		int pos = i;
		int cellLine = pos / nbOfColumns;
		int cellColumn = pos - cellLine * nbOfColumns;

		cells[i].column = cellColumn;
		cells[i].line = cellLine;
		cells[i].alpha = alpha[i];
	}

	superPXs[0].totalDistance = 500000 * numberOfBoxes;
	superPXs[0].nbOfCells = numberOfBoxes;

	for (size_t iterations = 0; iterations < 1000; iterations++)
	{
		for (size_t i = 0; i < numberOfSuperPx; i++)
		{
			superPXsNewPosition[i] = 0;

			int lowRadiusCol = superPXs[i].column - 2 * S;
			int highRadiusCol = superPXs[i].column + 2 * S;
			int lowRadiusLine = superPXs[i].line - 2 * S;
			int highRadiusLine = superPXs[i].line + 2 * S;

			if (lowRadiusCol < 0) lowRadiusCol = 0;
			if (highRadiusCol >= nbOfColumns) highRadiusCol = nbOfColumns - 1;
			if (lowRadiusLine < 0) lowRadiusLine = 0;
			if (highRadiusLine >= nbOfLines) highRadiusLine = nbOfLines - 1;


			for (size_t column = lowRadiusCol; column < highRadiusCol; column++)
			{
				for (size_t line = lowRadiusLine; line < highRadiusLine; line++)
				{
					int currentCellPos = line * nbOfColumns + column;
					double distanceWithCenter = distanceComputationForSuperPx(cells[currentCellPos], superPXs[i]);

					if (distanceWithCenter < cells[currentCellPos].distance) {

						superPXs[cells[currentCellPos].superPX].totalDistance -= cells[currentCellPos].distance;
						superPXs[cells[currentCellPos].superPX].nbOfCells--;

						superPXs[i].totalDistance += distanceWithCenter;
						superPXs[i].nbOfCells++;

						cells[currentCellPos].distance = distanceWithCenter;
						cells[currentCellPos].superPX = i;
					}
				}
			}

		}

		for (size_t i = 0; i < numberOfBoxes; i++)
		{
			cellAttributes currentCell = cells[i];
			int pos = currentCell.line * nbOfColumns + currentCell.column;
			superPXsNewPosition[currentCell.superPX] += pos;
		}

		for (size_t i = 0; i < numberOfSuperPx; i++)
		{
			superPXsNewPosition[i] /= superPXs[i].nbOfCells;
			int pos = superPXsNewPosition[i];

			superPXs[i].line = pos / nbOfColumns;
			superPXs[i].column = pos - superPXs[i].line * nbOfColumns;
		}
	}

	for (size_t i = 0; i < numberOfBoxes; i++)
	{
		cellAttributes currentCell = cells[i];
		if (!currentCell.isFlagged) {
			if (!superPXs[currentCell.superPX].isFlagged) {
				flagSuperPX(i);
				superPXs[currentCell.superPX].isFlagged = true;
			}

			else {
				int superPx = superPXs.size();
				superPXCenter newSuperPX;
				newSuperPX.isFlagged = true;
				superPXs.push_back(newSuperPX);
				generateNewSuperPX(superPx, i);
				numberOfSuperPx++;
			}

		}
	}

	ccScalarField* SuperPXsField = new ccScalarField("SuperPXsField");
	SuperPXsField->resize(nbOfPoints);

	for (size_t i = 0; i < nbOfPoints; i++)
	{
		SuperPXsField->at(i) = NAN;
	}

	for (size_t i = 0; i < numberOfBoxes; i++) {
		if (!isnan(cells[i].alpha))
			superPXs[cells[i].superPX].alpha += cells[i].alpha;
	}

	for (size_t i = 0; i < numberOfSuperPx; i++)
	{
		superPXs[i].alpha /= superPXs[i].nbOfCells;
	}

	for (size_t i = 0; i < numberOfBoxes; i++)
	{
		ScalarType color = superPXs[cells[i].superPX].alpha;
		for (size_t j = 0; j < grid[i].size(); j++)
		{
			SuperPXsField->at(grid[i].at(j)) = color;
		}
	}


	SuperPXsField->computeMinAndMax();

	SuperPXsField->resizeSafe(nbOfPoints, true, CCCoreLib::NAN_VALUE);
	cloud->addScalarField(SuperPXsField);
	cloud->setCurrentDisplayedScalarField(cloud->getScalarFieldIndexByName("SuperPXsField"));
	cloud->showSF(true);


	m_glWindow->redraw();

	delete[] cells;
	delete[] superPXsNewPosition;
}

void elevationModel::flagSuperPX(int position) {

	cells[position].isFlagged = true;

	int line = cells[position].line;
	int column = cells[position].column;

	for (int currentNeighbourColumn = -1; currentNeighbourColumn < 2; currentNeighbourColumn++) {
		for (int currentNeighbourLine = -1; currentNeighbourLine < 2; currentNeighbourLine++)
		{
			int currentNeighbourPos = position + currentNeighbourColumn + currentNeighbourLine * nbOfColumns;
			if (currentNeighbourPos != position && currentNeighbourPos > -1 && currentNeighbourPos < numberOfBoxes) {
				if (cells[currentNeighbourPos].superPX == cells[position].superPX && !cells[currentNeighbourPos].isFlagged) {
					flagSuperPX(currentNeighbourPos);
				}
			}
		}
	}
}

void elevationModel::fillNewSuperPX(int formerSuperPX, int superPX, int position) {

	cells[position].isFlagged = true;
	superPXs[formerSuperPX].nbOfCells--;
	superPXs[superPX].nbOfCells++;
	cells[position].superPX = superPX;

	for (int currentNeighbourColumn = -1; currentNeighbourColumn < 2; currentNeighbourColumn++) {
		for (int currentNeighbourLine = -1; currentNeighbourLine < 2; currentNeighbourLine++)
		{
			int currentNeighbourPos = position + currentNeighbourColumn + currentNeighbourLine * nbOfColumns;
			if (currentNeighbourPos != position && currentNeighbourPos > -1 && currentNeighbourPos < numberOfBoxes) {
				if (cells[currentNeighbourPos].superPX == formerSuperPX && !cells[currentNeighbourPos].isFlagged) {
					fillNewSuperPX(formerSuperPX, superPX, currentNeighbourPos);
				}
			}
		}
	}
}

void elevationModel::generateNewSuperPX(int superPX, int position) {
	int formerSuperPX = cells[position].superPX;
	fillNewSuperPX(formerSuperPX, superPX, position);
}

//Simple Euclidean distance computation method but with the alpha value taken in account, this is only used in the GenerateSuperPX method
double elevationModel::distanceComputationForSuperPx(cellAttributes cell, superPXCenter center) {
	double alphaImportance = alphaBOX->value() / 100.;
	return sqrt(
		(((center.line - cell.line) * (center.line - cell.line)) +
			((center.column - cell.column) * (center.column - cell.column))) * (1 - alphaImportance) +
		(center.alpha - cell.alpha) * (center.alpha - cell.alpha) * alphaImportance
	);
}

void elevationModel::on_superPXSlider_valueChanged() {
	generateSuperPixelsWithSLIC();
	SPXGenerated = true;
}



void elevationModel::LBP() {
	//LOCAL BINARY PATTERN

	if (!LBPGenerated) {

		double* cellsGrayValues = new double[numberOfBoxes];
		std::vector<int>* neighbours = new std::vector<int>[numberOfBoxes];

		// Defining the average gray value of each cell
		for (int currentCell = 0; currentCell < numberOfBoxes; currentCell++)
		{
			double currentCellGray = 0;
			int currentCellSize = grid[currentCell].size();
			for (int i = 0; i < currentCellSize; i++)
			{
				currentCellGray += grayValues[grid[currentCell].at(i)];
			}
			currentCellGray /= currentCellSize;
			cellsGrayValues[currentCell] = currentCellGray;
		}


		for (int currentCell = 0; currentCell < numberOfBoxes; currentCell++)
		{
			for (int currentNeighbourColumn = -1; currentNeighbourColumn < 2; currentNeighbourColumn++) {
				for (int currentNeighbourLine = -1; currentNeighbourLine < 2; currentNeighbourLine++)
				{
					int currentNeighbourPos = currentCell + currentNeighbourColumn + currentNeighbourLine * nbOfColumns;
					if (currentNeighbourPos > -1 && currentNeighbourPos < numberOfBoxes)
						neighbours[currentCell].push_back(currentNeighbourPos);
				}
			}
		}

		int* binaryValues = new int[numberOfBoxes];

		for (int currentCell = 0; currentCell < numberOfBoxes; currentCell++)
		{
			int nbOfneighbours = neighbours[currentCell].size();
			int binaryValue = 0;
			int currentCellGrayValue = cellsGrayValues[currentCell];

			for (size_t i = 0; i < nbOfneighbours; i++)
			{
				int currentNeighbourGrayValue = grayValues[neighbours[currentCell].at(i)];
				if (currentNeighbourGrayValue < currentCellGrayValue)
					binaryValue = binaryValue * 10;
				else binaryValue = binaryValue * 10 + 1;
			}
			localBinaryPattern[currentCell] = BinaryToDecimal(binaryValue);
		}

		//showing the result on the point cloud
		for (size_t currentCell = 0; currentCell < numberOfBoxes; currentCell++)
		{
			int currentCellNbOfPoints = grid[currentCell].size();
			ccColor::Rgb currentCellColor = ccColor::Rgb(localBinaryPattern[currentCell], localBinaryPattern[currentCell], localBinaryPattern[currentCell]);
			for (size_t i = 0; i < currentCellNbOfPoints; i++)
			{
				cloud->setPointColor(grid[currentCell].at(i), currentCellColor);
			}
		}
		m_glWindow->redraw();


		delete[] cellsGrayValues;
		delete[] neighbours;
		delete[] binaryValues;
	}

	//Writing the result
	if (!FilePath2->text().isEmpty()) {

		QFile file(FilePath2->text());
		if (file.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			QTextStream stream(&file);
			stream << cloud->getName() << "\n" << "case dimensions : " << floor(stepU * 100.) << "cm²" << "\n" << "centroid : x = " << centroid->x << " y = " << centroid->y << " z = " << centroid->z << "\n LBP\n";

			for (size_t currentLine = 0; currentLine < nbOfLines; currentLine++)
			{
				for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
				{
					stream << "C " << currentLine << "_" << currentColumn << ";";
				}
				stream << "\n";

				for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
				{
					stream << localBinaryPattern[currentLine * nbOfColumns + currentColumn] << "; ";
				}
				stream << "\n";
			}
		}
		file.close();
	}

	if (AvgCSV2->isChecked()) {
		QFile fileAVG(FilePath2->text().insert(FilePath2->text().size() - 4, "AVG"));
		QTextStream stream(&fileAVG);
		stream << cloud->getName() << "\n" << "case dimensions : " << floor(stepU * 100.) << "cm²" << "\n" << "centroid : x = " << centroid->x << " y = " << centroid->y << " z = " << centroid->z << "\n LBP AVG\n";

		if (fileAVG.open(QIODevice::WriteOnly | QIODevice::Text)) {
			for (size_t currentLine = 0; currentLine < nbOfLines; currentLine++)
			{
				double avgLBP = 0;

				for (size_t currentColumn = 0; currentColumn < nbOfColumns; currentColumn++)
				{
					avgLBP += localBinaryPattern[currentLine * nbOfColumns + currentColumn];
				}
				avgLBP /= (nbOfColumns);

				stream << "ligne " << currentLine << "; " << avgLBP << " ;" << "\n";
			}
			fileAVG.close();
		}
	}

	cloud->showSF(false);


	if (VariogramGenerated)
		showAlpha->setEnabled(true);

	Generate2->setText("Save");
	showGrid->setEnabled(true);
	showLBP->setEnabled(false);
	showColors->setEnabled(true);
	showGray->setEnabled(true);

	LBPGenerated = true;
	
	m_glWindow->redraw();
}

//Spin boxes 
void elevationModel::on_cmBOX_valueChanged(double arg1) {
	pxBOX->setValue(arg1 / 100);
}

void elevationModel::on_pxBOX_valueChanged(double arg1) {
	cmBOX->setValue(arg1 * 100);

	stepU = arg1;
	stepV = arg1;
	if (arg1 == 0)
		return;

	else {

		//makes the grid's squares length more adapted to the cloud point studied
		double roundedDown = floor(lengthU / arg1);
		double remainder = fmod(lengthU, arg1);
		stepU = arg1 + remainder / roundedDown;

		roundedDown = floor(lengthV / arg1);
		remainder = fmod(lengthV, arg1);
		stepV = arg1 + remainder / roundedDown;


		double truc = floor(lengthU / arg1);
		double bidule = lengthU / arg1 - roundedDown;
		stepU = arg1 * remainder / roundedDown + arg1;



		roundedDown = floor(lengthV / arg1);
		remainder = lengthV / arg1 - roundedDown;
		stepV = arg1 * remainder / roundedDown + arg1;


		generateGrid();

	}
}

//Buttons in the "show" box

void elevationModel::on_showColors_clicked() {
	cloud->showSF(false);

	for (size_t i = 0; i < nbOfPoints; i++)
	{
		cloud->setPointColor(i, originalColors[i]);
	}

	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (LBPGenerated)
		showLBP->setEnabled(true);
	if (stepU > 0 && stepV > 0)
		showGrid->setEnabled(true);
	if (SPXGenerated)
		showSuperPXs->setEnabled(true);

	showColors->setEnabled(false);
	showGray->setEnabled(true);



	m_glWindow->redraw();

}

void elevationModel::on_showGray_clicked() {
	cloud->showSF(false);

	for (size_t pointIndex = 0; pointIndex < nbOfPoints; pointIndex++)
		cloud->setPointColor(pointIndex, ccColor::Rgb(grayValues[pointIndex], grayValues[pointIndex], grayValues[pointIndex]));

	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (LBPGenerated)
		showLBP->setEnabled(true);
	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (LBPGenerated)
		showLBP->setEnabled(true);
	if (stepU > 0 && stepV > 0)
		showGrid->setEnabled(true);
	if (SPXGenerated)
		showSuperPXs->setEnabled(true);

	showColors->setEnabled(true);
	showGray->setEnabled(false);


	m_glWindow->redraw();
}

void elevationModel::on_showGrid_clicked() {
	cloud->showSF(false);

	for (size_t i = 0; i < nbOfPoints; i++)
	{
		cloud->setPointColor(i, gridColors[i]);
	}

	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (LBPGenerated)
		showLBP->setEnabled(true);
	if (SPXGenerated)
		showSuperPXs->setEnabled(true);

	cloud->showSF(false);
	showGrid->setEnabled(false);
	showColors->setEnabled(true);
	showGray->setEnabled(true);


	m_glWindow->redraw();
}

void elevationModel::on_showSuperPXs_clicked() {
	cloud->showSF(true);
	int spxIndex = cloud->getScalarFieldIndexByName("SuperPXsField");
	if (spxIndex != -1)
		cloud->setCurrentDisplayedScalarField(spxIndex);

	showSuperPXs->setEnabled(false);
	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (LBPGenerated)
		showLBP->setEnabled(true);

	showSuperPXs->setEnabled(false);
	showGrid->setEnabled(true);
	showColors->setEnabled(true);
	showGray->setEnabled(true);

	m_glWindow->redraw();
}

void elevationModel::on_showAlpha_clicked() {
	cloud->showSF(true);

	int alphaIndex = cloud->getScalarFieldIndexByName("alpha");

	if (alphaIndex != -1)
		cloud->setCurrentDisplayedScalarField(alphaIndex);

	if (LBPGenerated)
		showLBP->setEnabled(true);
	if (SPXGenerated)
		showSuperPXs->setEnabled(true);

	showGrid->setEnabled(true);
	showAlpha->setEnabled(false);
	showColors->setEnabled(true);
	showGray->setEnabled(true);

	m_glWindow->redraw();
}

void elevationModel::on_showLBP_clicked() {
	cloud->showSF(false);

	for (size_t currentCell = 0; currentCell < numberOfBoxes; currentCell++)
	{
		int currentCellNbOfPoints = grid[currentCell].size();
		ccColor::Rgb currentCellColor = ccColor::Rgb(localBinaryPattern[currentCell], localBinaryPattern[currentCell], localBinaryPattern[currentCell]);
		for (size_t i = 0; i < currentCellNbOfPoints; i++)
		{
			cloud->setPointColor(grid[currentCell].at(i), currentCellColor);
		}
	}

	if (VariogramGenerated)
		showAlpha->setEnabled(true);
	if (SPXGenerated)
		showSuperPXs->setEnabled(true);

	showGrid->setEnabled(true);
	showLBP->setEnabled(false);
	showColors->setEnabled(true);
	showGray->setEnabled(true);

	m_glWindow->redraw();
}

//Generate buttons

void elevationModel::on_Generate_clicked()
{
	Variogram();
}

void elevationModel::on_Generate2_clicked() {
	LBP();
}

//close button

void elevationModel::on_Close_clicked()
{
	close();
}
