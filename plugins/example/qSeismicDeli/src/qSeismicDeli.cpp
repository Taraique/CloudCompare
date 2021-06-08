//##########################################################################
//#                                                                        #
//#                CLOUDCOMPARE PLUGIN: qSeismicDeli                      #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                             COPYRIGHT: XXX                             #
//#                                                                        #
//##########################################################################

#include <QtGui>
#include <QSettings>
#include <QFileInfo>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QMainWindow>

#include <sstream>

#include "Math.h"
#include "CCCoreLib.h"
#include "SquareMatrix.h"
#include "ccPickingHub.h"

#include "qSeismicDeli.h"
#include "elevationModel.h"
#include "ActionA.h"

#include "ccPointCloud.h"
#include <ccDrawableObject.h>
#include <ccBasicTypes.h>
#include <ccBox.h>

#include <ccScalarField.h>
#include <ScalarField.h>
#include "Jacobi.h"

#include <vector>

qSeismicDeli::qSeismicDeli(QObject* parent)
	: QObject(parent)
	, ccStdPluginInterface(":/CC/plugin/qSeismicDeli/info.json")
	, m_elevationModelAction(nullptr)
	, m_generateEigenData(nullptr)
	, m_test(nullptr)
{
}

void qSeismicDeli::onNewSelection(const ccHObject::Container& selectedEntities)
{
	m_elevationModelAction->setEnabled(!selectedEntities.empty());
	m_generateEigenData->setEnabled(!selectedEntities.empty());
	m_test->setEnabled(true);
}

QList<QAction*> qSeismicDeli::getActions()
{
	QList<QAction*> group;

	if (!m_generateEigenData) {
		m_generateEigenData = new QAction("Generate Eigen Values and Vectors", this);
		m_generateEigenData->setToolTip(getDescription());
		m_generateEigenData->setIcon(getIcon());

		connect(m_generateEigenData, &QAction::triggered, this, &qSeismicDeli::GenerateEigens);
	}
	group.push_back(m_generateEigenData);

	if (!m_elevationModelAction)
	{

		m_elevationModelAction = new QAction("Generate Elevation Model", this);
		m_elevationModelAction->setToolTip(getDescription());
		m_elevationModelAction->setIcon(getIcon());

		connect(m_elevationModelAction, &QAction::triggered, this, &qSeismicDeli::ElevationModel);
	}
	group.push_back(m_elevationModelAction);

	if (!m_test)
	{

		m_test = new QAction("bouton de test", this);
		m_test->setToolTip(getDescription());
		m_test->setIcon(getIcon());

		connect(m_test, &QAction::triggered, this, &qSeismicDeli::test);
	}

	group.push_back(m_test);

	return group;
}

void qSeismicDeli::ElevationModel() {
	m_app->dispToConsole(QString::number(highestU - lowestU) + "  " + QString::number(highestV - lowestV));
	elevationModel eM(highestU, lowestU, highestV, lowestV, m_app);
	eM.exec();
}

void qSeismicDeli::test() {
	/*float size = 0.06;

	ccPointCloud* cloud = static_cast<ccPointCloud*>(m_app->getSelectedEntities().front());

	ccGLMatrix* boxPosition = new ccGLMatrix(cloud->getGLTransformation());

	double tr[3] = { x,y,0 };

	boxPosition->setTranslation(tr);
	ccBox* box = new ccBox(CCVector3(0.06, 0.06, 0.06), boxPosition, QString("boite"));

	ccBBox box = ccBBox(CCVector3(minX, 3.25, minY ), CCVector3(maxX, 3.6, maxY));
	box.setValidity(true);
	CC_DRAW_CONTEXT context;
	context.bbDefaultCol = ccColor::blue;
	context.display = m_app->getActiveGLWindow();
	box.draw(context, ccColor::orangeRGB);
	auto truc = cloud->crop(box,true);

	m_app->dispToConsole((m_app->getActiveGLWindow()->whatsThis()));

	m_app->dispToConsole(QString::number(truc->size()));
	m_app->redrawAll();*/
}


void qSeismicDeli::GenerateEigens() {
	ccPointCloud* currentCloud = static_cast<ccPointCloud*>(m_app->getSelectedEntities().front());

	m_app->dispToConsole(currentCloud->getName());

	CCVector3 centroid;
	CCVector3 standardisedAverage;

	float standardDeviation = 0;
	int currentCloudNumberOfPoints = currentCloud->size();

	CCVector3* standardisedPoints;
	standardisedPoints = new CCVector3[currentCloudNumberOfPoints];

	//centroid computation
	for (int i = 0; i < currentCloudNumberOfPoints; i++) {
		centroid.x += currentCloud->getPoint(i)->x;
		centroid.y += currentCloud->getPoint(i)->y;
		centroid.z += currentCloud->getPoint(i)->z;
	}

	centroid.x /= currentCloudNumberOfPoints;
	centroid.y /= currentCloudNumberOfPoints;
	centroid.z /= currentCloudNumberOfPoints;

	//standard deviation computation
	for (int i = 0; i < currentCloudNumberOfPoints; i++) {
		standardDeviation += pow(centroid.x - currentCloud->getPoint(i)->x, 2)
			+ pow(centroid.y - currentCloud->getPoint(i)->y, 2)
			+ pow(centroid.z - currentCloud->getPoint(i)->z, 2);
	}

	standardDeviation = sqrt(standardDeviation / currentCloudNumberOfPoints);



	//standardisation
	for (int i = 0; i < currentCloudNumberOfPoints; i++) {
		standardisedPoints[i].x = (currentCloud->getPoint(i)->x - centroid.x) / standardDeviation;
		standardisedPoints[i].y = (currentCloud->getPoint(i)->y - centroid.y) / standardDeviation;
		standardisedPoints[i].z = (currentCloud->getPoint(i)->z - centroid.z) / standardDeviation;
		standardisedAverage.x += standardisedPoints[i].x;
		standardisedAverage.y += standardisedPoints[i].y;
		standardisedAverage.z += standardisedPoints[i].z;
	}

	standardisedAverage.x /= currentCloudNumberOfPoints;
	standardisedAverage.y /= currentCloudNumberOfPoints;
	standardisedAverage.z /= currentCloudNumberOfPoints;

	//Covariance matrix processing

	//initialisation
	CCCoreLib::SquareMatrixTpl<ScalarType> covariances = CCCoreLib::SquareMatrixTpl<ScalarType>(3);

	for (int i = 0; i < currentCloudNumberOfPoints; i++) {

		//à optimiser

		//xx
		covariances.m_values[0][0] += pow(standardisedPoints[i].x - standardisedAverage.x, 2);
		//xy
		covariances.m_values[1][0] += (standardisedPoints[i].x - standardisedAverage.x) * (standardisedPoints[i].y - standardisedAverage.y);
		//xz
		covariances.m_values[2][0] += (standardisedPoints[i].x - standardisedAverage.x) * (standardisedPoints[i].z - standardisedAverage.z);
		//yx
		covariances.m_values[0][1] += (standardisedPoints[i].y - standardisedAverage.y) * (standardisedPoints[i].x - standardisedAverage.x);
		//yy
		covariances.m_values[1][1] += pow(standardisedPoints[i].y - standardisedAverage.y, 2);
		//yz
		covariances.m_values[2][1] += (standardisedPoints[i].y - standardisedAverage.y) * (standardisedPoints[i].z - standardisedAverage.z);
		//zx
		covariances.m_values[0][2] += (standardisedPoints[i].z - standardisedAverage.z) * (standardisedPoints[i].x - standardisedAverage.x);
		//zy
		covariances.m_values[1][2] += (standardisedPoints[i].z - standardisedAverage.z) * (standardisedPoints[i].y - standardisedAverage.y);
		//zz
		covariances.m_values[2][2] += pow(standardisedPoints[i].z - standardisedAverage.z, 2);
	}


	for (int i = 0; i < 3; i++) {
		covariances.m_values[i][0] /= currentCloudNumberOfPoints;
		covariances.m_values[i][1] /= currentCloudNumberOfPoints;
		covariances.m_values[i][2] /= currentCloudNumberOfPoints;
	}

	//diagonalization
	CCCoreLib::SquareMatrixTpl<ScalarType> eVectors = CCCoreLib::SquareMatrixTpl<ScalarType>(3);
	std::vector<ScalarType> eValues;


	CCCoreLib::Jacobi<ScalarType> jacobi;

	jacobi.ComputeEigenValuesAndVectors(covariances, eVectors, eValues, false, 50);
	jacobi.SortEigenValuesAndVectors(eVectors, eValues);

	eigenVectors = eVectors;
	eigenValues = eValues;

	//u<v<h

	ScalarType h;
	ScalarType v;
	ScalarType u;

	Vector3Tpl<ScalarType> direction;

	ScalarType v1[3];
	ScalarType v2[3];
	ScalarType v3[3];

	jacobi.GetEigenVector(eVectors, 0, v1);
	jacobi.GetEigenVector(eVectors, 1, v2);
	jacobi.GetEigenVector(eVectors, 2, v3);

	ccScalarField* scalarField = new ccScalarField("u");
	ccScalarField* scalarField2 = new ccScalarField("v");
	ccScalarField* scalarField3 = new ccScalarField("h");

	for (int i = 0; i < currentCloudNumberOfPoints; i++) {

		// A OPTI

		direction.x = currentCloud->getPoint(i)->x - centroid.x;
		direction.y = currentCloud->getPoint(i)->y - centroid.y;
		direction.z = currentCloud->getPoint(i)->z - centroid.z;

		u = direction.x * v1[0] + direction.y * v1[1] + direction.z * v1[2];
		v = direction.x * v2[0] + direction.y * v2[1] + direction.z * v2[2];
		h = direction.x * v3[0] + direction.y * v3[1] + direction.z * v3[2];

		if (u > highestU) {
			highestU = u;
		}
		else if (u < lowestU) {
			lowestU = u;
		}

		if (v > highestV) {
			highestV = v;
		}

		else if (v < lowestV) {
			lowestV = v;
		}

		scalarField->emplace_back(u);
		scalarField2->emplace_back(v);
		scalarField3->emplace_back(h);
	}

	scalarField->computeMinAndMax();
	scalarField->resizeSafe(currentCloudNumberOfPoints, true, CCCoreLib::NAN_VALUE);

	scalarField2->computeMinAndMax();
	scalarField2->resizeSafe(currentCloudNumberOfPoints, true, CCCoreLib::NAN_VALUE);

	scalarField3->computeMinAndMax();
	scalarField3->resizeSafe(currentCloudNumberOfPoints, true, CCCoreLib::NAN_VALUE);

	currentCloud->addScalarField(scalarField);
	currentCloud->addScalarField(scalarField2);
	currentCloud->addScalarField(scalarField3);
	delete[] standardisedPoints;
}

void qSeismicDeli::launchEM() {

}