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
#include <ccPolyline.h>

#include "ccPointCloud.h"
#include <ccDrawableObject.h>
#include <ccBasicTypes.h>
#include <ccBox.h>

#include <ccScalarField.h>
#include <ScalarField.h>
#include <Jacobi.h>

#include <vector>
#include <Neighbourhood.h>

qSeismicDeli::qSeismicDeli(QObject* parent)
	: QObject(parent)
	, ccStdPluginInterface(":/CC/plugin/qSeismicDeli/info.json")
	, m_elevationModelAction(nullptr)
	, m_generateEigenData(nullptr)
{
}

void qSeismicDeli::onNewSelection(const ccHObject::Container& selectedEntities)
{
	m_elevationModelAction->setEnabled(!selectedEntities.empty());
	m_generateEigenData->setEnabled(!selectedEntities.empty());
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

	

	return group;
}

void qSeismicDeli::ElevationModel() {
	GenerateEigens();
	elevationModel eM(highestU, lowestU, highestV, lowestV,centroid, m_app);
	eM.exec();
	m_app->refreshAll();
}

//Generate the Eigens Values/Vector of the Point Cloud and the PCA

void qSeismicDeli::GenerateEigens() {

	ccPointCloud* currentCloud = static_cast<ccPointCloud*>(m_app->getSelectedEntities().front());
	int currentCloudNumberOfPoints = currentCloud->size();

	{
		int uFieldIndex = currentCloud->getScalarFieldIndexByName("u");
		if (uFieldIndex != -1)
			currentCloud->deleteScalarField(uFieldIndex);

		int vFieldIndex = currentCloud->getScalarFieldIndexByName("v");
		if (vFieldIndex != -1)
			currentCloud->deleteScalarField(vFieldIndex);

		int hFieldIndex = currentCloud->getScalarFieldIndexByName("h");
		if (hFieldIndex != -1)
			currentCloud->deleteScalarField(hFieldIndex);
	}

	Vector3Tpl<ScalarType> direction;

	double v1[3];
	double v2[3];
	double v3[3];

	CCCoreLib::SquareMatrixd eigVectors;
	std::vector<double> eigValues;

	//using the Cloud Compare libraries to get the Eigen Vectors
	CCCoreLib::Neighbourhood neighbourhood(currentCloud);
	centroid = neighbourhood.getGravityCenter();
	CCCoreLib::Jacobi<double>::ComputeEigenValuesAndVectors(neighbourhood.computeCovarianceMatrix(), eigVectors, eigValues, false);
	CCCoreLib::Jacobi<double>::SortEigenValuesAndVectors(eigVectors, eigValues);

	CCCoreLib::Jacobi<double>::GetEigenVector(eigVectors, 0, v1);
	CCCoreLib::Jacobi<double>::GetEigenVector(eigVectors, 1, v2);
	CCCoreLib::Jacobi<double>::GetEigenVector(eigVectors, 2, v3);
	


	ccScalarField* scalarField = new ccScalarField("u");
	ccScalarField* scalarField2 = new ccScalarField("v");
	ccScalarField* scalarField3 = new ccScalarField("h");

	
	//generating Scalar Fields with the results
	for (int i = 0; i < currentCloudNumberOfPoints; i++) {
		ScalarType h;
		ScalarType v;
		ScalarType u;
		
		direction.x = currentCloud->getPoint(i)->x - centroid->x;
		direction.y = currentCloud->getPoint(i)->y - centroid->x;
		direction.z = currentCloud->getPoint(i)->z - centroid->x;

		u = direction.x * v1[0] + direction.y * v1[1] + direction.z * v1[2];
		v = direction.x * v2[0] + direction.y * v2[1] + direction.z * v2[2];
		h = direction.x * v3[0] + direction.y * v3[1] + direction.z * v3[2];
		scalarField->emplace_back(u);
		scalarField2->emplace_back(v);
		scalarField3->emplace_back(h);

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

}
