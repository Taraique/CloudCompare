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

#pragma once

#include "ccStdPluginInterface.h"
#include "elevationModel.h"

//! qSeismicDeli qCC plugin
/** Replace 'qSeismicDeli' by your own plugin class name throughout and then
		check 'qSeismicDeli.cpp' for more directions.

	Each plugin requires an info.json file to provide information about itself -
	the name, authors, maintainers, icon, etc..

	The one method you are required to implement is 'getActions'. This should
	return all actions (QAction objects) for the plugin. CloudCompare will
	automatically add these with their icons in the plugin toolbar and to the
	plugin menu. If	your plugin returns	several actions, CC will create a
	dedicated toolbar and a	sub-menu for your plugin. You are responsible for
	connecting these actions to	methods in your plugin.

	Use the ccStdPluginInterface::m_app variable for access to most of the CC
	components (database, 3D views, console, etc.) - see the ccMainAppInterface
	class in ccMainAppInterface.h.
**/

class qSeismicDeli : public QObject, public ccStdPluginInterface
{
	Q_OBJECT
		Q_INTERFACES(ccPluginInterface ccStdPluginInterface)

		// Replace "qSeismicDeli" by your plugin name (IID should be unique - let's hope your plugin name is unique ;)
	   // The info.json file provides information about the plugin to the loading system and
	   // it is displayed in the plugin information dialog.
		Q_PLUGIN_METADATA(IID "cccorp.cloudcompare.plugin.qSeismicDeli" FILE "../info.json")

public:
	explicit qSeismicDeli(QObject* parent = nullptr);
	~qSeismicDeli() override = default;

	// Inherited from ccStdPluginInterface
	void onNewSelection(const ccHObject::Container& selectedEntities) override;
	void ElevationModel();
	void GenerateEigens();
	QList<QAction*> getActions() override;

private:
	//! Default action
	/** You can add as many actions as you want in a plugin.
		Each action will correspond to an icon in the dedicated
		toolbar and an entry in the plugin menu.
	**/
	QAction* m_elevationModelAction;
	QAction* m_generateEigenData;

private:
	elevationModel* eM;
	bool isGenerated = true;
	CCCoreLib::SquareMatrixTpl<ScalarType> eigenVectors;
	std::vector<ScalarType> eigenValues;
	float highestU = -INFINITY;
	float highestV = -INFINITY;
	float lowestU = INFINITY;
	float lowestV = INFINITY;
	const CCVector3* centroid;



};
