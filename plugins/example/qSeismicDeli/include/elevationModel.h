#ifndef ELEVATIONMODEL_H
#define ELEVATIONMODEL_H

#include <QDialog>
#include "ui_elevationModel.h"

#include <iostream>
#include <vector>
#include <ccGLWidget.h>

class ccMainAppInterface;
class ccGLWindow;
class ccPointCloud;
class ccHObject;
class ccBox;


class elevationModel : public QDialog, public Ui::elevationModel
{
    Q_OBJECT

public:
    explicit elevationModel(float hU, float lU, float hV, float lV, const CCVector3* center, ccMainAppInterface* app = 0);
    virtual ~elevationModel();

private:

    ccGLWindow* m_glWindow;
    ccMainAppInterface* m_app;
    ccPointCloud* cloud;
    ccColor::Rgba* originalColors;
    ccColor::Rgba* gridColors;
    int numberOfBoxes;
    float highestU;
    float highestV;
    float lowestU;
    float lowestV;
    float lengthU;
    float lengthV;
    float* grayValues;
    CCVector3* centroid;

    //the results to put into files
    std::vector<int>* grid = new std::vector<int>[0];
   
    int* localBinaryPattern = new int[0];
    
    double stepU;
    double stepV;
    int nbOfColumns;
    int nbOfLines;
    int nbOfPoints;

    boolean VariogramGenerated;
    boolean LBPGenerated;
    boolean SPXGenerated; 

    //Variogram parameters
    int nbOfT;
    double* gammaResults = new double[0];
    float t;
    double* alpha = new double[0];

    //For the Variogram computation
    struct coupleOfPoints {
        int firstSpouse;
        int secondSpouse;
    };

    //For the generation of super Pixels
    struct superPXCenter {
        float column;
        float line;
        ScalarType alpha = 0;
        float totalDistance = 0;
        float nbOfCells = 0;
        float movementFromLastIteration = 0;
        bool isFlagged = false;
    };

    struct cellAttributes {
        int column;
        int line;
        ScalarType alpha;
        bool isFlagged = false;

        float distance;
        int superPX;
    };

    cellAttributes* cells;
    std::vector<superPXCenter> superPXs;


private :
    void generateGrid();
    double distance(const CCVector3* point1, const CCVector3* point2);
    double distanceComputationForSuperPx(cellAttributes cell, superPXCenter center);
    float findT(std::vector<int> box);
    int BinaryToDecimal(int binary);
    void Variogram();
    void LBP();
    void generateSuperPixelsWithSLIC();
    void flagSuperPX(int position);
    void fillNewSuperPX(int formerSuperPX, int superPX, int position);
    void generateNewSuperPX(int superPX, int position);
    
private slots:
    void on_cmBOX_valueChanged(double arg1);
    void on_pxBOX_valueChanged(double arg1);
    void on_superPXSlider_valueChanged();
    void on_browseToolButton_clicked();
    void on_browseToolButton2_clicked();
    void on_Generate_clicked();
    void on_Generate2_clicked();
    void on_showGray_clicked();
    void on_showLBP_clicked();
    void on_showSuperPXs_clicked();

    void on_showAlpha_clicked();
    void on_showGrid_clicked();
    void on_showColors_clicked();

    void on_Close_clicked();

};

#endif // ELEVATIONMODEL_H
