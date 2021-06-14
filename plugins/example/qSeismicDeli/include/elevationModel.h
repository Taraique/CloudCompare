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
    explicit elevationModel(float hU, float lU, float hV, float lV, ccMainAppInterface* app = 0);
    virtual ~elevationModel();
    double pixelValue = 0;
    double exampleValue = 0;

private:
    ccGLWindow* m_glWindow;
    ccMainAppInterface* m_app;
    //! Broom box
    ccBox* m_broomBox;
    //! Boxes container
    ccHObject* m_boxes;
    ccPointCloud* cloud;
    void Variogram();
    int numberOfBoxes;
    float highestU;
    float highestV;
    float lowestU;
    float lowestV;
    float lengthU;
    float lengthV;
    std::vector<int>* grid = new std::vector<int>[0];
    void generateGrid();
    double distance(const CCVector3* point1, const CCVector3* point2);
    float findT(std::vector<int> box);
    double step;
    int test;
    int nbOfColumns;
    int nbOfLines;

    struct coupleOfPoints {
        int firstSpouse;
        int secondSpouse;
    };

private slots:
    void on_cmBOX_valueChanged(double arg1);
    void on_pxBOX_valueChanged(double arg1);
    void on_browseToolButton_clicked();
    void on_Generate_clicked();
    void on_Close_clicked();
};

#endif // ELEVATIONMODEL_H
