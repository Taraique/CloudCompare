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
    void getUVH();
    //std::vector<CCVector3> clickedPoints;

    float highestU;
    float highestV;
    float lowestU;
    float lowestV;
    float lengthU;
    float lengthV;
    std::vector<const CCVector3 *>* grid = new std::vector<const CCVector3 *>[0];
    void generateGrid(double step);

private slots:
    void on_cmBOX_valueChanged(double arg1);
    void on_pxBOX_valueChanged(double arg1);
    void on_buttonBox_accepted();
    void on_buttonBox_rejected();
    void on_browseToolButton_clicked();
};

#endif // ELEVATIONMODEL_H
