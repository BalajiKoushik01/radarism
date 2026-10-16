#pragma once
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <vector>
#include "PlotData.h"

class PlotWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT
public:
    explicit PlotWidget(QWidget* parent = nullptr);
    void setData(PlotDataPtr data);

protected:
    void initializeGL() override;
    void paintGL() override;
    void wheelEvent(QWheelEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;

private:
    QRectF plotRect() const;
    void clampView();
    bool drawPlot(const QRectF& pr, double dpr, double& ylo, double& yhi);
    void drawVerts(const std::vector<float>& v, GLenum mode, float r, float g, float b, float a);

    PlotDataPtr data_;
    double x0_ = 0.0, x1_ = 1.0;   // visible window in matched-filter bins
    bool dragging_ = false;
    double lastX_ = 0.0;

    QOpenGLShaderProgram prog_;
    QOpenGLBuffer vbo_;
    QOpenGLVertexArrayObject vao_;
    std::vector<float> mn_, mx_, thr_, grid_, env_, mid_, line_;
};
