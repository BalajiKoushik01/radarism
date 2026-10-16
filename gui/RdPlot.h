#pragma once
#include <QOpenGLFunctions>
#include <QOpenGLWidget>
#include <QImage>
#include "CpiData.h"

class RdPlot : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT
public:
    explicit RdPlot(QWidget* parent = nullptr);
    ~RdPlot() override;

    void setData(radar::CpiResultPtr r, double rangeDb = 45.0);
    void setSelected(int row, int col);
    void setColourRange(double rangeDb);
    void resetView();

signals:
    void cellSelected(int row, int col);
    void cursorInfo(const QString& info);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void wheelEvent(QWheelEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;

private:
    void rebuildImage();
    void clampView();
    QRectF plotRect() const;
    bool cellAt(const QPointF& pt, int& row, int& col) const;
    void emitCursor(const QPointF& pt);

    radar::CpiResultPtr r_;
    QImage img_;
    GLuint texture_ = 0;
    bool dirty_ = true;
    double rangeDb_ = 45.0;
    double peakDb_ = 0.0;
    double vc0_ = 0.0, vc1_ = 1.0;
    double vr0_ = 0.0, vr1_ = 1.0;
    int selRow_ = -1, selCol_ = -1;
    bool dragging_ = false, moved_ = false;
    QPointF press_, last_;
};
