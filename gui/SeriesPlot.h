#pragma once
#include <QColor>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPointF>
#include <QString>
#include <memory>
#include <vector>
#include "radar/minmax_pyramid.hpp"

struct SeriesTrace {
    const radar::MinMaxPyramid* pyr = nullptr;   // must stay alive for as long as the plot's keep-alive handle does
    QColor color;                                // drawn fully opaque
    double yScale = 1.0;                         // displayed value = stored value * yScale (must be > 0)
    QString name;
};

// Min/max-decimated multi-trace plot. Sample i sits at display x = (xStart + i*dx) * xScale.
// The view is kept in DISPLAY units; every drawing path (grid, traces, crosshair, pan, zoom) goes through
// idxToX() / xToIdx(), so the axis and the data cannot disagree.
class SeriesPlot : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT
public:
    explicit SeriesPlot(QWidget* parent = nullptr);
    void setSeries(std::shared_ptr<const void> keepAlive, std::vector<SeriesTrace> traces, double xStart,
                   double dx, double xScale, const QString& xLabel, const QString& yLabel);
    void clearSeries();

protected:
    void initializeGL() override;
    void paintGL() override;
    void wheelEvent(QWheelEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void leaveEvent(QEvent* e) override;
    bool event(QEvent* e) override;

private:
    QRectF plotRect() const;
    double idxToX(double u) const { return (xStart_ + u * dx_) * xScale_; }
    double xToIdx(double x) const { return (x / xScale_ - xStart_) / dx_; }
    void resetView();
    void clampView();
    void zoomAt(double widgetX, double factor);
    bool drawTraces(const QRectF& pr, double dpr);
    void drawVerts(const std::vector<float>& v, GLenum mode, const QColor& c);
    double columnCentre(std::size_t c) const;   // sample-index position of query column c

    std::shared_ptr<const void> keep_;
    std::vector<SeriesTrace> traces_;
    std::size_t n_ = 0;
    double xStart_ = 0.0, dx_ = 1.0, xScale_ = 1.0;
    QString xLabel_, yLabel_;
    double viewXMin_ = 0.0, viewXMax_ = 1.0;     // visible window, display units

    // state of the last query (used by the overlay)
    std::size_t qi0_ = 0, qspan_ = 0, qnc_ = 0;
    double lastYlo_ = 0.0, lastYhi_ = 1.0;

    bool dragging_ = false, hover_ = false;
    QPointF last_, cursor_;

    QOpenGLShaderProgram prog_;
    QOpenGLBuffer vbo_;
    QOpenGLVertexArrayObject vao_;
    std::vector<std::vector<float>> mn_, mx_;
    std::vector<float> grid_, line_;
};
