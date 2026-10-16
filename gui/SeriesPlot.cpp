#include "SeriesPlot.h"
#include <QFontDatabase>
#include <QMouseEvent>
#include <QNativeGestureEvent>
#include <QPainter>
#include <QVector4D>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace {
constexpr int kLeft = 66, kRight = 16, kTop = 14, kBottom = 40;
const QColor kBg(11, 14, 20);
const QColor kGridMajor(58, 66, 80);
const QColor kGridMinor(28, 33, 42);
const QColor kFrame(110, 120, 135);
const QColor kText(170, 180, 195);

double niceStep(double span, int target) {
    if (!(span > 0.0) || target < 1) return 0.0;
    const double raw = span / target;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double norm = raw / mag;
    return (norm < 1.5 ? 1.0 : norm < 3.5 ? 2.0 : norm < 7.5 ? 5.0 : 10.0) * mag;
}

// multiples of `step` inside [lo, hi]; skipFifths drops every 5th (those are the majors when step is major/5)
std::vector<double> multiples(double lo, double hi, double step, bool skipFifths) {
    std::vector<double> v;
    if (!(step > 0.0)) return v;
    const long long k0 = static_cast<long long>(std::ceil(lo / step));
    const long long k1 = static_cast<long long>(std::floor(hi / step));
    for (long long k = k0; k <= k1; ++k) {
        if (skipFifths && k % 5 == 0) continue;
        v.push_back(static_cast<double>(k) * step);
    }
    return v;
}
}  // namespace

SeriesPlot::SeriesPlot(QWidget* parent) : QOpenGLWidget(parent), vbo_(QOpenGLBuffer::VertexBuffer) {
    setMinimumSize(360, 200);
    setMouseTracking(true);
}

QRectF SeriesPlot::plotRect() const {
    return QRectF(kLeft, kTop, width() - kLeft - kRight, height() - kTop - kBottom);
}

void SeriesPlot::setSeries(std::shared_ptr<const void> keepAlive, std::vector<SeriesTrace> traces, double xStart,
                           double dx, double xScale, const QString& xLabel, const QString& yLabel) {
    keep_ = std::move(keepAlive);
    traces_ = std::move(traces);
    xStart_ = xStart;
    dx_ = dx;
    xScale_ = xScale;
    xLabel_ = xLabel;
    yLabel_ = yLabel;
    n_ = (traces_.empty() || !traces_[0].pyr) ? 0 : traces_[0].pyr->size();
    if (n_ < 2 || !(dx_ > 0.0) || !(xScale_ > 0.0)) n_ = 0;
    resetView();
    update();
}

void SeriesPlot::clearSeries() {
    traces_.clear();
    keep_.reset();
    n_ = 0;
    update();
}

void SeriesPlot::resetView() {
    if (n_ == 0) { viewXMin_ = 0.0; viewXMax_ = 1.0; return; }
    viewXMin_ = idxToX(0.0);                              // first sample at the left edge
    viewXMax_ = idxToX(static_cast<double>(n_ - 1));      // last sample at the right edge
}

void SeriesPlot::clampView() {
    if (n_ == 0) return;
    const double lo = idxToX(0.0), hi = idxToX(static_cast<double>(n_ - 1)), full = hi - lo;
    const double minSpan = std::min(4.0 * dx_ * xScale_, full);   // never closer than ~4 samples
    const double span = std::clamp(viewXMax_ - viewXMin_, minSpan, full);
    viewXMin_ = std::clamp(viewXMin_, lo, hi - span);
    viewXMax_ = viewXMin_ + span;
}

void SeriesPlot::zoomAt(double widgetX, double factor) {
    if (n_ == 0) return;
    const QRectF pr = plotRect();
    const double frac = std::clamp((widgetX - pr.left()) / pr.width(), 0.0, 1.0);
    const double span = viewXMax_ - viewXMin_;
    const double anchor = viewXMin_ + frac * span;
    const double ns = span * factor;
    viewXMin_ = anchor - frac * ns;
    viewXMax_ = viewXMin_ + ns;
    clampView();
    update();
}

void SeriesPlot::initializeGL() {
    initializeOpenGLFunctions();
    prog_.addShaderFromSourceCode(QOpenGLShader::Vertex,
        "#version 330 core\nlayout(location=0) in vec2 pos;\nvoid main(){ gl_Position = vec4(pos, 0.0, 1.0); }");
    prog_.addShaderFromSourceCode(QOpenGLShader::Fragment,
        "#version 330 core\nuniform vec4 color;\nout vec4 frag;\nvoid main(){ frag = color; }");
    prog_.link();
    vao_.create();
    vbo_.create();
}

void SeriesPlot::drawVerts(const std::vector<float>& v, GLenum mode, const QColor& c) {
    if (v.size() < 4) return;
    vao_.bind();
    prog_.bind();
    prog_.setUniformValue("color", QVector4D(c.redF(), c.greenF(), c.blueF(), 1.0f));
    vbo_.bind();
    vbo_.allocate(v.data(), static_cast<int>(v.size() * sizeof(float)));
    prog_.enableAttributeArray(0);
    prog_.setAttributeBuffer(0, GL_FLOAT, 0, 2);
    glDrawArrays(mode, 0, static_cast<GLsizei>(v.size() / 2));
    vbo_.release();
    prog_.release();
    vao_.release();
}

// Column c of the last query covers samples [a, b); its position is the centre of those sample indices,
// (a + b - 1) / 2, so a single-sample column sits exactly on its sample.
double SeriesPlot::columnCentre(std::size_t c) const {
    const std::size_t a = qi0_ + qspan_ * c / qnc_;
    const std::size_t b = qi0_ + qspan_ * (c + 1) / qnc_;
    return 0.5 * (static_cast<double>(a) + static_cast<double>(b) - 1.0);
}

bool SeriesPlot::drawTraces(const QRectF& pr, double dpr) {
    if (n_ < 2 || !(viewXMax_ > viewXMin_)) return false;
    const double span = viewXMax_ - viewXMin_;
    const double u0 = xToIdx(viewXMin_), u1 = xToIdx(viewXMax_);
    const long long N = static_cast<long long>(n_);
    const long long lo = std::max<long long>(0, static_cast<long long>(std::floor(u0)) - 1);
    const long long hi = std::min<long long>(N, static_cast<long long>(std::ceil(u1)) + 2);
    if (hi - lo < 2) return false;
    const std::size_t i0 = static_cast<std::size_t>(lo), i1 = static_cast<std::size_t>(hi);
    const std::size_t want = std::max<std::size_t>(2, static_cast<std::size_t>(pr.width() * dpr));

    const std::size_t T = traces_.size();
    mn_.resize(T);
    mx_.resize(T);
    std::size_t nc = 0;
    double vlo = std::numeric_limits<double>::infinity(), vhi = -vlo;
    for (std::size_t t = 0; t < T; ++t) {
        if (!traces_[t].pyr) return false;
        nc = traces_[t].pyr->query(i0, i1, want, mn_[t], mx_[t]);   // one column per pixel (fewer when zoomed in)
        if (nc == 0) return false;
        const double s = traces_[t].yScale;
        vlo = std::min(vlo, s * static_cast<double>(*std::min_element(mn_[t].begin(), mn_[t].end())));
        vhi = std::max(vhi, s * static_cast<double>(*std::max_element(mx_[t].begin(), mx_[t].end())));
    }
    qi0_ = i0;
    qspan_ = i1 - i0;
    qnc_ = nc;
    if (!(vhi - vlo > 1e-12)) vhi = vlo + 1.0;
    const double pad = 0.06 * (vhi - vlo);
    const double ylo = vlo - pad, yhi = vhi + pad;
    lastYlo_ = ylo;
    lastYhi_ = yhi;

    auto ndcX = [&](double xd) { return static_cast<float>((xd - viewXMin_) / span * 2.0 - 1.0); };
    auto ndcY = [&](double v) { return static_cast<float>((v - ylo) / (yhi - ylo) * 2.0 - 1.0); };

    glViewport(static_cast<int>(pr.left() * dpr), static_cast<int>((height() - pr.bottom()) * dpr),
               static_cast<int>(pr.width() * dpr), static_cast<int>(pr.height() * dpr));

    // Grid: minor lines (major/5) then major lines, in display units
    const double xstep = niceStep(span, 8), ystep = niceStep(yhi - ylo, 6);
    grid_.clear();
    for (double xv : multiples(viewXMin_, viewXMax_, xstep / 5.0, true)) {
        const float x = ndcX(xv);
        grid_.insert(grid_.end(), {x, -1.0f, x, 1.0f});
    }
    for (double yv : multiples(ylo, yhi, ystep / 5.0, true)) {
        const float y = ndcY(yv);
        grid_.insert(grid_.end(), {-1.0f, y, 1.0f, y});
    }
    drawVerts(grid_, GL_LINES, kGridMinor);
    grid_.clear();
    for (double xv : multiples(viewXMin_, viewXMax_, xstep, false)) {
        const float x = ndcX(xv);
        grid_.insert(grid_.end(), {x, -1.0f, x, 1.0f});
    }
    for (double yv : multiples(ylo, yhi, ystep, false)) {
        const float y = ndcY(yv);
        grid_.insert(grid_.end(), {-1.0f, y, 1.0f, y});
    }
    drawVerts(grid_, GL_LINES, kGridMajor);

    // Traces: one connected polyline per trace. Decimated columns contribute (min, max) in whichever order
    // is closer to the previous point, so adjacent columns join without gaps; single-sample columns
    // (deep zoom) contribute one vertex, giving exact sample-to-sample lines.
    for (std::size_t t = 0; t < T; ++t) {
        const double s = traces_[t].yScale;
        line_.clear();
        line_.reserve(nc * 4);
        float lastY = 0.0f;
        bool have = false;
        for (std::size_t c = 0; c < nc; ++c) {
            const float x = ndcX(idxToX(columnCentre(c)));
            const float ylo_n = ndcY(s * static_cast<double>(mn_[t][c]));
            const float yhi_n = ndcY(s * static_cast<double>(mx_[t][c]));
            if (mn_[t][c] == mx_[t][c]) {
                line_.insert(line_.end(), {x, ylo_n});
                lastY = ylo_n;
            } else if (!have || std::fabs(lastY - ylo_n) <= std::fabs(lastY - yhi_n)) {
                line_.insert(line_.end(), {x, ylo_n, x, yhi_n});
                lastY = yhi_n;
            } else {
                line_.insert(line_.end(), {x, yhi_n, x, ylo_n});
                lastY = ylo_n;
            }
            have = true;
        }
        QColor col = traces_[t].color;
        col.setAlpha(255);
        drawVerts(line_, GL_LINE_STRIP, col);
    }
    return true;
}

void SeriesPlot::paintGL() {
    const double dpr = devicePixelRatioF();
    const int W = static_cast<int>(width() * dpr), H = static_cast<int>(height() * dpr);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_BLEND);
    glViewport(0, 0, W, H);
    glClearColor(static_cast<float>(kBg.redF()), static_cast<float>(kBg.greenF()), static_cast<float>(kBg.blueF()), 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    const QRectF pr = plotRect();
    bool drawn = false;
    if (n_ > 1 && !traces_.empty() && prog_.isLinked() && pr.width() > 20 && pr.height() > 20) drawn = drawTraces(pr, dpr);
    glViewport(0, 0, W, H);

    QPainter p(this);
    p.setRenderHint(QPainter::TextAntialiasing);
    QFont f = p.font();
    f.setPointSizeF(8.5);
    p.setFont(f);
    p.setPen(kFrame);
    p.drawRect(pr);
    if (!drawn) {
        p.setPen(kText);
        p.drawText(pr, Qt::AlignCenter, "No data");
        return;
    }

    const double span = viewXMax_ - viewXMin_;
    auto px = [&](double xd) { return pr.left() + (xd - viewXMin_) / span * pr.width(); };
    auto py = [&](double v) { return pr.bottom() - (v - lastYlo_) / (lastYhi_ - lastYlo_) * pr.height(); };

    p.setPen(kText);
    for (double xv : multiples(viewXMin_, viewXMax_, niceStep(span, 8), false))
        p.drawText(QRectF(px(xv) - 36, pr.bottom() + 4, 72, 16), Qt::AlignCenter, QString::number(xv, 'g', 7));
    for (double yv : multiples(lastYlo_, lastYhi_, niceStep(lastYhi_ - lastYlo_, 6), false))
        p.drawText(QRectF(2, py(yv) - 8, kLeft - 8, 16), Qt::AlignRight | Qt::AlignVCenter, QString::number(yv, 'g', 4));
    p.drawText(QRectF(pr.left(), height() - 18, pr.width(), 16), Qt::AlignCenter, xLabel_);
    p.save();
    p.translate(11, pr.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-120, -8, 240, 16), Qt::AlignCenter, yLabel_);
    p.restore();

    // Sample markers when each sample is several pixels wide (exact sample positions)
    const double pxPerSample = pr.width() / (span / (dx_ * xScale_));
    if (qnc_ == qspan_ && pxPerSample >= 9.0) {
        p.setClipRect(pr);
        for (std::size_t t = 0; t < traces_.size(); ++t) {
            QColor col = traces_[t].color;
            col.setAlpha(255);
            p.setPen(Qt::NoPen);
            p.setBrush(col);
            for (std::size_t c = 0; c < qnc_; ++c) {
                const double x = px(idxToX(columnCentre(c)));
                if (x < pr.left() - 4 || x > pr.right() + 4) continue;
                const double y = py(traces_[t].yScale * static_cast<double>(mn_[t][c]));
                p.drawRect(QRectF(x - 2.5, y - 2.5, 5.0, 5.0));
            }
        }
        p.setClipping(false);
        p.setBrush(Qt::NoBrush);
    }

    // Legend
    double lx = pr.right() - 6;
    for (auto it = traces_.rbegin(); it != traces_.rend(); ++it) {
        if (it->name.isEmpty()) continue;
        const int w = p.fontMetrics().horizontalAdvance(it->name) + 14;
        lx -= w;
        QColor col = it->color;
        col.setAlpha(255);
        p.setPen(col);
        p.drawText(QRectF(lx, pr.top() + 3, w, 14), Qt::AlignRight | Qt::AlignVCenter, it->name);
    }

    // Crosshair with live readout: cursor x/y plus each trace's value at the nearest sample
    if (hover_ && pr.contains(cursor_) && !dragging_) {
        p.setPen(QPen(QColor(205, 215, 230, 130), 1.0));
        p.drawLine(QPointF(cursor_.x(), pr.top()), QPointF(cursor_.x(), pr.bottom()));
        p.drawLine(QPointF(pr.left(), cursor_.y()), QPointF(pr.right(), cursor_.y()));

        const double xd = viewXMin_ + (cursor_.x() - pr.left()) / pr.width() * span;
        const double yv = lastYlo_ + (pr.bottom() - cursor_.y()) / pr.height() * (lastYhi_ - lastYlo_);
        const long long k = std::clamp<long long>(std::llround(xToIdx(xd)), 0, static_cast<long long>(n_) - 1);
        QStringList lines;
        lines << QString("x = %1   y = %2").arg(xd, 0, 'g', 8).arg(yv, 0, 'g', 5);
        QString vals = QString("sample %1").arg(k);
        std::vector<float> a, b;
        for (const auto& tr : traces_) {
            if (!tr.pyr) continue;
            tr.pyr->query(static_cast<std::size_t>(k), static_cast<std::size_t>(k) + 1, 1, a, b);
            if (!a.empty()) vals += QString("   %1 = %2").arg(tr.name.isEmpty() ? QStringLiteral("v") : tr.name).arg(static_cast<double>(b[0]) * tr.yScale, 0, 'g', 5);
        }
        lines << vals;
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPointSizeF(8.5);
        p.setFont(mono);
        const QFontMetrics fm(mono);
        int w = 0;
        for (const auto& l : lines) w = std::max(w, fm.horizontalAdvance(l));
        const QSizeF box(w + 12, 2 * fm.height() + 8);
        QPointF o(cursor_.x() + 14, cursor_.y() + 14);
        if (o.x() + box.width() > pr.right()) o.setX(cursor_.x() - 14 - box.width());
        if (o.y() + box.height() > pr.bottom()) o.setY(cursor_.y() - 14 - box.height());
        p.fillRect(QRectF(o, box), QColor(11, 14, 20, 225));
        p.setPen(kFrame);
        p.drawRect(QRectF(o, box));
        p.setPen(QColor(225, 232, 242));
        p.drawText(QPointF(o.x() + 6, o.y() + fm.ascent() + 4), lines[0]);
        p.drawText(QPointF(o.x() + 6, o.y() + 2 * fm.ascent() + 6), lines[1]);
    }
}

void SeriesPlot::wheelEvent(QWheelEvent* e) {
    if (n_ == 0) return;
    const QPoint ad = e->angleDelta();
    if (std::abs(ad.x()) > std::abs(ad.y())) {   // horizontal two-finger scroll pans
        const double span = viewXMax_ - viewXMin_;
        const double shift = -ad.x() * 0.0015 * span;
        viewXMin_ += shift;
        viewXMax_ += shift;
        clampView();
        update();
    } else {                                     // vertical scroll zooms around the cursor
        zoomAt(e->position().x(), std::exp(-ad.y() * 0.0015));
    }
    e->accept();
}

bool SeriesPlot::event(QEvent* ev) {
    if (ev->type() == QEvent::NativeGesture && n_ > 0) {
        auto* g = static_cast<QNativeGestureEvent*>(ev);
        if (g->gestureType() == Qt::ZoomNativeGesture) {   // delivered only where the platform supports pinch
            zoomAt(g->position().x(), std::exp(-g->value()));
            return true;
        }
    }
    return QOpenGLWidget::event(ev);
}

void SeriesPlot::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        dragging_ = true;
        last_ = e->position();
    }
}

void SeriesPlot::mouseMoveEvent(QMouseEvent* e) {
    cursor_ = e->position();
    hover_ = true;
    if (dragging_ && n_ > 0) {
        const QRectF pr = plotRect();
        const double shift = -(e->position().x() - last_.x()) / pr.width() * (viewXMax_ - viewXMin_);
        viewXMin_ += shift;
        viewXMax_ += shift;
        clampView();
        last_ = e->position();
    }
    update();
}

void SeriesPlot::mouseReleaseEvent(QMouseEvent*) {
    dragging_ = false;
    update();
}

void SeriesPlot::mouseDoubleClickEvent(QMouseEvent*) {
    if (n_ == 0) return;
    resetView();
    update();
}

void SeriesPlot::leaveEvent(QEvent*) {
    hover_ = false;
    update();
}
