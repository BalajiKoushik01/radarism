#include "PlotWidget.h"
#include <QMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QVector4D>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

namespace {
constexpr int kLeft = 62, kRight = 14, kTop = 14, kBottom = 38;

std::vector<double> niceTicks(double lo, double hi, int target) {
    std::vector<double> t;
    const double span = hi - lo;
    if (!(span > 0.0) || target < 1) return t;
    const double raw = span / target;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double norm = raw / mag;
    const double step = (norm < 1.5 ? 1.0 : norm < 3.5 ? 2.0 : norm < 7.5 ? 5.0 : 10.0) * mag;
    const long long k0 = static_cast<long long>(std::ceil(lo / step));
    const long long k1 = static_cast<long long>(std::floor(hi / step));
    for (long long k = k0; k <= k1; ++k) t.push_back(static_cast<double>(k) * step);
    return t;
}
}  // namespace

PlotWidget::PlotWidget(QWidget* parent) : QOpenGLWidget(parent), vbo_(QOpenGLBuffer::VertexBuffer) {
    setMinimumSize(520, 320);
}

QRectF PlotWidget::plotRect() const {
    return QRectF(kLeft, kTop, width() - kLeft - kRight, height() - kTop - kBottom);
}

void PlotWidget::setData(PlotDataPtr data) {
    data_ = std::move(data);
    x0_ = 0.0;
    x1_ = data_ ? static_cast<double>(data_->trace.size()) : 1.0;
    update();
}

void PlotWidget::clampView() {
    if (!data_) return;
    const double n = static_cast<double>(data_->trace.size());
    double span = std::clamp(x1_ - x0_, std::min(20.0, n), n);
    x0_ = std::clamp(x0_, 0.0, n - span);
    x1_ = x0_ + span;
}

void PlotWidget::initializeGL() {
    initializeOpenGLFunctions();
    prog_.addShaderFromSourceCode(QOpenGLShader::Vertex,
        "#version 330 core\nlayout(location=0) in vec2 pos;\nvoid main(){ gl_Position = vec4(pos, 0.0, 1.0); }");
    prog_.addShaderFromSourceCode(QOpenGLShader::Fragment,
        "#version 330 core\nuniform vec4 color;\nout vec4 frag;\nvoid main(){ frag = color; }");
    prog_.link();
    vao_.create();
    vbo_.create();
}

void PlotWidget::drawVerts(const std::vector<float>& v, GLenum mode, float r, float g, float b, float a) {
    if (v.size() < 4) return;
    vao_.bind();
    prog_.bind();
    prog_.setUniformValue("color", QVector4D(r, g, b, a));
    vbo_.bind();
    vbo_.allocate(v.data(), static_cast<int>(v.size() * sizeof(float)));
    prog_.enableAttributeArray(0);
    prog_.setAttributeBuffer(0, GL_FLOAT, 0, 2);
    glDrawArrays(mode, 0, static_cast<GLsizei>(v.size() / 2));
    vbo_.release();
    prog_.release();
    vao_.release();
}

bool PlotWidget::drawPlot(const QRectF& pr, double dpr, double& ylo, double& yhi) {
    const auto& res = *data_->result;
    const std::size_t n = data_->trace.size();
    if (n < 2 || !(x1_ > x0_)) return false;
    const double fs = res.params.wf.fs, span = x1_ - x0_;
    const auto i0 = static_cast<std::size_t>(std::max(0.0, std::floor(x0_)));
    const auto i1 = std::min(n, static_cast<std::size_t>(std::ceil(x1_)));
    const auto want = std::max<std::size_t>(2, static_cast<std::size_t>(pr.width() * dpr));
    const std::size_t nc = data_->trace.query(i0, i1, want, mn_, mx_);
    if (nc == 0) return false;

    const double wspan = static_cast<double>(i1 - i0);
    auto centerBin = [&](std::size_t c) {
        return static_cast<double>(i0) + wspan * (static_cast<double>(c) + 0.5) / static_cast<double>(nc);
    };
    auto ndcX = [&](double bin) { return static_cast<float>((bin - x0_) / span * 2.0 - 1.0); };

    // CFAR threshold in the same dB units, sampled once per column
    const bool have_thr = !res.cfar.threshold.empty();
    float vlo = *std::min_element(mn_.begin(), mn_.end());
    float vhi = *std::max_element(mx_.begin(), mx_.end());
    thr_.assign(nc, 0.0f);
    if (have_thr) {
        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t idx = std::min(n - 1, static_cast<std::size_t>(centerBin(c)));
            const double db = 10.0 * std::log10(std::max(static_cast<double>(res.cfar.threshold[idx]) / data_->ref, 1e-12));
            thr_[c] = static_cast<float>(db);
        }
        vlo = std::min(vlo, *std::min_element(thr_.begin(), thr_.end()));
        vhi = std::max(vhi, *std::max_element(thr_.begin(), thr_.end()));
    }
    ylo = std::min<double>(vlo, -3.0);
    yhi = std::max<double>(vhi, 3.0);
    const double pad = 0.06 * (yhi - ylo);
    ylo -= pad;
    yhi += pad;
    auto ndcY = [&](double v) { return static_cast<float>((v - ylo) / (yhi - ylo) * 2.0 - 1.0); };

    glViewport(static_cast<int>(pr.left() * dpr), static_cast<int>((height() - pr.bottom()) * dpr),
               static_cast<int>(pr.width() * dpr), static_cast<int>(pr.height() * dpr));

    // grid (x ticks in km, y ticks in dB)
    const double bins_per_km = 1e3 * 2.0 * fs / radar::kSpeedOfLight;
    grid_.clear();
    for (double km : niceTicks(x0_ / bins_per_km, x1_ / bins_per_km, 8)) {
        const float x = ndcX(km * bins_per_km);
        grid_.insert(grid_.end(), {x, -1.0f, x, 1.0f});
    }
    for (double db : niceTicks(ylo, yhi, 6)) {
        const float y = ndcY(db);
        grid_.insert(grid_.end(), {-1.0f, y, 1.0f, y});
    }
    drawVerts(grid_, GL_LINES, 0.22f, 0.24f, 0.28f, 1.0f);

    // min/max envelope + mid-line (the mid-line keeps the trace visible when zoomed in to single samples)
    env_.resize(nc * 4);
    mid_.resize(nc * 2);
    for (std::size_t c = 0; c < nc; ++c) {
        const float x = ndcX(centerBin(c));
        env_[4 * c + 0] = x; env_[4 * c + 1] = ndcY(mn_[c]);
        env_[4 * c + 2] = x; env_[4 * c + 3] = ndcY(mx_[c]);
        mid_[2 * c + 0] = x; mid_[2 * c + 1] = ndcY(0.5 * (static_cast<double>(mn_[c]) + mx_[c]));
    }
    drawVerts(env_, GL_LINES, 0.30f, 0.70f, 1.0f, 1.0f);
    drawVerts(mid_, GL_LINE_STRIP, 0.30f, 0.70f, 1.0f, 1.0f);

    if (have_thr) {
        line_.resize(nc * 2);
        for (std::size_t c = 0; c < nc; ++c) {
            line_[2 * c] = ndcX(centerBin(c));
            line_[2 * c + 1] = ndcY(thr_[c]);
        }
        drawVerts(line_, GL_LINE_STRIP, 1.0f, 0.62f, 0.15f, 1.0f);
    }

    // detections
    line_.clear();
    for (std::size_t d : res.cfar.detections) {
        const double b = static_cast<double>(d);
        if (b < x0_ || b > x1_) continue;
        const float x = ndcX(b);
        line_.insert(line_.end(), {x, -1.0f, x, 1.0f});
    }
    drawVerts(line_, GL_LINES, 1.0f, 0.30f, 0.30f, 1.0f);

    // expected target position (apparent bin, wraps modulo N): ticks at the top and bottom edges
    const double exp_bin = std::fmod(res.report.delay_samples, static_cast<double>(n));
    if (exp_bin >= x0_ && exp_bin <= x1_) {
        const float x = ndcX(exp_bin);
        line_.assign({x, -1.0f, x, -0.82f, x, 0.82f, x, 1.0f});
        drawVerts(line_, GL_LINES, 0.25f, 0.85f, 0.35f, 1.0f);
    }
    return true;
}

void PlotWidget::paintGL() {
    const double dpr = devicePixelRatioF();
    const int W = static_cast<int>(width() * dpr), H = static_cast<int>(height() * dpr);
    // QPainter may leave state changed between frames, so reset what we depend on
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_BLEND);
    glViewport(0, 0, W, H);
    glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    const QRectF pr = plotRect();
    double ylo = 0.0, yhi = 1.0;
    bool drawn = false;
    if (data_ && prog_.isLinked() && pr.width() > 20 && pr.height() > 20) drawn = drawPlot(pr, dpr, ylo, yhi);
    glViewport(0, 0, W, H);

    // Text and frame via QPainter (after all raw GL work is finished and released)
    QPainter p(this);
    p.setRenderHint(QPainter::TextAntialiasing);
    QFont f = p.font();
    f.setPointSizeF(8.5);
    p.setFont(f);
    p.setPen(QColor(150, 155, 165));
    p.drawRect(pr);
    if (!drawn) {
        p.drawText(pr, Qt::AlignCenter, data_ ? "No range profile (CW carries no delay information)" : "Set parameters and press Generate");
        return;
    }
    const double fs = data_->result->params.wf.fs;
    const double bins_per_km = 1e3 * 2.0 * fs / radar::kSpeedOfLight;
    const double span = x1_ - x0_;
    for (double km : niceTicks(x0_ / bins_per_km, x1_ / bins_per_km, 8)) {
        const double px = pr.left() + (km * bins_per_km - x0_) / span * pr.width();
        p.drawText(QRectF(px - 32, pr.bottom() + 3, 64, 16), Qt::AlignCenter, QString::number(km, 'g', 6));
    }
    for (double db : niceTicks(ylo, yhi, 6)) {
        const double py = pr.bottom() - (db - ylo) / (yhi - ylo) * pr.height();
        p.drawText(QRectF(24, py - 8, kLeft - 28, 16), Qt::AlignRight | Qt::AlignVCenter, QString::number(db, 'g', 4));
    }
    p.drawText(QRectF(pr.left(), height() - 18, pr.width(), 16), Qt::AlignCenter, "Range (km)");
    p.save();
    p.translate(10, pr.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-110, -8, 220, 16), Qt::AlignCenter, "dB re matched-filter noise floor");
    p.restore();

    // peak-detection markers (the strongest peak is labelled with its refined range)
    if (!data_->peaks.empty()) {
        p.setClipRect(pr);
        for (std::size_t i = 0; i < data_->peaks.size(); ++i) {
            const auto& pk = data_->peaks[i];
            if (pk.bin_refined < x0_ || pk.bin_refined > x1_) continue;
            const double x = pr.left() + (pk.bin_refined - x0_) / (x1_ - x0_) * pr.width();
            const double y = pr.bottom() - (pk.snr_db - ylo) / (yhi - ylo) * pr.height();
            const double sz = (i == 0) ? 7.0 : 5.0;
            p.setPen(Qt::NoPen);
            p.setBrush(i == 0 ? QColor(255, 235, 59) : QColor(255, 235, 59, 150));
            QPolygonF tri;
            tri << QPointF(x, y - 3.0) << QPointF(x - sz, y - 3.0 - 1.6 * sz) << QPointF(x + sz, y - 3.0 - 1.6 * sz);
            p.drawPolygon(tri);
            if (i == 0) {
                p.setPen(QColor(255, 235, 59));
                p.drawText(QRectF(x - 80, y - 3.0 - 1.6 * sz - 24, 160, 22), Qt::AlignCenter,
                           QString::number(pk.range_m / 1e3, 'f', 3) + " km");
            }
        }
        p.setClipping(false);
        p.setBrush(Qt::NoBrush);
    }

    // legend
    double lx = pr.right() - 4;
    const struct { const char* text; QColor c; } items[] = {
        {"expected", QColor(64, 217, 89)}, {"detection", QColor(255, 77, 77)},
        {"CFAR threshold", QColor(255, 158, 38)}, {"matched filter", QColor(77, 179, 255)}};
    for (const auto& it : items) {
        const int w = p.fontMetrics().horizontalAdvance(it.text) + 14;
        lx -= w;
        p.setPen(it.c);
        p.drawText(QRectF(lx, pr.top() + 3, w, 14), Qt::AlignRight | Qt::AlignVCenter, it.text);
    }
}

void PlotWidget::wheelEvent(QWheelEvent* e) {
    if (!data_) return;
    const QRectF pr = plotRect();
    const double frac = std::clamp((e->position().x() - pr.left()) / pr.width(), 0.0, 1.0);
    const double span = x1_ - x0_;
    const double bin = x0_ + frac * span;
    const double newSpan = span * std::exp(-e->angleDelta().y() * 0.0015);
    x0_ = bin - frac * newSpan;
    x1_ = x0_ + newSpan;
    clampView();
    update();
}

void PlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) { dragging_ = true; lastX_ = e->position().x(); }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* e) {
    if (!dragging_ || !data_) return;
    const QRectF pr = plotRect();
    const double dx = e->position().x() - lastX_;
    lastX_ = e->position().x();
    const double shift = -dx / pr.width() * (x1_ - x0_);
    x0_ += shift;
    x1_ += shift;
    clampView();
    update();
}

void PlotWidget::mouseReleaseEvent(QMouseEvent*) { dragging_ = false; }

void PlotWidget::mouseDoubleClickEvent(QMouseEvent*) {
    if (!data_) return;
    x0_ = 0.0;
    x1_ = static_cast<double>(data_->trace.size());
    update();
}
