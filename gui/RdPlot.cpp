#include "RdPlot.h"
#include <QColor>
#include <QDateTime>
#include <QDebug>
#include <QLinearGradient>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

namespace {
std::vector<double> niceTicks(double minVal, double maxVal, int maxTicks) {
    std::vector<double> ticks;
    double range = maxVal - minVal;
    if (range <= 0.0) return ticks;
    double approxStep = range / maxTicks;
    double mag = std::pow(10.0, std::floor(std::log10(approxStep)));
    double residual = approxStep / mag;
    double step;
    if (residual < 1.5) step = 1.0 * mag;
    else if (residual < 3.0) step = 2.0 * mag;
    else if (residual < 7.0) step = 5.0 * mag;
    else step = 10.0 * mag;

    double first = std::ceil(minVal / step) * step;
    for (double v = first; v <= maxVal + 1e-9 * step; v += step) {
        ticks.push_back(v);
    }
    return ticks;
}
} // namespace

RdPlot::RdPlot(QWidget* parent) : QOpenGLWidget(parent) {
    setMouseTracking(true);
}

RdPlot::~RdPlot() {
    makeCurrent();
    if (texture_) {
        glDeleteTextures(1, &texture_);
        texture_ = 0;
    }
    doneCurrent();
}

void RdPlot::setData(radar::CpiResultPtr r, double rangeDb) {
    r_ = r;
    rangeDb_ = rangeDb;
    selRow_ = -1;
    selCol_ = -1;
    resetView();
    dirty_ = true;
    update();
}

void RdPlot::setSelected(int row, int col) {
    selRow_ = row;
    selCol_ = col;
    update();
}

void RdPlot::setColourRange(double rangeDb) {
    rangeDb_ = std::max(rangeDb, 5.0);
    dirty_ = true;
    update();
}

void RdPlot::resetView() {
    if (!r_) return;
    vc0_ = 0.0;
    vc1_ = static_cast<double>(r_->gate_len);
    vr0_ = 0.0;
    vr1_ = static_cast<double>(r_->doppler_rows);
}

void RdPlot::clampView() {
    if (!r_) return;
    const double G = static_cast<double>(r_->gate_len), L = static_cast<double>(r_->doppler_rows);
    const double w = std::clamp(vc1_ - vc0_, std::min(4.0, G), G);
    const double h = std::clamp(vr1_ - vr0_, std::min(4.0, L), L);
    vc0_ = std::clamp(vc0_, 0.0, G - w);
    vc1_ = vc0_ + w;
    vr0_ = std::clamp(vr0_, 0.0, L - h);
    vr1_ = vr0_ + h;
}

void RdPlot::rebuildImage() {
    if (!r_) return;
    const std::size_t G = r_->gate_len, L = r_->doppler_rows;
    img_ = QImage(static_cast<int>(G), static_cast<int>(L), QImage::Format_RGBA8888);

    double maxVal = -1e30;
    for (float p : r_->rd_power) {
        if (p > maxVal) maxVal = p;
    }
    peakDb_ = 10.0 * std::log10(std::max(maxVal / r_->rd_noise_floor_w, 1e-14));
    const double minDb = peakDb_ - rangeDb_;

    for (std::size_t row = 0; row < L; ++row) {
        uchar* line = img_.scanLine(static_cast<int>(L - 1 - row));
        for (std::size_t col = 0; col < G; ++col) {
            const double p = static_cast<double>(r_->rd_power[row * G + col]);
            const double db = 10.0 * std::log10(std::max(p / r_->rd_noise_floor_w, 1e-14));
            const double norm = std::clamp((db - minDb) / rangeDb_, 0.0, 1.0);

            uchar red, green, blue;
            if (norm < 0.25) {
                double t = norm / 0.25;
                red = static_cast<uchar>(87 * t);
                green = static_cast<uchar>(16 * t);
                blue = static_cast<uchar>(110 * t);
            } else if (norm < 0.50) {
                double t = (norm - 0.25) / 0.25;
                red = static_cast<uchar>(87 + (188 - 87) * t);
                green = static_cast<uchar>(16 + (55 - 16) * t);
                blue = static_cast<uchar>(110 + (84 - 110) * t);
            } else if (norm < 0.75) {
                double t = (norm - 0.50) / 0.25;
                red = static_cast<uchar>(188 + (249 - 188) * t);
                green = static_cast<uchar>(55 + (142 - 55) * t);
                blue = static_cast<uchar>(84 + (9 - 84) * t);
            } else {
                double t = (norm - 0.75) / 0.25;
                red = static_cast<uchar>(249 + (252 - 249) * t);
                green = static_cast<uchar>(142 + (255 - 142) * t);
                blue = static_cast<uchar>(9 + (164 - 9) * t);
            }

            line[col * 4 + 0] = red;
            line[col * 4 + 1] = green;
            line[col * 4 + 2] = blue;
            line[col * 4 + 3] = 255;
        }
    }
    dirty_ = false;
}

void RdPlot::initializeGL() {
    initializeOpenGLFunctions();
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void RdPlot::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
}

void RdPlot::paintGL() {
    glClearColor(0.08f, 0.09f, 0.11f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    if (!r_) return;
    if (dirty_) rebuildImage();

    const double L = static_cast<double>(r_->doppler_rows);
    const QRectF pr = plotRect();

    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    const double u0 = vc0_ / static_cast<double>(r_->gate_len);
    const double u1 = vc1_ / static_cast<double>(r_->gate_len);
    const double v0 = vr0_ / L;
    const double v1 = vr1_ / L;

    const QRectF srcRect(u0 * img_.width(), (1.0 - v1) * img_.height(), (u1 - u0) * img_.width(), (v1 - v0) * img_.height());
    p.drawImage(pr, img_, srcRect);

    p.setPen(QColor(180, 190, 200));
    p.drawRect(pr);

    auto px = [&](double col) { return pr.left() + (col - vc0_) / (vc1_ - vc0_) * pr.width(); };
    auto py = [&](double row) { return pr.bottom() - (row - vr0_) / (vr1_ - vr0_) * pr.height(); };

    const double binM = radar::kSpeedOfLight / (2.0 * r_->sim.wf.fs);
    const double rMin = (r_->gate_start + vc0_) * binM / 1e3;
    const double rMax = (r_->gate_start + vc1_) * binM / 1e3;
    for (double rKm : niceTicks(rMin, rMax, 6)) {
        const double bin = rKm * 1e3 / binM - r_->gate_start;
        const double x = px(bin);
        if (x >= pr.left() && x <= pr.right()) {
            p.drawLine(QPointF(x, pr.bottom()), QPointF(x, pr.bottom() + 4));
            p.drawText(QRectF(x - 30, pr.bottom() + 5, 60, 15), Qt::AlignCenter, QString::number(rKm, 'f', 1));
        }
    }

    const double dv = r_->report.prf_hz / L * (radar::kSpeedOfLight / r_->sim.ch.fc) / 2.0;
    const double vMin = (vr0_ - L / 2.0) * dv;
    const double vMax = (vr1_ - L / 2.0) * dv;
    for (double vMps : niceTicks(vMin, vMax, 6)) {
        const double row = vMps / dv + L / 2.0;
        const double y = py(row);
        if (y >= pr.top() && y <= pr.bottom()) {
            p.drawLine(QPointF(pr.left() - 4, y), QPointF(pr.left(), y));
            p.drawText(QRectF(pr.left() - 55, y - 8, 50, 16), Qt::AlignRight | Qt::AlignVCenter, QString::number(vMps, 'f', 1));
        }
    }

    p.drawText(QRectF(pr.left(), height() - 18, pr.width(), 16), Qt::AlignCenter, "Range (km)");
    p.save();
    p.translate(10, pr.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-100, -8, 200, 16), Qt::AlignCenter, "Velocity (m/s)");
    p.restore();

    const double cbLeft = pr.right() + 12, cbWidth = 14;
    const QRectF cbRect(cbLeft, pr.top(), cbWidth, pr.height());
    QLinearGradient grad(cbRect.left(), cbRect.bottom(), cbRect.left(), cbRect.top());
    grad.setColorAt(0.00, QColor(0, 0, 4));
    grad.setColorAt(0.25, QColor(87, 16, 110));
    grad.setColorAt(0.50, QColor(188, 55, 84));
    grad.setColorAt(0.75, QColor(249, 142, 9));
    grad.setColorAt(1.00, QColor(252, 255, 164));
    p.fillRect(cbRect, grad);
    p.drawRect(cbRect);

    const double kRight = width() - 10;
    const double bottomDb = peakDb_ - rangeDb_;
    for (double db : niceTicks(bottomDb, peakDb_, 5)) {
        const double pyDb = pr.bottom() - (db - bottomDb) / rangeDb_ * pr.height();
        p.drawText(QRectF(cbLeft + cbWidth + 4, pyDb - 8, kRight - cbLeft - cbWidth - 4, 16),
                   Qt::AlignLeft | Qt::AlignVCenter, QString::number(db, 'f', 0) + " dB");
    }

    if (r_->sim.run_detection) {
        p.setPen(QPen(QColor(255, 60, 60), 1.5));
        p.setBrush(Qt::NoBrush);
        for (const auto& d : r_->detections) {
            const double c = static_cast<double>(d.gate_bin);
            const double r = static_cast<double>(d.doppler_bin);
            if (c < vc0_ || c > vc1_ || r < vr0_ || r > vr1_) continue;
            const double cx = px(c + 0.5), cy = py(r + 0.5);
            p.drawEllipse(QPointF(cx, cy), 6.0, 6.0);
        }
    }

    if (selRow_ >= 0 && selCol_ >= 0) {
        const double c = static_cast<double>(selCol_);
        const double r = static_cast<double>(selRow_);
        if (c >= vc0_ && c <= vc1_ && r >= vr0_ && r <= vr1_) {
            p.setPen(QPen(QColor(255, 255, 255), 2.0, Qt::DashLine));
            p.setBrush(Qt::NoBrush);
            const double cx = px(c + 0.5), cy = py(r + 0.5);
            p.drawRect(QRectF(cx - 5, cy - 5, 10, 10));
        }
    }
}

QRectF RdPlot::plotRect() const {
    return QRectF(60, 15, width() - 130, height() - 40);
}

bool RdPlot::cellAt(const QPointF& pt, int& row, int& col) const {
    const QRectF pr = plotRect();
    if (!r_ || !pr.contains(pt)) return false;
    const double u = (pt.x() - pr.left()) / pr.width();
    const double v = (pr.bottom() - pt.y()) / pr.height();
    col = static_cast<int>(std::floor(vc0_ + u * (vc1_ - vc0_)));
    row = static_cast<int>(std::floor(vr0_ + v * (vr1_ - vr0_)));
    col = std::clamp(col, 0, static_cast<int>(r_->gate_len) - 1);
    row = std::clamp(row, 0, static_cast<int>(r_->doppler_rows) - 1);
    return true;
}

void RdPlot::emitCursor(const QPointF& pt) {
    int row = 0, col = 0;
    if (!cellAt(pt, row, col)) {
        emit cursorInfo(QString());
        return;
    }
    const std::size_t G = r_->gate_len, L = r_->doppler_rows;
    const std::size_t absBin = (r_->gate_start + static_cast<std::size_t>(col)) % r_->n;
    const double rangeKm = radar::range_from_bin(static_cast<double>(absBin), r_->sim.wf.fs) / 1e3;
    const double fd = (static_cast<double>(row) - static_cast<double>(L / 2)) * r_->report.prf_hz / static_cast<double>(L);
    const double velMps = fd * (radar::kSpeedOfLight / r_->sim.ch.fc) / 2.0;
    const double pWatts = static_cast<double>(r_->rd_power[static_cast<std::size_t>(row) * G + static_cast<std::size_t>(col)]);
    const double snrDb = 10.0 * std::log10(pWatts / r_->rd_noise_floor_w);

    emit cursorInfo(QString("Range: %1 km | Velocity: %2 m/s (%3 Hz) | Power: %4 dB re floor")
                        .arg(rangeKm, 0, 'f', 3)
                        .arg(velMps, 0, 'f', 2)
                        .arg(fd, 0, 'f', 1)
                        .arg(snrDb, 0, 'f', 1));
}

void RdPlot::wheelEvent(QWheelEvent* e) {
    if (!r_) return;
    const QRectF pr = plotRect();
    const double u = std::clamp((e->position().x() - pr.left()) / pr.width(), 0.0, 1.0);
    const double v = std::clamp((pr.bottom() - e->position().y()) / pr.height(), 0.0, 1.0);
    const double cCur = vc0_ + u * (vc1_ - vc0_);
    const double rCur = vr0_ + v * (vr1_ - vr0_);
    const double factor = std::exp(-e->angleDelta().y() * 0.0015);

    const double newW = (vc1_ - vc0_) * factor;
    const double newH = (vr1_ - vr0_) * factor;
    vc0_ = cCur - u * newW;
    vc1_ = vc0_ + newW;
    vr0_ = rCur - v * newH;
    vr1_ = vr0_ + newH;

    clampView();
    update();
}

void RdPlot::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        dragging_ = true;
        moved_ = false;
        press_ = last_ = e->position();
    }
}

void RdPlot::mouseMoveEvent(QMouseEvent* e) {
    emitCursor(e->position());
    if (!dragging_ || !r_) return;
    const QRectF pr = plotRect();
    const QPointF pos = e->position();
    if ((pos - press_).manhattanLength() > 3) moved_ = true;

    const double dx = (pos.x() - last_.x()) / pr.width() * (vc1_ - vc0_);
    const double dy = (pos.y() - last_.y()) / pr.height() * (vr1_ - vr0_);
    last_ = pos;

    vc0_ -= dx; vc1_ -= dx;
    vr0_ += dy; vr1_ += dy;

    clampView();
    update();
}

void RdPlot::mouseReleaseEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        dragging_ = false;
        if (!moved_) {
            int row = 0, col = 0;
            if (cellAt(e->position(), row, col)) {
                setSelected(row, col);
                emit cellSelected(row, col);
            }
        }
    }
}

void RdPlot::mouseDoubleClickEvent(QMouseEvent*) {
    if (!r_) return;
    resetView();
    update();
}
