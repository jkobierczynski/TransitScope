// SPDX-License-Identifier: GPL-3.0-or-later
#include "PlotWidget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// "nice" tick step: 1, 2 or 5 x 10^n, about n_ticks ticks over the span
double niceStep(double span, int nTicks) {
    if (!(span > 0)) return 1;
    const double raw = span / nTicks;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / mag;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
}

QString tickLabel(double v, double step) {
    const int decimals = std::max(0, int(-std::floor(std::log10(step) + 1e-9)));
    if (std::fabs(v) < step * 1e-6) v = 0;
    return QString::number(v, 'f', std::min(decimals, 8));
}

}  // namespace

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
    setMouseTracking(true);
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
}

void PlotWidget::clear() {
    series_.clear();
    bands_.clear();
    hlines_.clear();
    message_.clear();
    zoomed_ = false;
    update();
}

void PlotWidget::addSeries(const Series& s) {
    series_.push_back(s);
    if (!zoomed_) resetZoom();
}

void PlotWidget::setLabels(const QString& t, const QString& x, const QString& y) {
    title_ = t;
    xlabel_ = x;
    ylabel_ = y;
    update();
}

void PlotWidget::setMessage(const QString& m) {
    message_ = m;
    update();
}

void PlotWidget::addVerticalBand(double a, double b, const QColor& c) {
    bands_.push_back({a, b, c});
    update();
}

void PlotWidget::addHorizontalLine(double y, const QColor& c) {
    hlines_.push_back({y, c});
    update();
}

void PlotWidget::setLog(bool lx, bool ly) {
    logX_ = lx;
    logY_ = ly;
    resetZoom();
}

void PlotWidget::linkX(PlotWidget* other) {
    connect(this, &PlotWidget::xRangeChanged, other, [other](double a, double b) { other->setXRange(a, b, false); });
    connect(other, &PlotWidget::xRangeChanged, this, [this](double a, double b) { setXRange(a, b, false); });
}

QRectF PlotWidget::dataBounds() const {
    double xa = std::numeric_limits<double>::infinity(), xb = -xa, ya = xa, yb = -xa;
    for (const auto& s : series_)
        for (int i = 0; i < s.x.size(); ++i) {
            const double x = tx(s.x[i]), y = ty(s.y[i]);
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            const double e = !logY_ && i < s.err.size() && std::isfinite(s.err[i]) ? s.err[i] : 0;
            xa = std::min(xa, x);
            xb = std::max(xb, x);
            ya = std::min(ya, y - e);
            yb = std::max(yb, y + e);
        }
    if (!std::isfinite(xa)) return {0, 0, 1, 1};
    return QRectF(QPointF(xa, ya), QPointF(xb, yb));
}

void PlotWidget::autoY() {
    double ya = std::numeric_limits<double>::infinity(), yb = -ya;
    for (const auto& s : series_)
        for (int i = 0; i < s.x.size(); ++i) {
            const double x = tx(s.x[i]), y = ty(s.y[i]);
            if (!std::isfinite(y) || !(x >= x0_ && x <= x1_)) continue;
            const double e = !logY_ && i < s.err.size() && std::isfinite(s.err[i]) ? s.err[i] : 0;
            ya = std::min(ya, y - e);
            yb = std::max(yb, y + e);
        }
    for (const auto& h : hlines_)
        if (std::isfinite(ty(h.y))) { ya = std::min(ya, ty(h.y)); yb = std::max(yb, ty(h.y)); }
    if (!std::isfinite(ya)) { ya = 0; yb = 1; }
    if (yb - ya < 1e-12) { ya -= 0.5; yb += 0.5; }
    const double pad = 0.06 * (yb - ya);
    y0_ = ya - pad;
    y1_ = yb + pad;
}

void PlotWidget::resetZoom() {
    const QRectF b = dataBounds();
    const double pad = 0.02 * std::max(b.width(), 1e-12);
    x0_ = b.left() - pad;
    x1_ = b.right() + pad;
    if (x1_ - x0_ < 1e-12) { x0_ -= 0.5; x1_ += 0.5; }
    zoomed_ = false;
    autoY();
    update();
}

void PlotWidget::setXRange(double a, double b, bool emitSignal) {
    if (!(b > a)) return;
    x0_ = a;
    x1_ = b;
    zoomed_ = true;
    autoY();
    update();
    if (emitSignal) emit xRangeChanged(a, b);
}

QRect PlotWidget::plotArea() const {
    const QFontMetrics fm(font());
    const int left = fm.horizontalAdvance("-0.00000") + fm.height() + 8;
    const int top = title_.isEmpty() ? 10 : fm.height() + 12;
    const int bottom = 2 * fm.height() + 12;
    return rect().adjusted(left, top, -12, -bottom);
}

void PlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    const QPalette pal = palette();
    const QColor fg = pal.color(QPalette::Text), grid = pal.color(QPalette::Mid);
    const QRect a = plotArea();
    const QFontMetrics fm(font());

    if (series_.isEmpty()) {
        p.setPen(pal.color(QPalette::PlaceholderText));
        p.drawText(rect(), Qt::AlignCenter | Qt::TextWordWrap, message_.isEmpty() ? tr("No data") : message_);
        return;
    }
    // X, Y: data value -> pixel; AX, AY: axis coordinate (log10 on a log axis) -> pixel
    auto AX = [&](double x) { return a.left() + (x - x0_) / (x1_ - x0_) * a.width(); };
    auto AY = [&](double y) { return a.bottom() - (y - y0_) / (y1_ - y0_) * a.height(); };
    auto X = [&](double x) { return AX(tx(x)); };
    auto Y = [&](double y) { return AY(ty(y)); };
    // ticks: axis coordinate and its label
    auto ticks = [](double lo, double hi, bool log, int n) {
        QVector<QPair<double, QString>> t;
        if (log) {
            const bool sub = hi - lo < 2.5;   // 1, 2, 5 per decade over short ranges, decades otherwise
            for (int d = int(std::floor(lo)); d <= int(std::ceil(hi)); ++d)
                for (double m : {1.0, 2.0, 5.0}) {
                    if (m != 1 && !sub) continue;
                    const double v = d + std::log10(m);
                    if (v < lo - 1e-9 || v > hi + 1e-9) continue;
                    const double val = m * std::pow(10.0, d);
                    t.push_back({v, QString::number(val, 'g', 6)});
                }
            return t;
        }
        const double s = niceStep(hi - lo, n);
        for (double v = std::ceil(lo / s) * s; v <= hi; v += s) t.push_back({v, tickLabel(v, s)});
        return t;
    };

    // bands, grid and ticks
    p.save();
    p.setClipRect(a);
    for (const auto& b : bands_) p.fillRect(QRectF(QPointF(X(b.x0), a.top()), QPointF(X(b.x1), a.bottom())), b.c);
    p.restore();
    const auto xt = ticks(x0_, x1_, logX_, std::max(3, a.width() / 90));
    const auto yt = ticks(y0_, y1_, logY_, std::max(3, a.height() / 45));
    p.setPen(QPen(grid, 0, Qt::DotLine));
    for (const auto& t : xt) p.drawLine(QPointF(AX(t.first), a.top()), QPointF(AX(t.first), a.bottom()));
    for (const auto& t : yt) p.drawLine(QPointF(a.left(), AY(t.first)), QPointF(a.right(), AY(t.first)));
    p.setPen(fg);
    for (const auto& t : xt) p.drawText(QPointF(AX(t.first) - fm.horizontalAdvance(t.second) / 2.0, a.bottom() + fm.ascent() + 4), t.second);
    for (const auto& t : yt) p.drawText(QPointF(a.left() - fm.horizontalAdvance(t.second) - 5, AY(t.first) + fm.ascent() / 2.0 - 1), t.second);
    p.drawRect(a);
    if (!title_.isEmpty()) {
        QFont bf = font();
        bf.setBold(true);
        p.setFont(bf);
        p.drawText(QRect(a.left(), 2, a.width(), fm.height() + 6), Qt::AlignHCenter | Qt::AlignVCenter, title_);
        p.setFont(font());
    }
    p.drawText(QRect(a.left(), a.bottom() + fm.height() + 6, a.width(), fm.height() + 4), Qt::AlignHCenter, xlabel_);
    p.save();
    p.translate(fm.height() / 2.0 + 2, a.center().y());
    p.rotate(-90);
    p.drawText(QRect(-a.height() / 2, -fm.height() / 2, a.height(), fm.height() + 2), Qt::AlignHCenter, ylabel_);
    p.restore();

    // data
    p.save();
    p.setClipRect(a.adjusted(1, 1, -1, -1));
    for (const auto& h : hlines_) {
        p.setPen(QPen(h.c, 1, Qt::DashLine));
        p.drawLine(QPointF(a.left(), Y(h.y)), QPointF(a.right(), Y(h.y)));
    }
    for (const auto& s : series_) {
        if (s.line) {
            QPainterPath path;
            bool started = false;
            for (int i = 0; i < s.x.size(); ++i) {
                if (!std::isfinite(ty(s.y[i])) || !std::isfinite(tx(s.x[i]))) { started = false; continue; }
                const QPointF pt(X(s.x[i]), Y(s.y[i]));
                if (started) path.lineTo(pt); else path.moveTo(pt);
                started = true;
            }
            p.setPen(QPen(s.color, s.size));
            p.setBrush(Qt::NoBrush);
            p.drawPath(path);
        } else {
            QColor ec = s.color;
            ec.setAlpha(90);
            p.setPen(QPen(ec, 1));
            for (int i = 0; i < s.err.size() && i < s.x.size(); ++i)
                if (!logY_ && std::isfinite(s.err[i]) && std::isfinite(s.y[i]))
                    p.drawLine(QPointF(X(s.x[i]), Y(s.y[i] - s.err[i])), QPointF(X(s.x[i]), Y(s.y[i] + s.err[i])));
            p.setPen(Qt::NoPen);
            p.setBrush(s.color);
            const double r = s.size / 2;
            for (int i = 0; i < s.x.size(); ++i)
                if (std::isfinite(ty(s.y[i])) && std::isfinite(tx(s.x[i]))) p.drawEllipse(QPointF(X(s.x[i]), Y(s.y[i])), r, r);
        }
    }
    p.restore();

    // legend
    int ly = a.top() + 6;
    for (const auto& s : series_) {
        if (s.name.isEmpty()) continue;
        const int w = fm.horizontalAdvance(s.name);
        const int lx = a.right() - w - 30;
        p.fillRect(QRect(lx - 4, ly - 2, w + 30, fm.height() + 4), QColor(pal.color(QPalette::Base).red(), pal.color(QPalette::Base).green(), pal.color(QPalette::Base).blue(), 200));
        p.setPen(QPen(s.color, s.line ? s.size : 1));
        p.setBrush(s.color);
        if (s.line) p.drawLine(lx, ly + fm.height() / 2, lx + 16, ly + fm.height() / 2);
        else p.drawEllipse(QPointF(lx + 8, ly + fm.height() / 2.0), 3, 3);
        p.setPen(fg);
        p.drawText(lx + 22, ly + fm.ascent(), s.name);
        ly += fm.height() + 4;
    }

    if (dragging_) {
        p.setPen(QPen(pal.color(QPalette::Highlight), 1, Qt::DashLine));
        QColor hc = pal.color(QPalette::Highlight);
        hc.setAlpha(40);
        p.setBrush(hc);
        p.drawRect(QRect(dragStart_, dragNow_).normalized());
    }
}

void PlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton && plotArea().contains(e->pos())) {
        dragging_ = true;
        dragStart_ = dragNow_ = e->pos();
    }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* e) {
    const QRect a = plotArea();
    if (dragging_) {
        dragNow_ = e->pos();
        update();
    }
    if (a.contains(e->pos()) && !series_.isEmpty()) {
        double x = x0_ + (e->pos().x() - a.left()) * (x1_ - x0_) / a.width();
        double y = y0_ + (a.bottom() - e->pos().y()) * (y1_ - y0_) / a.height();
        if (logX_) x = std::pow(10.0, x);
        if (logY_) y = std::pow(10.0, y);
        setToolTip(QString("%1, %2").arg(x, 0, 'g', 8).arg(y, 0, 'g', 6));
    }
}

void PlotWidget::mouseReleaseEvent(QMouseEvent* e) {
    if (!dragging_) return;
    dragging_ = false;
    const QRect a = plotArea();
    const QRect r = QRect(dragStart_, e->pos()).normalized();
    if (r.width() > 6) {
        const double xa = x0_ + (r.left() - a.left()) * (x1_ - x0_) / a.width();
        const double xb = x0_ + (r.right() - a.left()) * (x1_ - x0_) / a.width();
        setXRange(xa, xb);
    }
    update();
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent*) {
    resetZoom();
    emit xRangeChanged(x0_, x1_);
}
