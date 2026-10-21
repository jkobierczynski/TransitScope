// SPDX-License-Identifier: GPL-3.0-or-later
#include "ImageView.h"

#include <QContextMenuEvent>
#include <QGraphicsEllipseItem>
#include <QGraphicsPixmapItem>
#include <QGraphicsSimpleTextItem>
#include <QMenu>
#include <QScrollBar>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

#include "transit/stats.h"

ImageView::ImageView(QWidget* parent) : QGraphicsView(parent) {
    setScene(new QGraphicsScene(this));
    setBackgroundBrush(QColor(20, 20, 24));
    setDragMode(QGraphicsView::ScrollHandDrag);
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    setRenderHint(QPainter::Antialiasing);
    setMouseTracking(true);
    viewport()->setCursor(Qt::CrossCursor);
}

void ImageView::clearImage(const QString& message) {
    img_ = {};
    note_ = message;
    scene()->clear();
    pix_ = nullptr;
    markers_.clear();
    viewport()->update();
}

void ImageView::setImage(transit::Image img, const QString& note) {
    const bool firstImage = img_.empty() || img_.nx != img.nx || img_.ny != img.ny;
    img_ = std::move(img);
    note_ = note;
    render();
    redrawStars();
    if (firstImage) fitInView();
}

void ImageView::setStretch(double low, double high) {
    lowSigma_ = low;
    highSigma_ = high;
    render();
}

void ImageView::render() {
    if (img_.empty()) return;
    // background level and noise from a subsample: median and MAD
    std::vector<float> sample;
    const size_t step = std::max<size_t>(1, img_.px.size() / 200000);
    for (size_t i = 0; i < img_.px.size(); i += step)
        if (std::isfinite(img_.px[i])) sample.push_back(img_.px[i]);
    if (sample.empty()) return;
    const double med = transit::median(sample);
    for (auto& v : sample) v = static_cast<float>(std::fabs(v - med));
    const double sd = std::max(1e-6, 1.4826 * transit::median(sample));
    const double lo = med + lowSigma_ * sd, hi = med + highSigma_ * sd;
    // asinh stretch, so faint stars and bright cores are both visible
    const double beta = 8.0, norm = std::asinh(beta);
    QImage q(img_.nx, img_.ny, QImage::Format_Grayscale8);
    for (int y = 0; y < img_.ny; ++y) {
        uchar* row = q.scanLine(img_.ny - 1 - y);   // FITS orientation: row 0 at the bottom
        for (int x = 0; x < img_.nx; ++x) {
            double v = (img_.at(x, y) - lo) / (hi - lo);
            v = std::clamp(v, 0.0, 1.0);
            row[x] = static_cast<uchar>(255.0 * std::asinh(beta * v) / norm + 0.5);
        }
    }
    if (!pix_) {
        pix_ = scene()->addPixmap(QPixmap::fromImage(q));
        pix_->setZValue(0);
        scene()->setSceneRect(0, 0, img_.nx, img_.ny);
    } else {
        pix_->setPixmap(QPixmap::fromImage(q));
    }
}

void ImageView::fitInView() {
    if (img_.empty()) return;
    QGraphicsView::fitInView(QRectF(0, 0, img_.nx, img_.ny), Qt::KeepAspectRatio);
}

QPointF ImageView::toScene(double x, double y) const { return {x + 0.5, img_.ny - 0.5 - y}; }
QPointF ImageView::toImage(const QPointF& s) const { return {s.x() - 0.5, img_.ny - 0.5 - s.y()}; }

void ImageView::setStars(const transit::StarList& stars, size_t target, const QSet<size_t>& excluded) {
    stars_ = stars;
    target_ = target;
    excluded_ = excluded;
    if (selected_ >= long(stars_.size())) selected_ = -1;
    redrawStars();
}

void ImageView::redrawStars() {
    for (auto* m : markers_) {
        scene()->removeItem(m);
        delete m;
    }
    markers_.clear();
    if (img_.empty()) return;
    const double r = ap_.radii.empty() ? 10 : ap_.radii[ap_.radii.size() / 2];
    for (size_t i = 0; i < stars_.size(); ++i) {
        const QPointF c = toScene(stars_[i].first, stars_[i].second);
        QColor col = i == target_ ? QColor(255, 80, 80) : excluded_.contains(i) ? QColor(150, 150, 150) : QColor(80, 220, 120);
        QPen pen(col, 1.5);
        pen.setCosmetic(true);
        if (excluded_.contains(i)) pen.setStyle(Qt::DashLine);
        auto* e = scene()->addEllipse(c.x() - r, c.y() - r, 2 * r, 2 * r, pen);
        e->setZValue(2);
        markers_ << e;
        if (long(i) == selected_) {
            QPen ap(QColor(255, 210, 60), 1);
            ap.setCosmetic(true);
            ap.setStyle(Qt::DotLine);
            for (double rr : {ap_.ann_in, ap_.ann_out}) {
                auto* a = scene()->addEllipse(c.x() - rr, c.y() - rr, 2 * rr, 2 * rr, ap);
                a->setZValue(2);
                markers_ << a;
            }
        }
        auto* t = scene()->addSimpleText(i == target_ ? QString("%1 T").arg(i) : QString::number(i));
        t->setBrush(col);
        t->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        t->setPos(c.x() + r * 0.75, c.y() - r * 0.75);
        t->setZValue(3);
        markers_ << t;
    }
}

long ImageView::starAt(const QPointF& scene) const {
    const double r = ap_.radii.empty() ? 10 : ap_.radii[ap_.radii.size() / 2];
    const double pick = std::max(r, 8.0 / std::max(0.01, transform().m11()));   // at least 8 screen px
    long best = -1;
    double bd = pick * pick;
    for (size_t i = 0; i < stars_.size(); ++i) {
        const QPointF c = toScene(stars_[i].first, stars_[i].second);
        const double d = std::pow(c.x() - scene.x(), 2) + std::pow(c.y() - scene.y(), 2);
        if (d <= bd) { bd = d; best = long(i); }
    }
    return best;
}

void ImageView::wheelEvent(QWheelEvent* e) {
    const double f = std::pow(1.0015, e->angleDelta().y());
    scale(f, f);
}

void ImageView::mousePressEvent(QMouseEvent* e) {
    pressPos_ = e->pos();
    QGraphicsView::mousePressEvent(e);
}

void ImageView::mouseMoveEvent(QMouseEvent* e) {
    QGraphicsView::mouseMoveEvent(e);
    if (img_.empty()) return;
    const QPointF p = toImage(mapToScene(e->pos()));
    const int x = int(std::floor(p.x() + 0.5)), y = int(std::floor(p.y() + 0.5));
    if (x >= 0 && y >= 0 && x < img_.nx && y < img_.ny)
        emit statusMessage(QString("x %1  y %2  value %3").arg(p.x(), 0, 'f', 1).arg(p.y(), 0, 'f', 1).arg(img_.at(x, y), 0, 'f', 1));
}

void ImageView::mouseReleaseEvent(QMouseEvent* e) {
    QGraphicsView::mouseReleaseEvent(e);
    if (e->button() != Qt::LeftButton || img_.empty()) return;
    if ((e->pos() - pressPos_).manhattanLength() > 4) return;   // that was a drag
    const QPointF s = mapToScene(e->pos());
    const long hit = starAt(s);
    if (hit >= 0) {
        selected_ = hit;
        redrawStars();
        const auto& st = stars_[size_t(hit)];
        emit statusMessage(tr("star %1 at x %2, y %3").arg(hit).arg(st.first, 0, 'f', 2).arg(st.second, 0, 'f', 2));
        return;
    }
    const QPointF p = toImage(s);
    if (p.x() < 0 || p.y() < 0 || p.x() >= img_.nx || p.y() >= img_.ny) return;
    emit starAdded(p.x(), p.y());
}

void ImageView::contextMenuEvent(QContextMenuEvent* e) {
    if (img_.empty()) return;
    const long hit = starAt(mapToScene(e->pos()));
    QMenu menu(this);
    if (hit >= 0) {
        const size_t i = size_t(hit);
        selected_ = hit;
        redrawStars();
        menu.addSection(tr("Star %1").arg(i));
        auto* t = menu.addAction(tr("Set as target"));
        t->setEnabled(i != target_);
        auto* x = menu.addAction(excluded_.contains(i) ? tr("Allow as comparison star") : tr("Exclude from comparison stars"));
        x->setEnabled(i != target_);
        menu.addSeparator();
        auto* r = menu.addAction(tr("Remove star"));
        r->setEnabled(i != target_);
        if (i == target_) r->setToolTip(tr("Choose another target first"));
        QAction* a = menu.exec(e->globalPos());
        if (a == t) emit targetRequested(i);
        else if (a == x) emit excludeToggled(i);
        else if (a == r) emit starRemoveRequested(i);
    } else {
        auto* fit = menu.addAction(tr("Fit image in window"));
        if (menu.exec(e->globalPos()) == fit) fitInView();
    }
}

void ImageView::drawForeground(QPainter* p, const QRectF&) {
    if (note_.isEmpty()) return;
    p->save();
    p->resetTransform();
    p->setPen(img_.empty() ? QColor(200, 200, 200) : QColor(230, 230, 120));
    QFont f = font();
    p->setFont(f);
    const QRect area = viewport()->rect().adjusted(8, 6, -8, -6);
    p->drawText(area, img_.empty() ? Qt::AlignCenter | Qt::TextWordWrap : Qt::AlignLeft | Qt::AlignTop, note_);
    p->restore();
}
