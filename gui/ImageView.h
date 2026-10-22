// SPDX-License-Identifier: GPL-3.0-or-later
// The reference frame with the star list drawn on it. FITS orientation (row 0 at the bottom).
//   left click on the sky       add a star there (centroided)
//   left click on a star        select it
//   right click on a star       set as target / exclude from comparison / remove
//   wheel                       zoom around the cursor; drag pans
#pragma once
#include <QGraphicsView>
#include <QMap>
#include <QSet>

#include "transit/image.h"
#include "transit/photometry.h"
#include "transit/reduce.h"

class QGraphicsPixmapItem;

class ImageView : public QGraphicsView {
    Q_OBJECT
public:
    explicit ImageView(QWidget* parent = nullptr);

    /// The calibrated (or raw) reference frame; `note` is shown in the corner (file, calibration).
    void setImage(transit::Image img, const QString& note);
    void clearImage(const QString& message);
    const transit::Image& image() const { return img_; }

    void setStars(const transit::StarList& stars, size_t target, const QSet<size_t>& excluded);
    void setAperture(const transit::ApertureConfig& ap) { ap_ = ap; redrawStars(); }
    /// Comparison stars of the last light curve with their ensemble weights (empty: none shown).
    void setComparisons(const QMap<size_t, double>& weights) { comps_ = weights; redrawStars(); }
    /// Stretch: black and white points as multiples of the background noise above the median.
    void setStretch(double low_sigma, double high_sigma);
    void fitInView();

signals:
    void starAdded(double x, double y);
    void starRemoveRequested(size_t index);
    void targetRequested(size_t index);
    void excludeToggled(size_t index);
    void statusMessage(const QString& text);

protected:
    void wheelEvent(QWheelEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;
    void drawForeground(QPainter* p, const QRectF& rect) override;

private:
    void render();
    void redrawStars();
    QPointF toScene(double x, double y) const;        // image pixel -> scene (y flipped)
    QPointF toImage(const QPointF& scene) const;
    long starAt(const QPointF& scene) const;          // -1 if none within the marker radius

    transit::Image img_;
    QString note_;
    QGraphicsPixmapItem* pix_ = nullptr;
    QList<QGraphicsItem*> markers_;
    transit::StarList stars_;
    size_t target_ = 0;
    QSet<size_t> excluded_;
    QMap<size_t, double> comps_;
    long selected_ = -1;
    transit::ApertureConfig ap_;
    double lowSigma_ = -1.5, highSigma_ = 60;
    QPoint pressPos_;
};
