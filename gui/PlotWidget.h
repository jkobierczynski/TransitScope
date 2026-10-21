// SPDX-License-Identifier: GPL-3.0-or-later
// A small scatter/line plot drawn with QPainter: series of points (optionally with error bars) or
// lines, nice axis ticks, labels and a legend. Drag a rectangle to zoom, double-click to reset.
#pragma once
#include <QColor>
#include <QString>
#include <QVector>
#include <QWidget>

class PlotWidget : public QWidget {
    Q_OBJECT
public:
    struct Series {
        QString name;
        QVector<double> x, y, err;      // err empty: no error bars
        QColor color;
        bool line = false;              // false: points
        double size = 3;                // point diameter or line width, px
    };
    explicit PlotWidget(QWidget* parent = nullptr);

    void clear();
    void addSeries(const Series& s);
    void setLabels(const QString& title, const QString& x, const QString& y);
    void setMessage(const QString& msg);                  // shown when there is nothing to plot
    void addVerticalBand(double x0, double x1, const QColor& c);   // e.g. the transit window
    void addHorizontalLine(double y, const QColor& c);
    /// Shared x range with another plot (zooming one zooms the other).
    void linkX(PlotWidget* other);
    void setXRange(double x0, double x1, bool emitSignal = true);
    void resetZoom();

    QSize minimumSizeHint() const override { return {240, 160}; }

signals:
    void xRangeChanged(double x0, double x1);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;

private:
    QRectF dataBounds() const;
    QRect plotArea() const;
    void autoY();

    QVector<Series> series_;
    struct Band { double x0, x1; QColor c; };
    QVector<Band> bands_;
    struct HLine { double y; QColor c; };
    QVector<HLine> hlines_;
    QString title_, xlabel_, ylabel_, message_;
    double x0_ = 0, x1_ = 1, y0_ = 0, y1_ = 1;
    bool zoomed_ = false;
    QPoint dragStart_, dragNow_;
    bool dragging_ = false;
};
