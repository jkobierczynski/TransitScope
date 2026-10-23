// SPDX-License-Identifier: GPL-3.0-or-later
// Red-noise and systematics diagnostics of a fit (transit/diagnostics.h): binned RMS against the
// white-noise expectation (log-log), autocorrelation of the residuals, the correlation of the
// residuals with per-frame quantities, and the residuals plotted against a chosen quantity.
#pragma once
#include <QWidget>

#include <vector>

#include "transit/diagnostics.h"

class PlotWidget;
class QComboBox;
class QLabel;
class QTableWidget;

class DiagnosticsView : public QWidget {
    Q_OBJECT
public:
    explicit DiagnosticsView(QWidget* parent = nullptr);
    /// Computed from lc.csv, model_<mode>.csv and (if present) phot.csv in `outputDir`.
    void load(const QString& outputDir, const QString& mode, size_t target);
    const transit::Diagnostics& diagnostics() const { return d_; }

private:
    void showScatter(int quantity);
    QLabel* summary_;
    PlotWidget* bins_;
    PlotWidget* acf_;
    PlotWidget* scatter_;
    QTableWidget* table_;
    QComboBox* quantity_;
    transit::Diagnostics d_;
    std::vector<double> residual_;
    std::vector<transit::Quantity> quantities_;
};
