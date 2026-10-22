// SPDX-License-Identifier: GPL-3.0-or-later
// Light curve, fitted model, residuals and the fitted parameters, read from a project's output
// folder (lc.csv, model_<mode>.csv, fit_<mode>.json).
#pragma once
#include <QWidget>

#include <cmath>

class PlotWidget;
class QComboBox;
class QLabel;
class QPushButton;
class QTableWidget;

class ResultsView : public QWidget {
    Q_OBJECT
public:
    explicit ResultsView(QWidget* parent = nullptr);
    /// Re-read the files in `outputDir`; tc_pred from lc.csv's header is the time zero.
    /// t14_days: transit duration, shown as the predicted transit window.
    void load(const QString& outputDir, double t14_days);

signals:
    void exportRequested(int target);   // transit::ExportTarget

private:
    void refresh();
    QString dir_;
    double t14_ = NAN;
    QComboBox* mode_;
    QLabel* summary_;
    PlotWidget* lc_;
    PlotWidget* res_;
    QTableWidget* table_;
    QPushButton* exoclock_;
    QPushButton* etd_;
};
