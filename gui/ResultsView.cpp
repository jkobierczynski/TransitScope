// SPDX-License-Identifier: GPL-3.0-or-later
#include "ResultsView.h"

#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QPushButton>
#include <QSplitter>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>

#include <cmath>

#include "PlotWidget.h"

namespace {

struct Csv {
    QStringList header;
    QVector<QStringList> rows;
    QMap<QString, QString> comments;   // "# key=value" header lines
    int col(const QString& n) const { return int(header.indexOf(n)); }
};

Csv readCsv(const QString& path) {
    Csv c;
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return c;
    QTextStream in(&f);
    while (!in.atEnd()) {
        const QString line = in.readLine();
        if (line.startsWith('#')) {
            for (const QString& part : line.mid(1).split(' ', Qt::SkipEmptyParts)) {
                const int eq = part.indexOf('=');
                if (eq > 0) c.comments[part.left(eq)] = part.mid(eq + 1);
            }
            continue;
        }
        if (line.isEmpty()) continue;
        if (c.header.isEmpty()) c.header = line.split(',');
        else c.rows.push_back(line.split(','));
    }
    return c;
}

double num(const QString& s) {
    bool ok = false;
    const double v = s.toDouble(&ok);
    return ok ? v : NAN;
}

}  // namespace

ResultsView::ResultsView(QWidget* parent) : QWidget(parent) {
    auto* top = new QHBoxLayout;
    top->addWidget(new QLabel(tr("Fit:")));
    mode_ = new QComboBox;
    mode_->addItems({"constrained", "free"});
    top->addWidget(mode_);
    summary_ = new QLabel;
    summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    top->addWidget(summary_, 1);
    exoclock_ = new QPushButton(tr("Export for ExoClock…"));
    etd_ = new QPushButton(tr("Export for ETD…"));
    exoclock_->setToolTip(tr("Light curve as ExoClock's upload form reads it: JD_UTC at exposure start, flux"));
    etd_->setToolTip(tr("Light curve as ETD's upload reads it: geocentric JD_UTC at mid-exposure, flux"));
    exoclock_->setEnabled(false);
    etd_->setEnabled(false);
    top->addWidget(exoclock_);
    top->addWidget(etd_);
    connect(exoclock_, &QPushButton::clicked, this, [this] { emit exportRequested(0); });
    connect(etd_, &QPushButton::clicked, this, [this] { emit exportRequested(1); });

    lc_ = new PlotWidget;
    res_ = new PlotWidget;
    lc_->linkX(res_);
    table_ = new QTableWidget(0, 3);
    table_->setHorizontalHeaderLabels({tr("Parameter"), tr("Value"), tr("± (16–84 %)")});
    table_->verticalHeader()->hide();
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setStretchLastSection(true);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);

    auto* plots = new QSplitter(Qt::Vertical);
    plots->addWidget(lc_);
    plots->addWidget(res_);
    plots->setStretchFactor(0, 3);
    plots->setStretchFactor(1, 1);
    auto* split = new QSplitter(Qt::Horizontal);
    split->addWidget(plots);
    split->addWidget(table_);
    split->setStretchFactor(0, 3);
    split->setStretchFactor(1, 1);
    split->setSizes({860, 400});
    table_->setMinimumWidth(380);

    auto* lay = new QVBoxLayout(this);
    lay->addLayout(top);
    lay->addWidget(split, 1);
    connect(mode_, &QComboBox::currentIndexChanged, this, [this] { refresh(); });
    lc_->setMessage(tr("No results yet: run the project."));
    res_->setMessage(tr(" "));
}

void ResultsView::load(const QString& dir, double t14) {
    dir_ = dir;
    t14_ = t14;
    refresh();
}

void ResultsView::refresh() {
    lc_->clear();
    res_->clear();
    table_->setRowCount(0);
    summary_->clear();
    const QDir d(dir_);
    const Csv lc = readCsv(d.filePath("lc.csv"));
    exoclock_->setEnabled(!lc.rows.isEmpty());
    etd_->setEnabled(!lc.rows.isEmpty());
    if (lc.rows.isEmpty()) {
        lc_->setMessage(tr("No light curve yet (lc.csv in %1).\nRun the project up to the lightcurve stage.").arg(dir_));
        res_->setMessage(" ");
        return;
    }
    const double tc = num(lc.comments.value("tc_pred"));
    const double t0 = std::isfinite(tc) ? tc : num(lc.rows[0][lc.col("bjd_tdb")]);
    const QString mode = mode_->currentText();
    const Csv model = readCsv(d.filePath("model_" + mode + ".csv"));
    QFile jf(d.filePath("fit_" + mode + ".json"));
    QJsonObject fit;
    if (jf.open(QIODevice::ReadOnly)) fit = QJsonDocument::fromJson(jf.readAll()).object();

    const QString xl = tr("hours from predicted mid-transit (BJD_TDB %1)").arg(t0, 0, 'f', 5);
    if (!model.rows.isEmpty()) {
        // detrended flux and model
        PlotWidget::Series pts{tr("detrended flux"), {}, {}, {}, QColor(40, 110, 200), false, 3};
        PlotWidget::Series mod{tr("model"), {}, {}, {}, QColor(220, 60, 40), true, 2};
        PlotWidget::Series res{QString(), {}, {}, {}, QColor(40, 110, 200), false, 3};
        const int ct = model.col("bjd_tdb"), cf = model.col("flux_detrended"), cm = model.col("model_transit"), cr = model.col("residual");
        for (const auto& r : model.rows) {
            const double h = (num(r[ct]) - t0) * 24;
            pts.x << h; pts.y << num(r[cf]);
            mod.x << h; mod.y << num(r[cm]);
            res.x << h; res.y << num(r[cr]) * 1000;
        }
        lc_->addSeries(pts);
        lc_->addSeries(mod);
        lc_->setLabels(tr("Light curve, %1 fit").arg(mode), QString(), tr("relative flux (detrended)"));
        res_->addSeries(res);
        res_->addHorizontalLine(0, QColor(220, 60, 40));
        res_->setLabels(QString(), xl, tr("residual (ppt)"));
    } else {
        // the light curve alone: good frames, and flagged ones in grey
        PlotWidget::Series good{tr("good frames"), {}, {}, {}, QColor(40, 110, 200), false, 3};
        PlotWidget::Series bad{tr("flagged"), {}, {}, {}, QColor(160, 160, 160), false, 3};
        const int ct = lc.col("bjd_tdb"), cf = lc.col("flux"), cs = lc.col("sigma"), cg = lc.col("good");
        for (const auto& r : lc.rows) {
            auto& s = r[cg] == "1" ? good : bad;
            s.x << (num(r[ct]) - t0) * 24;
            s.y << num(r[cf]);
            s.err << num(r[cs]);
        }
        lc_->addSeries(good);
        if (!bad.x.isEmpty()) lc_->addSeries(bad);
        lc_->setLabels(tr("Light curve (not fitted yet)"), xl, tr("relative flux"));
        res_->setMessage(tr("Residuals appear after the fit stage."));
    }
    // the predicted transit window
    if (std::isfinite(tc)) {
        const double half = 0.5 * t14_;
        if (std::isfinite(half)) lc_->addVerticalBand(-half * 24, half * 24, QColor(255, 200, 0, 40));
    }
    if (fit.isEmpty()) {
        summary_->setText(tr("Comparison stars: %1   aperture %2 px").arg(lc.comments.value("comparison_stars")).arg(lc.comments.value("aperture_radius_px")));
        return;
    }
    // parameter table
    struct Row { const char* key; QString label; double scale; int dec; };
    const Row rows[] = {{"Tc", tr("Tc (BJD_TDB)"), 1, 6},       {"rp", tr("Rp/Rs"), 1, 4},
                        {"a", tr("a/Rs"), 1, 3},                {"inc", tr("inclination (°)"), 1, 2},
                        {"u1", tr("u1"), 1, 3},                 {"u2", tr("u2"), 1, 3},
                        {"c0", tr("baseline c0"), 1, 5},        {"c1", tr("baseline c1 (/day)"), 1, 5},
                        {"c2", tr("baseline c2 (/day²)"), 1, 4}, {"lnf", tr("ln noise scale"), 1, 3}};
    for (const auto& r : rows) {
        const QJsonObject o = fit.value(r.key).toObject();
        if (o.isEmpty()) continue;
        const int i = table_->rowCount();
        table_->insertRow(i);
        const double med = o.value("med").toDouble(), lo = o.value("lo").toDouble(), hi = o.value("hi").toDouble();
        table_->setItem(i, 0, new QTableWidgetItem(r.label));
        table_->setItem(i, 1, new QTableWidgetItem(QString::number(med * r.scale, 'f', r.dec)));
        table_->setItem(i, 2, new QTableWidgetItem(QString("+%1 / −%2").arg(hi, 0, 'f', r.dec).arg(lo, 0, 'f', r.dec)));
    }
    auto extra = [&](const QString& label, const QString& v) {
        const int i = table_->rowCount();
        table_->insertRow(i);
        table_->setItem(i, 0, new QTableWidgetItem(label));
        table_->setItem(i, 1, new QTableWidgetItem(v));
    };
    extra(tr("O–C (min)"), QString("%1 ± %2").arg(fit.value("O_C_min").toDouble(), 0, 'f', 2).arg(fit.value("Tc_err_inflated_min").toDouble(), 0, 'f', 2));
    extra(tr("depth (%)"), QString::number(fit.value("depth_pct").toDouble(), 'f', 3));
    extra(tr("impact parameter b"), QString::number(fit.value("b").toDouble(), 'f', 3));
    extra(tr("red-noise β"), QString::number(fit.value("beta").toDouble(), 'f', 2));
    extra(tr("rms residual (ppt)"), QString::number(fit.value("rms_resid").toDouble() * 1000, 'f', 2));
    extra(tr("acceptance"), QString::number(fit.value("acceptance").toDouble(), 'f', 3));
    extra(tr("max autocorrelation (steps)"), QString::number(fit.value("autocorr_max").toDouble(), 'f', 0));
    const QJsonObject tcO = fit.value("Tc").toObject();
    summary_->setText(tr("Tc = %1 ± %2 BJD_TDB   O–C %3 min   Rp/Rs = %4")
                          .arg(tcO.value("med").toDouble(), 0, 'f', 5)
                          .arg(0.5 * (tcO.value("lo").toDouble() + tcO.value("hi").toDouble()), 0, 'f', 5)
                          .arg(fit.value("O_C_min").toDouble(), 0, 'f', 2)
                          .arg(fit.value("rp").toObject().value("med").toDouble(), 0, 'f', 4));
}
