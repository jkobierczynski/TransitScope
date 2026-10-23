// SPDX-License-Identifier: GPL-3.0-or-later
#include "DiagnosticsView.h"

#include <QComboBox>
#include <QDir>
#include <QFileInfo>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>

#include <cmath>

#include "PlotWidget.h"
#include "transit/fit.h"
#include "transit/phot_table.h"

DiagnosticsView::DiagnosticsView(QWidget* parent) : QWidget(parent) {
    summary_ = new QLabel;
    summary_->setWordWrap(true);
    summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    bins_ = new PlotWidget;
    bins_->setLog(true, true);
    acf_ = new PlotWidget;
    scatter_ = new PlotWidget;
    table_ = new QTableWidget(0, 5);
    table_->setHorizontalHeaderLabels({tr("Quantity"), tr("r"), tr("limit"), tr("effect (ppt)"), tr("")});
    table_->horizontalHeaderItem(2)->setToolTip(tr("|r| above this is significant (3 σ with N/β² effective points)"));
    table_->horizontalHeaderItem(3)->setToolTip(tr("change of the residual across the 5–95 % range of the quantity"));
    table_->verticalHeader()->hide();
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setStretchLastSection(true);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    quantity_ = new QComboBox;

    auto* scatterBox = new QWidget;
    auto* sl = new QVBoxLayout(scatterBox);
    sl->setContentsMargins(0, 0, 0, 0);
    auto* qrow = new QHBoxLayout;
    qrow->addWidget(new QLabel(tr("Residual against:")));
    qrow->addWidget(quantity_, 1);
    sl->addLayout(qrow);
    sl->addWidget(scatter_, 1);

    auto* topRow = new QSplitter(Qt::Horizontal);
    topRow->addWidget(bins_);
    topRow->addWidget(acf_);
    auto* bottomRow = new QSplitter(Qt::Horizontal);
    bottomRow->addWidget(table_);
    bottomRow->addWidget(scatterBox);
    bottomRow->setStretchFactor(0, 2);
    bottomRow->setStretchFactor(1, 3);
    auto* split = new QSplitter(Qt::Vertical);
    split->addWidget(topRow);
    split->addWidget(bottomRow);

    auto* lay = new QVBoxLayout(this);
    lay->addWidget(summary_);
    lay->addWidget(split, 1);

    connect(quantity_, &QComboBox::currentIndexChanged, this, [this](int i) { showScatter(i); });
    connect(table_, &QTableWidget::currentCellChanged, this, [this](int row, int, int, int) {
        if (row >= 0) quantity_->setCurrentIndex(row);
    });
    bins_->setMessage(tr("Diagnostics appear after the fit stage."));
    acf_->setMessage(" ");
    scatter_->setMessage(" ");
}

void DiagnosticsView::load(const QString& dir, const QString& mode, size_t target) {
    bins_->clear();
    acf_->clear();
    scatter_->clear();
    table_->setRowCount(0);
    quantity_->blockSignals(true);
    quantity_->clear();
    quantity_->blockSignals(false);
    summary_->clear();
    d_ = {};
    residual_.clear();
    quantities_.clear();
    const QDir d(dir);
    const QString lc = d.filePath("lc.csv"), model = d.filePath("model_" + mode + ".csv"), phot = d.filePath("phot.csv");
    if (!QFileInfo::exists(lc) || !QFileInfo::exists(model)) {
        bins_->setMessage(tr("Diagnostics appear after the fit stage (no model_%1.csv yet).").arg(mode));
        acf_->setMessage(" ");
        scatter_->setMessage(" ");
        return;
    }
    try {
        const transit::FitData data = transit::load_fit_data(lc.toStdString());
        residual_ = transit::load_residuals(model.toStdString());
        if (residual_.size() != data.t.size())
            throw std::runtime_error(tr("the fit is older than the light curve; run the fit stage again").toStdString());
        if (QFileInfo::exists(phot)) {
            const transit::PhotTable t = transit::load_phot_csv(phot.toStdString());
            std::vector<size_t> comps;
            double radius = NAN;
            transit::read_lc_header(lc.toStdString(), comps, radius);
            for (size_t f : data.frame)
                if (f >= t.nf) throw std::runtime_error("lc.csv does not belong to phot.csv");
            quantities_ = transit::frame_quantities(t, data.frame, target, comps, radius);
        }
        d_ = transit::compute_diagnostics(residual_, quantities_);
    } catch (const std::exception& e) {
        bins_->setMessage(tr("No diagnostics: %1").arg(QString::fromStdString(e.what())));
        acf_->setMessage(" ");
        scatter_->setMessage(" ");
        residual_.clear();
        quantities_.clear();
        return;
    }

    // summary
    QString html = tr("<b>%1 fit:</b> residual rms %2 ppt, red-noise β = %3 (largest %4), %5 effective of %6 points.")
                       .arg(mode)
                       .arg(d_.rms_ppt, 0, 'f', 2)
                       .arg(d_.beta, 0, 'f', 2)
                       .arg(d_.beta_max, 0, 'f', 2)
                       .arg(d_.n_eff, 0, 'f', 0)
                       .arg(d_.n);
    if (d_.warnings.empty())
        html += tr(" <span style='color:#2a8a2a'>No significant red noise or systematics.</span>");
    for (const auto& w : d_.warnings)
        html += "<br><span style='color:#c05000'>⚠ " + QString::fromStdString(w).toHtmlEscaped() + "</span>";
    summary_->setText(html);

    // binned RMS, log-log, against the white-noise expectation
    PlotWidget::Series meas{tr("measured"), {}, {}, {}, QColor(40, 110, 200), false, 6};
    PlotWidget::Series measLine{QString(), {}, {}, {}, QColor(40, 110, 200), true, 1};
    PlotWidget::Series white{tr("white noise"), {}, {}, {}, QColor(220, 60, 40), true, 2};
    for (const auto& b : d_.bins) {
        meas.x << double(b.n); meas.y << b.rms * 1000;
        measLine.x << double(b.n); measLine.y << b.rms * 1000;
        white.x << double(b.n); white.y << b.white * 1000;
    }
    bins_->addSeries(white);
    bins_->addSeries(measLine);
    bins_->addSeries(meas);
    bins_->setLabels(tr("Binned residual rms (β = %1)").arg(d_.beta, 0, 'f', 2), tr("points per bin"), tr("rms (ppt)"));

    // autocorrelation with the ±3/√N band
    PlotWidget::Series acf{tr("autocorrelation"), {}, {}, {}, QColor(40, 110, 200), false, 6};
    PlotWidget::Series acfLine{QString(), {}, {}, {}, QColor(40, 110, 200), true, 1};
    for (size_t k = 0; k < d_.acf.size(); ++k) {
        acf.x << double(k + 1); acf.y << d_.acf[k];
        acfLine.x << double(k + 1); acfLine.y << d_.acf[k];
    }
    acf_->addHorizontalLine(0, palette().color(QPalette::Mid));
    acf_->addHorizontalLine(d_.acf_threshold, QColor(220, 60, 40));
    acf_->addHorizontalLine(-d_.acf_threshold, QColor(220, 60, 40));
    acf_->addSeries(acfLine);
    acf_->addSeries(acf);
    acf_->setLabels(tr("Autocorrelation of the residuals"), tr("lag (points)"), tr("correlation"));

    // correlations
    for (const auto& c : d_.correlations) {
        const int i = table_->rowCount();
        table_->insertRow(i);
        const bool ok = std::isfinite(c.r);
        auto item = [&](const QString& s) {
            auto* it = new QTableWidgetItem(s);
            if (c.flagged) it->setForeground(QColor(192, 80, 0));
            return it;
        };
        QString name = QString::fromStdString(c.name);
        if (!c.unit.empty()) name += " (" + QString::fromStdString(c.unit) + ")";
        table_->setItem(i, 0, item(name));
        table_->setItem(i, 1, item(ok ? QString::number(c.r, 'f', 2) : tr("constant")));
        table_->setItem(i, 2, item(QString::number(c.threshold, 'f', 2)));
        table_->setItem(i, 3, item(ok ? QString::number(c.span_ppt, 'f', 2) : QString()));
        table_->setItem(i, 4, item(c.flagged ? tr("⚠ significant") : QString()));
        quantity_->addItem(name);
    }
    if (quantities_.empty()) {
        scatter_->setMessage(tr("No phot.csv in the output folder: correlations with frame quantities need it."));
        return;
    }
    // start with the most significant correlation
    int best = 0;
    double bestScore = -1;
    for (size_t i = 0; i < d_.correlations.size(); ++i) {
        const auto& c = d_.correlations[i];
        const double s = std::isfinite(c.r) ? std::fabs(c.r) / c.threshold : -1;
        if (s > bestScore) { bestScore = s; best = int(i); }
    }
    quantity_->setCurrentIndex(best);
    showScatter(best);
}

void DiagnosticsView::showScatter(int qi) {
    scatter_->clear();
    if (qi < 0 || size_t(qi) >= quantities_.size()) return;
    const auto& q = quantities_[size_t(qi)];
    const auto& c = d_.correlations[size_t(qi)];
    PlotWidget::Series pts{QString(), {}, {}, {}, QColor(40, 110, 200, 150), false, 3};
    double xa = INFINITY, xb = -INFINITY, sx = 0, sy = 0;
    size_t n = 0;
    for (size_t i = 0; i < residual_.size() && i < q.values.size(); ++i) {
        if (!std::isfinite(q.values[i])) continue;
        pts.x << q.values[i];
        pts.y << residual_[i] * 1000;
        xa = std::min(xa, q.values[i]);
        xb = std::max(xb, q.values[i]);
        sx += q.values[i];
        sy += residual_[i] * 1000;
        ++n;
    }
    if (n == 0) {
        scatter_->setMessage(tr("%1 is not recorded for these frames.").arg(QString::fromStdString(q.name)));
        return;
    }
    scatter_->addSeries(pts);
    if (std::isfinite(c.r)) {   // least-squares line through the means
        const double mx = sx / double(n), my = sy / double(n);
        PlotWidget::Series fit{tr("slope %1 ppt per unit").arg(c.slope_ppt, 0, 'g', 3), {xa, xb},
                               {my + c.slope_ppt * (xa - mx), my + c.slope_ppt * (xb - mx)}, {},
                               c.flagged ? QColor(220, 60, 40) : QColor(120, 120, 120), true, 2};
        scatter_->addSeries(fit);
    }
    QString xl = QString::fromStdString(q.name);
    if (!q.unit.empty()) xl += " (" + QString::fromStdString(q.unit) + ")";
    scatter_->setLabels(tr("Residual vs %1: r = %2").arg(QString::fromStdString(q.name)).arg(c.r, 0, 'f', 2), xl, tr("residual (ppt)"));
}
