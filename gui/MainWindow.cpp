// SPDX-License-Identifier: GPL-3.0-or-later
#include "MainWindow.h"

#include <QAction>
#include <QApplication>
#include <QCloseEvent>
#include <QComboBox>
#include <QDateTime>
#include <QDesktopServices>
#include <QUrl>
#include <QDir>
#include <QDockWidget>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QFutureWatcher>
#include <QLabel>
#include <QListWidget>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QScreen>
#include <QSettings>
#include <QStatusBar>
#include <QTabWidget>
#include <QThread>
#include <QTimer>
#include <QToolBar>
#include <QtConcurrent/QtConcurrentRun>

#include <chrono>
#include <cstdio>
#include <cmath>
#include <functional>
#include <streambuf>

#include "ImageView.h"
#include "ResultsView.h"
#include "SettingsEditor.h"
#include "transit/export.h"
#include "transit/image_reader.h"

using namespace transit;

namespace {

const Stage kStages[] = {Stage::Masters, Stage::Reduce, Stage::LightCurve, Stage::LimbDarkening, Stage::Fit};

// std::ostream sink that hands each finished line to a callback (run_project's log)
class LineBuf : public std::streambuf {
public:
    explicit LineBuf(std::function<void(const std::string&)> f) : f_(std::move(f)) {}
    ~LineBuf() override { if (!line_.empty()) f_(line_); }
protected:
    int overflow(int c) override {
        if (c == '\n') { f_(line_); line_.clear(); }
        else if (c != EOF) line_ += char(c);
        return c;
    }
private:
    std::function<void(const std::string&)> f_;
    std::string line_;
};

struct LoadedFrame {
    Image img;
    QString note, error;
};

QString describeStage(Stage s) {
    switch (s) {
        case Stage::Masters: return QObject::tr("master frames");
        case Stage::Reduce: return QObject::tr("reduce");
        case Stage::LightCurve: return QObject::tr("light curve");
        case Stage::LimbDarkening: return QObject::tr("limb darkening");
        case Stage::Fit: return QObject::tr("fit");
    }
    return {};
}

}  // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowIcon(QIcon(":/transitscope.png"));
    settings_ = new SettingsEditor;
    image_ = new ImageView;
    results_ = new ResultsView;
    tabs_ = new QTabWidget;
    tabs_->addTab(image_, tr("Reference frame"));
    tabs_->addTab(results_, tr("Results"));
    setCentralWidget(tabs_);

    auto* settingsDock = new QDockWidget(tr("Settings"), this);
    settingsDock->setObjectName("settings");
    settingsDock->setWidget(settings_);
    addDockWidget(Qt::LeftDockWidgetArea, settingsDock);

    problems_ = new QListWidget;
    auto* problemsDock = new QDockWidget(tr("Problems"), this);
    problemsDock->setObjectName("problems");
    problemsDock->setWidget(problems_);
    addDockWidget(Qt::BottomDockWidgetArea, problemsDock);

    log_ = new QPlainTextEdit;
    log_->setReadOnly(true);
    log_->setMaximumBlockCount(5000);
    log_->setFont(QFont("monospace"));
    auto* logDock = new QDockWidget(tr("Log"), this);
    logDock->setObjectName("log");
    logDock->setWidget(log_);
    addDockWidget(Qt::BottomDockWidgetArea, logDock);
    tabifyDockWidget(problemsDock, logDock);
    problemsDock->raise();

    // menus
    // menu entries as addAction(text, receiver, slot) + setShortcut: the form every Qt 6 version has
    auto* file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&New project"), this, &MainWindow::newProject)->setShortcut(QKeySequence::New);
    file->addAction(tr("&Open…"), this, &MainWindow::openDialog)->setShortcut(QKeySequence::Open);
    recentMenu_ = file->addMenu(tr("Open &recent"));
    file->addAction(tr("&Save"), this, &MainWindow::save)->setShortcut(QKeySequence::Save);
    file->addAction(tr("Save &as…"), this, &MainWindow::saveAs)->setShortcut(QKeySequence::SaveAs);
    file->addSeparator();
    file->addAction(tr("&Quit"), this, &QWidget::close)->setShortcut(QKeySequence::Quit);
    auto* stars = menuBar()->addMenu(tr("&Stars"));
    stars->addAction(tr("&Detect stars on the reference frame"), this, &MainWindow::detectStars);
    stars->addAction(tr("&Fit image in window"), image_, &ImageView::fitInView);
    auto* view = menuBar()->addMenu(tr("&View"));
    view->addAction(settingsDock->toggleViewAction());
    view->addAction(problemsDock->toggleViewAction());
    view->addAction(logDock->toggleViewAction());
    auto* help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&About TransitScope"), this, [this] {
        QMessageBox::about(this, tr("About TransitScope"),
                           tr("<b>TransitScope</b> %1<br>Exoplanet transit photometry: calibration, aperture photometry, "
                              "BJD_TDB timing, comparison-star ensemble, limb darkening and MCMC transit fit."
                              "<p>GPL-3.0-or-later. Uses Qt %2, CFITSIO, ERFA and toml++.")
                               .arg(QApplication::applicationVersion(), qVersion()));
    });

    // run toolbar
    auto* tb = addToolBar(tr("Run"));
    tb->setObjectName("run");
    runAct_ = tb->addAction(QIcon::fromTheme("media-playback-start"), tr("Run"));
    runAct_->setToolTip(tr("Run the stages from … to … (Ctrl+R)"));
    runAct_->setShortcut(QKeySequence("Ctrl+R"));
    tb->addWidget(new QLabel(tr("  from ")));
    fromStage_ = new QComboBox;
    toStage_ = new QComboBox;
    for (Stage s : kStages) {
        fromStage_->addItem(describeStage(s));
        toStage_->addItem(describeStage(s));
    }
    toStage_->setCurrentIndex(toStage_->count() - 1);
    tb->addWidget(fromStage_);
    tb->addWidget(new QLabel(tr(" to ")));
    tb->addWidget(toStage_);
    cancelAct_ = tb->addAction(QIcon::fromTheme("process-stop"), tr("Cancel"));
    cancelAct_->setEnabled(false);
    tb->addSeparator();
    stageLabel_ = new QLabel;
    stageLabel_->setMinimumWidth(220);
    tb->addWidget(stageLabel_);
    progress_ = new QProgressBar;
    progress_->setMaximumWidth(260);
    progress_->setTextVisible(true);
    progressAct_ = tb->addWidget(progress_);
    progressAct_->setVisible(false);
    connect(runAct_, &QAction::triggered, this, [this] {
        startRun(kStages[fromStage_->currentIndex()], kStages[toStage_->currentIndex()]);
    });
    connect(cancelAct_, &QAction::triggered, this, &MainWindow::cancelRun);
    connect(fromStage_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (toStage_->currentIndex() < i) toStage_->setCurrentIndex(i);
        revalidate();
    });
    connect(toStage_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (fromStage_->currentIndex() > i) fromStage_->setCurrentIndex(i);
        revalidate();
    });

    revalidateTimer_ = new QTimer(this);
    revalidateTimer_->setSingleShot(true);
    revalidateTimer_->setInterval(250);
    connect(revalidateTimer_, &QTimer::timeout, this, &MainWindow::revalidate);
    connect(settings_, &SettingsEditor::changed, this, [this] {
        setModified(true);
        revalidateTimer_->start();
    });

    connect(image_, &ImageView::statusMessage, this, [this](const QString& m) { statusBar()->showMessage(m, 4000); });
    connect(results_, &ResultsView::exportRequested, this, [this](int t) { exportLightCurve(t, QString()); });
    connect(image_, &ImageView::starAdded, this, &MainWindow::addStar);
    connect(image_, &ImageView::starRemoveRequested, this, &MainWindow::removeStar);
    connect(image_, &ImageView::targetRequested, this, &MainWindow::setTarget);
    connect(image_, &ImageView::excludeToggled, this, &MainWindow::toggleExclude);

    QSettings st;
    restoreGeometry(st.value("window/geometry").toByteArray());
    restoreState(st.value("window/state").toByteArray());
    if (!st.contains("window/state")) {
        resize(1400, 900);
        // dock sizes only take effect once the window is laid out
        QTimer::singleShot(0, this, [this, settingsDock, problemsDock] {
            resizeDocks({settingsDock}, {430}, Qt::Horizontal);
            resizeDocks({problemsDock}, {150}, Qt::Vertical);
        });
    }
    updateRecent(QString());
    newProject();
}

MainWindow::~MainWindow() {
    if (runThread_) {
        cancel_->store(true);
        runThread_->wait();
    }
}

// ---------------------------------------------------------------------------------------------
// project file

QString MainWindow::baseDir() const { return path_.isEmpty() ? QDir::homePath() : QFileInfo(path_).absolutePath(); }

void MainWindow::setPath(const QString& path) {
    path_ = path;
    settings_->setBaseDir(baseDir());
    setModified(modified_);
}

void MainWindow::setModified(bool m) {
    modified_ = m;
    const QString name = path_.isEmpty() ? tr("untitled") : QFileInfo(path_).fileName();
    setWindowTitle(QString("%1%2 — TransitScope").arg(name, m ? "*" : ""));
}

void MainWindow::loadSettingsFromProject(const Project& p) {
    StarList stars = p.reduce.stars;
    settings_->load(QString::fromStdString(project_to_toml(p)), stars);
}

void MainWindow::newProject() {
    if (!maybeSave()) return;
    const Project p = parse_project(project_template(), QDir::homePath().toStdString());
    setPath(QString());
    loadSettingsFromProject(p);
    log_->clear();
    imageKey_.clear();
    image_->clearImage(tr("New project: an example to edit.\n\nSet the folders under [frames], the target and the "
                          "ephemeris, then save it (File → Save as) next to your data."));
    setModified(false);
    revalidate();
}

bool MainWindow::openProject(const QString& file) {
    try {
        const Project p = load_project(file.toStdString());
        setPath(QFileInfo(file).absoluteFilePath());
        loadSettingsFromProject(p);
        setModified(false);
        log_->clear();
        imageKey_.clear();
        updateRecent(path_);
        revalidate();
        return true;
    } catch (const std::exception& e) {
        QMessageBox::warning(this, tr("Open project"), QString::fromStdString(e.what()));
        return false;
    }
}

void MainWindow::openDialog() {
    if (!maybeSave()) return;
    const QString f = QFileDialog::getOpenFileName(this, tr("Open project"), baseDir(), tr("TransitScope projects (*.toml);;All files (*)"));
    if (!f.isEmpty()) openProject(f);
}

bool MainWindow::save() {
    if (path_.isEmpty()) return saveAs();
    revalidate();
    if (!project_) {
        QMessageBox::warning(this, tr("Save"), tr("The settings have problems (see Problems); fix them before saving."));
        return false;
    }
    try {
        save_project(path_.toStdString(), *project_);
    } catch (const std::exception& e) {
        QMessageBox::warning(this, tr("Save"), QString::fromStdString(e.what()));
        return false;
    }
    setModified(false);
    updateRecent(path_);
    statusBar()->showMessage(tr("Saved %1").arg(path_), 3000);
    return true;
}

bool MainWindow::saveAs() {
    QString f = QFileDialog::getSaveFileName(this, tr("Save project"), path_.isEmpty() ? QDir::homePath() + "/project.toml" : path_,
                                             tr("TransitScope projects (*.toml)"));
    if (f.isEmpty()) return false;
    if (!f.endsWith(".toml")) f += ".toml";
    return saveAsPath(f);
}

bool MainWindow::saveAsPath(const QString& file) {
    // Relative paths were relative to the old location (the home folder for a new project): make
    // them absolute so they keep pointing at the same folders and files. A new project's output
    // folder stays relative: it then goes next to the project file.
    const QString oldBase = baseDir();
    const int n = settings_->makePathsAbsolute(oldBase, path_.isEmpty() ? QStringList{"output_dir"} : QStringList{});
    if (n) appendLog(tr("%n relative path(s) made absolute for the new project location.", nullptr, n));
    setPath(QFileInfo(file).absoluteFilePath());
    return save();
}

bool MainWindow::maybeSave() {
    if (!modified_) return true;
    const auto r = QMessageBox::question(this, tr("Unsaved changes"), tr("Save the changes to the project?"),
                                         QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    if (r == QMessageBox::Save) return save();
    return r == QMessageBox::Discard;
}

void MainWindow::updateRecent(const QString& path) {
    QSettings st;
    QStringList recent = st.value("recent").toStringList();
    if (!path.isEmpty()) {
        recent.removeAll(path);
        recent.prepend(path);
        while (recent.size() > 10) recent.removeLast();
        st.setValue("recent", recent);
    }
    recentMenu_->clear();
    for (const QString& r : recent)
        recentMenu_->addAction(r, this, [this, r] { if (maybeSave()) openProject(r); });
    recentMenu_->setEnabled(!recent.isEmpty());
}

void MainWindow::closeEvent(QCloseEvent* e) {
    if (runThread_) {
        if (QMessageBox::question(this, tr("Run in progress"), tr("A run is in progress. Cancel it and quit?")) != QMessageBox::Yes) {
            e->ignore();
            return;
        }
        cancel_->store(true);
        runThread_->wait();
    }
    if (!maybeSave()) {
        e->ignore();
        return;
    }
    QSettings st;
    st.setValue("window/geometry", saveGeometry());
    st.setValue("window/state", saveState());
    e->accept();
}

// ---------------------------------------------------------------------------------------------
// validation

void MainWindow::revalidate() {
    revalidateTimer_->stop();
    problems_->clear();
    try {
        project_ = parse_project(settings_->toToml().toStdString(), baseDir().toStdString());
    } catch (const std::exception& e) {
        project_.reset();
        const QStringList lines = QString::fromStdString(e.what()).split('\n');
        for (const QString& l : lines.mid(lines.size() > 1 ? 1 : 0)) {
            auto* it = new QListWidgetItem(QIcon::fromTheme("dialog-error"), l.trimmed());
            it->setForeground(QColor(200, 40, 40));
            problems_->addItem(it);
        }
        runAct_->setEnabled(false);
        return;
    }
    RunOptions opt;
    opt.from = kStages[fromStage_->currentIndex()];
    opt.to = kStages[toStage_->currentIndex()];
    const auto bad = check_project_inputs(*project_, opt);
    for (const auto& b : bad) problems_->addItem(new QListWidgetItem(QIcon::fromTheme("dialog-warning"), QString::fromStdString(b)));
    if (bad.empty()) {
        auto* it = new QListWidgetItem(QIcon::fromTheme("dialog-ok"), tr("Ready to run %1 → %2").arg(fromStage_->currentText(), toStage_->currentText()));
        it->setForeground(QColor(30, 140, 60));
        problems_->addItem(it);
    }
    runAct_->setEnabled(!runThread_ && bad.empty());
    image_->setAperture(project_->reduce.aperture);
    syncStarsToImage();
    maybeReloadImage();
    reloadResults();
}

// ---------------------------------------------------------------------------------------------
// reference frame

void MainWindow::maybeReloadImage() {
    if (!project_) return;
    const Project& p = *project_;
    const auto& f = p.frames;
    auto stamp = [&](const std::string& path) {
        const QFileInfo fi(QString::fromStdString(path));
        return fi.exists() ? fi.lastModified().toString(Qt::ISODateWithMs) : QString("-");
    };
    const std::string dark = f.master_dark.empty() ? p.output("master_dark.fits") : p.resolve(f.master_dark);
    const std::string flat = f.master_flat.empty() ? p.output("master_flat.fits") : p.resolve(f.master_flat);
    const std::string bias = f.master_bias.empty() ? (f.bias_dir.empty() ? std::string() : p.output("master_bias.fits")) : p.resolve(f.master_bias);
    const QString key = QStringList{QString::fromStdString(f.backend), QString::fromStdString(p.resolve(f.lights_dir)),
                                    QString::fromStdString(f.lights_prefix), QString::number(p.reduce.ref_index),
                                    stamp(dark), stamp(flat), stamp(bias)}.join('|');
    if (key == imageKey_) return;
    imageKey_ = key;
    const int ticket = ++imageLoads_;
    const float hot = f.hot_adu;
    auto* watcher = new QFutureWatcher<LoadedFrame>(this);
    connect(watcher, &QFutureWatcher<LoadedFrame>::finished, this, [this, watcher, ticket] {
        watcher->deleteLater();
        if (ticket != imageLoads_) return;   // a newer load was started
        LoadedFrame r = watcher->result();
        if (!r.error.isEmpty()) image_->clearImage(r.error);
        else image_->setImage(std::move(r.img), r.note);
        syncStarsToImage();
    });
    image_->clearImage(tr("Loading the reference frame…"));
    watcher->setFuture(QtConcurrent::run([p, dark, flat, bias, hot]() -> LoadedFrame {
        LoadedFrame r;
        try {
            const auto& fr = p.frames;
            const auto lights = find_light_frames(p.resolve(fr.lights_dir), fr.lights_prefix);
            if (lights.empty()) {
                r.error = QObject::tr("No light frames named %1* in\n%2").arg(QString::fromStdString(fr.lights_prefix),
                                                                             QString::fromStdString(p.resolve(fr.lights_dir)));
                return r;
            }
            const size_t idx = std::min(p.reduce.ref_index, lights.size() - 1);
            auto rd = make_reader(fr.backend);
            rd->open(lights[idx]);
            Image raw = rd->pixels(0);
            const QString name = QFileInfo(QString::fromStdString(lights[idx])).fileName();
            if (QFileInfo::exists(QString::fromStdString(dark)) && QFileInfo::exists(QString::fromStdString(flat))) {
                const Calibration cal = load_calibration(fr.backend, dark, flat,
                                                         QFileInfo::exists(QString::fromStdString(bias)) ? bias : std::string(), hot);
                r.img = cal.apply(raw);
                r.note = QObject::tr("frame %1 (%2), calibrated").arg(idx).arg(name);
            } else {
                r.img = std::move(raw);
                r.note = QObject::tr("frame %1 (%2), raw: master frames not built yet").arg(idx).arg(name);
            }
        } catch (const std::exception& e) {
            r.error = QString::fromStdString(e.what());
        }
        return r;
    }));
}

StarList MainWindow::currentStars() const {
    if (!project_) return settings_->stars();
    if (!project_->stars_file.empty()) {
        try {
            return read_star_list(project_->resolve(project_->stars_file));
        } catch (...) {
            return {};
        }
    }
    return project_->reduce.stars;
}

void MainWindow::syncStarsToImage() {
    if (!project_) return;
    QSet<size_t> ex;
    for (size_t i : project_->lightcurve.exclude_stars) ex.insert(i);
    image_->setStars(currentStars(), project_->reduce.target_star, ex);
}

void MainWindow::editStars(const StarList& stars) {
    // an edited list is kept in the project file itself
    if (!settings_->value("reduce", "stars_file").isEmpty()) {
        settings_->setValue("reduce", "stars_file", QString());
        appendLog(tr("The star list is now kept in the project file (stars_file cleared)."));
    }
    settings_->setStars(stars);
    revalidate();
}

QList<size_t> MainWindow::indexList(const QString& section, const QString& key) const {
    QList<size_t> v;
    for (const QString& s : settings_->value(section, key).split(QRegularExpression("[,\\s]+"), Qt::SkipEmptyParts)) v << s.toULongLong();
    return v;
}

void MainWindow::setIndexList(const QString& section, const QString& key, const QList<size_t>& v) {
    QStringList s;
    for (size_t i : v) s << QString::number(i);
    settings_->setValue(section, key, s.join(", "));
}

void MainWindow::addStar(double x, double y) {
    if (!project_ || image_->image().empty()) return;
    double cx = x, cy = y;
    try {
        centroid(image_->image(), cx, cy, project_->reduce.aperture);
    } catch (...) {
    }
    if (!std::isfinite(cx) || !std::isfinite(cy) || std::hypot(cx - x, cy - y) > project_->reduce.aperture.centroid_radius) {
        cx = x;
        cy = y;
    }
    StarList s = currentStars();
    for (size_t i = 0; i < s.size(); ++i)
        if (std::hypot(s[i].first - cx, s[i].second - cy) < 3) {
            statusBar()->showMessage(tr("Star %1 is already there").arg(i), 3000);
            return;
        }
    s.emplace_back(cx, cy);
    statusBar()->showMessage(tr("Added star %1 at x %2, y %3").arg(s.size() - 1).arg(cx, 0, 'f', 2).arg(cy, 0, 'f', 2), 4000);
    editStars(s);
}

void MainWindow::removeStar(size_t k) {
    if (!project_) return;
    StarList s = currentStars();
    if (k >= s.size() || k == project_->reduce.target_star) return;
    s.erase(s.begin() + long(k));
    // star indices after k move down by one
    auto remap = [&](const QString& sec, const QString& key) {
        QList<size_t> out;
        for (size_t i : indexList(sec, key)) {
            if (i == k) continue;
            out << (i > k ? i - 1 : i);
        }
        setIndexList(sec, key, out);
    };
    remap("lightcurve", "exclude_stars");
    remap("lightcurve", "comparison_stars");
    const size_t t = project_->reduce.target_star;
    if (t > k) settings_->setValue("reduce", "target_star", QString::number(t - 1));
    editStars(s);
}

void MainWindow::setTarget(size_t i) {
    settings_->setValue("reduce", "target_star", QString::number(i));
    QList<size_t> ex = indexList("lightcurve", "exclude_stars");
    ex.removeAll(i);
    setIndexList("lightcurve", "exclude_stars", ex);
    QList<size_t> comps = indexList("lightcurve", "comparison_stars");
    comps.removeAll(i);
    setIndexList("lightcurve", "comparison_stars", comps);
    revalidate();
}

void MainWindow::toggleExclude(size_t i) {
    QList<size_t> ex = indexList("lightcurve", "exclude_stars");
    if (ex.contains(i)) ex.removeAll(i);
    else {
        ex << i;
        std::sort(ex.begin(), ex.end());
    }
    setIndexList("lightcurve", "exclude_stars", ex);
    revalidate();
}

void MainWindow::detectStars() {
    if (!project_ || image_->image().empty()) return;
    const Image& img = image_->image();
    const double m = project_->reduce.detect_margin;
    StarList found;
    for (const auto& d : find_stars(img))
        if (d.x > m && d.x < img.nx - m && d.y > m && d.y < img.ny - m) found.emplace_back(d.x, d.y);
    if (found.empty()) {
        QMessageBox::information(this, tr("Detect stars"), tr("No stars found."));
        return;
    }
    const StarList old = currentStars();
    if (!old.empty() &&
        QMessageBox::question(this, tr("Detect stars"),
                              tr("Replace the %1 stars in the list by the %2 detected ones?\n"
                                 "The excluded and chosen comparison stars are cleared; the target becomes the "
                                 "detected star nearest to the current one.").arg(old.size()).arg(found.size())) != QMessageBox::Yes)
        return;
    size_t target = 0;
    if (project_->reduce.target_star < old.size()) {
        const auto& t = old[project_->reduce.target_star];
        double best = 1e300;
        for (size_t i = 0; i < found.size(); ++i) {
            const double d = std::hypot(found[i].first - t.first, found[i].second - t.second);
            if (d < best) { best = d; target = i; }
        }
    }
    settings_->setValue("reduce", "target_star", QString::number(target));
    settings_->setValue("lightcurve", "exclude_stars", QString());
    settings_->setValue("lightcurve", "comparison_stars", QString());
    editStars(found);
}

// ---------------------------------------------------------------------------------------------
// runs

void MainWindow::appendLog(const QString& line) { log_->appendPlainText(line); }

void MainWindow::startRun(Stage from, Stage to) {
    if (runThread_) return;
    revalidate();
    if (!project_) return;
    if ((path_.isEmpty() || modified_) && !save()) return;
    const Project p = *project_;
    RunOptions opt;
    opt.from = from;
    opt.to = to;
    cancel_ = std::make_shared<std::atomic<bool>>(false);
    auto cancel = cancel_;
    MainWindow* self = this;
    fitModeIndex_ = 0;
    lastDone_ = 0;
    lastStage_ = -1;
    log_->appendPlainText(QString("— %1: %2 → %3 —").arg(QDateTime::currentDateTime().toString("HH:mm:ss"),
                                                        describeStage(from), describeStage(to)));
    for (QObject* o : findChildren<QDockWidget*>())
        if (o->objectName() == "log") static_cast<QDockWidget*>(o)->raise();

    runThread_ = QThread::create([p, opt, cancel, self] {
        LineBuf buf([self](const std::string& l) {
            QMetaObject::invokeMethod(self, [self, l = QString::fromStdString(l)] { self->appendLog(l); }, Qt::QueuedConnection);
        });
        std::ostream os(&buf);
        RunOptions o = opt;
        o.log = &os;
        auto last = std::chrono::steady_clock::now() - std::chrono::seconds(1);
        int lastStage = -1;
        StageProgress prog = [&](Stage s, size_t done, size_t total) {
            if (cancel->load()) return false;
            const auto now = std::chrono::steady_clock::now();
            if (int(s) != lastStage || done == total || done <= 1 || now - last > std::chrono::milliseconds(60)) {
                last = now;
                lastStage = int(s);
                QMetaObject::invokeMethod(self, [self, s, done, total] { self->showProgress(int(s), done, total); }, Qt::QueuedConnection);
            }
            return true;
        };
        QString error;
        bool ok = false;
        try {
            run_project(p, o, prog);
            ok = true;
        } catch (const Cancelled&) {
            error = QObject::tr("cancelled");
        } catch (const std::exception& e) {
            error = QString::fromStdString(e.what());
        }
        os.flush();
        QMetaObject::invokeMethod(self, [self, ok, error] { self->runFinished(ok, error); }, Qt::QueuedConnection);
    });
    connect(runThread_, &QThread::finished, this, [this] {
        runThread_->deleteLater();
        runThread_ = nullptr;
        cancelAct_->setEnabled(false);
        revalidate();
    });
    runAct_->setEnabled(false);
    cancelAct_->setEnabled(true);
    progress_->setRange(0, 0);
    progressAct_->setVisible(true);
    stageLabel_->setText(tr("starting…"));
    runThread_->start();
}

void MainWindow::cancelRun() {
    if (!runThread_) return;
    cancel_->store(true);
    stageLabel_->setText(tr("cancelling…"));
}

void MainWindow::showProgress(int stage, qulonglong done, qulonglong total) {
    const Stage s = Stage(stage);
    if (s == Stage::Fit && stage == lastStage_ && done < lastDone_) ++fitModeIndex_;
    lastStage_ = stage;
    lastDone_ = done;
    QString label = describeStage(s);
    if (s == Stage::Fit && project_ && fitModeIndex_ < int(project_->fit_modes.size()))
        label += project_->fit_modes[size_t(fitModeIndex_)] == FitConfig::Mode::Free ? tr(" (free)") : tr(" (constrained)");
    stageLabel_->setText(label);
    progress_->setRange(0, int(std::min<qulonglong>(total, 1u << 30)));
    progress_->setValue(int(done));
    progress_->setFormat(QString("%v / %m"));
}

void MainWindow::runFinished(bool ok, const QString& message) {
    progressAct_->setVisible(false);
    stageLabel_->setText(ok ? tr("finished") : message == tr("cancelled") ? tr("cancelled") : tr("failed"));
    if (!ok && message != tr("cancelled")) {
        appendLog(tr("ERROR: %1").arg(message));
        QMessageBox::warning(this, tr("Run failed"), message);
    } else if (!ok) {
        appendLog(tr("Cancelled."));
    } else {
        appendLog(tr("Done."));
        tabs_->setCurrentWidget(results_);
    }
    imageKey_.clear();   // masters may have changed
    reloadResults();
    emit runDone(ok);
}

bool MainWindow::exportLightCurve(int target, const QString& fileForTest) {
    if (!project_) return false;
    const ExportTarget t = ExportTarget(target);
    const QString name = export_target_name(t);
    QString file = fileForTest;
    if (file.isEmpty()) {
        file = QFileDialog::getSaveFileName(this, tr("Export for %1").arg(name),
                                            QString::fromStdString(default_export_path(*project_, t)), tr("Text files (*.txt)"));
        if (file.isEmpty()) return false;
    }
    ExportResult r;
    try {
        r = export_light_curve(*project_, t, file.toStdString());
    } catch (const std::exception& e) {
        if (fileForTest.isEmpty()) QMessageBox::warning(this, tr("Export for %1").arg(name), QString::fromStdString(e.what()));
        return false;
    }
    appendLog(tr("%1 export: %2 points -> %3").arg(name).arg(r.points).arg(QString::fromStdString(r.path)));
    if (!fileForTest.isEmpty()) return true;
    QMessageBox box(this);
    box.setWindowTitle(tr("Export for %1").arg(name));
    box.setIcon(QMessageBox::Information);
    box.setText(tr("Wrote %1 points to\n%2").arg(r.points).arg(QString::fromStdString(r.path)));
    box.setInformativeText(QString::fromStdString(r.instructions));
    box.setTextInteractionFlags(Qt::TextSelectableByMouse);
    QPushButton* open = box.addButton(tr("Open the upload page"), QMessageBox::ActionRole);
    box.addButton(QMessageBox::Close);
    box.exec();
    if (box.clickedButton() == open)
        QDesktopServices::openUrl(QUrl(t == ExportTarget::ExoClock ? "https://www.exoclock.space/upload/" : "https://var.astro.cz/ETD"));
    return true;
}

void MainWindow::reloadResults() {
    if (!project_) return;
    const QString out = QString::fromStdString(project_->resolve(project_->output_dir));
    results_->load(out, project_->lightcurve.eph.T14, project_->lightcurve.target);
    updateComparisonMarks(out);
}

void MainWindow::updateComparisonMarks(const QString& outDir) {
    // comparison stars and weights from lc.csv's header; star indices refer to the star list of the
    // run that wrote phot.csv, so they are shown only while the list has the same number of stars
    QMap<size_t, double> marks;
    QFile lc(QDir(outDir).filePath("lc.csv"));
    QStringList comps, weights;
    if (lc.open(QIODevice::ReadOnly | QIODevice::Text)) {
        while (!lc.atEnd()) {
            const QString line = QString::fromUtf8(lc.readLine()).trimmed();
            if (!line.startsWith('#')) break;
            if (line.startsWith("# comparison_stars=")) comps = line.mid(19).split(';', Qt::SkipEmptyParts);
            if (line.startsWith("# weights=")) weights = line.mid(10).split(';', Qt::SkipEmptyParts);
        }
    }
    long nsPhot = -1;
    QFile ph(QDir(outDir).filePath("phot.csv"));
    if (ph.open(QIODevice::ReadOnly) && ph.size() > 0) {
        ph.seek(std::max<qint64>(0, ph.size() - 4096));
        const QStringList lines = QString::fromUtf8(ph.readAll()).split('\n', Qt::SkipEmptyParts);
        if (!lines.isEmpty()) {
            const QStringList f = lines.last().split(',');
            if (f.size() > 7) nsPhot = f[7].toLong() + 1;
        }
    }
    const bool sameList = nsPhot == long(currentStars().size());
    if (sameList && comps.size() == weights.size())
        for (int i = 0; i < comps.size(); ++i) marks[comps[i].toULongLong()] = weights[i].toDouble();
    image_->setComparisons(marks);
    if (!comps.isEmpty() && !sameList)
        statusBar()->showMessage(tr("Comparison stars of the last run not shown: the star list changed since."), 6000);
}

void MainWindow::runForTest(const QString& screenshot, double cancelAfterSeconds) {
    connect(this, &MainWindow::runDone, this, [this, screenshot] {
        QTimer::singleShot(800, this, [this, screenshot] {
            grab().save(screenshot);
            std::fprintf(stderr, "%s\n", qPrintable(log_->toPlainText()));
            QApplication::exit(stageLabel_->text() == tr("finished") ? 0 : 3);
        });
    });
    QTimer::singleShot(500, this, [this, cancelAfterSeconds, screenshot] {
        runAct_->trigger();
        if (cancelAfterSeconds > 0)
            QTimer::singleShot(int(cancelAfterSeconds * 1000), this, [this, screenshot] {
                QString running = screenshot;
                running.insert(running.lastIndexOf('.'), "-running");
                grab().save(running);
                cancelRun();
            });
    });
}

void MainWindow::screenshotAndQuit(const QString& file, int tab) {
    tabs_->setCurrentIndex(std::min(tab, 1));   // 2: the results' diagnostics
    results_->showTab(tab == 2 ? 1 : 0);
    auto* timer = new QTimer(this);
    auto start = std::make_shared<QDateTime>(QDateTime::currentDateTime());
    connect(timer, &QTimer::timeout, this, [this, timer, file, start] {
        const bool imageReady = !image_->image().empty() || start->msecsTo(QDateTime::currentDateTime()) > 20000;
        if (!imageReady || start->msecsTo(QDateTime::currentDateTime()) < 1500) return;
        timer->stop();
        grab().save(file);
        QApplication::quit();
    });
    timer->start(200);
}
