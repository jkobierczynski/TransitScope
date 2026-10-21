// SPDX-License-Identifier: GPL-3.0-or-later
// TransitScope main window: the project's settings (left), the reference frame with the star list
// and the results (tabs), problems and the run log (bottom). Runs happen in a worker thread.
#pragma once
#include <QMainWindow>
#include <QSet>

#include <atomic>
#include <memory>
#include <optional>

#include "transit/project.h"

class ImageView;
class QAction;
class QComboBox;
class QLabel;
class QListWidget;
class QMenu;
class QPlainTextEdit;
class QProgressBar;
class QTabWidget;
class QThread;
class QTimer;
class ResultsView;
class SettingsEditor;

class MainWindow : public QMainWindow {
    Q_OBJECT
    friend class GuiTest;
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;
    bool openProject(const QString& path);
    /// For tests and documentation: write a screenshot of the window once the reference frame is
    /// shown (or after a timeout), with tab `tab` in front, then quit.
    void screenshotAndQuit(const QString& file, int tab);
    /// For tests: run with the toolbar's stages, optionally cancel after some seconds, save a
    /// screenshot when the run has ended, print the log and quit (exit status 0 if it finished).
    void runForTest(const QString& screenshot, double cancelAfterSeconds);

signals:
    void runDone(bool ok);

protected:
    void closeEvent(QCloseEvent* e) override;

private:
    // project file
    void newProject();
    void openDialog();
    bool save();
    bool saveAs();
    bool maybeSave();
    void setPath(const QString& path);
    void setModified(bool m);
    void updateRecent(const QString& path);
    QString baseDir() const;

    // settings -> project
    void revalidate();
    void loadSettingsFromProject(const transit::Project& p);

    // reference frame and stars
    void maybeReloadImage();
    void syncStarsToImage();
    transit::StarList currentStars() const;
    void editStars(const transit::StarList& stars);
    void addStar(double x, double y);
    void removeStar(size_t i);
    void setTarget(size_t i);
    void toggleExclude(size_t i);
    void detectStars();
    QList<size_t> indexList(const QString& section, const QString& key) const;
    void setIndexList(const QString& section, const QString& key, const QList<size_t>& v);

    // runs
    void startRun(transit::Stage from, transit::Stage to);
    void cancelRun();
    void runFinished(bool ok, const QString& message);
    void appendLog(const QString& line);
    void showProgress(int stage, qulonglong done, qulonglong total);
    void reloadResults();

    SettingsEditor* settings_;
    ImageView* image_;
    ResultsView* results_;
    QTabWidget* tabs_;
    QListWidget* problems_;
    QPlainTextEdit* log_;
    QProgressBar* progress_;
    QLabel* stageLabel_;
    QComboBox* fromStage_;
    QComboBox* toStage_;
    QAction* runAct_;
    QAction* cancelAct_;
    QAction* progressAct_;
    QMenu* recentMenu_;
    QTimer* revalidateTimer_;

    QString path_;
    bool modified_ = false;
    std::optional<transit::Project> project_;   // the last valid state of the settings
    QString imageKey_;                          // what the shown reference frame was loaded from
    int imageLoads_ = 0;                        // ignore results of superseded loads

    QThread* runThread_ = nullptr;
    std::shared_ptr<std::atomic<bool>> cancel_;
    int fitModeIndex_ = 0;
    qulonglong lastDone_ = 0;
    int lastStage_ = -1;
};
