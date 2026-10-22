// SPDX-License-Identifier: GPL-3.0-or-later
// GUI checks (run offscreen: QT_QPA_PLATFORM=offscreen):
//   - a project opens, validates and shows its reference frame with the star list
//   - clicking on a star-free spot adds a centroided star; clicking near a star does not duplicate it
//   - set target / exclude / remove keep the indices in exclude_stars, comparison_stars and
//     target_star consistent
//   - editing a value in the settings tree revalidates (an invalid value shows a problem)
//   - save writes a project that loads back with the edited star list
// usage: gui_test PROJECT.toml   (a project whose reference frame and star list exist)
#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QListWidget>
#include <QTemporaryDir>
#include <QTest>

#include <cstdio>

#include "../../gui/ImageView.h"
#include "../../gui/MainWindow.h"
#include "../../gui/SettingsEditor.h"

static int failures = 0;
static void check(bool ok, const QString& what) {
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", qPrintable(what));
    if (!ok) ++failures;
}

class GuiTest {
public:
    static int run(const QString& project) {
        QTemporaryDir tmp;
        // work on a copy, so saving does not touch the original
        const QString copy = tmp.filePath("p.toml");
        QFile::copy(project, copy);
        MainWindow w;
        w.show();
        check(w.openProject(copy), "project opens");
        check(w.project_.has_value(), "settings are valid");
        (void)QTest::qWaitFor([&] { return !w.image_->image().empty(); }, 20000);
        check(!w.image_->image().empty(), "reference frame shown");
        const auto stars0 = w.currentStars();
        check(stars0.size() == 37, QString("37 stars from the stars file (%1)").arg(stars0.size()));

        // click on the sky far from any star: adds one (the editor keeps it inline from now on)
        const auto& img = w.image_->image();
        double bx = -1, by = -1;
        for (int y = 100; y < img.ny - 100 && bx < 0; y += 37)
            for (int x = 100; x < img.nx - 100; x += 41) {
                bool clear = true;
                for (const auto& s : stars0) clear &= std::hypot(s.first - x, s.second - y) > 60;
                if (clear) { bx = x; by = y; break; }
            }
        w.image_->fitInView();
        QApplication::processEvents();
        const QPoint vp = w.image_->mapFromScene(QPointF(bx + 0.5, img.ny - 0.5 - by));
        QTest::mouseClick(w.image_->viewport(), Qt::LeftButton, {}, vp);
        QApplication::processEvents();
        auto stars1 = w.currentStars();
        check(stars1.size() == 38, "a click on the sky adds a star");
        check(w.settings_->value("reduce", "stars_file").isEmpty(), "the edited list moves into the project (stars_file cleared)");
        // a click on an existing star selects it instead
        const QPoint onStar = w.image_->mapFromScene(QPointF(stars0[5].first + 0.5, img.ny - 0.5 - stars0[5].second));
        QTest::mouseClick(w.image_->viewport(), Qt::LeftButton, {}, onStar);
        QApplication::processEvents();
        check(w.currentStars().size() == 38, "a click on a star does not add a duplicate");

        // index bookkeeping
        w.settings_->setValue("lightcurve", "exclude_stars", "10, 32, 36");
        w.settings_->setValue("lightcurve", "comparison_stars", "1, 2, 33");
        w.revalidate();
        emit w.image_->starRemoveRequested(20);
        QApplication::processEvents();
        check(w.settings_->value("lightcurve", "exclude_stars") == "10, 31, 35", "removing star 20 shifts exclude_stars");
        check(w.settings_->value("lightcurve", "comparison_stars") == "1, 2, 32", "... and comparison_stars");
        check(w.currentStars().size() == 37, "the star is gone");
        emit w.image_->targetRequested(1);
        QApplication::processEvents();
        check(w.settings_->value("reduce", "target_star") == "1", "set target");
        check(w.settings_->value("lightcurve", "comparison_stars") == "2, 32", "the new target is no comparison star");
        emit w.image_->excludeToggled(4);
        QApplication::processEvents();
        check(w.settings_->value("lightcurve", "exclude_stars") == "4, 10, 31, 35", "exclude toggles on (sorted)");
        emit w.image_->excludeToggled(4);
        QApplication::processEvents();
        check(w.settings_->value("lightcurve", "exclude_stars") == "10, 31, 35", "... and off");
        emit w.image_->starRemoveRequested(0);
        QApplication::processEvents();
        check(w.settings_->value("reduce", "target_star") == "0" && w.currentStars().size() == 36,
              "removing a star before the target moves the target index");
        check(w.settings_->value("lightcurve", "exclude_stars") == "9, 30, 34", "... and the excluded indices");

        // export (needs the original project's lc.csv and phot.csv)
        {
            const transit::Project orig = transit::load_project(project.toStdString());
            w.settings_->setValue("project", "output_dir", QString::fromStdString(orig.resolve(orig.output_dir)));
            w.revalidate();
        }
        const QString exo = QDir(QString::fromStdString(w.project_->resolve(w.project_->output_dir))).exists("lc.csv")
                                ? QString() : QStringLiteral("skip");
        if (exo.isEmpty()) {
            const QString f = tmp.filePath("exo.txt");
            check(w.exportLightCurve(0, f), "export for ExoClock");
            QFile ef(f);
            ef.open(QIODevice::ReadOnly);
            const QString text = QString::fromUtf8(ef.readAll());
            check(text.startsWith("#Planet: ") && text.contains("#Time format: JD_UTC\n") && text.count('\n') > 100,
                  "the ExoClock file has the header and the points");
        }

        // validation of edits
        w.settings_->setValue("site", "lat", "north");
        w.revalidate();
        bool flagged = false;
        for (int i = 0; i < w.problems_->count(); ++i) flagged |= w.problems_->item(i)->text().contains("[site] lat must be a number");
        check(flagged && !w.project_.has_value(), "an invalid value is reported in Problems");
        w.settings_->setValue("site", "lat", "50.93");
        w.revalidate();
        check(w.project_.has_value(), "fixing it makes the settings valid again");

        // save and load back
        const auto starsBefore = w.currentStars();
        check(w.save(), "save");
        const transit::Project back = transit::load_project(copy.toStdString());
        check(back.reduce.stars == starsBefore && back.reduce.target_star == 0 &&
                  back.lightcurve.exclude_stars == std::vector<size_t>{9, 30, 34},
              "the saved project holds the edited stars, target and exclusions");
        // save as elsewhere: relative paths keep pointing at the same place
        const QString oldBase = QFileInfo(copy).absolutePath();
        w.settings_->setValue("reduce", "stars_file", "");
        w.settings_->setValue("frames", "master_flat", "sub/../masters/flat.fits");
        QDir().mkpath(tmp.filePath("elsewhere/deeper"));
        const QString moved = tmp.filePath("elsewhere/deeper/moved.toml");
        w.saveAsPath(moved);
        check(w.settings_->value("frames", "master_flat") == QDir::cleanPath(oldBase + "/masters/flat.fits"),
              "save as: a relative path becomes absolute against the old location (" + w.settings_->value("frames", "master_flat") + ")");
        const transit::Project mp = transit::load_project(moved.toStdString());
        check(QString::fromStdString(mp.resolve(mp.frames.master_flat)) == QDir::cleanPath(oldBase + "/masters/flat.fits"),
              "... and resolves to the same file from the new location");
        w.modified_ = false;
        std::printf("%d failure(s)\n", failures);
        return failures ? 1 : 0;
    }
};

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    QApplication::setOrganizationName("transitscope-test");
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s PROJECT.toml\n", argv[0]);
        return 2;
    }
    return GuiTest::run(QString::fromLocal8Bit(argv[1]));
}
