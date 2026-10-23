// SPDX-License-Identifier: GPL-3.0-or-later
// transitscope-gui [project.toml] [--screenshot FILE [--tab N]]
#include <QApplication>
#include <QCommandLineParser>

#include <clocale>

#include "MainWindow.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    // QApplication adopts the user's locale for the C library; with a decimal comma (nl_BE, fr_BE,
    // de_DE...) strtod/stod would stop at the '.' of every number in the CSV files and FITS headers.
    std::setlocale(LC_NUMERIC, "C");
    QApplication::setApplicationName("TransitScope");
    QApplication::setOrganizationName("TransitScope");
    QApplication::setApplicationVersion(TRANSITSCOPE_VERSION);

    QCommandLineParser cli;
    cli.setApplicationDescription("Exoplanet transit photometry");
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument("project", "Project file to open (.toml)");
    QCommandLineOption shot("screenshot", "Save a screenshot of the window to FILE and quit (testing).", "FILE");
    QCommandLineOption tab("tab", "Tab in front for --screenshot (0 reference frame, 1 results).", "N", "0");
    QCommandLineOption run("test-run", "Run the project, save a screenshot to FILE when done and quit (testing).", "FILE");
    QCommandLineOption cancel("test-cancel-after", "With --test-run: cancel after SECONDS.", "SECONDS", "0");
    cli.addOption(shot);
    cli.addOption(tab);
    cli.addOption(run);
    cli.addOption(cancel);
    cli.process(app);

    MainWindow w;
    w.show();
    if (!cli.positionalArguments().isEmpty()) w.openProject(cli.positionalArguments().first());
    if (cli.isSet(shot)) w.screenshotAndQuit(cli.value(shot), cli.value(tab).toInt());
    else if (cli.isSet(run)) w.runForTest(cli.value(run), cli.value(cancel).toDouble());
    return app.exec();
}
