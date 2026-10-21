// SPDX-License-Identifier: GPL-3.0-or-later
#include "SettingsEditor.h"

#include <QDir>
#include <QFileDialog>
#include <QHeaderView>
#include <QRegularExpression>
#include <QToolButton>

#include <charconv>
#include <cmath>

#define TOML_EXCEPTIONS 1
#include <toml.hpp>

#include "transit/project.h"

namespace {

// split "value   # comment" at the first '#' outside a string
std::pair<QString, QString> splitComment(const QString& s) {
    bool inStr = false;
    QChar quote;
    for (int i = 0; i < s.size(); ++i) {
        const QChar c = s[i];
        if (inStr) {
            if (c == '\\') { ++i; continue; }
            if (c == quote) inStr = false;
        } else if (c == '"' || c == '\'') {
            inStr = true;
            quote = c;
        } else if (c == '#') {
            return {s.left(i).trimmed(), s.mid(i + 1).trimmed()};
        }
    }
    return {s.trimmed(), QString()};
}

QString quoted(const QString& s) {
    QString o = "\"";
    for (QChar c : s) {
        if (c == '"') o += "\\\"";
        else if (c == '\\') o += "\\\\";
        else if (c == '\n') o += "\\n";
        else if (c == '\t') o += "\\t";
        else o += c;
    }
    return o + "\"";
}

// a TOML number as typed, or the text as a string so the project reader names the problem
QString numberOrQuoted(const QString& v) {
    try {
        const toml::table t = toml::parse("x = " + v.toStdString());
        if (const toml::node* n = t.get("x"); n && n->is_number()) return v;
    } catch (...) {
    }
    return quoted(v);
}

bool isNumber(const QString& s) {
    bool ok = false;
    s.trimmed().toDouble(&ok);
    return ok;
}

QString numberText(double v) {
    char buf[64];
    auto r = std::to_chars(buf, buf + sizeof buf, v);
    return QString::fromLatin1(buf, int(r.ptr - buf));
}

QString nodeText(const toml::node& n) {
    if (auto s = n.value<std::string>(); s && n.is_string()) return QString::fromStdString(*s);
    if (n.is_integer()) return QString::number(*n.value<int64_t>());
    if (n.is_floating_point()) return numberText(*n.value<double>());
    if (const toml::array* a = n.as_array()) {
        QStringList parts;
        for (const auto& el : *a) parts << nodeText(el);
        return parts.join(", ");
    }
    if (n.is_boolean()) return *n.value<bool>() ? "true" : "false";
    return {};
}

}  // namespace

SettingsEditor::SettingsEditor(QWidget* parent) : QTreeWidget(parent) {
    setColumnCount(3);
    setHeaderLabels({tr("Setting"), tr("Value"), QString()});
    header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    header()->setSectionResizeMode(1, QHeaderView::Stretch);
    header()->setSectionResizeMode(2, QHeaderView::Fixed);
    header()->resizeSection(2, 28);
    header()->setStretchLastSection(false);
    setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed | QAbstractItemView::SelectedClicked);
    setAlternatingRowColors(true);
    buildSchema();
    connect(this, &QTreeWidget::itemChanged, this, [this](QTreeWidgetItem* it, int col) {
        if (loading_ || col != 1 || !it->parent()) return;
        emit changed();
    });
}

void SettingsEditor::buildSchema() {
    // keys, kinds and comments from the commented template
    const QString tpl = QString::fromStdString(transit::project_template());
    const toml::table parsed = toml::parse(transit::project_template());
    static const QRegularExpression sectionRe(R"(^\s*\[([A-Za-z_]+)\])");
    static const QRegularExpression keyRe(R"(^\s*(#\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$)");
    QString section;
    for (const QString& line : tpl.split('\n')) {
        if (auto m = sectionRe.match(line); m.hasMatch()) {
            section = m.captured(1);
            auto* top = new QTreeWidgetItem(this, {section});
            QFont f = top->font(0);
            f.setBold(true);
            top->setFont(0, f);
            top->setFirstColumnSpanned(true);
            top->setFlags(Qt::ItemIsEnabled);
            top->setExpanded(true);
            sections_[section] = top;
            continue;
        }
        auto m = keyRe.match(line);
        if (section.isEmpty() || !m.hasMatch()) continue;
        Key k;
        k.section = section;
        k.name = m.captured(2);
        const auto [valueText, comment] = splitComment(m.captured(3));
        k.comment = comment;
        if (m.capturedLength(1)) k.comment = (comment.isEmpty() ? QString() : comment + "\n") + tr("Optional: commented out in the template.");
        if (find(section, k.name)) continue;
        // kind from the example value
        toml::table one;
        try {
            one = toml::parse("x = " + valueText.toStdString());
        } catch (...) {
            continue;
        }
        const toml::node& n = *one.get("x");
        if (k.name == "stars") k.kind = Kind::Stars;
        else if (k.name == "ra" || k.name == "dec") k.kind = Kind::Angle;
        else if (n.is_string()) k.kind = Kind::String;
        else if (n.is_integer()) k.kind = Kind::Integer;
        else if (n.is_floating_point()) k.kind = Kind::Number;
        else if (const toml::array* a = n.as_array()) {
            if (!a->empty() && (*a)[0].is_string()) k.kind = Kind::Strings;
            else if (k.name.endsWith("stars")) k.kind = Kind::Integers;
            else k.kind = Kind::Numbers;
        } else continue;
        // integers written as floats in the template stay floats; arrays named *_stars hold indices
        k.isPath = k.name.endsWith("_dir") || k.name.startsWith("master_") || k.name == "stars_file" || k.name == "output_dir";
        k.isDir = k.name.endsWith("_dir");

        auto* it = new QTreeWidgetItem(sections_[section], {k.name, QString()});
        it->setToolTip(0, k.comment);
        it->setToolTip(1, k.comment);
        if (k.kind == Kind::Stars) {
            it->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
            it->setForeground(1, palette().brush(QPalette::Disabled, QPalette::Text));
        } else {
            it->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
        }
        k.item = it;
        keys_.push_back(k);
        if (k.isPath) {
            auto* b = new QToolButton(this);
            b->setText("…");
            b->setToolTip(k.isDir ? tr("Choose a folder") : tr("Choose a file"));
            const int idx = keys_.size() - 1;
            connect(b, &QToolButton::clicked, this, [this, idx] { browse(keys_[idx]); });
            setItemWidget(it, 2, b);
        }
    }
    (void)parsed;
}

SettingsEditor::Key* SettingsEditor::find(const QString& section, const QString& key) {
    for (auto& k : keys_) if (k.section == section && k.name == key) return &k;
    return nullptr;
}
const SettingsEditor::Key* SettingsEditor::find(const QString& section, const QString& key) const {
    for (const auto& k : keys_) if (k.section == section && k.name == key) return &k;
    return nullptr;
}

QString SettingsEditor::starsText() const {
    if (stars_.empty() && !value("reduce", "stars_file").isEmpty()) return tr("from stars_file");
    return stars_.empty() ? tr("none (detected on the reference frame)")
                          : tr("%n star(s): edit on the reference frame", nullptr, int(stars_.size()));
}

void SettingsEditor::load(const QString& tomlText, const transit::StarList& stars) {
    loading_ = true;
    stars_ = stars;
    toml::table t;
    try {
        t = toml::parse(tomlText.toStdString());
    } catch (...) {
    }
    for (auto& k : keys_) {
        QString text;
        if (k.kind == Kind::Stars) {
            text = starsText();
        } else if (const toml::node* sec = t.get(k.section.toStdString())) {
            if (const toml::table* st = sec->as_table())
                if (const toml::node* n = st->get(k.name.toStdString())) text = nodeText(*n);
        }
        k.item->setText(1, text);
    }
    loading_ = false;
}

QString SettingsEditor::value(const QString& section, const QString& key) const {
    const Key* k = find(section, key);
    return k ? k->item->text(1).trimmed() : QString();
}

void SettingsEditor::setValue(const QString& section, const QString& key, const QString& text) {
    if (Key* k = find(section, key)) {
        if (k->item->text(1) == text) return;
        loading_ = true;
        k->item->setText(1, text);
        loading_ = false;
        emit changed();
    }
}

void SettingsEditor::setStars(const transit::StarList& s) {
    stars_ = s;
    if (Key* k = find("reduce", "stars")) {
        loading_ = true;
        k->item->setText(1, starsText());
        loading_ = false;
    }
    emit changed();
}

QString SettingsEditor::toToml() const {
    QString out;
    QString section;
    for (const auto& k : keys_) {
        if (k.section != section) {
            section = k.section;
            out += "\n[" + section + "]\n";
        }
        if (k.kind == Kind::Stars) {
            if (stars_.empty()) continue;
            out += "stars = [\n";
            for (const auto& [x, y] : stars_) out += "  [" + numberText(x) + ", " + numberText(y) + "],\n";
            out += "]\n";
            continue;
        }
        const QString v = k.item->text(1).trimmed();
        if (v.isEmpty() && k.kind != Kind::Integers && k.kind != Kind::Numbers) continue;
        QString val;
        switch (k.kind) {
            case Kind::String: val = quoted(v); break;
            case Kind::Angle: val = isNumber(v) ? v : quoted(v); break;
            case Kind::Number:
            case Kind::Integer: val = numberOrQuoted(v); break;
            case Kind::Numbers:
            case Kind::Integers: {
                if (v.isEmpty()) {
                    // an empty list is a value of its own for these ([] = automatic / all)
                    if (k.name == "radii" && k.section == "reduce") continue;
                    val = "[]";
                    break;
                }
                QStringList parts;
                for (const QString& p : v.split(QRegularExpression("[,\\s]+"), Qt::SkipEmptyParts)) parts << numberOrQuoted(p);
                val = "[" + parts.join(", ") + "]";
                break;
            }
            case Kind::Strings: {
                QStringList parts;
                for (const QString& p : v.split(',', Qt::SkipEmptyParts)) parts << quoted(p.trimmed());
                val = "[" + parts.join(", ") + "]";
                break;
            }
            case Kind::Stars: break;
        }
        out += k.name + " = " + val + "\n";
    }
    return out;
}

void SettingsEditor::browse(Key& k) {
    const QString current = k.item->text(1).trimmed();
    QString start = current.isEmpty() ? baseDir_ : QDir(baseDir_).absoluteFilePath(current);
    if (start.startsWith("~/")) start = QDir::homePath() + start.mid(1);
    QString chosen = k.isDir ? QFileDialog::getExistingDirectory(this, k.name, start)
                             : QFileDialog::getOpenFileName(this, k.name, start);
    if (chosen.isEmpty()) return;
    // keep paths inside the project's folder relative, so the project can move with its data
    if (!baseDir_.isEmpty()) {
        const QString rel = QDir(baseDir_).relativeFilePath(chosen);
        if (!rel.startsWith("../../")) chosen = rel;
    }
    setValue(k.section, k.name, chosen);
}
