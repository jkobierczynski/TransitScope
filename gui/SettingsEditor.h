// SPDX-License-Identifier: GPL-3.0-or-later
// Tree of every project setting: [section] -> key = value, with the template's comments as
// tooltips. The list of keys and their types comes from project_template(), so the editor follows
// the file format without a second copy of it. Values are edited as plain text (strings without
// quotes, arrays as "8, 10, 12"); toToml() turns the tree back into a project file.
#pragma once
#include <QMap>
#include <QString>
#include <QTreeWidget>

#include "transit/reduce.h"

class SettingsEditor : public QTreeWidget {
    Q_OBJECT
public:
    enum class Kind { String, Angle, Number, Integer, Numbers, Integers, Strings, Stars };

    explicit SettingsEditor(QWidget* parent = nullptr);

    /// Fill from a project's TOML text (project_to_toml); keys absent there are left empty.
    void load(const QString& tomlText, const transit::StarList& stars);
    /// The tree as project file text; empty values are left out (their defaults apply).
    QString toToml() const;

    QString value(const QString& section, const QString& key) const;
    void setValue(const QString& section, const QString& key, const QString& text);   // emits changed()
    const transit::StarList& stars() const { return stars_; }
    void setStars(const transit::StarList& s);                                         // emits changed()

    /// Base directory for the browse buttons of path keys.
    void setBaseDir(const QString& dir) { baseDir_ = dir; }

signals:
    void changed();

private:
    struct Key {
        QString section, name, comment;
        Kind kind;
        bool isPath = false, isDir = false;
        QTreeWidgetItem* item = nullptr;
    };
    void buildSchema();
    void browse(Key& k);
    Key* find(const QString& section, const QString& key);
    const Key* find(const QString& section, const QString& key) const;
    QString starsText() const;

    QList<Key> keys_;
    QMap<QString, QTreeWidgetItem*> sections_;
    transit::StarList stars_;
    QString baseDir_;
    bool loading_ = false;
};
