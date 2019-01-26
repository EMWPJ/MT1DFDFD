#-------------------------------------------------
#
# Project created by QtCreator 2017-06-21T15:18:18
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MT1DFDFD
TEMPLATE = app


SOURCES += main.cpp\
        widget.cpp \
    mt1d.cpp

HEADERS  += widget.h \
    mt1d.h

FORMS    += widget.ui

CONFIG += C++11
