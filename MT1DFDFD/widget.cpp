﻿#include "widget.h"
#include "ui_widget.h"

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    QVector < double > dz;
    QVector < double > sigma;
    QVector < double > fre;
    dz.clear();
    sigma.clear();
    fre.clear();
    dz<<50000<<100000<<200000<<400000<<800000
     <<7500<<10000<<15000<<20000<<30000
    <<5000<<5000<<5000<<5000<<5000
    <<3000<<3000<<3000<<3000<<3000
    <<2000<<2000<<2000<<2000<<2000
    <<1500<<1500<<1500<<1500<<1500
    <<1000<<1000<<1000<<1000<<1000
    <<750<<750<<750<<750<<750
    <<500<<500<<500<<500<<500
    <<300<<300<<300<<300<<300
    <<200<<200<<200<<200<<200
    <<150<<150<<150<<150<<150
    <<100<<100<<100<<100<<100
    <<75<<75<<75<<75<<75
    <<50<<50<<50<<50<<50
    <<30<<30<<30<<30<<30
    <<20<<20<<20<<20<<20
    <<15<<15<<15<<15<<15
    <<10<<10<<10<<10<<10
    <<5<<5<<5<<5<<5
    <<5<<5<<5<<5<<5
    <<10<<10<<10<<10<<10
    <<15<<15<<15<<15<<15
    <<20<<20<<20<<20<<20
    <<30<<30<<30<<30<<30
    <<50<<50<<50<<50<<50
    <<75<<75<<75<<75<<75
    <<100<<100<<100<<100<<100
    <<150<<150<<150<<150<<150
    <<200<<200<<200<<200<<200
    <<300<<300<<300<<300<<300
    <<500<<500<<500<<500<<500
    <<750<<750<<750<<750<<750
    <<1000<<1000<<1000<<1000<<1000
    <<1500<<1500<<1500<<1500<<1500
    <<2000<<2000<<2000<<2000<<2000
    <<3000<<3000<<3000<<3000<<3000
    <<5000<<5000<<5000<<5000<<5000
    <<7500<<10000<<15000<<20000<<30000
    <<50000<<100000<<200000<<400000<<800000;
    sigma<<0<<0<<0<<0<<0
        <<0<<0<<0<<0<<0
       <<0<<0<<0<<0<<0
      <<0<<0<<0<<0<<0
     <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0<<0<<0<<0<<0
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.1<<0.1<<0.1<<0.1<<0.1
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001
    <<0.001<<0.001<<0.001<<0.001<<0.001;
    for(int i=0;i<601;i++)
    {
        fre<<pow(10,(double)i/100)*0.001;
    }
    MT1D  mt(dz,sigma,fre);
    mt.Forward();
    int fcount = fre.count();
    for(int i=0;i<fcount;i++)
    {
        qDebug()<<(double)fre[i]<<(double)mt.Rhoxy.at(i)<<(double)mt.Rhoyx.at(i);
    }
}

Widget::~Widget()
{
    delete ui;
}
