#ifndef MT1D_H
#define MT1D_H
#define Mu (4*M_PI*pow(10,-7))
#define Ep (pow(10,-9)/36.0/M_PI)
#include <QVector>
#include "qmath.h"
#include <ccomplex>
#include <QDebug>

using cmplx = std::complex<double>;
class MT1D
{
public:
    QVector < double > dz;
    QVector < double > sigma;
    QVector < double > fre;
    QVector < QVector < cmplx > > Ex;
    QVector < QVector < cmplx > > Hx;
    QVector < cmplx > Zxy;
    QVector < cmplx > Zyx;
    QVector < double > Rhoxy;
    QVector < double > Rhoyx;
    QVector < double > Phasexy;
    QVector < double > Phaseyx;
public:
    MT1D();
    MT1D(QVector < double > _dz,
         QVector < double > _sigma,
         QVector < double > _fre);
    void setModel(QVector < double > _dz,
                  QVector < double > _sigma,
                  QVector < double > _fre);
    void setDz(QVector < double > _dz);
    void setSigma(QVector < double > _sigma);
    void setFre(QVector < double > _fre);
    void Forward();
    QVector < cmplx > catchingMethod(QVector < QVector < cmplx > > _a,
                                     QVector < cmplx > _b);
    void calZxyZyx();
    void calRhoPhase();

};

#endif // MT1D_H
