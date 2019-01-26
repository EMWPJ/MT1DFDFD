#include "mt1d.h"

MT1D::MT1D()
{

}

MT1D::MT1D(QVector < double > _dz,
     QVector < double > _sigma,
     QVector < double > _fre)
{
    setModel(_dz,_sigma,_fre);
}

void MT1D::setModel(QVector < double > _dz,
              QVector < double > _sigma,
              QVector < double > _fre)
{
    setDz(_dz);
    setSigma(_sigma);
    setFre(_fre);
    if(dz.count()<=3
            ||sigma.count()!=dz.count()
            ||fre.count()<1)
    {
        return;
    }
    Forward();
    calZxyZyx();
    calRhoPhase();
}

void MT1D::setDz(QVector < double > _dz)
{
    dz = _dz;
}

void MT1D::setSigma(QVector < double > _sigma)
{
    sigma = _sigma;
}

void MT1D::setFre(QVector < double > _fre)
{
    fre = _fre;
}

void MT1D::Forward()
{
    Ex.clear();
    Hx.clear();
    QVector < QVector < cmplx > > a;
    QVector < cmplx > b;
    int ncount = dz.count()+1;
    a.resize(ncount);
    b.resize(ncount);
    int fcount = fre.count();
    /* TE */
    a[0].clear();
    a[0]<<cmplx(0 , 0)<<cmplx(1.0E10 , 0)<<cmplx(0 , 0);
    b[0] = cmplx(1.0E10 , 0);
    for(int j=1;j<ncount-1;j++)
    {
        a[j].clear();
        a[j].resize(3);
        double dzj = (dz[j] + dz[j-1])/2.0;
        a[j][0] = 1/dzj/dz[j-1];
        a[j][2] = 1/dzj/dz[j];
        b[j] = cmplx(0 , 0);
    }
    a[ncount-1].resize(3);
    a[ncount-1][0] = cmplx(0 , 0);
    a[ncount-1][1] = cmplx(1.0E10 , 0);
    a[ncount-1][2] = cmplx(0 , 0);
    b[ncount-1] = cmplx(0 , 0);
    for(int i=0;i<fcount;i++)
    {
        double omega = fre.at(i) * 2.0 * M_PI;
        for(int j=1;j<ncount-1;j++)
        {
            double sigmaj = (sigma[j-1]*dz[j-1]+sigma[j]*dz[j])/(dz[j-1] + dz[j]);
            a[j][1] = cmplx(0.0,1.0) * omega * Mu * sigmaj - a[j][0] - a[j][2];
        }
        Ex<<catchingMethod(a,b);
    }
    /* TM */
    QVector < double > rho;
    rho.clear();
    for(int i=0;i<ncount-2;i++)
    {
        if(fabs(sigma[i])<1.0E-15)
        {
            rho<<1.0E20;
        }
        else
        {
            rho<<1.0/sigma[i];
        }
    }
    rho<<rho.last();
    a[0].clear();
    a[0]<<cmplx(0 , 0)<<cmplx(1.0E10 , 0)<<cmplx(0 , 0);
    b[0] = cmplx(1E10 , 0);
    for(int j=1;j<ncount-1;j++)
    {
        a[j].clear();
        a[j].resize(3);
        double dzj = (dz[j] + dz[j-1])/2.0;
        a[j][0] = rho[j-1]/dzj/dz[j-1];
        a[j][2] = rho[j]/dzj/dz[j];
        b[j] = cmplx(0 , 0);
    }
    a[ncount-1].resize(3);
    a[ncount-1][0] = cmplx(0 , 0);
    a[ncount-1][1] = cmplx(1.0E10 , 0);
    a[ncount-1][2] = cmplx(0 , 0);
    b[ncount-1] = cmplx(0 , 0);
    for(int i=0;i<fcount;i++)
    {
        double omega = fre.at(i) * 2.0 * M_PI;
        for(int j=1;j<ncount-1;j++)
        {
            a[j][1] = cmplx(0.0,1.0) * omega * Mu - a[j][0] - a[j][2];
        }
        Hx<<catchingMethod(a,b);
    }
}


void MT1D::calZxyZyx()
{
    int fcount = fre.count();
    int ncount = sigma.count();
    int k;
    for(int j=0;j<ncount-1;j++)
    {
        if(fabs(sigma[j])<1.0E-15&&fabs(sigma[j+1])>1.0E-15)
        {
            k = j + 1;
            break;
        }
    }
    for(int i=0;i<fcount;i++)
    {
        double omega = fre.at(i) * 2.0 * M_PI;
        cmplx Hs = sigma[k]*dz[k]/2*(0.75*Ex[i][k]+0.25*Ex[i][k+1])+(Ex[i][k+1]-Ex[i][k])/(cmplx(0.0,1.0)*omega*Mu*dz[k]);
        Zxy<<Ex[i][k]/Hs;
        cmplx Es = (Hx[i][k]-Hx[i][k+1])/sigma[k]/dz[k]-2.0*cmplx(0.0,1.0)*omega*Mu*(0.75*Hx[i][k]+0.25*Hx[i][k+1])/dz[k];
        Zyx<<Es/Hx[i][k];
    }
}

void MT1D::calRhoPhase()
{
    Rhoxy.clear();
    Rhoyx.clear();
    Phasexy.clear();
    Phaseyx.clear();
    int fcount = fre.count();
    for(int i=0;i<fcount;i++)
    {
        double omega = fre.at(i) * 2.0 * M_PI;
        Rhoxy<<abs(Zxy[i])*abs(Zxy[i])/omega/Mu;
        Rhoyx<<abs(Zyx[i])*abs(Zyx[i])/omega/Mu;
        Phasexy<<arg(Zxy[i])*180.0/M_PI;
        Phaseyx<<arg(Zyx[i])*180.0/M_PI;
        qDebug()<<fre[i]<<Rhoxy.last()<<Phasexy.last()<<Rhoyx.last()<<Phaseyx.last();
    }
}

QVector < cmplx > MT1D::catchingMethod(QVector < QVector < cmplx > > _a,
                                 QVector < cmplx > _b)
{
    int size = _b.count();
    QVector < cmplx > alpha,beta,gamma,x,y;
    alpha.resize(size);
    beta.resize(size);
    gamma.resize(size);
    x.resize(size);
    y.resize(size);
    alpha[0] = _a[0][1];
    beta[0] = _a[0][2]/alpha[0];
    for(int i=1;i<size-1;i++)
    {
        gamma[i] = _a[i][0];
        alpha[i] = _a[i][1] - gamma[i]*beta[i-1];
        beta[i] = _a[i][2] / alpha[i];
    }
    gamma[size-1] = _a[size-1][0];
    alpha[size-1] = _a[size-1][1] - gamma[size-1]*beta[size-2];
    y[0] = _b[0]/alpha[0];
    for(int i=1;i<size;i++)
    {
        y[i] = (_b[i]-gamma[i]*y[i-1])/alpha[i];
    }
    x[size-1]=y[size-1];
    for(int i=size-2;i>=0;i--)
    {
        x[i] = y[i] - beta[i]*x[i+1];
    }
    return x;
}
