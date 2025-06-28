// January 16, 2016
#include "ch_sym_DQPM.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_deriv.h>
#include <vector>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <fstream>

const double integration_accuracy = 1e-3;
const double integration_accuracy_a = 1e-6;
const int numb = 20000;

const double Nc = 3.0;
const double c2=0.5517, c1=.830185185; //Corrected on August 30, 2015
const double fpi=92;



double dPq_dqrds(double qr, double qi, double s);
double dPq_dqids(double qr, double qi, double s);

//const double m_l = 1, m_h = 1;
//const double m_l = 3.5, m_h = 95;
//const double m_l = 3.5, m_h = 95;
const double m_l = 3.5, m_h = 95;
//const double m_l = 0.0, m_h = 0.0;
//const double m_l = 50, m_h = 95;

double qr_p = 0.011091, qi_p=0.0;
double sl_p = 1.6572;
double sh_p = 40.361468;

double T, mu, mua;
double T0=270; // test August 30, 2015

double g  = 4.0;
double& Y=g; //Synonim
double M=Y*fpi/2.0;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%Parmeters

double h_u     = pow(121.716,3);
double h_s     = pow(384.424,3);
double c       = 4557.82;
double lambda  = 18.2552+0.0396*pow(Y, 4);
double m2		   = +pow(537.606,2) - pow(11.2915,2)*pow(Y, 4);

struct parameters
{
    double* qr;
    double* qi;
    double* s;
};

double SB_dVdm(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(m*(-((f*(Power(E,8*MPi*qi) + 2*f + 3*Power(E,4*MPi*qi)*Power(f,2) + 2*Power(E,2*MPi*qi)*(1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr)))) - (fa*(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + 2*fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr)))))/(Eq*T)
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*Y;
}



double Usigma(double sl, double sh, double qr, double qi)
{
    double out =  m2* ( 2.0*x2(sl) + x2(sh) )
                  + lambda * (2.0*pow(sl,4) + pow(sh,4))
                  - 2.0*c*sl*sl*sh
                  - 2.0*sl*(h_u  )
                  - 1.0*sh*(h_s  )
                  ;
    return out;
}

double dUsigma_dsl(double sl, double sh, double qr, double qi)
{
    double out =
        m2       *      4*sl
        + lambda *      8*pow(sl,3)
        - c      *      4*sl*sh
        - 2      *      (h_u  )
        ;
    return out;
}


double dUsigma_dsldsl(double sl, double sh)
{
// August 28, 2015
    double out =
        m2       *      4
        + lambda *      24*pow(sl,2)
        - c      *      4*sh
        ;

    return out;
}


double dUsigma_dsh(double sl, double sh, double qr, double qi)
{
    double out =
        m2       *      2*sh
        + lambda *      4*pow(sh,3)
        - c      *      2*x2(sl)
        - 1      *      (h_s  )
        ;
    return out;
}


double dUsigma_dshdsh(double sl, double sh)
{
    double out =
        m2       *      2
        + lambda *      12*pow(sh,2);
    return out;
}



double dUsigma_dshdsl(double sl, double sh)
{
    double out =
        -c*4*sl;
    return out;
}



double dUsigma_dsldqi(double qr, double qi)
{
    return 0.0;
}


double dUsigma_dshdqi(double qr, double qi)
{
    return 0.0;
}

double dUsigma_dsldqr(double qr, double qi)
{
    return 0.0;
}


double dUsigma_dshdqr(double qr, double qi)
{
    return 0.0;
}


double d1func()
{
    return 2.0*M_PI*M_PI/15.0*c1*pow(T0/T,2);
}

double d2func()
{
    return 2.0*M_PI*M_PI/3.0*(1-c2*pow(T0/T,2));
}

double gluonpotential_dqr(double qr, double qi, double sl, double sh)
{

    double d1=d1func();
    double d2=d2func();
    double glue=
      8*(d1*(-1 + 3*qr) - 3*d2*(-1 + 2*qr)*(9*Power(qi,2) + qr - 3*Power(qr,2)))
        ;

    double out = glue;
    return out;
}

double gluonpotential_dqi(double qr, double qi, double sl, double sh)
{

    double d1=d1func();
    double d2=d2func();
    double glue=
-72*qi*(d1 + d2*(1 - 18*Power(qi,2) - 6*qr + 6*Power(qr,2)))
        ;
    double out = glue;
    return out;
}


double gluonpotential_dqrdqr(double qr, double qi, double sl, double sh)
{
// to update
    double d1=d1func();
    double d2=d2func();
    double glue=
24*(d1 + d2*(1 - 18*Power(qi,2) - 10*qr + 18*Power(qr,2)))
        ;
    double out = glue;
    return out;
}


double gluonpotential_dqidqi(double qr, double qi, double sl, double sh)
{
// to update
    double d1=d1func();
    double d2=d2func();
    double glue=
-72*(d1 + d2*(1 - 54*Power(qi,2) - 6*qr + 6*Power(qr,2)))
        ;
    double out = glue;
    return out;
}

double pressure(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s ;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
log(((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr)))/Power(E,4*MPi*qi)) + log(((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr)))/Power(E,4*MPi*qi))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return result;
}


double nb(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s ;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
((f*(Power(E,8*MPi*qi) + 2*f + 3*Power(E,4*MPi*qi)*Power(f,2) + 2*Power(E,2*MPi*qi)*(1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))) - (fa*(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + 2*fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))))/T
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return result;
}


double dPq_dsds(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
-(((Eq*Power(f,2)*Power(m,2)*Power(Power(E,8*MPi*qi) + 2*f + 3*Power(E,4*MPi*qi)*Power(f,2) + 2*Power(E,2*MPi*qi)*(1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr),2))/(Power(1 + Power(E,4*MPi*qi)*f,2)*Power(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr),2)) + (Eq*Power(fa,2)*Power(m,2)*Power(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + 2*fa)*Cos(2*MPi*qr),2))/(Power(Power(E,4*MPi*qi) + fa,2)*Power(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr),2)) + (f*(3*Power(E,4*MPi*qi)*Power(f,2)*(-3*Eq*Power(m,2) + Power(Eq,2)*T - Power(m,2)*T) + 2*f*(-2*Eq*Power(m,2) + Power(Eq,2)*T - Power(m,2)*T) + Power(E,8*MPi*qi)*(-(Eq*Power(m,2)) + Power(Eq,2)*T - Power(m,2)*T) + 2*Power(E,2*MPi*qi)*(-(Eq*(1 + 4*Power(E,4*MPi*qi)*f)*Power(m,2)) + Power(Eq,2)*(1 + 2*Power(E,4*MPi*qi)*f)*T - (1 + 2*Power(E,4*MPi*qi)*f)*Power(m,2)*T)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))) + (fa*(-(Eq*(1 + 4*Power(E,8*MPi*qi)*fa + 9*Power(E,4*MPi*qi)*Power(fa,2))*Power(m,2)) + Power(Eq,2)*(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2))*T - (1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2))*Power(m,2)*T + 2*Power(E,2*MPi*qi)*(2*fa*(-2*Eq*Power(m,2) + Power(Eq,2)*T - Power(m,2)*T) + Power(E,4*MPi*qi)*(-(Eq*Power(m,2)) + Power(Eq,2)*T - Power(m,2)*T))*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))))/(Power(Eq,3)*Power(T,2)))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*x2(g);
}



double dPq_dsdsds(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(2*Power((-3*Power(f,3)*m)/(Eq*T) - (f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (2*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T),3))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3) + (2*Power((-3*Power(fa,3)*m)/(Eq*T) - (2*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T),3))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3) - (3*((-3*Power(f,3)*m)/(Eq*T) - (f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (2*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*((9*Power(f,3)*Power(m,2))/(Power(Eq,2)*Power(T,2)) - (3*Power(f,3))/(Eq*T) + (3*Power(f,3)*Power(m,2))/(Power(Eq,3)*T) + (f*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) + (f*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) + (4*Power(f,2)*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (2*Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T) + (2*Power(f,2)*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T)))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) - (3*((-3*Power(fa,3)*m)/(Eq*T) - (2*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*((9*Power(fa,3)*Power(m,2))/(Power(Eq,2)*Power(T,2)) - (3*Power(fa,3))/(Eq*T) + (3*Power(fa,3)*Power(m,2))/(Power(Eq,3)*T) + (4*Power(fa,2)*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (2*Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) + (2*Power(fa,2)*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) + (fa*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T) + (fa*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T)))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) + ((-27*Power(f,3)*Power(m,3))/(Power(Eq,3)*Power(T,3)) + (27*Power(f,3)*m)/(Power(Eq,2)*Power(T,2)) - (27*Power(f,3)*Power(m,3))/(Power(Eq,4)*Power(T,2)) + (9*Power(f,3)*m)/(Power(Eq,3)*T) - (9*Power(f,3)*Power(m,3))/(Power(Eq,5)*T) - (f*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*Power(T,3)) + (3*f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (3*f*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,4)*Power(T,2)) + (3*f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) - (3*f*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,5)*T) - (8*Power(f,2)*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*Power(T,3)) + (12*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (12*Power(f,2)*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,4)*Power(T,2)) + (6*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T) - (6*Power(f,2)*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,5)*T))/(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr))) + ((-27*Power(fa,3)*Power(m,3))/(Power(Eq,3)*Power(T,3)) + (27*Power(fa,3)*m)/(Power(Eq,2)*Power(T,2)) - (27*Power(fa,3)*Power(m,3))/(Power(Eq,4)*Power(T,2)) + (9*Power(fa,3)*m)/(Power(Eq,3)*T) - (9*Power(fa,3)*Power(m,3))/(Power(Eq,5)*T) - (8*Power(fa,2)*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*Power(T,3)) + (12*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (12*Power(fa,2)*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,4)*Power(T,2)) + (6*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) - (6*Power(fa,2)*Power(m,3)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,5)*T) - (fa*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*Power(T,3)) + (3*fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (3*fa*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,4)*Power(T,2)) + (3*fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T) - (3*fa*Power(m,3)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,5)*T))/(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*x2(g)*g;
}


double dPq_dsdsdqr(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(2*Power((-3*Power(f,3)*m)/(Eq*T) - (f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (2*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T),2)*((-4*f*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi) - 4*Power(E,2*MPi*qi)*Power(f,2)*MPi*Sin(2*MPi*qr)))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3) - (((9*Power(f,3)*Power(m,2))/(Power(Eq,2)*Power(T,2)) - (3*Power(f,3))/(Eq*T) + (3*Power(f,3)*Power(m,2))/(Power(Eq,3)*T) + (f*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) + (f*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) + (4*Power(f,2)*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (2*Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T) + (2*Power(f,2)*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T))*((-4*f*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi) - 4*Power(E,2*MPi*qi)*Power(f,2)*MPi*Sin(2*MPi*qr)))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) + (2*Power((-3*Power(fa,3)*m)/(Eq*T) - (2*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T),2)*(-4*Power(E,2*MPi*qi)*fa*MPi*Sin(2*MPi*qr) - (4*Power(fa,2)*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi)))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3) - (((9*Power(fa,3)*Power(m,2))/(Power(Eq,2)*Power(T,2)) - (3*Power(fa,3))/(Eq*T) + (3*Power(fa,3)*Power(m,2))/(Power(Eq,3)*T) + (4*Power(fa,2)*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,2)*Power(T,2)) - (2*Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) + (2*Power(fa,2)*Power(m,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Power(Eq,3)*T) + (fa*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,2)*Power(T,2)) - (fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T) + (fa*Power(m,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Power(Eq,3)*T))*(-4*Power(E,2*MPi*qi)*fa*MPi*Sin(2*MPi*qr) - (4*Power(fa,2)*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi)))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) - (2*((-3*Power(f,3)*m)/(Eq*T) - (f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (2*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*((4*f*m*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T) + (8*Power(E,2*MPi*qi)*Power(f,2)*m*MPi*Sin(2*MPi*qr))/(Eq*T)))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) - (2*((-3*Power(fa,3)*m)/(Eq*T) - (2*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*((4*Power(E,2*MPi*qi)*fa*m*MPi*Sin(2*MPi*qr))/(Eq*T) + (8*Power(fa,2)*m*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T)))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) + ((-4*f*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Power(Eq,2)*Power(T,2)) - (16*Power(E,2*MPi*qi)*Power(f,2)*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(Eq,2)*Power(T,2)) + (4*f*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T) + (8*Power(E,2*MPi*qi)*Power(f,2)*MPi*Sin(2*MPi*qr))/(Eq*T) - (4*f*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Power(Eq,3)*T) - (8*Power(E,2*MPi*qi)*Power(f,2)*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(Eq,3)*T))/(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr))) + ((-4*Power(E,2*MPi*qi)*fa*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(Eq,2)*Power(T,2)) - (16*Power(fa,2)*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Power(Eq,2)*Power(T,2)) + (4*Power(E,2*MPi*qi)*fa*MPi*Sin(2*MPi*qr))/(Eq*T) + (8*Power(fa,2)*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T) - (4*Power(E,2*MPi*qi)*fa*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(Eq,3)*T) - (8*Power(fa,2)*Power(m,2)*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Power(Eq,3)*T))/(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*x2(g);
}


double dPq_dsdqrdqr(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
((8*f*m*Power(MPi,2)*Cos(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T) + (16*Power(E,2*MPi*qi)*Power(f,2)*m*Power(MPi,2)*Cos(2*MPi*qr))/(Eq*T))/(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr))) + ((8*Power(E,2*MPi*qi)*fa*m*Power(MPi,2)*Cos(2*MPi*qr))/(Eq*T) + (16*Power(fa,2)*m*Power(MPi,2)*Cos(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T))/(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr))) - (2*((-4*f*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi) - 4*Power(E,2*MPi*qi)*Power(f,2)*MPi*Sin(2*MPi*qr))*((4*f*m*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T) + (8*Power(E,2*MPi*qi)*Power(f,2)*m*MPi*Sin(2*MPi*qr))/(Eq*T)))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) - (2*(-4*Power(E,2*MPi*qi)*fa*MPi*Sin(2*MPi*qr) - (4*Power(fa,2)*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi))*((4*Power(E,2*MPi*qi)*fa*m*MPi*Sin(2*MPi*qr))/(Eq*T) + (8*Power(fa,2)*m*MPi*Sin(2*MPi*qr))/(Power(E,2*MPi*qi)*Eq*T)))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2) + ((-3*Power(f,3)*m)/(Eq*T) - (f*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (2*Power(f,2)*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*(-(((-8*f*Power(MPi,2)*Cos(2*MPi*qr))/Power(E,2*MPi*qi) - 8*Power(E,2*MPi*qi)*Power(f,2)*Power(MPi,2)*Cos(2*MPi*qr))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2)) + (2*Power((-4*f*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi) - 4*Power(E,2*MPi*qi)*Power(f,2)*MPi*Sin(2*MPi*qr),2))/Power(1 + Power(f,3) + f*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + Power(f,2)*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3)) + ((-3*Power(fa,3)*m)/(Eq*T) - (2*Power(fa,2)*m*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)))/(Eq*T) - (fa*m*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)))/(Eq*T))*(-((-8*Power(E,2*MPi*qi)*fa*Power(MPi,2)*Cos(2*MPi*qr) - (8*Power(fa,2)*Power(MPi,2)*Cos(2*MPi*qr))/Power(E,2*MPi*qi))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),2)) + (2*Power(-4*Power(E,2*MPi*qi)*fa*MPi*Sin(2*MPi*qr) - (4*Power(fa,2)*MPi*Sin(2*MPi*qr))/Power(E,2*MPi*qi),2))/Power(1 + Power(fa,3) + Power(fa,2)*(Power(E,4*MPi*qi) + (2*Cos(2*MPi*qr))/Power(E,2*MPi*qi)) + fa*(Power(E,-4*MPi*qi) + 2*Power(E,2*MPi*qi)*Cos(2*MPi*qr)),3))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*(g);
}

double dPq_dqrdqr(double qr, double qi, double s) //Per flavor
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
8*Power(E,2*MPi*qi)*Power(MPi,2)*(Cos(2*MPi*qr)*(-(f/(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))) - fa/(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))) + 2*Power(E,2*MPi*qi)*(-(Power(f,2)/Power(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr),2)) - Power(fa,2)/Power(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr),2))*Power(Sin(2*MPi*qr),2))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result;
}



double dPq_dqrds(double qr, double qi, double s) //Per flavor
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(-4*Power(E,2*MPi*qi)*(-1 + f*fa)*m*MPi*(8*Power(E,6*MPi*qi)*f*fa*(1 + f*fa)*Cos(2*MPi*qr) + (f + Power(E,4*MPi*qi)*fa)*(-Power(f,2) - Power(E,8*MPi*qi)*Power(fa,2) + Power(E,4*MPi*qi)*(1 + 6*f*fa + Power(f,2)*Power(fa,2)) + 2*Power(E,4*MPi*qi)*f*fa*Cos(4*MPi*qr)))*Sin(2*MPi*qr))/(Eq*T*Power(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr),2)*Power(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr),2))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*g;
}


double dPq_dqids(double qr, double qi, double s) //Per flavor
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);
            double out =
                pow(p,2) * (
(4*m*MPi*(-((f*(Power(E,8*MPi*qi) - 2*f + Power(E,2*MPi*qi)*(-1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr)))) + (Power(f,2)*(Power(E,8*MPi*qi) - f + Power(E,2*MPi*qi)*(-1 + Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr))*(Power(E,8*MPi*qi) + 2*f + 3*Power(E,4*MPi*qi)*Power(f,2) + 2*Power(E,2*MPi*qi)*(1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/(Power(1 + Power(E,4*MPi*qi)*f,2)*Power(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr),2)) + (fa*(1 - 2*Power(E,8*MPi*qi)*fa - Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) - 2*fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))) + (Power(fa,2)*(-1 + Power(E,8*MPi*qi)*fa + Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) - fa)*Cos(2*MPi*qr))*(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + 2*fa)*Cos(2*MPi*qr)))/(Power(Power(E,4*MPi*qi) + fa,2)*Power(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr),2))))/(Eq*T)
                );//pow(T,3);
            if (fabs(out)<1e-10) return 0.0;
            return out;//*pow(T,3);
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper = 10000;
    gsl_set_error_handler_off();
    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0,0.00001, numb, GSL_INTEG_GAUSS21,  w, &result, &error);
    gsl_integration_workspace_free(w);
    gsl_set_error_handler (NULL);
    return  result*g;
}



double dOmegadsldsl(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 * T ;
    double quark  = factor_light*dPq_dsds(qr, qi, sl);
    double meson = dUsigma_dsldsl(sl, sh);
    double vacuum = -2.0*Nc/(8*x2(M_PI)) * pow(g,4) * x2(sl) * ( 7 + 12*log(g*sl/M));

    double SB = 	- m_l / Y *  factor_light * dPq_dsdsds(qr, qi, sl) ;
    return quark + meson + vacuum + SB;
}



double dOmegadshdsh(double qr, double qi, double sl, double sh)
{
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 * T ;
    double meson = dUsigma_dshdsh(sl, sh);
    double quark  = factor_heavy*dPq_dsds(qr, qi, sh);
    double vacuum = -Nc/(8*x2(M_PI)) * pow(g,4) * x2(sh) * ( 7 + 12*log(g*sh/M));
    double SB = 	- m_h / Y *  factor_heavy * dPq_dsdsds(qr, qi, sh) ;
    return quark + meson+ vacuum + SB;
}



double dOmegadshdsl(double qr, double qi, double sl, double sh)
{
    double meson = dUsigma_dshdsl(sl, sh);
    return meson;
}


double dOmegadsldqr(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 * T ;
    double quark  = factor_light*dPq_dqrds(qr, qi, sl);
    double meson = dUsigma_dsldqr(qr,qi);
    double SB = 	- m_l / Y *  factor_light * dPq_dsdsdqr(qr, qi, sl) ;
    return quark + meson + SB;
}


double dOmegadshdqr(double qr, double qi, double sl, double sh)
{
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 * T ;
    double meson = dUsigma_dshdqr(qr,qi);
    double quark  = factor_heavy*dPq_dqrds(qr, qi, sh);
    double SB = 	- m_h / Y *  factor_heavy * dPq_dsdsdqr(qr, qi, sh) ;
    return quark + meson + SB;
}


double dOmegadqrdqr(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 * T ;
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 * T ;
    double gluon = gluonpotential_dqrdqr(qr, qi, sl, sh);
    double quark  = factor_heavy*dPq_dqrdqr(qr, qi, sh)  + factor_light*dPq_dqrdqr(qr, qi, sl)  ;
    double SB =
        - m_l / Y *  factor_light * dPq_dsdqrdqr(qr, qi, sl)
        - m_h / Y *  factor_heavy * dPq_dsdqrdqr(qr, qi, sh) ;
    return quark + gluon*x2(T)*x2(T)+SB;
}




double EoM_s(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(m*(-((f*(Power(E,8*MPi*qi) + 2*f + 3*Power(E,4*MPi*qi)*Power(f,2) + 2*Power(E,2*MPi*qi)*(1 + 2*Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr)))) - (fa*(1 + 2*Power(E,8*MPi*qi)*fa + 3*Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + 2*fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr)))))/(Eq*T)
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result*g;
}



double EoM_qr(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
(-4*Power(E,2*MPi*qi)*MPi*((f + Power(E,4*MPi*qi)*fa)*(1 + f*fa) + 4*Power(E,2*MPi*qi)*f*fa*Cos(2*MPi*qr))*Sin(2*MPi*qr))/((Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr)))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;

    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result;
}


double EoM_qi(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);

            return
                pow(p,2) * (
4*MPi*((f*(Power(E,8*MPi*qi) - f + Power(E,2*MPi*qi)*(-1 + Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))) + (fa*(-1 + Power(E,8*MPi*qi)*fa + Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) - fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,0.000001, numb, GSL_INTEG_GAUSS15,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result;
}





double pressure_final(double qr, double qi, double sl, double sh)
{
// updated August 28, 2015
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 * T ;
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 * T ;
    double norm = 1.0/pow(T,4);
    double Us = Usigma(sl, sh, qr, qi)*norm;
    double thermal = factor_light*pressure(qr, qi, sl)*norm +  factor_heavy*pressure(qr, qi, sh)*norm;
    double m_light = g * sl;
    double m_heavy = g * sh;
    double vacuum_light = - Nc * 2 /(8 * M_PI * M_PI) * 	pow(m_light, 4) * log(m_light/M) * norm;
    double vacuum_heavy = - Nc * 1 /(8 * M_PI * M_PI) * 	pow(m_heavy, 4) * log(m_heavy/M) * norm;
    
	//= -2.0*Nc/(8*x2(M_PI)) * pow(g,4) * x2(sl) * ( 7 + 12*log(g*sl/M));
	double d1=d1func();
    double d2=d2func();
    double Ug =
        (
-4*d1*(9*Power(qi,2) + (2 - 3*qr)*qr) + 4*d2*(81*Power(qi,4) - 9*Power(qi,2)*(1 - 6*qr + 6*Power(qr,2)) + Power(qr,2)*(3 - 10*qr + 9*Power(qr,2)))
        );
    return thermal + vacuum_light + vacuum_heavy + Us + Ug
           - m_l/Y*factor_light*norm*EoM_s(qr, qi, sl)
           - m_h/Y*factor_heavy*norm*EoM_s(qr, qi, sh);


    ;
}



double dPdqi2(double qr, double qi, double s)
{
    class uint
    {
    public:
        static double func(double p, void* params)
        {
            parameters* P = (parameters*) params;
            double qr = *(P->qr);
            double qi = *(P->qi);
            double s = *(P->s);
            double m = g*s;
            double Eq = sqrt(p*p+x2(m));
            double f = exp(-(Eq-mu)/T);
            double fa = exp(-(Eq-mua)/T);
			/*cout << " " <<   pow(p,2) * (
#include "dmu4.dat"
                ) << "\n";*/
            return
                pow(p,2) * (
8*Power(MPi,2)*((-2*Power(f,2)*Power(Power(E,8*MPi*qi) - f + Power(E,2*MPi*qi)*(-1 + Power(E,4*MPi*qi)*f)*Cos(2*MPi*qr),2))/(Power(1 + Power(E,4*MPi*qi)*f,2)*Power(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr),2)) + (f*(2*(Power(E,8*MPi*qi) + f) + (Power(E,2*MPi*qi) + Power(E,6*MPi*qi)*f)*Cos(2*MPi*qr)))/((1 + Power(E,4*MPi*qi)*f)*(Power(E,4*MPi*qi) + Power(f,2) + 2*Power(E,2*MPi*qi)*f*Cos(2*MPi*qr))) - (2*Power(fa,2)*Power(-1 + Power(E,8*MPi*qi)*fa + Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) - fa)*Cos(2*MPi*qr),2))/(Power(Power(E,4*MPi*qi) + fa,2)*Power(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr),2)) + (fa*(2 + 2*Power(E,8*MPi*qi)*fa + Power(E,2*MPi*qi)*(Power(E,4*MPi*qi) + fa)*Cos(2*MPi*qr)))/((Power(E,4*MPi*qi) + fa)*(1 + Power(E,4*MPi*qi)*Power(fa,2) + 2*Power(E,2*MPi*qi)*fa*Cos(2*MPi*qr))))
                );
        }
    } func;

    parameters to_fun;
    to_fun.qr = &qr;
    to_fun.qi = &qi;
    to_fun.s = &s;

    gsl_function F;
    F.function = &func.func;
    F.params = &to_fun;
    double result,error;
    size_t num;
    double upper =  10000;

    gsl_integration_workspace * w   = gsl_integration_workspace_alloc (numb);
    int status = gsl_integration_qag(&F, 0.0 , upper ,0.0,1e-3, numb, GSL_INTEG_GAUSS21,  w, &result, &error);

    gsl_integration_workspace_free(w);
    return  result;
}

double nb_final(double qr, double qi, double sl, double sh)
{
// to update
    double factor_light =  1.0/(2.0*x2(M_PI))  * 2 * 2  ;
    double factor_heavy =  1.0/(2.0*x2(M_PI))  * 2 * 1  ;
    double norm = 1.0/pow(T,3);
    double thermal = factor_light*nb(qr, qi, sl)*norm +  factor_heavy*nb(qr, qi, sh)*norm;
    return thermal;
}



double pressure_final_vacuum(void)
{
    double norm = 1.0/pow(T,4);
    double sl = 46.0029;
    double sh = 76.0827;
    double qr = 0.3317;
    double qi = 0.0;
    double Us = Usigma(sl, sh, qr, qi)*norm;
    double m_light = g * sl;
    double m_heavy = g * sh;
    double vacuum_light = - Nc * 2 /(8 * M_PI * M_PI) * 	pow(m_light, 4) * log(m_light/M) * norm;
    double vacuum_heavy = - Nc * 1 /(8 * M_PI * M_PI) * 	pow(m_heavy, 4) * log(m_heavy/M) * norm;
    double d1=d1func();
    double d2=d2func();
    double Ug =
        (
-4*d1*(9*Power(qi,2) + (2 - 3*qr)*qr) + 4*d2*(81*Power(qi,4) - 9*Power(qi,2)*(1 - 6*qr + 6*Power(qr,2)) + Power(qr,2)*(3 - 10*qr + 9*Power(qr,2)))
        );
    return vacuum_light + vacuum_heavy + Us + Ug;
}


double SigmaMass(double qr, double qi, double sl, double sh)
{
    double out = ( dOmegadsldsl(qr,qi,sl,sh)  +  dOmegadshdsh(qr,qi,sl,sh) + 2.0*dOmegadshdsl(qr,qi,sl,sh) )/6.0;
    ;
    return out/fpi/fpi;
}


double PionMass(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T  ;
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;


    double H_u = h_u + m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);

    double out =
        (H_u  )/(2.0*sl);
    return out/fpi/fpi;
}


double KaonMass(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 *T ;

    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;


    double H_u = h_u + m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);
    double H_s = h_s + m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh);
    double out =
        (H_u +H_s )/(2.0*(sl+sh));
    return out/fpi/fpi;
}


double KappaMass(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;

    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;


    double H_u = h_u + m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);
    double H_s = h_s + m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh);

    double out =
        (H_u -H_s )/(2.0*(sl-sh));
    return out/fpi/fpi;
}


double EtaMass(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;

    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;


    double H_u = h_u + m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);
    double H_s = h_s + m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh);


    double sqrtdiff = sqrt(pow(0.5*H_u/sl - 0.5*H_s/sh + 2 * c * sh * (1.0-0.5*pow(sl/sh,2) ),2) + 8*c*c*sl*sl);
    double sum = (0.5*H_u/sl + 0.5*H_s/sh + 2 * c * sh * (1.0+0.5*pow(sl/sh,2) ));
    double out = 0.5*(sum-sqrtdiff);
    return out/fpi/fpi;
}

double EtaPrimeMass(double qr, double qi, double sl, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 *T  ;

    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 *T  ;


    double H_u = h_u + m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);
    double H_s = h_s + m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh);

    double sqrtdiff = sqrt(pow(0.5*H_u/sl - 0.5*H_s/sh + 2 * c * sh * (1.0-0.5*pow(sl/sh,2) ),2) + 8*c*c*sl*sl);
    double sum = (0.5*H_u/sl + 0.5*H_s/sh + 2 * c * sh * (1.0+0.5*pow(sl/sh,2) ));
    double out = 0.5*(sum+sqrtdiff);
    return out/fpi/fpi;
}

int Interface (const gsl_vector * v, void *params,
               gsl_vector * f)
{
    double qr = (gsl_vector_get(v, 2)) ;
    double qi = (gsl_vector_get(v, 3));
    double sl = fabs(gsl_vector_get(v, 0));
    double sh = fabs(gsl_vector_get(v, 1));

    //cout << " " << sl <<" \n"	<< flush;

    //if (qr<0) qr=0;

    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 * T ;
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 * T ;


    gsl_vector_set (f, 1, factor_light * EoM_qr(qr,qi,sl)/pow(T,4) +
                    factor_heavy * EoM_qr(qr,qi,sh)/pow(T,4) +
                    gluonpotential_dqr(qr, qi,sl,sh)
                    - m_l / Y *  factor_light * dPq_dqrds(qr, qi, sl)/pow(T,4)
                    - m_h / Y *  factor_heavy * dPq_dqrds(qr, qi, sh)/pow(T,4)
                   ) ;

    gsl_vector_set (f, 2, factor_light * EoM_qi(qr,qi,sl)/pow(T,4) +
                    factor_heavy * EoM_qi(qr,qi,sh)/pow(T,4) +
                    gluonpotential_dqi(qr, qi,sl,sh)
                    - m_l / Y *  factor_light * dPq_dqids(qr, qi, sl)/pow(T,4)
                    - m_h / Y *  factor_heavy * dPq_dqids(qr, qi, sh)/pow(T,4)
                   ) ;

    gsl_vector_set (f, 0, factor_light * EoM_s(qr,qi,sl)/pow(T,4) -
                    (Power(g,4)*Nc*2*Power(sl,3)*(1 + 4*log((g*sl)/M)))/(8.*Power(Pi,2))/pow(T,4) // vacuum term light Nf=2
                    + dUsigma_dsl(sl,sh,qr,qi)/pow(T,4)
                    - m_l / Y *  factor_light * dPq_dsds(qr, qi, sl)/pow(T,4)
                   ) ;

    gsl_vector_set (f, 3, factor_heavy * EoM_s(qr,qi,sh)/pow(T,4) -
                    (Power(g,4)*Nc*1*Power(sh,3)*(1 + 4*log((g*sh)/M)))/(8.*Power(Pi,2))/pow(T,4) //vacuum term heavy Nf=1
                    + dUsigma_dsh(sl,sh,qr,qi)/pow(T,4)
                    - m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh)/pow(T,4)
                   ) ;

    return GSL_SUCCESS;

}


double DU2ds(double qr, double qi, double sl_i, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 *T ;
    double step = sl_i*0.01;
    double sl = sl_i+step;
    double up=factor_light * EoM_s(qr,qi,sl) -
              (Power(g,4)*Nc*2*Power(sl,3)*(1 + 4*log((g*sl)/M)))/(8.*Power(Pi,2)) // vacuum term light Nf=2
              + dUsigma_dsl(sl,sh,qr,qi)
              - m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);

    sl = sl_i-step;
    double down=factor_light * EoM_s(qr,qi,sl) -
                (Power(g,4)*Nc*2*Power(sl,3)*(1 + 4*log((g*sl)/M)))/(8.*Power(Pi,2)) // vacuum term light Nf=2
                + dUsigma_dsl(sl,sh,qr,qi)
                - m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);
    return 0.5*(up-down)/step;
}


double DU2dh(double qr, double qi, double sl, double sh_i)
{
    double factor_heavy =  - 1.0/(2.0*x2(M_PI))  * 2 * 1  *T ;
    double step = sh_i*0.01;
    double sh = sh_i+step;
    double up= factor_heavy * EoM_s(qr,qi,sh) -
               (Power(g,4)*Nc*1*Power(sh,3)*(1 + 4*log((g*sh)/M)))/(8.*Power(Pi,2)) //vacuum term heavy Nf=1
               +  dUsigma_dsh(sl,sh,qr,qi)
               -  m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh)
               ;

    sh = sh_i-step;
    double down=factor_heavy * EoM_s(qr,qi,sh) -
                (Power(g,4)*Nc*1*Power(sh,3)*(1 + 4*log((g*sh)/M)))/(8.*Power(Pi,2)) //vacuum term heavy Nf=1
                + dUsigma_dsh(sl,sh,qr,qi)
                - m_h / Y *  factor_heavy * dPq_dsds(qr, qi, sh);
    return 0.5*(up-down)/step;
}



double DU2dhds(double qr, double qi, double sl, double sh)
{
    return -4.0*c*sl;
}

double a0Mass(double qr, double qi, double sl_i, double sh)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 2 *T ;
    double step = sl_i*0.01;
    double sl = sl_i+step;
    double up=factor_light * EoM_s(qr,qi,sl) -
              (Power(g,4)*Nc*2*Power(sl,3)*(1 + 4*log((g*sl)/M)))/(8.*Power(Pi,2)) // vacuum term light Nf=2
              - m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);

    sl = sl_i-step;
    double down=factor_light * EoM_s(qr,qi,sl) -
                (Power(g,4)*Nc*2*Power(sl,3)*(1 + 4*log((g*sl)/M)))/(8.*Power(Pi,2)) // vacuum term light Nf=2
                - m_l / Y *  factor_light * dPq_dsds(qr, qi, sl);


    double out =  0.25*0.5*(up-down)/step + m2 + 6. * lambda*sl_i*sl_i + c*sh;

    return out/fpi/fpi;
}

void chis(double qr, double qi, double sl, double sh, double& chi_pipi, double& chi_a0a0)
{
    double factor_light =  - 1.0/(2.0*x2(M_PI))  * 2 * 1 *T ;

	double pipi = (factor_light * EoM_s(qr,qi,sl) / (2.0*sl)  - (Power(g,4)*Nc*Power(sl,2)*(1 + 4*log((g*sl)/M)))/(16.*Power(Pi,2)) );
	double a0a0 = (factor_light * dPq_dsds(qr, qi, sl)   - (Power(g,4)*Nc*Power(sl,2)*(7 + 12*log((g*sl)/M)))/(16.*Power(Pi,2)) );

	double pi_0 = (-c*sh + 2*lambda*sl*sl + m2); 
	double a0_0 = (c*sh + 6*lambda*sl*sl + m2); 

	double piressumed = 1.0/(PionMass(qr, qi, sl, sh)*fpi*fpi);
	double a0ressumed = 1.0/(a0Mass(qr, qi, sl, sh)*fpi*fpi);
; 
	
	//chi_pipi = pipi + 3.0*pipi*piressumed*pipi*g*g;
	chi_pipi = pipi*pi_0*piressumed;
	//chi_a0a0 = a0a0 + 3.0*a0a0*a0ressumed*a0a0*g*g;
	chi_a0a0 = a0a0*a0_0*a0ressumed;
}


int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
    fprintf(stdout, "iter = %3u x = %f %f  "
            "f(x) = % .3e % .3e   \n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1)
           );
    return 1;
}




double solver (void* p)
{

    double qr=qr_p, qi = qi_p, sl=sl_p, sh=sh_p;

    {

        const gsl_multiroot_fsolver_type *Tt;
        gsl_multiroot_fsolver *solv;

        int status;
        size_t i, iter = 0;

        const size_t n = 4;
        gsl_multiroot_function f = {&Interface, n, NULL};

        double x_init[4] = {sl, sh,  qr, qi};
        gsl_vector *x = gsl_vector_alloc (n);

        gsl_vector_set (x, 0, x_init[0]);
        gsl_vector_set (x, 1, x_init[1]);
        gsl_vector_set (x, 2, x_init[2]);
        gsl_vector_set (x, 3, x_init[3]);

        Tt = gsl_multiroot_fsolver_hybrids;
        solv = gsl_multiroot_fsolver_alloc (Tt, n);
        gsl_multiroot_fsolver_set (solv, &f, x);

        //print_state (iter, solv);

        do
        {
            iter++;
            status = gsl_multiroot_fsolver_iterate (solv);

            // print_state (iter, solv);

            if (status)   /* check if solver is stuck */
                break;

            status =
                gsl_multiroot_test_residual (solv->f, 1e-12);
        }
        while (status == GSL_CONTINUE && iter < 2000);

        //printf ("status = %s\n", gsl_strerror (status));

        sl = fabs(gsl_vector_get (solv->x, 0) );
        sh = fabs(gsl_vector_get (solv->x, 1) );
        qr = (gsl_vector_get (solv->x, 2) );
        qi = (gsl_vector_get (solv->x, 3) );

        qr_p=qr;
        qi_p=qi;
        sl_p=sl;
        sh_p=sh;

        //printf("%f %f %f %f\n", T, l , ls, sigma);
        fflush(stdout);
        gsl_multiroot_fsolver_free (solv);
        gsl_vector_free (x);
    }
    return 1.0;
}



double Lfunc(double qr, double qi)
{
    return 1.0/3.0 * ( exp(4.0*M_PI*qi) + 2.0  * exp(-2.0*M_PI*qi)*cos(2.0*M_PI*qr) );
}

double Ldegfunc(double qr, double qi)
{
    return 1.0/3.0 * ( exp(-4.0*M_PI*qi) + 2.0  * exp(2.0*M_PI*qi)*cos(2.0*M_PI*qr) );
}

int main(void)
{
    FILE * pOp;
    FILE * pTher;
    FILE * pMasses;
    FILE * pSusc;
    /*pOp = fopen ("data/OP_u.dat","w");
    pTher = fopen ("data/test_u.dat","w");
    pMasses = fopen ("data/masses.dat","w");
    */
	pSusc = fopen ("data/chi_pipi.dat","w");
    //for simplicity
    Y=5;
	g=Y;
    lambda  = 18.255188335605816+.03962588514828971*pow(Y, 4);
    m2		   = +pow(537.606,2) - pow(11.2915,2)*pow(Y, 4);
    M=Y*fpi/2.0;
    //for simplicity

    mu=0;
    mua=-mu;
    gsl_set_error_handler_off();
    for(T=500; T>50; T-=1)
    {
        if(T>3000) return 1;
        mu=T*0.0;
        mua=-mu;

        solver(NULL);
        cout <<  T << " " << sl_p << " " << sh_p<< " " << qr_p << " "<< qi_p<< "\n";

        double chi_pipi; 
        double chi_a0a0;

		chis(qr_p, qi_p, sl_p, sh_p, chi_pipi, chi_a0a0); 

		
		double kappa = pow(g,4)*3/16/M_PI/M_PI; 
		double vacuum_diff = 6.0 * kappa * sl_p*sl_p + 8.0 * kappa * sl_p*sl_p*log(g*sl_p/M); 


		fprintf(pSusc,"%.25f ", T);
		fprintf(pSusc,"%.25f ", chi_pipi );
		fprintf(pSusc,"%.25f ", chi_a0a0 );
		fprintf(pSusc,"%.25f ", vacuum_diff );
		fprintf(pSusc,"\n");
	}
}
