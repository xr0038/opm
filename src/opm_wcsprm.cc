/**
 * @file opm_wcsprm.cc
 * @brief Optimistic Position Matching Algorithm
 *
 * @author Ryou Ohsawa
 * @date 2015
 */
#include "opm.h"
#include <algorithm>

namespace opm {
  namespace {
    std::array<double, 4>
    first_solver
    (const std::array<double, 6> matelem,
     const std::array<double, 6> vecelem)
    {
      const double &x_ = matelem[0], &y_ = matelem[1];
      const double &xx = matelem[2];
      const double &yy = matelem[4], &N  = matelem[5];
      const double &xa = vecelem[0], &ya = vecelem[1], &a_ = vecelem[2];
      const double &xd = vecelem[3], &yd = vecelem[4], &d_ = vecelem[5];
      
      double D =
        std::sqrt((x_*a_+y_*d_-N*(xa+yd))*(x_*a_+y_*d_-N*(xa+yd))
                  +(x_*d_-y_*a_-N*(xd-ya))*(x_*d_-y_*a_-N*(xd-ya)));
      double L = (x_*x_+y_*y_-N*(xx+yy)+D)/N;
      double c_ = (x_*a_+y_*d_-N*(xa+yd))/(x_*x_+y_*y_-N*(xx+yy+L));
      double s_ = (x_*d_-y_*a_-N*(xd-ya))/(x_*x_+y_*y_-N*(xx+yy+L));
      double Tx = (a_-x_*c_+y_*s_)/N;
      double Ty = (d_-y_*c_-x_*s_)/N;
      return std::array<double, 4>{c_,s_,Tx,Ty};
    }
  }
  conversion
  solve_wcsprm(wcsprm* pwcs, const opm::matched_list matched)
  {
    const auto &s = matched.cbegin();
    const auto &e = matched.cend();

    if (pwcs == NULL) {
      fprintf(stderr, "error: wcsprm is not defined.\n");
      exit(1);
    }

    if (pwcs->crval == NULL) {
      fprintf(stderr, "error: CRVALn is not defined.\n");
      exit(1);
    }
    //double A = pwcs->crval[0];
    //double D = pwcs->crval[1];

    double xa = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*p.second.xi(); });
    double ya = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*p.second.xi(); });
    double a_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.second.xi(); });
    double xd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*p.second.eta(); });
    double yd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*p.second.eta(); });
    double d_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.second.eta(); });

    double x_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x; });
    double y_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y; });
    double xx = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*p.first.x; });
    double xy = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*p.first.y; });
    double yy = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*p.first.y; });
    double N = matched.size();

    auto res = opm::first_solver({x_,y_,xx,xy,yy,N},{xa,ya,a_,xd,yd,d_});
    double &&CRPIX1  =  res[0]*res[2] + res[1]*res[3];
    double &&CRPIX2  = -res[1]*res[2] + res[0]*res[3];
    pwcs->crpix[0] = CRPIX1;
    pwcs->crpix[1] = CRPIX2;
    pwcs->pc[0] =  res[0];
    pwcs->pc[1] = -res[1];
    pwcs->pc[2] =  res[1];
    pwcs->pc[3] =  res[0];
    if (wcsset(pwcs) !=0) {
      fprintf(stderr,"error: wcsset failed.\n");
      exit(1);
    };
    return opm::conversion({res[0],-res[1],res[2],res[3],res[4],res[5]});
  }
}
