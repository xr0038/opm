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
      /*
       * x = cc*x'-ss*y' + Tx
       * y = ss*x'+cc*y' + Ty
       * cc*cc + ss*ss = 1
       */
      const double &x_ = matelem[0], &y_ = matelem[1];
      const double &xx = matelem[2], &yy = matelem[4], &N  = matelem[5];
      const double &xa = vecelem[0], &ya = vecelem[1], &a_ = vecelem[2];
      const double &xd = vecelem[3], &yd = vecelem[4], &d_ = vecelem[5];
      
      double D =
        std::sqrt((x_*a_+y_*d_-N*(xa+yd))*(x_*a_+y_*d_-N*(xa+yd))
                  +(x_*d_-y_*a_-N*(xd-ya))*(x_*d_-y_*a_-N*(xd-ya)));
      double L =  (x_*x_+y_*y_-N*(xx+yy)+D)/N;
      double cc = (x_*a_+y_*d_-N*(xa+yd))/(x_*x_+y_*y_-N*(xx+yy-L));
      double ss = (x_*d_-y_*a_-N*(xd-ya))/(x_*x_+y_*y_-N*(xx+yy-L));
      double Tx = (a_-x_*cc+y_*ss)/N;
      double Ty = (d_-y_*cc-x_*ss)/N;
      return std::array<double, 4>{cc,ss,Tx,Ty};
    }
  }
  conversion
  update_wcsprm(wcsprm* pwcs, const opm::matched_list matched)
  {
    const auto &s = matched.cbegin();
    const auto &e = matched.cend();

    const auto &dx = pwcs->cdelt[0];
    const auto &dy = pwcs->cdelt[1];

    if (pwcs == NULL) {
      fprintf(stderr, "error: wcsprm is not defined.\n");
      exit(1);
    }

    if (pwcs->crval == NULL) {
      fprintf(stderr, "error: CRVALn is not defined.\n");
      exit(1);
    }

    double xa = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*(p.second.xi()/dx); });
    double ya = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*(p.second.xi()/dx); });
    double a_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + (p.second.xi()/dx); });
    double xd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*(p.second.eta()/dy); });
    double yd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*(p.second.eta()/dy); });
    double d_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + (p.second.eta()/dy); });

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
#ifdef DEBUG
    fprintf(stderr,"(c,s,x,y) = ( %lf %lf %lf %lf )\n",
            res[0], res[1], res[2], res[3]);
#endif
    double &cc = res[0];
    double &ss = res[1];
    double &Tx = res[2];
    double &Ty = res[3];
    double &&CRPIX1 = - (  cc*Tx + ss*Ty );
    double &&CRPIX2 = - (- ss*Tx + cc*Ty );
    pwcs->crpix[0] = CRPIX1;
    pwcs->crpix[1] = CRPIX2;
    pwcs->pc[0] =  cc;
    pwcs->pc[1] = -ss;
    pwcs->pc[2] =  ss;
    pwcs->pc[3] =  cc;
    if (wcsset(pwcs) !=0) {
    wcsprt(pwcs);
      fprintf(stderr,"error: wcsset failed.\n");
      exit(1);
    };

    return opm::conversion({res[0],-res[1],res[1],res[0],res[2],res[3]});
  }
}
