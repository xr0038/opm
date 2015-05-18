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
    constexpr double epsilon = 0.2;

    opm::conversion
    rigid_solver
    (wcsprm* pwcs,
     const std::array<double, 6> &&matelem,
     const std::array<double, 6> &&vecelem)
    {
      /*
       * x = cc*x'-ss*y' + Tx
       * y = ss*x'+cc*y' + Ty
       * cc*cc + ss*ss = 1
       */
      const auto &da = pwcs->cdelt[0];
      const auto &dd = pwcs->cdelt[1];

      const double &x_ = matelem[0], &y_ = matelem[1];
      const double &xx = matelem[2], &yy = matelem[4], &N  = matelem[5];
      const double &&xa = vecelem[0]/da;
      const double &&ya = vecelem[1]/da;
      const double &&a_ = vecelem[2]/da;
      const double &&xd = vecelem[3]/dd;
      const double &&yd = vecelem[4]/dd;
      const double &&d_ = vecelem[5]/dd;

      double D =
        std::sqrt((x_*a_+y_*d_-N*(xa+yd))*(x_*a_+y_*d_-N*(xa+yd))
                  +(x_*d_-y_*a_-N*(xd-ya))*(x_*d_-y_*a_-N*(xd-ya)));
      double L =  (x_*x_+y_*y_-N*(xx+yy)+D)/N;
      double cc = (x_*a_+y_*d_-N*(xa+yd))/(x_*x_+y_*y_-N*(xx+yy-L));
      double ss = (x_*d_-y_*a_-N*(xd-ya))/(x_*x_+y_*y_-N*(xx+yy-L));
      double Tx = (a_-x_*cc+y_*ss)/N;
      double Ty = (d_-y_*cc-x_*ss)/N;

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

      return opm::conversion({cc,-ss,ss,cc,Tx,Ty});
    }


    inline void
    array_merge
    (std::array<double, 4> &X,
     const std::array<double, 4> &d,
     const double eps)
    { for (size_t i=0; i<4; i++) X[i] = (1.0-eps)*X[i] + eps*d[i]; }
    inline double
    array_dot
    (const std::array<double, 4> &a, const std::array<double, 4> &b)
    {
      double s(0);
      for (size_t i=0; i<4; i++) s += a[i]*b[i];
      return s;
    }
    inline double
    array_norm
    (const std::array<double, 4> &a)
    { return opm::array_dot(a,a); }
    inline double
    array_distance
    (const std::array<double, 4> &a, const std::array<double, 4> &b)
    {
      std::array<double, 4> c;
      for (size_t i=0; i<4; i++) c[i] = a[i]-b[i];
      return array_norm(c);
    }
    std::array<double, 4>
    elastic_sub_solver
    (const std::array<double, 4> &cn,
     const std::array<double, 3> &A,
     const std::array<double, 4> &b)
    {
      std::array<double, 4> c0;
      c0[0] = ( A[2]*b[0]-A[1]*b[1]);
      c0[1] = (-A[1]*b[0]+A[0]*b[1]);
      c0[2] = ( A[2]*b[2]-A[1]*b[3]);
      c0[3] = (-A[1]*b[2]+A[0]*b[3]);

      std::array<double, 4> cx = {cn[2], cn[3], cn[0], cn[1]};

      double x = opm::array_dot(c0,cx)/opm::array_norm(cx);
      for (size_t i=0; i<4; i++)
        c0[i] = c0[i] - x*cx[i];

      return c0;
    }

    std::array<double, 6>
    elastic_solver
    (wcsprm* pwcs,
     const std::array<double, 6> &&matelem,
     const std::array<double, 6> &&vecelem)
    {
      /*
       * x = c11*x'+c12*y' + Tx
       * y = c21*x'+c22*y' + Ty
       * c11*c21 + c12*c22 = 0
       */
      const double &x_ = matelem[0], &y_ = matelem[1], &xx = matelem[2];
      const double &xy = matelem[3], &yy = matelem[4], &N  = matelem[5];
      const double &xa = vecelem[0], &ya = vecelem[1], &a_ = vecelem[2];
      const double &xd = vecelem[3], &yd = vecelem[4], &d_ = vecelem[5];

      const double &&Axx = xx - x_*x_/N;
      const double &&Axy = xy - x_*y_/N;
      const double &&Ayy = yy - y_*y_/N;
      const double det = (Axx*Ayy - Axy*Axy);
      const std::array<double, 3> A = {Axx/det, Axy/det, Ayy/det};

      const double &&Bax = xa - x_*a_/N;
      const double &&Bay = ya - y_*a_/N;
      const double &&Bdx = xd - x_*d_/N;
      const double &&Bdy = yd - y_*d_/N;
      const std::array<double, 4> b = {Bax, Bay, Bdx, Bdy};

      std::array<double, 4> c =
        {pwcs->cdelt[0]*pwcs->pc[0], pwcs->cdelt[0]*pwcs->pc[1],
         pwcs->cdelt[1]*pwcs->pc[2], pwcs->cdelt[1]*pwcs->pc[3]};

      for (size_t i=0; i<30; i++) {
        std::array<double, 4> cs = elastic_sub_solver(c, A, b);
        opm::array_merge(c,cs,epsilon);
        if (opm::array_distance(c,cs) < 1e-17) break;
      }

      double &c11 = c[0];
      double &c12 = c[1];
      double &c21 = c[2];
      double &c22 = c[3];
      double Tx = (a_-x_*c11-y_*c12)/N;
      double Ty = (d_-x_*c21-y_*c22)/N;

      double cdet = c11*c22 - c12*c21;
      double &&CRPIX1 = - (   c22*Tx - c12*Ty ) / cdet;
      double &&CRPIX2 = - ( - c21*Tx + c11*Ty ) / cdet;
      double cd1 =  -std::sqrt(c11*c11+c12*c12);
      double cd2 =   std::sqrt(c21*c21+c22*c22);
      pwcs->crpix[0] = CRPIX1;
      pwcs->crpix[1] = CRPIX2;
      pwcs->cdelt[0] = cd1;
      pwcs->cdelt[1] = cd2;
      pwcs->pc[0] =  c11/cd1;
      pwcs->pc[1] =  c12/cd1;
      pwcs->pc[2] =  c21/cd2;
      pwcs->pc[3] =  c22/cd2;
      if (wcsset(pwcs) !=0) {
        wcsprt(pwcs);
        fprintf(stderr,"error: wcsset failed.\n");
        exit(1);
      };
      wcsprt(pwcs);
      return std::array<double, 6>{c11,c12,c21,c22,Tx,Ty};
    }
  }

  conversion
  update_wcsprm(wcsprm* pwcs, const opm::matched_list matched)
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

    double xa = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*(p.second.xi()); });
    double ya = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*(p.second.xi()); });
    double a_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + (p.second.xi()); });
    double xd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.x*(p.second.eta()); });
    double yd = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + p.first.y*(p.second.eta()); });
    double d_ = std::accumulate
      (s,e,0.,[&](double s, opm::matched_star p)
       { return s + (p.second.eta()); });

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

    auto res =
      opm::elastic_solver(pwcs,{x_,y_,xx,xy,yy,N},{xa,ya,a_,xd,yd,d_});

    return res;
  }
}
