/**
 * @file opm_final.cc
 * @brief Optimistic Position Matching Algorithm
 *
 * @author Ryou Ohsawa
 * @date 2015
 */
#include "opm.h"

namespace opm {
  namespace {
    inline opm::object
    convert_coordinates(const opm::conversion &coeff,
                        const opm::object &t)
    {
      double &&x = coeff[0]*t.x + coeff[1]*t.y + coeff[4];
      double &&y = coeff[2]*t.x + coeff[3]*t.y + coeff[5];
      return opm::object(x,y,t.m);
    }

    inline opm::objectlist
    convert_starlist(const opm::conversion &coeff,
                     const opm::objectlist &obj)
    {
      opm::objectlist retval;
      for (auto &p : obj)
        retval.push_back(convert_coordinates(coeff,p));
      return retval;
    }

    inline opm::referencelist
    sift_by_declination
    (const opm::object &T, const double d, const opm::referencelist &src)
    {
      opm::reference Tpd(T.x,T.y+d,T.m);
      opm::reference Tmd(T.x,T.y-d,T.m);
      auto lb = 
        std::lower_bound(src.begin(),src.end(),Tmd,
                         [&](opm::reference l, opm::reference r )
                         { return l.py() < r.py(); });
      auto ub 
        = std::upper_bound(src.begin(),src.end(),Tpd,
                           [&](opm::reference l, opm::reference r )
                           { return l.py() < r.py(); });
      opm::referencelist retval(lb,ub);
      return retval;
    }
    
    inline bool
    neighbor
    (const opm::object &T, const opm::reference &R,
     const double &dx, const double &dy)
    {
      return std::abs(T.x-R.px()) < dx && std::abs(T.y-R.py()) < dy;
    }

    inline opm::referencelist
    sift(const opm::object &T, const double dx, const double dy,
         const opm::referencelist &src)
    {
      opm::referencelist retval;
      opm::referencelist &&list = opm::sift_by_declination(T, dy, src);
      for (const auto &t: list) {
        if (neighbor(T,t,dx,dy)) retval.push_back(t);
      }
      return retval;
    }
  }

  opm::matched_list
  final::match(const opm::conversion coeff,
               const opm::objectlist &obj, const double dx, const double dy,
               const opm::referencelist &src)
  {
    opm::objectlist &&converted = opm::convert_starlist(coeff, obj);

    opm::matched_list retval;
    size_t Nobj = converted.size();

    for (size_t i=0; i < Nobj; i++) {
      auto &p = obj[i];
      auto &q = converted[i];
      opm::referencelist &&sifted = opm::sift(q, dx, dy, src);
      if (sifted.size() != 1) continue;
      retval.push_back(std::make_pair(p, sifted[0]));
    }
    return retval;
  }
}
