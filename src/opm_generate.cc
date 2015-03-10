/**
 * @file opm_generate.cc
 * @brief Optimistic Position Matching Algorithm
 *
 * @author Ryou Ohsawa
 * @date 2015
 */
#include "opm.h"

namespace opm
{
  namespace {
    inline double
    sqdist(const opm::xym &obj1, const opm::xy &obj2)
    {
      return obj1.x * obj1.x + obj1.y * obj1.y
        + obj2.x * obj2.x + obj2.y * obj2.y
        - 2.0 * (obj1 * obj2);
    }
    inline int16_t
    check_order(const double &a, const double &b, const double &c)
    {
      if (a >= b && a >= c) {
        if (b >= c) return 0; else return 1;
      } else if (b >= c && b >= a) {
        if (c >= a) return 2; else return 3;
      } else if (c >= a && c >= b) {
        if (a >= b) return 4; else return 5;
      }
      return -1;
    }
    inline double
    d(const opm::xym A, const opm::xym B, const opm::xym C)
    { return (A-C) * (B-C); }
    inline double
    r(const double a, const double c)
    { return a/c; }
  }

  objectlist
  generate_objectlist
  (const std::vector<opm::xym> &obj)
  {
    opm::objectlist retval;
    for (const auto &p : obj) retval.push_back(opm::object(p));
    std::sort(retval.begin(), retval.end(),
              [&](const opm::object l, const opm::object r)
              { return l.y < r.y; });
    return retval;
  }

  referencelist
  generate_referencelist
  (const std::vector<opm::xym> &obj,
   wcsprm* wcs_)
  {
    opm::referencelist retval;
    for (const auto &p : obj) 
      retval.push_back(opm::reference(p,wcs_));
    std::sort(retval.begin(), retval.end(),
              [&](const opm::reference l, const opm::reference r)
              { return l.eta() < r.eta(); });
    return retval;
  }

  triangles
  generate_triangles
  (const opm::objectlist::const_iterator &it_begin,
   const opm::objectlist::const_iterator &it_end)
  {
    std::vector<const opm::object*> tmp_list;
    const size_t n = (size_t)(it_end-it_begin);

    for (auto it = it_begin; it != it_end; it++)
      tmp_list.push_back(&(*it));

    std::sort(tmp_list.begin(), tmp_list.end(),
              [&](const opm::object* pA, const opm::object* pB)
              { return pA->m < pB->m; });
    auto it = tmp_list.begin();

    if (n < 3) {
      fprintf(stderr,"error: number of stars not enough.\n");
      exit(1);
    }
    triangles return_list;

    for (size_t i=0; i<n; i++) {
      for (size_t j=i+1; j<n; j++) {
        for (size_t k=j+1; k<n; k++) {
          const opm::object* pA = *(it+i);
          const opm::object* pB = *(it+j);
          const opm::object* pC = *(it+k);
          double c = sqdist(**(it+i),**(it+j));
          double a = sqdist(**(it+j),**(it+k));
          double b = sqdist(**(it+k),**(it+i));
          switch (check_order(a,b,c)) {
          case 0: /* a > b > c */
            return_list
              .push_back(opm::triangle(d(*pA,*pB,*pC),r(a,c),0,pA,pB,pC));
            break;
          case 1: /* a > c > b */
            return_list
              .push_back(opm::triangle(d(*pA,*pC,*pB),r(a,b),1,pA,pC,pB));
            break;
          case 2: /* b > c > a */
            return_list
              .push_back(opm::triangle(d(*pB,*pC,*pA),r(b,a),2,pB,pC,pA));
            break;
          case 3: /* b > a > c */
            return_list
              .push_back(opm::triangle(d(*pB,*pA,*pC),r(b,c),3,pB,pA,pC));
            break;
          case 4: /* c > a > b */
            return_list
              .push_back(opm::triangle(d(*pC,*pA,*pB),r(c,b),4,pC,pA,pB));
            break;
          case 5: /* c > b > a */
            return_list
              .push_back(opm::triangle(d(*pC,*pB,*pA),r(c,a),5,pC,pB,pA));
            break;
          default:
            fprintf(stderr,"error: \n");
            exit(1);
            break;
          } /* switch (check_order(a,b,c)) */
        } /* for-loop k */
      } /* for-loop j */
    } /* for-loop i */

    sort(return_list.begin(),return_list.end(),
         [&](opm::triangle a, opm::triangle b){ return a.x*a.y > b.x*b.y; });
    return return_list;
  }


  database
  generate_database
  (const opm::referencelist::const_iterator &it_begin,
   const opm::referencelist::const_iterator &it_end)
  {
    std::vector<const opm::reference*> tmp_list;
    const size_t n = (size_t)(it_end-it_begin);


    for (auto it = it_begin; it != it_end; it++)
      tmp_list.push_back(&(*it));

    std::sort(tmp_list.begin(), tmp_list.end(),
              [&](const opm::reference* pA, const opm::reference* pB)
              { return pA->m < pB->m; });
    auto it = tmp_list.begin();

    if (n < 3) {
      fprintf(stderr,"error: number of stars not enough.\n");
      exit(1);
    }
    database return_database;

    for (size_t i=0; i<n; i++) {
      for (size_t j=i+1; j<n; j++) {
        for (size_t k=j+1; k<n; k++) {
          const opm::reference* pA = *(it+i);
          const opm::reference* pB = *(it+j);
          const opm::reference* pC = *(it+k);
          const opm::xym tA = {pA->xi(),pA->eta()};
          const opm::xym tB = {pB->xi(),pB->eta()};
          const opm::xym tC = {pC->xi(),pC->eta()};
          double c = sqdist(tA,tB);
          double a = sqdist(tB,tC);
          double b = sqdist(tC,tA);
          switch (check_order(a,b,c)) {
          case 0: /* a > b > c */
            return_database[0]
              .push_back(opm::triangle(d(tA,tB,tC),r(a,c),0,pA,pB,pC));
            break;
          case 1: /* a > c > b */
            return_database[1]
              .push_back(opm::triangle(d(tA,tC,tB),r(a,b),1,pA,pC,pB));
            break;
          case 2: /* b > c > a */
            return_database[2]
              .push_back(opm::triangle(d(tB,tC,tA),r(b,a),2,pB,pC,pA));
            break;
          case 3: /* b > a > c */
            return_database[3]
              .push_back(opm::triangle(d(tB,tA,tC),r(b,c),3,pB,pA,pC));
            break;
          case 4: /* c > a > b */
            return_database[4]
              .push_back(opm::triangle(d(tC,tA,tB),r(c,b),4,pC,pA,pB));
            break;
          case 5: /* c > b > a */
            return_database[5]
              .push_back(opm::triangle(d(tC,tB,tA),r(c,a),5,pC,pB,pA));
            break;
          default:
            fprintf(stderr,"error: \n");
            exit(1);
            break;
          } /* switch (check_order(a,b,c)) */
        } /* for-loop k */
      } /* for-loop j */
    } /* for-loop i */

    for (auto &triangle_array : return_database)
      sort(triangle_array.begin(),triangle_array.end(),
           [&](opm::triangle a, opm::triangle b){ return a > b; });
    return return_database;
  }

}
