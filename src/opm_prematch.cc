/**
 * @file opm_prematch.cc
 * @brief Optimistic Position Matching Algorithm
 *
 * @author Ryou Ohsawa
 * @date 2015
 */
#include "opm.h"

namespace opm {

  namespace {
    
    /**
     * @brief triangle list から @c y 座標で @c d 近傍にある triangle を検索
     * @note @c d は現在の @c y 座標に対する割合
     * @param T 近傍を検索する triangle インスタンス
     * @param d 検索する @c y 座標方向の距離 (T.y * d) を定める [0..1]
     * @param src triangle インスタンスを検索する triangle list
     * @return 近傍に存在する triangle インスタンスのリスト
     */
    inline triangles
    sift_by_triangleY
    (const opm::triangle &T, const double d, const opm::triangles &src)
    {
      opm::xy dy(0,T.y*d);
      auto lb = std::lower_bound(src.begin(),src.end(),T+dy,
                                 [&](opm::triangle a, opm::triangle b)
                                 { return a > b; });
      auto ub = std::upper_bound(src.begin(),src.end(),T-dy,
                                 [&](opm::triangle a, opm::triangle b)
                                 { return a > b; });

      opm::triangles retval(lb,ub);
      return retval;
    }
    /**
     * @brief triangle list から @c y 座標で @c d 近傍にある triangle を検索
     * @note @c d は現在の @c y 座標に対する割合
     * @param T 近傍を検索する triangle インスタンス
     * @param d 検索する @c y 座標方向の距離 (T.y * d) を定める [0..1]
     * @param src triangle インスタンスを検索する triangle database
     * @return 近傍に存在する triangle インスタンスのリスト
     */
    inline triangles
    sift_by_triangleY
    (const opm::triangle &T, const double d, const opm::database &src)
    {
      const opm::triangles &list = src[T.n];
      return opm::sift_by_triangleY(T, d, list);
    }

    /**
     * @brief triangle list から近傍にあるインスタンスをふるいにかける
     * @note @c dx, @c dy は @c x, @c y 座標に対する割合
     * @param T 近傍を検索する triangle インスタンス
     * @param dx 検索する @c x 座標方向の距離 (T.x * dx) を定める [0..1]
     * @param dy 検索する @c y 座標方向の距離 (T.y * dy) を定める [0..1]
     * @param src triangle インスタンスを検索する triangle list
     * @return 近傍に存在する triangle インスタンスのリスト
     */
    inline triangles
    sift
    (const opm::triangle &T, const double dx, const double dy,
     const opm::triangles &src)
    {
      opm::triangles retval;
      opm::triangles &&list = opm::sift_by_triangleY(T, dy, src);
      for (const auto &t: list) {
        if (T.neighbor(t,dx,dy)) retval.push_back(t);
      }
      return retval;
    }
    /**
     * @brief triangle list から近傍にあるインスタンスをふるいにかける
     * @note @c dx, @c dy は @c x, @c y 座標に対する割合
     * @param T 近傍を検索する triangle インスタンス
     * @param dx 検索する @c x 座標方向の距離 (T.x * dx) を定める [0..1]
     * @param dy 検索する @c y 座標方向の距離 (T.y * dy) を定める [0..1]
     * @param src triangle インスタンスを検索する triangle database
     * @return 近傍に存在する triangle インスタンスのリスト
     */
    inline triangles
    sift
    (const opm::triangle &T, const double dx, const double dy,
     const opm::database &src)
    {
      const opm::triangles &list = src[T.n];
      return opm::sift(T, dx, dy, list);
    }

    /**
     * @brief 2 つの triangle インスタンスから変換を推定する
     * @param ref 基準となる triangle インスタンス
     * @param obj 変換元となる triangle インスタンス
     * @return @c obj -> @c ref への変換係数
     */
    opm::conversion
    estimate_conversion
    (const opm::triangle &ref, const opm::triangle &obj);

    /**
     * @brief 与えられた conversion が妥当なものかチェックする
     *
     * plate scale が適切に求まっており画像に大きな歪みがなければ
     * R_11 ~ R_22, R12 ~ -R_21, R_11*R_22-R_12*R_21 ~ 1
     * @param obj チェックする conversion インスタンス
     * @param tol conversion の値をチェックするときのしきい値
     * @return conversion がまともなものであれば @c true
     */
    bool
    check_conversion(const opm::conversion &obj, const double tol);

    opm::conversion
    estimate_conversion
    (const opm::triangle &ref, const opm::triangle &obj)
    {
      const double &x0 = obj.v[0]->x, &y0 = obj.v[0]->y;
      const double &x1 = obj.v[1]->x, &y1 = obj.v[1]->y;
      const double &x2 = obj.v[2]->x, &y2 = obj.v[2]->y;
      const double &&x0_ref = ((opm::reference*)ref.v[0])->xi();
      const double &&y0_ref = ((opm::reference*)ref.v[0])->eta();
      const double &&x1_ref = ((opm::reference*)ref.v[1])->xi();
      const double &&y1_ref = ((opm::reference*)ref.v[1])->eta();
      const double &&x2_ref = ((opm::reference*)ref.v[2])->xi();
      const double &&y2_ref = ((opm::reference*)ref.v[2])->eta();

      double D = x0*(y2-y1) + y0*(x1-x2) + y1*x2 - x1*y2;
      double R11 = ( (y2-y1)*x0_ref + (y0-y2)*x1_ref + (y1-y0)*x2_ref ) / D;
      double R12 = ( (x1-x2)*x0_ref + (x2-x0)*x1_ref + (x0-x1)*x2_ref ) / D;
      double R21 = ( (y2-y1)*y0_ref + (y0-y2)*y1_ref + (y1-y0)*y2_ref ) / D;
      double R22 = ( (x1-x2)*y0_ref + (x2-x0)*y1_ref + (x0-x1)*y2_ref ) / D;
      double Tx  = ( (x2*y1-x1*y2)*x0_ref + 
                     (x0*y2-x2*y0)*x1_ref + (x1*y0-x0*y1)*x2_ref ) / D;
      double Ty  = ( (x2*y1-x1*y2)*y0_ref + 
                     (x0*y2-x2*y0)*y1_ref + (x1*y0-x0*y1)*y2_ref ) / D;

      return conversion({R11,R12,R21,R22,Tx,Ty});
    }

    bool
    check_conversion(const opm::conversion &obj, const double tol)
    {
      double &&cond1 = std::abs(obj[0]-obj[3]);
      double &&cond2 = std::abs(obj[1]+obj[2]);
      double &&cond3 = std::abs(obj[0]*obj[3]-obj[1]*obj[2]-1.0);
#ifdef DEBUG
      if (cond1 < tol && cond2 < tol && cond3 < tol) {
        fprintf(stderr,"(R11,R12,R21,R22,Tx,Ty)"
                " = (%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf)\n",
                obj[0],obj[1],obj[2],obj[3],obj[4],obj[5]);
        fprintf(stderr, "|R11 - R22|         = %lf\n", cond1);
        fprintf(stderr, "|R12 + R21|         = %lf\n", cond2);
        fprintf(stderr, "|R11xR22-R12xR21-1| = %lf\n", cond3);
      }
#endif
      return (cond1 < tol && cond2 < tol && cond3 < tol)? true : false;
    }
  }

  std::vector<opm::conversion>
  pre::match(const opm::triangle &T, const double dx, const double dy,
             const double tol, const opm::database &src)
  {
    std::vector<opm::conversion> retval;
    opm::triangles sifted = opm::sift(T, dx, dy, src);

    for (auto &r : sifted) {
      opm::conversion c = opm::estimate_conversion(r, T);
      if (opm::check_conversion(c, tol)) {
        retval.push_back(c);
      }
    }
    return retval;
  }
}
