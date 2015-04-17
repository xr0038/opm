/**
 * @file opm.h
 * @brief Optimistic Position Matching Algorithm
 *
 * @author Ryou Ohsawa
 * @date 2015
 */
#ifndef OPM_H
#define OPM_H

#include <iostream>
#include <vector>
#include <cstdio>
#include <array>
#include <algorithm>
#include <initializer_list>
#include <wcslib/wcs.h>

/**
 * @brief Optimistic Position Matching Algorithm
 * @note ここでは OPM_A (Tabur 2007) を定義する
 */
namespace opm
{
  /**
   * @brief 二次元空間での座標
   * @note @c x, @c y 座標はそれぞれ実数 (@c double) 値をとる
   */
  class xy {
  public:
    double x; /**< @brief 1 次元目の座標 */
    double y; /**< @brief 2 次元目の座標 */

    /** @brief デフォルトコンストラクタ / 座標は @c (0,0) */
    xy() : x(0), y(0) {}

    /** @brief 座標を引数にとるコンストラクタ */
    xy(double x_, double y_): x(x_), y(y_) {}
    
    /** @brief 座標を initializer_list で引数にとるコンストラクタ */
    xy(std::initializer_list<double> li)
    { 
      auto ptr = begin(li);
      x = *ptr;
      y = *++ptr;
    }

    /** @brief 内積を計算する */
    double
    operator*(const opm::xy &obj) const
    { return x * obj.x + y * obj.y; }

    /** @brief 座標のベクトル和を計算 */
    opm::xy
    operator+(const opm::xy &obj) const
    {
      opm::xy tmp(*this); 
      tmp.x += obj.x; tmp.y += obj.y;
      return tmp;
    }
    /** @brief 座標のベクトル差分を計算 */
    opm::xy
    operator-(const opm::xy &obj) const
    {
      opm::xy tmp(*this); 
      tmp.x -= obj.x; tmp.y -= obj.y;
      return tmp;
    }

    /** @brief @c x 座標で比較/@c x 座標が同じ時は @c y 座標で比較 */
    virtual bool
    operator<(const opm::xy &obj) const
    { return (x==obj.x)?(y < obj.y):(x < obj.x); }
    /** @brief @c x 座標で比較/@c x 座標が同じ時は @c y 座標で比較 */
    virtual bool
    operator>(const opm::xy &obj) const
    { return (x==obj.x)?(y > obj.y):(x > obj.x); }
    /** @brief どちらの座標も等しい時にのみ @c true */
    virtual bool
    operator==(const opm::xy &obj) const
    { return (x==obj.x) && (y==obj.y); }
  };

  /**
   * @brief 二次元空間での座標 + 等級
   * @note @c x, @c y 座標はそれぞれ整数 (@c size_t) 値をとる
   * @note @c value 値は @c double とする
   */
  class xym : public xy {
  public:
    double m;
    /** @breif デフォルトコンストラクタ */
    xym() : opm::xy(), m(0) {}
    /** @brief コピーコンストラクタ */
    xym(const xym& obj) : opm::xy(obj.x,obj.y), m(obj.m) {}
    /** @brief 座標と等級を明示するコンストラクタ */
    xym(double x_, double y_, double m_) : opm::xy(x_,y_), m(m_) {}
    /** @brief 初期化リストを用いたコンストラクタ */
    xym(std::initializer_list<double> li)
    { auto ptr = begin(li); x = *ptr; y = *++ptr; m = *++ptr; }

    /** @brief インスタンスのメンバを出力する debug 用関数 */
    void
    dump_data() const
    { printf("x:%8.3f, y:%8.3f, m:%12.5f\n", x, y, m); }

    /** @brief magnitude @c m で比較 */
    bool
    operator<(const opm::xym &obj) const
    { return m < obj.m; }
    /** @brief magnitude @c m で比較 */
    bool
    operator>(const opm::xym &obj) const
    { return m > obj.m; }
    /** @brief すべてが等しい場合のみ @c true */
    bool
    operator==(const opm::xym &obj) const
    { return (x==obj.x) && (y==obj.y) && (m==obj.m); }

    /** @brief @c xy とのベクトル和を定義 */
    opm::xym
    operator+(const opm::xy &obj) const
    { 
      opm::xym tmp(*this);
      tmp.x += obj.x; tmp.y += obj.y;
      return tmp;
    }
    /** @brief @c xy とのベクトル差を定義 */
    opm::xym
    operator-(const opm::xy &obj) const
    { 
      opm::xym tmp(*this);
      tmp.x -= obj.x; tmp.y -= obj.y;
      return tmp;
    }
  };

  /**
   * @brief fits 画像上の天体を表すクラス
   */
  class object : public xym {
  public:
    object(const xym& obj) : opm::xym(obj) {}
    object(const double x_, const double y_, const double m_)
      : opm::xym(x_,y_,m_) {}
    /** @brief @c xy とのベクトル和を定義 */
    opm::object
    operator+(const opm::xy &obj) const
    { 
      opm::object tmp(*this);
      tmp.x += obj.x; tmp.y += obj.y;
      return tmp;
    }
    /** @brief @c xy とのベクトル差を定義 */
    opm::object
    operator-(const opm::xy &obj) const
    { 
      opm::object tmp(*this);
      tmp.x -= obj.x; tmp.y -= obj.y;
      return tmp;
    }
  };


  /**
   * @brief カタログ上の天体を表すクラス
   *
   */
  class reference : public xym {
  public:
    /** @brief 座標と等級だけ指定するコンストラクタ */
    reference(double x_, double y_, double m_) : xym(x_,y_,m_)
    {
      init_celestial(NULL);
    }
    /** @brief 全メンバ変数を指定するコンストラクタ */
    reference(double x_, double y_, double m_, 
              wcsprm* pwcs)
      : xym(x_,y_,m_)
    {
      init_celestial(pwcs);
    }
    /** @brief コピーコンストラクタ */
    reference(const reference& obj)
      : xym(obj.x,obj.y,obj.m)
    {
      px_ = obj.px_;
      py_ = obj.py_;
      xi_  = obj.xi_;
      eta_ = obj.eta_;
    }
    /** @brief @c xym をコピーするコンストラクタ */
    reference(const xym& obj,
              wcsprm* pwcs)
      : opm::xym(obj)
    {
      init_celestial(pwcs);
    }
    double
    px() const
    { return px_; }
    double
    py() const
    { return py_; }
    double
    xi() const
    { return xi_; }
    double
    eta() const
    { return eta_; }
  private:
    double px_;
    double py_;
    double xi_;
    double eta_;
    void
    init_celestial(wcsprm* pwcs)
    {
      if (pwcs == NULL) {
        px_ = x; py_ = y;
        xi_ = x; eta_ = y;
      } else {
        int status;
        double imgcrd[2], pixcrd[2], phi[1], theta[1];
        double world[2] = {this->x, this->y};
        wcss2p(pwcs, 1, 2, world, phi, theta, imgcrd, pixcrd, &status);
        px_ = pixcrd[0];
        py_ = pixcrd[1];
        xi_  = imgcrd[0];
        eta_ = imgcrd[1];
      }
    }
  };

  /** @brief define objectlist by vector<object> */
  typedef std::vector<opm::object> objectlist;
  /** @brief define objectlist by vector<object> */
  typedef std::vector<opm::reference> referencelist;

  /**
   * @brief 3 天体の作る三角形を表現するクラス
   *
   * メンバ変数として @c x, @c y, @c n, @c v を持つ．
   * @c (x,y) は triangle space での座標である．
   * @c n は天体の明るさと辺の長さの対応関係である．
   * @c v は三角形を構成する頂点が対応する辺の長さ順に含まれる．
   */
  class triangle : public opm::xy {
  public:
    /** @brief 天体の明るさと辺の長さの対応関係 */
    size_t n;
    /** @brief 三角形を構成する頂点(ポインタ) */
    std::array<const opm::xym*,3> v;
    /**
     * @brief コンストラクタ
     * @param x_ 最短辺に対応する頂点からみたベクトルの内積
     * @param y_ 三角形の最長辺と最短辺の比 (a/c)
     * @param n_ 天体の明るさと辺の長さの対応
     * @param A 最も長い辺 (a) に対応する頂点
     * @param B 長さが真ん中の辺 (b) に対応する頂点
     * @param C 最も短い辺 (c) に対応する頂点
     */
    triangle(double x_, double y_, size_t n_,
             const opm::xym* pA, const opm::xym* pB, const opm::xym* pC)
      : opm::xy(x_,y_), n(n_), v({pA,pB,pC}) {}

    /**
     * @brief 2 つのインスタンスがマンハッタン距離 @c d で近傍にいるか調べる
     * @param obj 調べるインスタンス
     * @param dx @c x 方向の近傍距離 (x*dx) を定義する値
     * @param dy @c y 方向の近傍距離 (y*dy) を定義する値
     * @return @c obj1 と @c obj2 が @c d で近傍なら @c true
     */
    bool
    neighbor
    (const opm::triangle &obj, const double dx, const double dy) const
    { return (std::abs(x-obj.x) < x*dx) && (std::abs(y-obj.y) < y*dy); }

    /** @brief @c y 座標でソートする */
    bool
    operator<(const triangle &obj) const
    { return (y==obj.y)? x < obj.x : y < obj.y; }

    /** @brief @c y 座標でソートする */
    bool
    operator>(const triangle &obj) const
    { return (y==obj.y)? x > obj.x : y > obj.y; }

    opm::triangle
    operator+(const opm::triangle &obj) const
    { 
      opm::triangle tmp(*this);
      tmp.x += obj.x; tmp.y += obj.y;
      return tmp;
    }
    opm::triangle
    operator-(const opm::triangle &obj) const
    { 
      opm::triangle tmp(*this);
      tmp.x -= obj.x; tmp.y -= obj.y;
      return tmp;
    }
    opm::triangle
    operator+(const opm::xy &obj) const
    { 
      opm::triangle tmp(*this);
      tmp.x += obj.x; tmp.y += obj.y;
      return tmp;
    }
    opm::triangle
    operator-(const opm::xy &obj) const
    { 
      opm::triangle tmp(*this);
      tmp.x -= obj.x; tmp.y -= obj.y;
      return tmp;
    }
  };

  /** @brief define triangles by vector<triangle> */
  typedef std::vector<opm::triangle> triangles;
  /** @brief define database by array<triangles,6> */
  typedef std::array< triangles, 6 > database;

  /** @brief define matched_star by pair<object,reference> */
  typedef std::pair<opm::object, opm::reference> matched_star;
  /** @brief define matched_list by vector<matched_star> */
  typedef std::vector<matched_star> matched_list;

  /**
   * @brief @c xym のベクタから @c ojectlist を生成する
   */
  opm::objectlist
  generate_objectlist(const std::vector<opm::xym> &obj);
  /**
   * @brief @c xym のベクタから @c referencelist を生成する
   */
  opm::referencelist
  generate_referencelist(const std::vector<opm::xym> &obj,
                         wcsprm* pwcs);
  /**
   * @brief @c objectlist から triangle list を生成する
   * @note 生成後に @c objectlist をソートしてはいけない
   * @param it_begin triangle list を生成する @c objectlist の先頭イテレータ
   * @param it_end triangle list を生成する @c objectlist の末尾イテレータ
   * @return 生成した triangle list
   */
  opm::triangles
  generate_triangles
  (const opm::objectlist::const_iterator &it_begin,
   const opm::objectlist::const_iterator &it_end);
  /**
   * @brief @c objectlist から triangle list を生成する
   * @note 生成後に @c objectlist をソートしてはいけない
   * @param ol triangle list を生成する @c objectlist
   * @return 生成した triangle list
   */
  inline opm::triangles
  generate_triangles(const opm::objectlist &ol)
  { return generate_triangles(ol.cbegin(), ol.cend()); }
  /**
   * @brief @c objectlist から triangle list を生成する
   * @note 生成後に @c objectlist をソートしてはいけない
   * @param ol triangle list を生成する @c objectlist
   * @param s triangle list を生成する @c object の先頭位置
   * @param n triangle list を生成する @c object の総数
   * @return 生成した triangle list
   */
  inline opm::triangles
  generate_triangles(const opm::objectlist &ol,
                     const size_t s, const size_t n)
  { return generate_triangles(ol.cbegin()+s, ol.cbegin()+s+n); }

  /**
   * @brief @c referencelist から triangle database を生成する
   * @note 生成後に @c referencelist をソートしてはいけない
   * @param it_begin database を生成する @c referencelist の先頭イテレータ
   * @param it_end database を生成する @c referencelist の末尾イテレータ
   * @return 生成した triangle database
   */
  opm::database
  generate_database
  (const opm::referencelist::const_iterator &it_begin,
   const opm::referencelist::const_iterator &it_end);
  /**
   * @brief @c referencelist から triangle database を生成する
   * @note 生成後に @c referencelist をソートしてはいけない
   * @param rl database を生成する @c referencelist
   * @return 生成した triangle database
   */
  inline opm::database
  generate_database(const opm::referencelist &rl)
  { return generate_database(rl.cbegin(), rl.cend()); }
  /**
   * @brief @c referencelist から triangle database を生成する
   * @note 生成後に @c referencelist をソートしてはいけない
   * @param rl database を生成する @c referencelist
   * @param s triangle list を生成する @c object の先頭位置
   * @param n triangle list を生成する @c object の総数
   * @return 生成した triangle database
   */
  inline opm::database
  generate_database(const opm::referencelist &rl,
                    const size_t s, const size_t n)
  { return generate_database(rl.cbegin()+s, rl.cbegin()+s+n); }

  /**
   * @brief 座標の線形変換を定義するセット
   *
   * 変換は以下の式で定義される
   * @verbatim
   x' = R_11*x + R_12*y + T_x,
   y' = R_21*x + R_22*y + T_y
   * @endverbatim
   * メンバは {R_11, R_12, R_21, R_22, T_x, T_y} の順
   */
  typedef std::array<double, 6> conversion;
  
  /**
   * @brief prematch に使用する関数および構造体など
   * 
   * prematch では triangle space 上でマッチングを行なう．
   * 観測した位置からリファレンス画像への変換式のラフな推定を行なう．
   */
  namespace pre {
    /** @brief @c check_conversion で使うデフォルトの tolerance 値 */
    constexpr double tol_conversion = 0.01;

    /**
     * @brief triangle space 上でマッチングして変換係数を推定する
     * @param T 近傍を検索する triangle インスタンス
     * @param dx 検索する @c x 座標方向の距離 (T.x * dx) を定める [0..1]
     * @param dy 検索する @c y 座標方向の距離 (T.y * dy) を定める [0..1]
     * @param tol conversion のチェックで使用する tolerance 値
     * @param src triangle インスタンスを検索する triangle database
     * @return prematch で得られた obj -> ref 変換の係数のリスト
     */
    std::vector<opm::conversion>
    match(const opm::triangle &T,
          const double dx, const double dy, const double tol,
          const opm::database &src);
    inline std::vector<opm::conversion>
    match(const opm::triangle &T,
          const double dx, const double dy,
          const opm::database &src)
    { return opm::pre::
        match(T, dx, dy, opm::pre::tol_conversion, src); }
  }

  /**
   * @brief final match に使用する関数および構造体など
   * 
   * final match では天球面上でマッチングを行なう．
   */
  namespace final {
    opm::matched_list
    match(const opm::conversion coeff,
          const opm::objectlist &obj, const double dx, const double dy,
          const opm::referencelist &src);
    inline opm::matched_list
    match(const opm::conversion coeff,
          const opm::objectlist &obj, const double d,
          const opm::referencelist &src)
    { return opm::final::match(coeff, obj, d, d, src); }
  }

  /**
   * @brief matched list から wcs 情報をアップデートする
   * 
   * アップデートされる情報は
   * @li CRPIXn - レファレンスピクセルの位置
   * @li PCi_j  - ピクセル座標から中間座標への変換行列
   * @param pwcs 更新すべき wcsprm 構造体へのポインタ
   * @param matched マッチングした天体のリスト
   * @return 推定した変換係数リスト
   */
  conversion
  update_wcsprm(wcsprm* pwcs, const opm::matched_list matched);

}

#endif /* OPM_H */
