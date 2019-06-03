#include "qcdloop/qcdloop.h"
#include "qcdloop/types.h"
#include <stdexcept>
#include <iostream>

using namespace ql;

extern "C" {

  complex cln(complex const& x, double const& isig)
  {
    TadPole<complex,double,double> loc_td;
    return loc_td.cLn(x, isig);
  }

  bool qlzero(double const& x)
  {
    TadPole<complex,double,double> loc_td;
    return loc_td.iszero(x);
  }

  bool qlnonzero(double const& x)
  {
    TadPole<complex,double,double> loc_td;
    return !loc_td.iszero(x);
  }

  complex qli1(double const& m1, double const& mu2, int const& ep)
  {
    TadPole<complex,double,double> loc_td;
    vector<complex> loc_r(3);
    vector<double> loc_mI1(1);
    try
    {
      loc_mI1[0] = m1;
      loc_td.integral(loc_r, mu2, loc_mI1);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return loc_r[abs(ep)];
  }

  complex qli1c(complex const& m1, double const& mu2, int const& ep)
  {
    TadPole<complex,complex,double> loc_td;
    vector<complex> loc_r(3);
    vector<complex> loc_mI1(1);
    try
    {
      loc_mI1[0] = m1;
      loc_td.integral(loc_r, mu2, loc_mI1);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return loc_r[abs(ep)];
  }

  complex qli2(double const& p1, double const& m1, double const& m2,
      double const& mu2, int const& ep)
  {
    Bubble<complex,double,double> loc_bb;
    vector<complex> loc_r(3);
    vector<double> loc_mI2(2);
    vector<double> loc_pI2(1);

    try
    {
      loc_mI2[0] = m1;
      loc_mI2[1] = m2;
      loc_pI2[0] = p1;
      loc_bb.integral(loc_r, mu2, loc_mI2, loc_pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }    
    return loc_r[abs(ep)];
  }

  complex qli2c(double const& p1, complex const& m1, complex const& m2,
      double const& mu2, int const& ep)
  {
    Bubble<complex,complex,double> loc_bb;
    vector<complex> loc_r(3);
    vector<complex> loc_mI2(2);
    vector<double> loc_pI2(1);

    try
    {
      loc_mI2[0] = m1;
      loc_mI2[1] = m2;
      loc_pI2[0] = p1;
      loc_bb.integral(loc_r, mu2, loc_mI2, loc_pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return loc_r[abs(ep)];
  }

  complex qli3(double const& p1, double const& p2, double const& p3, double const& m1, double const& m2, double const& m3, double const& mu2, int const& ep)
  {
    Triangle<complex,double,double> loc_tr;
    vector<complex> loc_r(3);
    vector<double> loc_mI3(3);
    vector<double> loc_pI3(3);

    try
    {
      loc_mI3[0] = m1;
      loc_mI3[1] = m2;
      loc_mI3[2] = m3;
      loc_pI3[0] = p1;
      loc_pI3[1] = p2;
      loc_pI3[2] = p3;
      loc_tr.integral(loc_r, mu2, loc_mI3, loc_pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return loc_r[abs(ep)];
  }

  complex qli3c(double const& p1, double const& p2, double const& p3, complex const& m1, complex const& m2, complex const& m3, double const& mu2, int const& ep)
  {
    Triangle<complex,complex,double> loc_tr;
    vector<complex> loc_r(3);
    vector<complex> loc_mI3(3);
    vector<double> loc_pI3(3);

    try
    {
      loc_mI3[0] = m1;
      loc_mI3[1] = m2;
      loc_mI3[2] = m3;
      loc_pI3[0] = p1;
      loc_pI3[1] = p2;
      loc_pI3[2] = p3;
      loc_tr.integral(loc_r, mu2, loc_mI3, loc_pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return loc_r[abs(ep)];
  }

  complex qli4(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, double const& m1, double const& m2, double const& m3, double const& m4, double const& mu2, int const& ep)
  {
    Box<complex,double,double> loc_bo;
    vector<complex> loc_r(3);
    vector<double> loc_mI4(4);
    vector<double> loc_pI4(6);

    try
    {
      loc_mI4[0] = m1;
      loc_mI4[1] = m2;
      loc_mI4[2] = m3;
      loc_mI4[3] = m4;
      loc_pI4[0] = p1;
      loc_pI4[1] = p2;
      loc_pI4[2] = p3;
      loc_pI4[3] = p4;
      loc_pI4[4] = s12;
      loc_pI4[5] = s23;
      loc_bo.integral(loc_r, mu2, loc_mI4, loc_pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    
    return loc_r[abs(ep)];

}

  complex qli4c(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& mu2, int const& ep)
  {
    Box<complex,complex,double> loc_bo;
    vector<complex> loc_r(3);
    vector<complex> loc_mI4(4);
    vector<double> loc_pI4(6);

    try
    {
      loc_mI4[0] = m1;
      loc_mI4[1] = m2;
      loc_mI4[2] = m3;
      loc_mI4[3] = m4;
      loc_pI4[0] = p1;
      loc_pI4[1] = p2;
      loc_pI4[2] = p3;
      loc_pI4[3] = p4;
      loc_pI4[4] = s12;
      loc_pI4[5] = s23;
      loc_bo.integral(loc_r, mu2, loc_mI4, loc_pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }

    return loc_r[abs(ep)];

}

  void qli1q(qdouble const& m1, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    TadPole<qcomplex,qdouble,qdouble> loc_tdq;
    vector<qcomplex> loc_rq(3);
    vector<qdouble> loc_mI1q(1);
    try
    {
      loc_mI1q[0] = m1;
      loc_tdq.integral(loc_rq, mu2, loc_mI1q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

  void qli2q(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Bubble<qcomplex,qdouble,qdouble> loc_bbq;
    vector<qcomplex> loc_rq(3);
    vector<qdouble> loc_mI2q(2);
    vector<qdouble> loc_pI2q(1);

    try
    {
      loc_mI2q[0] = m1;
      loc_mI2q[1] = m2;
      loc_pI2q[0] = p1;
      loc_bbq.integral(loc_rq, mu2, loc_mI2q, loc_pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret =loc_rq[abs(ep)];
  }

  void qli3q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Triangle<qcomplex,qdouble,qdouble> loc_trq;
    vector<qcomplex> loc_rq(3);
    vector<qdouble> loc_mI3q(3);
    vector<qdouble> loc_pI3q(3);

    try
    {
      loc_mI3q[0] = m1;
      loc_mI3q[1] = m2;
      loc_mI3q[2] = m3;
      loc_pI3q[0] = p1;
      loc_pI3q[1] = p2;
      loc_pI3q[2] = p3;
      loc_trq.integral(loc_rq, mu2, loc_mI3q, loc_pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

  void qli4q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Box<qcomplex,qdouble,qdouble> loc_boq;
    vector<qcomplex> loc_rq(3);
    vector<qdouble> loc_mI4q(4);
    vector<qdouble> loc_pI4q(6);

    try
    {
      loc_mI4q[0] = m1;
      loc_mI4q[1] = m2;
      loc_mI4q[2] = m3;
      loc_mI4q[3] = m4;
      loc_pI4q[0] = p1;
      loc_pI4q[1] = p2;
      loc_pI4q[2] = p3;
      loc_pI4q[3] = p4;
      loc_pI4q[4] = s12;
      loc_pI4q[5] = s23;
      loc_boq.integral(loc_rq, mu2, loc_mI4q, loc_pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

  void qli1qc(qcomplex const& m1, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    TadPole<qcomplex,qcomplex,qdouble> loc_tdq;
    vector<qcomplex> loc_rq(3);
    vector<qcomplex> loc_mI1q(1);
    try
    {
      loc_mI1q[0] = m1;
      loc_tdq.integral(loc_rq, mu2, loc_mI1q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

  void qli2qc(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Bubble<qcomplex,qcomplex,qdouble> loc_bbq;
    vector<qcomplex> loc_rq(3);
    vector<qcomplex> loc_mI2q(2);
    vector<qdouble> loc_pI2q(1);

    try
    {
      loc_mI2q[0] = m1;
      loc_mI2q[1] = m2;
      loc_pI2q[0] = p1;
      loc_bbq.integral(loc_rq, mu2, loc_mI2q, loc_pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret =loc_rq[abs(ep)];
  }

  void qli3qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Triangle<qcomplex,qcomplex,qdouble> loc_trq;
    vector<qcomplex> loc_rq(3);
    vector<qcomplex> loc_mI3q(3);
    vector<qdouble> loc_pI3q(3);

    try
    {
      loc_mI3q[0] = m1;
      loc_mI3q[1] = m2;
      loc_mI3q[2] = m3;
      loc_pI3q[0] = p1;
      loc_pI3q[1] = p2;
      loc_pI3q[2] = p3;
      loc_trq.integral(loc_rq, mu2, loc_mI3q, loc_pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

  void qli4qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& mu2, int const& ep, qcomplex& ret)
  {
    Box<qcomplex,qcomplex,qdouble> loc_boq;
    vector<qcomplex> loc_rq(3);
    vector<qcomplex> loc_mI4q(4);
    vector<qdouble> loc_pI4q(6);

    try
    {
      loc_mI4q[0] = m1;
      loc_mI4q[1] = m2;
      loc_mI4q[2] = m3;
      loc_mI4q[3] = m4;
      loc_pI4q[0] = p1;
      loc_pI4q[1] = p2;
      loc_pI4q[2] = p3;
      loc_pI4q[3] = p4;
      loc_pI4q[4] = s12;
      loc_pI4q[5] = s23;
      loc_boq.integral(loc_rq, mu2, loc_mI4q, loc_pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    ret = loc_rq[abs(ep)];
  }

}
