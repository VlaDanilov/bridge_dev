//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/wrapper.h"
#include "qcdloop/qcdloop.h"
#include <stdexcept>
#include <iostream>
#include <omp.h>

extern"C" {

// result container
extern vector<complex> r;
#pragma omp threadprivate(r)
vector<complex> r(3);

extern vector<qcomplex> rq;
#pragma omp threadprivate(rq)
vector<qcomplex> rq(3);

// topologies
extern TadPole<complex,double,double> td;
#pragma omp threadprivate(td)
TadPole<complex,double,double> td;

extern TadPole<complex,complex,double> tdc;
#pragma omp threadprivate(tdc)
TadPole<complex,complex,double> tdc;


extern TadPole<qcomplex,qdouble,qdouble> tdq;
#pragma omp threadprivate(tdq)
TadPole<qcomplex,qdouble,qdouble> tdq;

extern TadPole<qcomplex,qcomplex,qdouble> tdcq;
#pragma omp threadprivate(tdcq)
TadPole<qcomplex,qcomplex,qdouble> tdcq;

extern Bubble<complex,double,double> bb;
#pragma omp threadprivate(bb)
Bubble<complex,double,double> bb;

extern Bubble<complex,complex,double> bbc;
#pragma omp threadprivate(bbc)
Bubble<complex,complex,double> bbc;

extern Bubble<qcomplex,qdouble,qdouble> bbq;
#pragma omp threadprivate(bbq)
Bubble<qcomplex,qdouble,qdouble> bbq;

extern Bubble<qcomplex,qcomplex,qdouble> bbcq;
#pragma omp threadprivate(bbcq)
Bubble<qcomplex,qcomplex,qdouble> bbcq;

extern Triangle<complex,double,double> tr;
#pragma omp threadprivate(tr)
Triangle<complex,double,double> tr;

extern Triangle<complex,complex,double> trc;
#pragma omp threadprivate(trc)
Triangle<complex,complex,double> trc;

extern Triangle<qcomplex,qdouble,qdouble> trq;
#pragma omp threadprivate(trq)
Triangle<qcomplex,qdouble,qdouble> trq;

extern Triangle<qcomplex,qcomplex,qdouble> trcq;
#pragma omp threadprivate(trcq)
Triangle<qcomplex,qcomplex,qdouble> trcq;

extern Box<complex,double,double> bo;
#pragma omp threadprivate(bo)
Box<complex,double,double> bo;

extern Box<complex,complex,double> boc;
#pragma omp threadprivate(boc)
Box<complex,complex,double> boc;

extern Box<qcomplex,qdouble,qdouble> boq;
#pragma omp threadprivate(boq)
Box<qcomplex,qdouble,qdouble> boq;

extern Box<qcomplex,qcomplex,qdouble> bocq;
#pragma omp threadprivate(bocq)
Box<qcomplex,qcomplex,qdouble> bocq;

// global vectors for improving parsing speed
extern vector<double>   mI1;
#pragma omp threadprivate(mI1)
vector<double>   mI1(1);

extern vector<complex>  mI1c;
#pragma omp threadprivate(mI1c)
vector<complex>  mI1c(1);

extern vector<qdouble>  mI1q;
#pragma omp threadprivate(mI1q)
vector<qdouble>  mI1q(1);

extern vector<qcomplex> mI1cq;
#pragma omp threadprivate(mI1cq)
vector<qcomplex> mI1cq(1);

extern vector<double>   mI2;
#pragma omp threadprivate(mI2)
vector<double>   mI2(2);

extern vector<complex>  mI2c;
#pragma omp threadprivate(mI2c)
vector<complex>  mI2c(2);

extern vector<qdouble>  mI2q;
#pragma omp threadprivate(mI2q)
vector<qdouble>  mI2q(2);

extern vector<qcomplex> mI2cq;
#pragma omp threadprivate(mI2cq)
vector<qcomplex> mI2cq(2);

extern vector<double>   mI3;
#pragma omp threadprivate(mI3)
vector<double>   mI3(3);

extern vector<complex>  mI3c;
#pragma omp threadprivate(mI3c)
vector<complex>  mI3c(3);

extern vector<qdouble>  mI3q;
#pragma omp threadprivate(mI3q)
vector<qdouble>  mI3q(3);

extern vector<qcomplex> mI3cq;
#pragma omp threadprivate(mI3cq)
vector<qcomplex> mI3cq(3);

extern vector<double>   mI4;
#pragma omp threadprivate(mI4)
vector<double>   mI4(4);

extern vector<complex>  mI4c;
#pragma omp threadprivate(mI4c)
vector<complex>  mI4c(4);

extern vector<qdouble>  mI4q;
#pragma omp threadprivate(mI4q)
vector<qdouble>  mI4q(4);

extern vector<qcomplex> mI4cq;
#pragma omp threadprivate(mI4cq)
vector<qcomplex> mI4cq(4);

extern vector<double>   pI2;
#pragma omp threadprivate(pI2)
vector<double>   pI2(1);

extern vector<qdouble>  pI2q;
#pragma omp threadprivate(pI2q)
vector<qdouble>  pI2q(1);

extern vector<double>   pI3;
#pragma omp threadprivate(pI3)
vector<double>   pI3(3);

extern vector<qdouble>  pI3q;
#pragma omp threadprivate(pI3q)
vector<qdouble>  pI3q(3);

extern vector<double>   pI4;
#pragma omp threadprivate(pI4)
vector<double>   pI4(6);

extern vector<qdouble>  pI4q;
#pragma omp threadprivate(pI4q)
vector<qdouble>  pI4q(6);

void qlcachesize(const int &size)
{
  td.setCacheSize(size);
  tdc.setCacheSize(size);
  tdq.setCacheSize(size);
  tdcq.setCacheSize(size);

  bb.setCacheSize(size);
  bbc.setCacheSize(size);
  bbq.setCacheSize(size);
  bbcq.setCacheSize(size);

  tr.setCacheSize(size);
  trc.setCacheSize(size);
  trq.setCacheSize(size);
  trcq.setCacheSize(size);

  bo.setCacheSize(size);
  boc.setCacheSize(size);
  boq.setCacheSize(size);
  bocq.setCacheSize(size);
}

bool qlzero(double const& x)
{
  return td.iszero(x);
}

bool qlzeroq(qdouble const& x)
{
  return tdq.iszero(x);
}

bool qlnonzero(double const& x)
{
  return !td.iszero(x);
}

bool qlnonzeroq(qdouble const& x)
{
  return !tdq.iszero(x);
}

complex cln(complex const& x, double const& isig)
{
  return td.cLn(x, isig);
}

qcomplex clnq(qcomplex const& x, qdouble const& isig)
{
  return tdq.cLn(x, isig);
}

void qltadpole(complex (&out)[3], double const& mu2, double const& m1)
{
  try
  {
    mI1[0] = m1;
    td.integral(r, mu2, mI1);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolec(complex (&out)[3], double const& mu2, complex const& m1)
{
  try
  {
    mI1c[0] = m1;
    tdc.integral(r, mu2, mI1c);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpoleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1)
{
  try
  {
    mI1q[0] = m1;
    tdq.integral(rq, mu2, mI1q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1)
{
  try
  {
    mI1cq[0] = m1;
    tdcq.integral(rq, mu2, mI1cq);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubble(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& p1)
{
  try
  {
    mI2[0] = m1;
    mI2[1] = m2;
    pI2[0] = p1;
    bb.integral(r, mu2, mI2, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblec(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, double const& p1)
{
  try
  {
    mI2c[0] = m1;
    mI2c[1] = m2;
    pI2[0]  = p1;
    bbc.integral(r, mu2, mI2c, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubbleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& p1)
{
  try
  {
    mI2q[0] = m1;
    mI2q[1] = m2;
    pI2q[0]  = p1;
    bbq.integral(rq, mu2, mI2q, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qdouble const& p1)
{
  try
  {
    mI2cq[0] = m1;
    mI2cq[1] = m2;
    pI2q[0]  = p1;
    bbcq.integral(rq, mu2, mI2cq, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangle(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3[0] = m1;
    mI3[1] = m2;
    mI3[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    tr.integral(r, mu2, mI3, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltrianglec(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3c[0] = m1;
    mI3c[1] = m2;
    mI3c[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    trc.integral(r, mu2, mI3c, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3q[0] = m1;
    mI3q[1] = m2;
    mI3q[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trq.integral(rq, mu2, mI3q, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}
void qltrianglecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3cq[0] = m1;
    mI3cq[1] = m2;
    mI3cq[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trcq.integral(rq, mu2, mI3cq, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbox(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4[0] = m1;
    mI4[1] = m2;
    mI4[2] = m3;
    mI4[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    bo.integral(r, mu2, mI4, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxc(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4c[0] = m1;
    mI4c[1] = m2;
    mI4c[2] = m3;
    mI4c[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    boc.integral(r, mu2, mI4c, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4q[0] = m1;
    mI4q[1] = m2;
    mI4q[2] = m3;
    mI4q[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    boq.integral(rq, mu2, mI4q, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxcq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4cq[0] = m1;
    mI4cq[1] = m2;
    mI4cq[2] = m3;
    mI4cq[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    bocq.integral(rq, mu2, mI4cq, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

#ifdef QL_NAMES

  // for backward compatibility.
  void qlinit()
  {
    std::cout << ql::yellow << "[QCDLoop warning]: this wrapper is not thread-safe." << std::endl;
    std::cout << "[QCDLoop suggestion]: consider developing object-oriented code." << ql::def << std::endl;
  }

  complex qli1(double const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1[0] = m1;
      td.integral(r, mu2, mI1);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli1c(complex const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1c[0] = m1;
      tdc.integral(r, mu2, mI1c);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli1q(qdouble const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1q[0] = m1;
      tdq.integral(rq, mu2, mI1q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli1qc(qcomplex const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1cq[0] = m1;
      tdcq.integral(rq, mu2, mI1cq);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli2(double const& p1, double const& m1, double const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2[0] = m1;
      mI2[1] = m2;
      pI2[0] = p1;
      bb.integral(r, mu2, mI2, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }    
    return r[abs(ep)];
  }

  complex qli2c(double const& p1, complex const& m1, complex const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2c[0] = m1;
      mI2c[1] = m2;
      pI2[0] = p1;
      bbc.integral(r, mu2, mI2c, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli2q(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2q[0] = m1;
      mI2q[1] = m2;
      pI2q[0] = p1;
      bbq.integral(rq, mu2, mI2q, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli2qc(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2cq[0] = m1;
      mI2cq[1] = m2;
      pI2q[0] = p1;
      bbcq.integral(rq, mu2, mI2cq, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli3(double const& p1, double const& p2, double const& p3, double const& m1, double const& m2, double const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3[0] = m1;
      mI3[1] = m2;
      mI3[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      tr.integral(r, mu2, mI3, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli3c(double const& p1, double const& p2, double const& p3, complex const& m1, complex const& m2, complex const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3c[0] = m1;
      mI3c[1] = m2;
      mI3c[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      trc.integral(r, mu2, mI3c, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli3q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3q[0] = m1;
      mI3q[1] = m2;
      mI3q[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trq.integral(rq, mu2, mI3q, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli3qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3cq[0] = m1;
      mI3cq[1] = m2;
      mI3cq[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trcq.integral(rq, mu2, mI3cq, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli4(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, double const& m1, double const& m2, double const& m3, double const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4[0] = m1;
      mI4[1] = m2;
      mI4[2] = m3;
      mI4[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      bo.integral(r, mu2, mI4, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }

    return r[abs(ep)];
  }

  complex qli4c(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4c[0] = m1;
      mI4c[1] = m2;
      mI4c[2] = m3;
      mI4c[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      boc.integral(r, mu2, mI4c, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli4q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4q[0] = m1;
      mI4q[1] = m2;
      mI4q[2] = m3;
      mI4q[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      boq.integral(rq, mu2, mI4q, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli4qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4cq[0] = m1;
      mI4cq[1] = m2;
      mI4cq[2] = m3;
      mI4cq[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      bocq.integral(rq, mu2, mI4cq, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

#endif

}
