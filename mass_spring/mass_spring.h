#ifndef MASS_SPRING_H
#define MASS_SPRING_H



#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;

#include <vector.hpp>
using namespace ngbla;



template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
};

class Connection
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

ostream & operator<< (ostream & ost, const Connection & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  Connection connections[2];
};

template <int D>
class MassSpringSystem
{
  std::vector<Vec<D>> fixes;
  std::vector<Mass<D>> masses;
  std::vector<Spring> springs;
  Vec<D> gravity=0.0;
public:
  void SetGravity (Vec<D> _gravity) { gravity = _gravity; }
  Vec<D> Gravity() const { return gravity; }
  
  Connection AddFix (Vec<D> p)
  {
    fixes.push_back(p);
    return { Connection::FIX, fixes.size()-1 };
  }

  Connection AddMass (Mass<D> m)
  {
    masses.push_back (m);
    return { Connection::MASS, masses.size()-1 };
  }
  
  size_t AddSpring (double length, double stiffness, Connection c1, Connection c2)
  {
    springs.push_back (Spring{length, stiffness, { c1, c2 } });
    return springs.size()-1;
  }

  

  auto & Fixes() const { return fixes; } 
  auto & Masses() const { return masses; } 
  auto & Springs() const { return springs; } 
};

template <int D>
ostream & operator<< (ostream & ost, const MassSpringSystem<D> mss)
{
  ost << "fixes:" << endl;
  for (auto f : mss.Fixes())
    ost << f << endl;

  ost << "masses: " << endl;
  for (auto m : mss.Masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << endl;

  ost << "springs: " << endl;
  for (auto sp : mss.Springs())
    ost << "length = " << sp.length << "stiffness = " << sp.stiffness
        << ", C1 = " << sp.connections[0] << ", C2 = " << sp.connections[1] << endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  const MassSpringSystem<D> & mss;
public:
  MSS_Function (const MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t DimX() const { return D*mss.Masses().size(); }
  virtual size_t DimF() const { return D*mss.Masses().size(); }
  
  virtual void Evaluate (VectorView<double> x, VectorView<double> f) const
  {
    f = 0.0;
    
    auto xmat = x.AsMatrix(mss.Masses().size(), 2);
    auto fmat = f.AsMatrix(mss.Masses().size(), 2);
    
    for (size_t i = 0; i < mss.Masses().size(); i++)
      fmat.Row(i) = mss.Masses()[i].mass*mss.Gravity();
    
    for (auto spring : mss.Springs())
      {
        auto [c1,c2] = spring.connections;
        Vec<D> p1, p2;
        if (c1.type == Connection::FIX)
          p1 = mss.Fixes()[c1.nr];
        else
          p1 = xmat.Row(c1.nr);
        if (c2.type == Connection::FIX)
          p2 = mss.Fixes()[c2.nr];
        else
          p2 = xmat.Row(c2.nr);

        double force = spring.stiffness * (L2Norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/L2Norm(p1-p2) * (p2-p1);
        if (c1.type == Connection::MASS)
          fmat.Row(c1.nr) += force*dir12;
        if (c2.type == Connection::MASS)
          fmat.Row(c2.nr) -= force*dir12;
      }

    for (size_t i = 0; i < mss.Masses().size(); i++)
      fmat.Row(i) /= mss.Masses()[i].mass;
  }
  
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    // TODO: exact differentiation
    double eps = 1e-6;
    Vector<> xl(DimX()), xr(DimX()), fl(DimF()), fr(DimF());
    for (size_t i = 0; i < DimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        Evaluate (xl, fl);
        Evaluate (xr, fr);
        df.Col(i) = 1/(2*eps) * (fr-fl);
      }
  }
  
};

#endif
