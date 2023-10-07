#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  double dt = tend/100;
  Vector<> x { 1, };
  auto rhs = make_shared<RHS>();
  
  SolveODE_Verlet(tend, dt, x, rhs,
                  [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << endl; });
}
