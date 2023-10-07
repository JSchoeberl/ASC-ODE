#include "mass_spring.h"


int main()
{
  MassSpringSystem<2> mss;
  mss.SetGravity( {0,-9.81} );
  auto fA = mss.AddFix( { 0.0, 0.0 } );
  auto mA = mss.AddMass( { 1, { 1.0, 0.0 } } );
  mss.AddSpring (1, 10, fA, mA );

  auto mB = mss.AddMass( { 1, { 2.0, 0.0 } } );
  mss.AddSpring (1, 20, mA, mB);
  
  cout << "mss: " << endl << mss << endl;

  auto mss_func = make_shared<MSS_Function<2>> (mss);



  double tend = 10;
  double dt = tend/1000;
  
  Vector<> x(2*mss.Masses().size());
  for (size_t i = 0; i < mss.Masses().size(); i++)
    x.Range(2*i,2*i+2) = mss.Masses()[i].pos;
  
  SolveODE_Verlet(tend, dt, x, mss_func,
                  [](double t, VectorView<double> x) { cout << "t = " << t
                                                            << ", x = " << Vec<4>(x) << endl; });
}
