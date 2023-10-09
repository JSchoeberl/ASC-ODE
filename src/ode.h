#ifndef ODE_h
#define ODE_h

#include <functional>
#include <calcinverse.hpp>


namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < 8; i++)
      {
        func->Evaluate(x, res);
        // cout << "|res| = " << L2Norm(res) << endl;
        func->EvaluateDeriv(x, fprime);
        CalcInverse(fprime);
        x -= fprime*res;
      }
  }


  
  void SolveODE_IE(double tend, double dt,
                   VectorView<double> x, shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xold = make_shared<ConstantFunction>(x);
    auto xnew = make_shared<IdenticFunction>(x.Size());
    auto equ = xnew-xold + dt * rhs;

    double t = 0;
    while (t < tend)
      {
        NewtonSolver (equ, x);
        xold->Set(x);
        if (callback) callback(t, x);
        t += dt;
      }
  }

  void SolveODE_Verlet(double tend, double dt,
                       VectorView<double> x, VectorView<double> dx,
                       shared_ptr<NonlinearFunction> rhs,   // x->f(x)
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xold = make_shared<ConstantFunction>(x);
    auto xoldold = make_shared<ConstantFunction>(Vector<double>(x-dt*dx));    
    auto xnew = make_shared<IdenticFunction>(x.Size());
    auto rhsold = make_shared<ConstantFunction>(x);        
    auto equ = xnew-2*xold+xoldold - dt*dt * rhsold;
    
    double t = 0;
    while (t < tend)
      {
        rhs->Evaluate(xold->Get(), rhsold->Get());

        NewtonSolver (equ, x);
        xoldold->Set(xold->Get());
        xold->Set(x);
        if (callback) callback(t, x);
        t += dt;
      }
    dx = 1/dt * (xold->Get()-xoldold->Get());
  }

  void SolveODE_Shake(double tend, double dt,
                      VectorView<double> xlam, size_t nx,
                      shared_ptr<NonlinearFunction> dLagrangef,
                      shared_ptr<NonlinearFunction> dLagrangeg,                      
                      std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xnew = make_shared<IdenticFunction>(xlam.Size());
    auto xold = make_shared<ConstantFunction>(xlam);
    auto xoldold = make_shared<ConstantFunction>(xlam);    
    auto dLagrangef_old = make_shared<ConstantFunction>(xlam);        
    
    auto ddx = make_shared<ProjectFunction> (xnew-2*xold+xoldold, 0, nx);
    auto equ = ddx - dt*dt*(dLagrangeg+dLagrangef_old);
    double t = 0;
    while (t < tend)
      {
        dLagrangef->Evaluate(xold->Get(), dLagrangef_old->Get());        
        NewtonSolver (equ, xlam);
        xoldold->Set(xold->Get());        
        xold->Set(xlam);
        if (callback) callback(t, xlam);
        t += dt;
      }
  }
  
}


#endif
