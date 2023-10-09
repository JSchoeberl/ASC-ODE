#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mass_spring.h"

namespace py = pybind11;


PYBIND11_MODULE(mass_spring, m) {
    m.doc() = "mass-spring-system simulator"; 

    py::class_<Mass<2>> (m, "Mass2d")
      ;
      
    m.def("Mass", [](double m, std::array<double,2> p)
    {
      return Mass<2>{m, { p[0], p[1] }};
    });


    

    
    py::class_<Mass<3>> (m, "Mass3d")
      .def_property_readonly("mass",
                             [](Mass<3> & m) { return m.mass; })
      .def_property_readonly("pos",
                             [](Mass<3> & m) { return m.pos; });
    ;
    
    m.def("Mass", [](double m, std::array<double,3> p)
    {
      return Mass<3>{m, { p[0], p[1], p[2] }};
    });

    
    py::class_<Fix<3>> (m, "Fix3d")
      .def_property_readonly("pos",
                             [](Fix<3> & f) { return f.pos; });
    

    m.def("Fix", [](std::array<double,3> p)
    {
      return Fix<3>{ { p[0], p[1], p[2] } };
    });

    py::class_<Connector> (m, "Connector");

    py::class_<Spring> (m, "Spring")
      .def(py::init<double, double, std::array<Connector,2>>())
      ;

    
    
    py::class_<MassSpringSystem<2>> (m, "MassSpringSystem2d")
      .def(py::init<>())
      .def("Add", [](MassSpringSystem<2> & mss, Mass<2> m) { return mss.AddMass(m); })
      ;
      
        
    py::class_<MassSpringSystem<3>> (m, "MassSpringSystem3d")
      .def(py::init<>())
      .def("__str__", [](MassSpringSystem<3> & mss) {
        stringstream sstr;
        sstr << mss;
        return sstr.str();
      })
      .def_property("gravity", [](MassSpringSystem<3> & mss) { return mss.Gravity(); },
                    [](MassSpringSystem<3> & mss, std::array<double,3> g) { mss.SetGravity(Vec<3>{g[0],g[1],g[2]}); })
      .def("Add", [](MassSpringSystem<3> & mss, Mass<3> m) { return mss.AddMass(m); })
      .def("Add", [](MassSpringSystem<3> & mss, Fix<3> f) { return mss.AddFix(f); })
      .def("Add", [](MassSpringSystem<3> & mss, Spring s) { return mss.AddSpring(s); })            
      .def_property_readonly("masses", [](MassSpringSystem<3> & mss) { return mss.Masses(); })
      .def_property_readonly("fixes", [](MassSpringSystem<3> & mss) { return mss.Fixes(); })      

      .def("GetState", [] (MassSpringSystem<3> & mss) {
        Vector<> x(3*mss.Masses().size());
        Vector<> dx(3*mss.Masses().size());
        mss.GetState (x, dx);
        return x;
      })
      ;
    

    m.def("Simulate", [](MassSpringSystem<3> & mss, double tend, size_t steps) {
      Vector<> x(3*mss.Masses().size());
      Vector<> dx(3*mss.Masses().size());
      mss.GetState (x, dx);
      
      auto mss_func = make_shared<MSS_Function<3>> (mss);
      
      SolveODE_Verlet(tend, tend/steps, x, dx, mss_func);
        /*
                      [](double t, VectorView<double> x) { cout << "t = " << t
                                                                << ", x = " << x << endl; });
        */
      mss.SetState (x, dx);
    });
      
    
}
