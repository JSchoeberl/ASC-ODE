#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mass_spring.h"

namespace py = pybind11;


PYBIND11_MODULE(mass_spring, m) {
    m.doc() = "mass-spring-system simulator"; 

    py::class_<Mass<2>> (m, "Mass2d");
      ;
      
    m.def("Mass", [](double m, std::array<double,2> p)
    {
      return Mass<2>{m, { p[0], p[1] }};
    });

    py::class_<Mass<3>> (m, "Mass3d");
      ;
      
    m.def("Mass", [](double m, std::array<double,3> p)
    {
      return Mass<3>{m, { p[0], p[1], p[2] }};
    });

    
    py::class_<MassSpringSystem<2>> (m, "MassSpringSystem2d")
      .def(py::init<>())
      .def("Add", [](MassSpringSystem<2> & mss, Mass<2> m) { mss.AddMass(m); })
      ;
      
    
    /*
    py::class_<Vector<double>> (m, "Vector")
      .def(py::init<size_t>(),
           py::arg("size"), "create vector of given size")
      .def("__len__", &Vector<double>::Size,
           "return size of vector")
      
      .def("__setitem__", [](Vector<double> & self, size_t i, double v) { self(i) = v; })
      .def("__getitem__", [](Vector<double> & self, size_t i) { return self(i); })
      
      .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
      {
        size_t start, stop, step, n;
        if (!inds.compute(self.Size(), &start, &stop, &step, &n))
          throw py::error_already_set();
        self.Range(start, stop).Slice(0,step) = val;
      })
      
      .def("__add__", [](Vector<double> & self, Vector<double> & other)
      { return Vector<double> (self+other); })

      .def("__rmul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })
      
      .def("__str__", [](const Vector<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })
    ;
    */
}
