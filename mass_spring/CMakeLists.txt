add_executable (mass_spring mass_spring.cc)

# include_directories(../src)

# target_sources (mass_spring PUBLIC ../src/ode.h ../src/nonlinfunc.h)
# target_sources (mass_spring PUBLIC ode.h)
target_link_libraries (mass_spring PUBLIC ngbla)

