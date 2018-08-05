add_executable( solution-example-moving main.cc )
target_compile_definitions( solution-example-moving PRIVATE ALUGRID_SIMPLEX;GRIDDIM=2;WORLDDIM=3;WANT_ISTL=1;POLORDER=1)
add_dune_mpi_flags( solution-example-moving )
target_link_libraries( solution-example-moving  )
