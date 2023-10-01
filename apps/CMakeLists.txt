foreach(_target IN ITEMS
    fenics-dg
    lorenzo
    myPF
    arc-length
    Bayat-pf
    Bayat2
    Nguyen
    CZMs
    CZMs2
    deal
    Adaptivity
    )
  add_subdirectory(${_target})
  install(TARGETS ${_target})
endforeach()
