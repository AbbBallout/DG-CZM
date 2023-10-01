cmake_minimum_required(VERSION 3.8)


# Assume that deal.II_DIR won't be provided on the command line in favor
# of DEAL_II_DIR and can thus be used to detect whether CMake was
# already run.
if(DEFINED deal.II_DIR)
  set(_had_deal_ii TRUE)
endif()
find_package(deal.II 9.4.0 QUIET REQUIRED
  HINTS "${deal.II_DIR}" "${DEAL_II_DIR}" "$ENV{DEAL_II_DIR}")
deal_ii_initialize_cached_variables()

project(CZM_DG)

add_subdirectory(apps)


