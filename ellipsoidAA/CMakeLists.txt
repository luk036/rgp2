cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Debug )
add_definitions ( -Wall )
set ( LIBRARY_OUTPUT_PATH ../lib )
include_directories( ../ ../aaf )
set ( ellipsoidAA_HDRS
  ellipsoid.h
  gp_base.h
  gp_prob1.h
  posynomial.h
  monomial.h
  rgp_base.h
  )

set ( ellipsoidAA_SRCS
  gp_prob1.cpp
  rgp_base.cpp
  rgp_prob1.cpp
  )

add_library ( ellipsoidAA ${ellipsoidAA_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )

link_directories( ../lib )

add_executable ( cvx_test1 cvx_test1.cpp ${ellipsoidAA_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
add_executable ( cvx_test2 cvx_test2.cpp ${ellipsoidAA_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
add_executable ( gp_test gp_test.cpp ${ellipsoidAA_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
add_executable ( rgp_test rgp_test.cpp ${ellipsoidAA_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( cvx_test1 aaf )
target_link_libraries ( cvx_test2 aaf )
target_link_libraries ( gp_test ellipsoidAA aaf )
target_link_libraries ( rgp_test ellipsoidAA aaf )
