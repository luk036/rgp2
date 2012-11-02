cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Debug )
add_definitions ( -Wall )
set ( LIBRARY_OUTPUT_PATH ../lib )
set ( aaf_HDRS
	aa.h
	aa_aaf.h
	aa_interval.h
	aa_rounding.h
	aa_util.h
	)

set ( aaf_SRCS
	aa_aafapprox.cpp
	aa_aafarithm.cpp
	aa_aafcommon.cpp
	aa_aaftrigo.cpp
	aa_interval.cpp
	aa_rounding.cpp
	)

add_library ( aaf ${aaf_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( aaf  )
