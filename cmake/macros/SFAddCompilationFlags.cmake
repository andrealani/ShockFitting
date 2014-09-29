##############################################################################
# this macro separates the sources form the headers
##############################################################################
MACRO( SF_ADD_C_FLAGS m_c_flags )

  IF ( NOT DEFINED N_CFLAG )
    SET( N_CFLAG 0 )
  ENDIF()

 	MATH ( EXPR N_CFLAG '${N_CFLAG}+1'  )    

  CHECK_C_COMPILER_FLAG ( ${m_c_flags} C_FLAG_TEST_${N_CFLAG} )

#  MESSAGE ( STATUS "FLAG [${m_c_flags}] is [${C_FLAG_TEST_${N_CFLAG}}] " )
  IF ( C_FLAG_TEST_${N_CFLAG} )
	  SET ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${m_c_flags}" )
    LOG ( "C FLAG [${m_c_flags}] added" )
  ELSE  ()
    LOG ( "C FLAG [${m_c_flags}] skipped" )
  ENDIF()

ENDMACRO( SF_ADD_C_FLAGS )
##############################################################################

##############################################################################
MACRO( SF_ADD_CXX_FLAGS m_cxx_flags )

  IF ( NOT DEFINED N_CXXFLAG )
    SET( N_CXXFLAG 0 )
  ENDIF()

 	MATH ( EXPR N_CXXFLAG '${N_CXXFLAG}+1'  )    

  CHECK_CXX_COMPILER_FLAG ( ${m_cxx_flags} CXX_FLAG_TEST_${N_CXXFLAG} )

#  MESSAGE ( STATUS "FLAG [${m_cxx_flags}] is [${CXX_FLAG_TEST_${N_CXXFLAG}}] " )
  IF ( CXX_FLAG_TEST_${N_CXXFLAG} )
	  SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${m_cxx_flags}" )
    LOG ( "C++ FLAG [${m_cxx_flags}] added" )
  ELSE  ()
    LOG ( "C++ FLAG [${m_cxx_flags}] skipped" )
  ENDIF()

ENDMACRO( SF_ADD_CXX_FLAGS )
##############################################################################

##############################################################################
MACRO( SF_ADD_C_FLAGS_SIGNAL_ERROR m_c_flags )

  SF_ADD_C_FLAGS ( ${m_c_flags}  )
  IF ( NOT C_FLAG_TEST_${N_CFLAG} )
    MESSAGE(FATAL_ERROR "C compiler [${CMAKE_C_COMPILER}] cannot used requested C flag [${m_c_flags}]")
  ENDIF ()

ENDMACRO( SF_ADD_C_FLAGS_SIGNAL_ERROR )
##############################################################################

##############################################################################
MACRO( SF_ADD_CXX_FLAGS_SIGNAL_ERROR m_cxx_flags )

  SF_ADD_CXX_FLAGS ( ${m_cxx_flags}  )
  IF ( NOT CXX_FLAG_TEST_${N_CXXFLAG} )
    MESSAGE(FATAL_ERROR "C++ compiler [${CMAKE_CXX_COMPILER}] cannot used requested C++ flag [${m_cxx_flags}]")
  ENDIF ()

ENDMACRO( SF_ADD_CXX_FLAGS_SIGNAL_ERROR )
##############################################################################

