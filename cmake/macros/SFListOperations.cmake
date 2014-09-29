##############################################################################
# this macro adds to a cached list if element not yet present
##############################################################################
MACRO( SF_CACHE_LIST_REMOVE THELIST THEVAR )
  IF ( DEFINED ${THELIST} )
    LIST( REMOVE_ITEM ${THELIST} ${THEVAR} )
  ENDIF ( DEFINED ${THELIST} )
ENDMACRO( SF_CACHE_LIST_REMOVE THELIST THEVAR )
##############################################################################

##############################################################################
# this macro adds to a cached list if element not yet present
##############################################################################
MACRO( SF_CACHE_LIST_APPEND THELIST THEVAR )
  SF_CACHE_LIST_REMOVE(${THELIST} ${THEVAR})
  SET ( ${THELIST} ${${THELIST}} ${THEVAR} CACHE INTERNAL "" FORCE )
ENDMACRO( SF_CACHE_LIST_APPEND THELIST THEVAR )
##############################################################################