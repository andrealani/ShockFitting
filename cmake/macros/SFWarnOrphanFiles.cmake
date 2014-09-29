##############################################################################
# macro for warning about missing files in the listing of files
# usefull for headers  that are forgotten and give errors compiling
# from installed sources.
# it warns about files that are present in the directory but are neither
# in the headers or the sources of all the libraries and applications
# defined in this directory
##############################################################################
MACRO( SF_WARN_ORPHAN_FILES )

  # first find all the files in the directory
  FOREACH( CFEXT ${SF_FILE_EXTENSIONS} )
    FILE ( GLOB THISEXT_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.${CFEXT})
    LIST ( LENGTH THISEXT_FILES NELEMS )
#     SF_DEBUG_VAR(NELEMS)
    IF ( NELEMS GREATER 0 )
      LIST ( APPEND alldirfiles "${THISEXT_FILES}")
    ENDIF ( NELEMS GREATER 0 )
  ENDFOREACH(CFEXT)

  # remove all found files from orphan list to avoid duplicates
  FOREACH( AFILE ${alldirfiles} )
    SF_CACHE_LIST_REMOVE ( SF_ORPHAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE} )
  ENDFOREACH(AFILE)

  # remove files present in optional list
  FOREACH( AFILE ${OPTIONAL_dirfiles} )
    LIST ( REMOVE_ITEM alldirfiles ${AFILE})
  ENDFOREACH(AFILE)

  # remove files present in libs
  FOREACH( LOCALLIB ${SF_LOCAL_LIBNAMES} )
    FOREACH( AFILE ${${LOCALLIB}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALLIB ${SF_LOCAL_LIBNAMES} )

  # remove files present in apps
  FOREACH( LOCALAPP ${SF_LOCAL_APPNAMES} )
    FOREACH( AFILE ${${LOCALAPP}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALAPP ${SF_LOCAL_APPNAMES} )

  # remove files present in tests
  FOREACH( LOCALTEST ${SF_LOCAL_TESTNAMES} )
    FOREACH( AFILE ${${LOCALTEST}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALTEST ${SF_LOCAL_TESTNAMES} )

  # warn about the other files
  FOREACH( AFILE ${alldirfiles} )
    # temporarily ignore Test files while the unit tests are not added
    IF(${AFILE} MATCHES "^Test" )
    ELSE(${AFILE} MATCHES "^Test" )
      LOGFILE(" +++ WARNING : orphan file ${AFILE}")
      SF_CACHE_LIST_APPEND ( SF_ORPHAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE} )
    ENDIF(${AFILE} MATCHES "^Test" )
  ENDFOREACH(AFILE)

ENDMACRO( SF_WARN_ORPHAN_FILES )
##############################################################################
