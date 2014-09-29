##############################################################################
# kernel library macro
##############################################################################
MACRO( SF_ADD_KERNEL_LIBRARY LIBNAME )

  # declare this library as part of the kernel
  SET ( ${LIBNAME}_kernellib ON )
  
  SF_ADD_LIBRARY( ${LIBNAME} )

ENDMACRO( SF_ADD_KERNEL_LIBRARY )
##############################################################################
