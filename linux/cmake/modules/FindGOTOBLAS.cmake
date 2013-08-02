# - Try to find GOTOBLAS
#
#  User should define GOTOBLAS_LIBRARY_PATH and GFORTRAN_LIBRARY_PATH
#
#  Once done this will define 
#  GOTOBLAS_FOUND - System has LAMA_BLAS
#  GOTOBLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

## Search GOTOBLAS library

find_library ( GOTOBLAS_LIBRARY goto2 HINTS ${GOTOBLAS_LIBRARY_PATH} )
if ( EXISTS ${GOTOBLAS_LIBRARY} )
    set ( GOTOBLAS_LIBRARIES ${GOTOBLAS_LIBRARY} )
    
    ## Search for gfortran. Required by GOTOBLAS

    find_library( GFORTRAN_LIBRARY gfortran HINTS ${GFORTRAN_LIBRARY_PATH} )
    
    if( EXISTS ${GFORTRAN_LIBRARY} )
        list ( APPEND GOTOBLAS_LIBRARIES ${GFORTRAN_LIBRARY} )
    else( EXISTS ${GFORTRAN_LIBRARY} )
        message ( STATUS "WARNING Library gfortran not found. Required by GOTO BLAS. Please define GFORTRAN_LIBRARY_PATH." )
    endif( EXISTS ${GFORTRAN_LIBRARY} )
    
    ## Search for pthreads. Required by GOTOBLAS
    
    set ( CMAKE_THREAD_PREFER_PTHREAD TRUE )
    
    find_package( Threads QUIET )
    
    if ( Threads_FOUND )
        if ( CMAKE_USE_PTHREADS_INIT )
            # TODO is this correct???
            list ( APPEND GOTOBLAS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} )
            # add_definitions(${CMAKE_THREAD_LIBS_INIT}) 
        else ( CMAKE_USE_PTHREADS_INIT )
            message ( STATUS "WARNING No pthreads found. Required by GOTOBLAS." )    
        endif ( CMAKE_USE_PTHREADS_INIT )
    else ( Threads_FOUND )
        message ( STATUS "WARNING No threads found.Required by GOTOBLAS." )
    endif ( Threads_FOUND )
    
elseif ( NOT DEFINED GOTOBLAS_LIBRARY_PATH )
    message ( STATUS "WARNING GOTOBLAS not found. Please define GOTOBLAS_LIBRARY_PATH." )
else ( EXISTS ${GOTOBLAS_LIBRARY} )
    message ( STATUS "WARNING GOTOBLAS not found. GOTOBLAS_LIBRARY_PATH=${GOTOBLAS_LIBRARY_PATH} directory does not exist." )
endif ( EXISTS ${GOTOBLAS_LIBRARY} )

## Module footer

include( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set GOTOBLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( GOTOBLAS DEFAULT_MSG GOTOBLAS_LIBRARIES )

if( GOTOBLAS_FOUND )
   add_definitions( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
endif()

mark_as_advanced( GOTOBLAS_LIBRARIES GOTOBLAS_LIBRARY )