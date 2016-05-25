#Author: Giovanni Balduzzi

function(create_ispc_library SOURCES)
  # create a library named ispc containing the .ispc source files in $SOURCES 
  # and the generated .ispc.h headers
  # Note: if sources is a list it needs to be passed inside quotation marks
  list(LENGTH SOURCES length)
  message("Adding " ${length} " files to ispc library")
  
  set(ISPC_COMPILER ispc CACHE STRING "ispc compiler")
  set(ISPCFLAGS -O2 --cpu=corei7-avx CACHE STRING "ispc compiler flags")
  set(OBJ_LOCATION ${PROJECT_BINARY_DIR}/ispc/obj)
  set(INC_LOCATION ${PROJECT_BINARY_DIR}/ispc/include)
  file(MAKE_DIRECTORY ${OBJ_LOCATION})
  file(MAKE_DIRECTORY ${INC_LOCATION})
  
  foreach(SOURCE  ${SOURCES})
    get_filename_component(NAME ${SOURCE} NAME)
    set(OBJECT ${OBJ_LOCATION}/${NAME}.o)
    set(HEADER ${INC_LOCATION}/${NAME}.h)
    add_custom_command( 
        OUTPUT ${OBJECT}
        COMMAND ${ISPC_COMPILER} ${ISPCFLAGS} -o ${OBJECT} ${SOURCE} -h ${HEADER}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
set_source_files_properties( ${OBJECT} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE )
set_source_files_properties( ${HEADER} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE )
list(APPEND OBJECTS ${OBJECT})
  endforeach(SOURCE)

  add_library(ispc STATIC  ${OBJECTS})
  target_include_directories(ispc PUBLIC ${INC_LOCATION})
  set_target_properties(ispc PROPERTIES
    LINKER_LANGUAGE CXX
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/ispc
    )
endfunction()
