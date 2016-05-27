set(M4DIR ${PROJECT_BINARY_DIR}/m4out )
file(MAKE_DIRECTORY ${M4DIR})

function(m4generate SOURCE)
  get_filename_component(FULLNAME ${SOURCE} NAME)
  #get_filename_component(EXT ${SOURCE} extension)
  #assert(${extension} STREQ ".m4")
  string(REPLACE ".m4" "" NAME ${FULLNAME})
  set(OUTPUT ${M4DIR}/${NAME})
  message("m4 output: " ${OUTPUT})
    add_custom_command( 
        OUTPUT ${OUTPUT}
        COMMAND m4  ${SOURCE} > ${OUTPUT}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${SOURCE}
	VERBATIM
        )
set_source_files_properties( ${OUTPUT} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE )
endfunction()
