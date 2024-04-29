set(LIBRARY_NAME
        variational)

set(SOURCES
	   src/SimpleAdam.cxx)

add_library(${LIBRARY_NAME} ${SOURCES})

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(${LIBRARY_NAME} PUBLIC eigen)
