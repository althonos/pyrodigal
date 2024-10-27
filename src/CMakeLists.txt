file(COPY ../vendor/Prodigal DESTINATION .)
file(COPY ./Prodigal DESTINATION .)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/Prodigal ${CMAKE_CURRENT_BINARY_DIR}/Prodigal/build)
add_subdirectory("pyrodigal")
