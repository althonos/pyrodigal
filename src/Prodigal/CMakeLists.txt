set(PRODIGAL_SOURCES
  bitmap.c
  dprog.c
  gene.c
  metagenomic.c
  node.c
  sequence.c
)

foreach(i RANGE 0 49)
  add_custom_command(
    OUTPUT training${i}.c
    COMMENT
      "Extracting initialize_metagenome_${i} from training.c"
    COMMAND
      ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/split_training.py 
        -i ${CMAKE_CURRENT_LIST_DIR}/training.c 
        -o ${CMAKE_CURRENT_BINARY_DIR}/training${i}.c
        --index ${i}
    MAIN_DEPENDENCY
      ${CMAKE_CURRENT_LIST_DIR}/training.c 
    VERBATIM
  )
  set(PRODIGAL_SOURCES ${PRODIGAL_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/training${i}.c)
endforeach()

add_library(prodigal ${PRODIGAL_SOURCES})
target_include_directories(prodigal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
