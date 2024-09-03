add_custom_command(
  OUTPUT training.c
  COMMENT
    "Rewriting training.c into sparse initializer"
  COMMAND
    Python::Interpreter ${CMAKE_CURRENT_LIST_DIR}/../scripts/split_training.py 
      -i ${CMAKE_CURRENT_LIST_DIR}/training.c 
      -o ${CMAKE_CURRENT_BINARY_DIR}/training.c
  MAIN_DEPENDENCY
    ${CMAKE_CURRENT_LIST_DIR}/training.c 
  VERBATIM
)

add_library(prodigal
  bitmap.c
  dprog.c
  gene.c
  metagenomic.c
  node.c
  sequence.c
  ${CMAKE_CURRENT_BINARY_DIR}/training.c
)