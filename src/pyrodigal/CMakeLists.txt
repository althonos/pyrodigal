add_subdirectory(impl)

set(LIB_LINKS impl)
if(HAVE_AVX2)
    set(LIB_LINKS ${LIB_LINKS} implavx)
endif()
if(HAVE_SSE2)
    set(LIB_LINKS ${LIB_LINKS} implsse)
endif()
if(HAVE_NEON)
    set(LIB_LINKS ${LIB_LINKS} implneon)
endif()

cython_extension(lib LINKS prodigal ${LIB_LINKS})

