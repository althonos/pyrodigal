impl = []

impl += static_library(
    'impl_generic',
    'generic.c',
    include_directories : prodigal_inc,
    extra_files : ['generic.h', 'generic.pxd'],
)

if mmx_supported
    impl += static_library(
        'impl_mmx', 
        'mmx.c', 
        include_directories : prodigal_inc, 
        c_args : mmx_flags,
        extra_files : ['mmx.h', 'mmx.pxd'],
    )
endif

if sse2_supported
    impl += static_library(
        'impl_sse2', 
        'sse.c', 
        include_directories : prodigal_inc, 
        c_args : sse2_flags,
        extra_files : ['sse.h', 'sse.pxd'],
    )
endif

if avx2_supported
    impl += static_library(
        'impl_avx', 
        'avx.c', 
        include_directories : prodigal_inc, 
        c_args : avx2_flags,
        extra_files : ['avx.h', 'avx.pxd'],
    )
endif

if avx512_supported
    impl += static_library(
        'impl_avx512', 
        'avx512.c', 
        include_directories : prodigal_inc, 
        c_args : avx512_flags,
        extra_files : ['avx512.h', 'avx512.pxd'],
    )
endif

if neon_supported
    impl += static_library(
        'impl_neon', 
        'neon.c', 
        include_directories : prodigal_inc, 
        c_args : neon_flags,
        extra_files : ['neon.h', 'neon.pxd'],
    )
endif

