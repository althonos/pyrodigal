python_sources = [
    '__init__.py',
    '__main__.py',
    '_version.py',
    'cli.py',
    'lib.pxd',
    'lib.pyi',
]

c_headers = [
    '_connection.h',
    '_sequence.h',
    '_translation.h',
]

cython_headers = [
    '_connection.pxd',
]

sources = python_sources + c_headers + cython_headers
python.install_sources(sources, subdir : 'pyrodigal')

subdir('impl')
subdir('tests')

sys_version_info_major = run_command([python, '-c', 'import sys; print(sys.version_info[0])' ], check : true).stdout().strip()
sys_version_info_minor = run_command([python, '-c', 'import sys; print(sys.version_info[1])' ], check : true).stdout().strip()
sys_implementation_name = run_command([python, '-c', 'import sys; print(sys.implementation.name)' ], check : true).stdout().strip()

cython_args = [
    '--module-name=pyrodigal.lib',
    '-EMMX_BUILD_SUPPORT=' + mmx_supported.to_string('True', 'False'),
    '-ESSE2_BUILD_SUPPORT=' + sse2_supported.to_string('True', 'False'),
    '-EAVX2_BUILD_SUPPORT=' + avx2_supported.to_string('True', 'False'),
    '-EAVX512_BUILD_SUPPORT=' + avx512_supported.to_string('True', 'False'),
    '-ENEON_BUILD_SUPPORT=' + neon_supported.to_string('True', 'False'),
    '-ESYS_VERSION_INFO_MAJOR=' + sys_version_info_major,
    '-ESYS_VERSION_INFO_MINOR=' + sys_version_info_minor,
    '-ESYS_IMPLEMENTATION_NAME=' + sys_implementation_name,
    '-ETARGET_CPU=' + target_machine.cpu_family(),
    '-Xcdivision=True',
    '-Xnonecheck=False',
]

if get_option('buildtype') == 'debug'
    cython_args += [
        '-Xcdivision_warnings=True',
        '-Xwarn.undeclared=True',
        '-Xwarn.unreachable=True',
        '-Xwarn.maybe_uninitialized=True',
        '-Xwarn.unused=True',
        '-Xwarn.unused_arg=True',
        '-Xwarn.unused_result=True',
        '-Xwarn.multiple_declarators=True',
    ]
else
    cython_args += [
        '-Xboundscheck=False',
        '-Xwraparound=False',
    ]
endif

if has_getid
    c_args = ['-DHAS_PYINTERPRETERSTATE_GETID']
else
    c_args = []
endif

python.extension_module('lib',
    ['lib.pyx'],
    include_directories: prodigal_inc,
    dependencies : py_dep,
    link_with : [prodigal, impl],
    install: true,
    subdir : 'pyrodigal',
    cython_args : cython_args,
    c_args : c_args,
    extra_files : sources,
)


