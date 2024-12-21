project(
  'champ',
  'cpp',
  default_options: [
    'cpp_std=c++17',
  ],
)

include_dir = include_directories('./include')

#py = import('python').find_installation('python3', required: true)
#py_dep = py.dependency(required: true)
py_dep = dependency('python3-embed')
assert(meson.get_compiler('cpp').has_header('Python.h', dependencies: py_dep))

pybind11_dep = subproject('pybind11').get_variable('pybind11_dep')
assert(meson.get_compiler('cpp').has_header('pybind11/pybind11.h', dependencies: pybind11_dep))

source_files = files(
  'atom.cpp',
  'champ.cpp',
)

#champ_lib = static_library(
#  'champ',
#  source_files,
#  include_directories: include_dir,
#  dependencies: [py_dep, pybind11_dep],
#  link_args : ['-Wl,-V'],
#)

subdir('tests')
