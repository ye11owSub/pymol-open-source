catch_dep = subproject('catch2').get_variable('catch2_with_main_dep')

test_files = source_files + files(
  'atom_test.cpp'
)

champ_test = executable(
  'champ_test',
  test_files,
  include_directories: include_dir,
  dependencies: [py_dep, pybind11_dep, catch_dep],
  #  link_with: champ_lib,
)

test('champ_test', champ_test)
