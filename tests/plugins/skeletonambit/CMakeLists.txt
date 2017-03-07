# check out and build fresh plugin
# creates a dir skeleton w/i skeleton then copies files up to satisfy requirements of both CMake and Psi4 plugins
# symlinks input.dat back to source dir to satisfy TestingMacro

add_skeleton_plugin("skeletonambit" "ambit" "psi;plug")
