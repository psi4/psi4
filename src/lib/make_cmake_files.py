from glob import glob
import os
exit() # Not working right now....
libs = []
for dirname, dirnames, filenames in os.walk('.'):
    lib =  os.path.basename(dirname)
    print lib
    if not lib.startswith("lib"):
        continue
    # These are handled differently
    if lib in ["libint", "libderiv"]:
        continue
    libs.append(lib)
    # Remove the lib part now
    name = lib[3:]
    outfile = open("%s/CMakeLists.txt" % lib, "w")
    files = []
    for f in filenames:
        if f.endswith(".cc"):
            files.append(f)
    #outfile.write("${SRC} = %s\n" % (" ".join(files)))
    #outfile.write("add_library(%s ${SRC})\n" % name)
    outfile.close()
outfile = open("CMakeLists.txt", w)
for l in libs:
    continue
    outfile.write("add_subdirectory(%s)\n" % l)
outfile.close()

