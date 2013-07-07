from glob import glob
import os
excludes = ["libchkpt/exist_add_prefix.cc", "libchkpt/vib_freqs.cc"]
libs = []
for d in [name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name))]:
    if not d.startswith("lib"):
        continue
    if d in ["libint", "libderiv"]:
        continue
    print "Processing", d
    libs.append(d)
    files = glob("%s/*.cc" % d)
    # Remove the lib part now
    libname = d[3:]
    filenames = []
    for f in files:
        if f in excludes:
            continue
        filenames.append(os.path.basename(f))

    # These are handled differently
    outfile = open("%s/CMakeLists.txt" % d, "w")
    outfile.write("set(SRC %s)\n" % (" ".join(filenames)))
    outfile.write("add_library(%s ${SRC})\n" % libname)
    outfile.close()
outfile = open("CMakeLists.txt", "w")
for l in libs:
    outfile.write("add_subdirectory(%s)\n" % l)
outfile.write("add_subdirectory(libint)\n")
outfile.write("add_subdirectory(libderiv)\n")
outfile.close()

