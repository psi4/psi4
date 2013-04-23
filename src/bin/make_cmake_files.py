from glob import glob
import os
exclude = ["cceom/read_guess.cc", "detci/calc_pt_block.cc", "detci/graphset.cc"]
libs = []
for d in [name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name))]:
    print "Processing", d
    if d in ["attic", "psi4"]:
        continue
    libs.append(d)
    files = glob("%s/*.cc" % d)
    # Remove the lib part now
    filenames = []
    for f in files:
        if f in exclude:
            continue
        filenames.append(os.path.basename(f))

    # These are handled differently
    outfile = open("%s/CMakeLists.txt" % d, "w")
    outfile.write("set(SRC %s)\n" % (" ".join(filenames)))
    outfile.write("add_library(%s ${SRC})\n" % d)
    outfile.close()
outfile = open("CMakeLists.txt", "w")
for l in libs:
    outfile.write("add_subdirectory(%s)\n" % l)
outfile.write("add_subdirectory(psi4)\n")
outfile.close()

