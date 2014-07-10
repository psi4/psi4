np = None

def set_np(module):
    global np
    np = module

# generate nPt random points within the grid
def generate_random_points(shape, nPt):
    pt = None
    if False:
        pt = []
        for k in range(nPt):
            x = shape[0] * np.random.random_sample()
            y = shape[1] * np.random.random_sample()
            pt.append((x, y))
    else:
        pt = np.random.random_sample((nPt,2))
        pt[:,0] *= shape[0]
        pt[:,1] *= shape[1]
    return pt

def dist((x, y), (x0, y0), expo):
    """Returns the Euclidean distance between two (2D) points."""
    return ((x-x0)**expo + (y-y0)**expo)**(1.0/expo)

def dist_array(shape, pt, expo):
    """Returns an array whose elements are the distance (in units of
    rows and columns) to a given point (x0, y0)."""
    x0 = pt[0]
    y0 = pt[1]
    f = lambda x, y: dist((x, y), (x0, y0), expo)
    return np.fromfunction(f, shape, dtype=np.float32)

def distmap(shape, nPt, expo):
    pt = generate_random_points(shape, nPt)
    # distMin[i,j] is the distance from pixel (i,j) to the closest 'pt'.
    distMin = dist_array(shape, pt[0], expo)
    for i in range(1, nPt):
        distMin = np.minimum(dist_array(shape, pt[i], expo), distMin)
    pxlMin = 0
    pxlMax = np.maximum.reduce(distMin)
    # scale and offset distances to lie between 0 and 1
    distMinScaled = (distMin - pxlMin) / (pxlMax - pxlMin)
    # (gray) pixel values lie between 0 and 255
    pxls = 255 * distMinScaled
    return pxls
    
def main():
    # python libs
    from getopt import getopt
    import sys

    # 3rd party libs
    #from PIL import Image
    import Image

    # command-line options
    nPt = 8 # default
    expo = 2 # default distance exponent (Euclidean)
    nGrid = 512 # image size (in each dimension)
    seed = 0 # default seed
    useGain = False
    help = False
    imgOutput = True
    rawOutput = False
    img_name = "distmap.png"
    raw_name = "distmap.dat"
    usage = """
    distmap.py [options]

      -h        help (prints this message)
      -g        use gain instead of numpy
      -e expo   set exponent (default: 2 (Euclidean))
      -o out    output filename (default: distmap.png)
      -r raw    raw output filename (default: distmap.dat)
      -p points number of points (default: 8)
      -s seed   for random number generation
      -w width  image size e.g. Width X Width square
      -d        whether to enable profiling

    """

    (optsvals, args) = getopt(sys.argv[1:], 'hge:o:r:p:s:w:')
    for (opt, val) in optsvals:
        if opt == '-h':
            help = True
        elif opt == '-g':
            useGain = True
        elif opt == '-e':
            expo = float(val)
        elif opt == '-o':
            img_name = str(val)
            imgOutput = True
        elif opt == '-r':
            raw_name = str(val)
            rawOutput = True
        elif opt == '-p':
            nPt = int(val)
        elif opt == '-s':
            seed = int(val)
        elif opt == '-w':
            nGrid = int(val)

    if useGain:
        import ga.gain as np
        from ga.gain import me
        set_np(np)
        if help:
            if not me():
                print usage
            sys.exit()
    else:
        import numpy as np
        set_np(np)
        if help:
            print usage
            sys.exit()

    shape = (nGrid,nGrid)
    # seed the random number generator
    np.random.seed(seed)
    # generate the distmap
    pxls = distmap(shape, nPt, expo)
    # write pxls to file
    if (imgOutput or rawOutput) and useGain:
        if me() == 0:
            pxls = pxls.get()
        else:
            return
    if imgOutput:
        im_r = Image.frombuffer("F", pxls.shape, pxls, "raw", "F", 0, 1)
        im_r = im_r.convert("L")
        im_r.save(img_name)
    if rawOutput:
        f = open(raw_name, "wb")
        f.write(buffer(pxls.astype(np.float32)))
        f.close()

if __name__ == '__main__':
    profile = False
    if profile:
        import cProfile
        if useGain:
            cProfile.run("main()", "distmap.prof" + str(me()))
        else:
            cProfile.run("main()", "distmap.prof")
    else:
        main()
