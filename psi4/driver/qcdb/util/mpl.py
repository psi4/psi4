
def plot_coord(ref, cand=None, orig=None, comment=None):
    """Display target geometry `ref` as black dots in 3D plot. If present, also
    plot candidate geometry `cand` as red dots and starting geometry `orig` as
    pale blue dots. Plot has text `comment`. For assessing alignment, red and
    black should overlap and pale blue shows where red started.
    
    """
    import numpy as np
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D

    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    bound = max(np.amax(ref), -1 * np.amin(ref))
    ax.scatter(ref[:, 0], ref[:, 1], ref[:, 2], c='k', label='goal')
    if cand is not None:
        ax.scatter(cand[:, 0], cand[:, 1], cand[:, 2], c='r', label='post-align')
    if orig is not None:
        ax.scatter(orig[:, 0], orig[:, 1], orig[:, 2], c='lightsteelblue', label='pre-align')

    if comment is not None:
        ax.text2D(0.05, 0.95, comment, transform=ax.transAxes)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(-bound, bound)
    ax.set_ylim(-bound, bound)
    ax.set_zlim(-bound, bound)
    ax.legend()
    
    pyplot.show()
