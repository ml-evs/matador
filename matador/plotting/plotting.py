# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some useful auxiliary routines for use in
the other plotting functions, and some generic routines for plotting
simple (x, y) line data (e.g. pair distribution functions), or lists of
(x, y) data against some third parameter (e.g. powder x-ray spectrum vs
voltage).

"""

SAVE_EXTS = ['pdf', 'png', 'svg']
MATADOR_STYLE = '/'.join(__file__.split('/')[:-1]) + '/../config/matador.mplstyle'


def set_style(style=None):
    """ Set the matplotlib style for all future plots, manually. This
    will conflict with the context manager used by the `plotting_function`
    wrapper.

    """
    import matplotlib.pyplot as plt
    if style is None or style == 'matador':
        style = MATADOR_STYLE
    if not isinstance(style, list):
        style = [style]
    # apply multiple compound styles, if present
    for styles in style:
        plt.style.use(styles)


def plotting_function(function):
    """ Wrapper for plotting functions to safely fail on X-forwarding
    errors and handle the plot style context manager.
    """

    from functools import wraps
    from matador.utils.print_utils import print_warning, print_failure
    from matador.config import load_custom_settings

    @wraps(function)
    def wrapped_plot_function(*args, **kwargs):
        """ Wrap and return the plotting function. """
        saving = False
        result = None

        # if we're going to be saving a figure, switch to Agg to avoid X-forwarding
        try:
            for arg in args:
                if arg.savefig:
                    import matplotlib
                    # don't warn as backend might have been set externally by e.g. Jupyter
                    matplotlib.use('Agg', force=False)
                    saving = True
                    break
        except AttributeError:
            pass
        if not saving:
            if any(kwargs.get(ext) for ext in SAVE_EXTS):
                import matplotlib
                matplotlib.use('Agg', force=False)
                saving = True

        settings = load_custom_settings(kwargs.get('config_fname'), quiet=True, no_quickstart=True)
        try:
            style = settings.get('plotting', {}).get('default_style')
            if kwargs.get('style'):
                style = kwargs['style']
            if style is not None and not isinstance(style, list):
                style = [style]
            if style is None:
                style = ['matador']
            if 'matador' in style:
                for ind, styles in enumerate(style):
                    if styles == 'matador':
                        style[ind] = MATADOR_STYLE

            # now actually call the function
            set_style(style)
            result = function(*args, **kwargs)

        except Exception as exc:
            if 'TclError' not in type(exc).__name__:
                raise exc

            print_failure('Caught exception: {}'.format(type(exc).__name__))
            print_warning('Error message was: {}'.format(exc))
            print_warning('This is probably an X-forwarding error')
            print_failure('Skipping plot...')

        return result

    return wrapped_plot_function


def get_linear_cmap(colours, num_colours=100, list_only=False):
    """ Create a linear colormap from a list of colours.

    Parameters:
        colours (:obj:`list` of :obj:`str`): list of fractional RGB/hex
            values of colours

    Keyword arguments:
        num_colours (int): number of colours in resulting cmap
        list_only (bool): return only a list of colours

    Returns:
        :obj:`matplotlib.colors.LinearSegmentedColormap` or :obj:`list`:
            returns list of colours if `list_only` is True, otherwise
            :obj:`matplotlib.colors.LinearSegmentedColormap`.

    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    colours = [to_rgb(colour) for colour in colours]

    uniq_colours = []
    _colours = [tuple(colour) for colour in colours]
    for colour in _colours:
        if colour not in uniq_colours:
            uniq_colours.append(colour)
    _colours = uniq_colours
    linear_cmap = []
    repeat = int(num_colours / len(_colours))
    for ind, colour in enumerate(_colours):
        if ind == len(_colours) - 1:
            break
        diff = np.asarray(_colours[ind + 1]) - np.asarray(_colours[ind])
        diff_norm = diff / repeat
        for i in range(repeat):
            linear_cmap.append(np.asarray(colour) + i * diff_norm)
    if list_only:
        return linear_cmap

    return LinearSegmentedColormap.from_list('linear_cmap', linear_cmap, N=num_colours)


class XYvsZPlot:
    """ This class wraps plotting (x, y) lines against a third
    variable.

    """
    def __init__(self, xys, zs, y_scale=1.0, **kwargs):
        """ Construct plot from data.

        Parameters:
            xys (:obj:`list` of :obj:`list` or numpy.ndarray): list or
                array of data to be plotted. For N lines of M samples,
                this can be provided as an (N, M, 2) or (2, M, N) array,
                or corresponding list/sublist format.
            zs (:obj:`list`): third parameter to plot lines against. The
                y-values are rescaled relative to the maximum across all
                lines so that no lines overlap (this can be overridden
                using offset_factor keyword).

        Keyword arguments:
            y_scale (float): controls the scale factor between the
                arbitrary y-scale and the z-scale.

        """

        import numpy as np

        self.plot_kwargs = kwargs

        _xys = np.asarray(xys)
        shape = np.shape(_xys)
        if shape[0] != 2 and shape[-1] != 2:
            raise RuntimeError('Data of shape {} is not compatible with XYvsZPlot.'
                               .format(shape))
        if shape[0] == 2:
            _xys = _xys.T

        self._xs = _xys[:, :, 0]
        self._ys = _xys[:, :, 1]
        self._zs = np.asarray(zs).flatten()

        if len(self._zs) != np.shape(self._xs)[0]:
            raise RuntimeError('x/y and z data do not match in shape!')

        self._y_scale = y_scale

        self.plot(**self.plot_kwargs)

    @property
    def y_scale(self):
        return self._y_scale

    @y_scale.setter
    def y_scale(self, value):
        """ Reset the y_scale and replot. """
        self._y_scale = value
        self.plot(**self.plot_kwargs)

    def get_plot(self):
        return self.fig, self.ax

    @plotting_function
    def plot(self, *args, **kwargs):
        """ Actually plot the data and optionally save it. """

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(len(self._xs)):
            ax.plot(self._xs[i, :], self._y_scale * self._ys[i, :] + self._zs[i])

        self.fig = fig
        self.ax = ax
        plt.show()
