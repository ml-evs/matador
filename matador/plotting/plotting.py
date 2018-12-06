# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some useful auxiliary routines
for use in the other plotting functions.

"""


def plotting_function(function):
    """ Wrapper for plotting functions to safely fail on X-forwarding
    errors.
    """

    from functools import wraps
    from matador.utils.print_utils import print_warning, print_failure
    from matador.config import load_custom_settings

    @wraps(function)
    def wrapped_plot_function(*args, **kwargs):
        """ Wrap and return the plotting function. """
        from tkinter import TclError
        result = None
        # if we're going to be saving a figure, switch to Agg to avoid X-forwarding
        saving = False

        try:
            for arg in args:
                if arg.savefig:
                    import matplotlib
                    matplotlib.use('Agg')
                    saving = True
                    break
        except AttributeError:
            pass
        if not saving:
            if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
                import matplotlib
                matplotlib.use('Agg')
                saving = True

        settings = load_custom_settings(kwargs.get('config_fname'), quiet=True, override=kwargs.get('override'))
        try:
            import matplotlib.pyplot as plt
            style = settings.get('plotting', {}).get('default_style')
            if style is None or style == 'matador':
                style = '/'.join(__file__.split('/')[:-1]) + '/../config/matador.mplstyle'
            plt.style.use(style)
            if kwargs.get('debug'):
                print('Using style {}'.format(style))
                print(plt.rcParams)
            result = function(*args, **kwargs)
        except TclError as exc:
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
