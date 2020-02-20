""" This submodule implements functions useful for plotting
the results of infrared and Raman spectroscopy calculations.

"""

from matador.plotting import plotting_function


def read_ir_file(seed):
    """ This function should take the seed filename and
    read in the IR data, probably into a dictionary.

    Parameters:
        seed (str): the filename to read.

    Returns:
        dict: containing IR data

    """

    raise NotImplementedError


# this decorator will automatically set the style of the plot
# to match the user's config/matador defaults, and safely handles
# e.g. X-forwarding and file writing
@plotting_function
def plot_ir_spectrum(seed):
    """ This function plots the IR spectrum found in the given file.
    
    """

    raise NotImplementedError

    # read IR output
    read_ir_file()

    # plot the data
    # fig, ax = plt.subplots()
    # ax.plot(energies, intensity)
    # return ax
