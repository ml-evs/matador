""" This submodule implements functions useful for plotting
the results of infrared and Raman spectroscopy calculations.

"""

from matador.plotting.plotting import plotting_function


def read_ir_file(seed):
    """ This function should take the seed filename and
    read in the IR data, probably into a dictionary.

    Parameters:
        seed (str): the filename to read.

    Returns:
        dict: containing IR data

    """

    if seed.endswith('.phonon'):
         seed = seed.replace('.phonon', '')
    with open(seed + '.phonon', 'r') as f:
        flines = f.readlines()
    for ind, line in enumerate(flines):
        if 'Number of branches' in line:
            no_branches = int(line.split()[3])
        if 'END header' in line:
            header_end = ind

    # contains IR and Raman data
    ir_dat = dict()
    for line in flines[header_end+2:header_end+no_branches+2]:
        line_split = line.split()
        if len(line_split) < 2:
            print_failure('Failed to read IR data. Check .phonon file.')
            exit()
        wavenumber = float(line_split[1])
        ir = float(line_split[2])
        # checks for raman data
        if len(line_split) > 3:
            raman = float(line_split[3])
            if wavenumber in ir_dat:
                ir_dat[wavenumber] = [ir_dat[wavenumber][0]+ir, ir_dat[wavenumber][1]+raman]
            else:
                ir_dat[wavenumber] = [ir, raman]
        else:
            if wavenumber in ir_dat:
                ir_dat[wavenumber] = [ir_dat[wavenumber][0]+ir]
            else:
                ir_dat[wavenumber] = [ir]

    return(ir_dat)

# this decorator will automatically set the style of the plot
# to match the user's config/matador defaults, and safely handles
# e.g. X-forwarding and file writing
@plotting_function
def plot_ir_spectrum(seed, ir_ss, ir_bs):
    """ This function plots the IR spectrum found in the given file.

        Parameters:
        bin_width (float): The width of the wavenumber bin in cm^-1. Determined by step_size and bin_scaler

    Keyword Arguments:
        step_size (float): change in wavenumber between points on x-axis,
        bin_scaler (float): scaled to the step_size (1 is minimum vale); used for gaussian broadening,
        
        step_size: change in wavenumber between points on x-axis,
        bin_scaler: scaled to the step_size (1 is minimum vale); used for gaussian broadening,
    """
    import matplotlib.pyplot as plt

    step_size = ir_ss
    bin_scaler = ir_bs
    bin_width = step_size * bin_scaler

    # read IR output
    ir_dat = read_ir_file(seed)
    wavenumbers_castep = list(ir_dat.keys())
    if len(ir_dat[wavenumbers_castep[0]]) == 2:
        raman = True
    else:
        raman = False
    min_wavenumber = int(wavenumbers_castep[0])-(2*step_size)
    max_wavenumber = int(wavenumbers_castep[len(wavenumbers_castep)-1])+(2*step_size)
    wavenumbers_plot = []
    for i in range(min_wavenumber, max_wavenumber, step_size):
        wavenumbers_plot.append(i)

    ir_plot = []
    raman_plot = []
    for i in wavenumbers_plot:
        ir_intensity = 0
        raman_intensity = 0
        for j in wavenumbers_castep:
            if j >= (i-(0.5*bin_width)) and j < (i+(0.5*bin_width)):
                ir_intensity = ir_intensity + ir_dat[j][0]
                if raman == True:
                    raman_intensity = raman_intensity + ir_dat[j][1]
                pass
        ir_plot.append(ir_intensity)
        if raman == True:
            raman_plot.append(raman_intensity)

    ir_max = 0
    raman_max = 0
    for i in wavenumbers_castep:
        if ir_dat[i][0] > ir_max:
            ir_max = ir_dat[i][0]
        if raman == True:
            if ir_dat[i][1] > raman_max:
                raman_max = ir_dat[i][1]

    fig, ax1 = plt.subplots(figsize = (14,8))
    plt.title(seed)
    ax1.plot(wavenumbers_plot, ir_plot, color='#EE3425')
    ax1.set_xlabel('Wavenumbers (cm$^{-1}$)')
    ax1.set_ylabel('IR intensities  ((D/A)**2/amu)', color='#EE3425')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()

    if raman == True:
        ax1.set_ylim(2.1*ir_max,-(ir_max*0.05))
        ax2 = ax1.twinx()
        ax2.plot(wavenumbers_plot, raman_plot, color='#236DE8')
        ax2.set_ylabel('Raman activies (A**4 amu**(-1))', color='#236DE8')
        ax2.set_ylim(-(raman_max*0.05), 2.1*raman_max)

    plt.show()
    #return ax
