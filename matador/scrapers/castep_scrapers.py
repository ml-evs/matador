# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the scraper functions for CASTEP-related
inputs and outputs.

"""


from collections import defaultdict
from os import stat
from os.path import isfile
from time import strptime
from pwd import getpwuid
import glob
import gzip
from matador.utils.cell_utils import abc2cart, calc_mp_spacing, cart2volume, wrap_frac_coords
from matador.utils.chem_utils import get_stoich
from matador.scrapers.utils import DFTError, CalculationError, scraper_function


@scraper_function
def res2dict(seed, db=True, **kwargs):
    """ Extract available information from .res file; preferably
    used in conjunction with cell or param file.

    Parameters:
        seed (str/list): filename or list of filenames of res file(s)
            (with or without extension).

    Keyword arguments:
        db (bool): whether to fail if unable to scrape energies.

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped data and True,
            if not, then an error string and False.

    """
    res = dict()
    # read .res file into array
    if seed.endswith('.res'):
        seed = seed.replace('.res', '')
    with open(seed + '.res', 'r') as f:
        flines = f.readlines()
    # add .res to source
    res['source'] = []
    res['source'].append(seed + '.res')
    # grab file owner username
    res['user'] = getpwuid(stat(seed + '.res').st_uid).pw_name
    # try to grab ICSD CollCode
    if 'CollCode' in seed:
        res['icsd'] = int(seed.split('CollCode')[-1].split('-')[0].split('_')[0].split('.')[0])
    if '-mp-' in seed:
        res['mp-id'] = int(seed.split('-mp-')[-1].split('-')[0].split('_')[0].split('.')[0])
    # alias special lines in res file
    titl = ''
    cell = ''
    remark = ''
    for line in flines:
        if 'TITL' in line and db:
            # if not db, then don't read title
            titl = line.split()
            if len(titl) != 12:
                raise RuntimeError('missing some TITL info')
        elif 'CELL' in line:
            cell = line.split()
        elif 'REM' in line:
            remark = line.split()
    if cell == '':
        raise RuntimeError('missing CELL info')
    elif titl == '' and db:
        raise RuntimeError('missing TITL')
    if db:
        res['pressure'] = float(titl[2])
        res['cell_volume'] = float(titl[3])
        res['enthalpy'] = float(titl[4])
        res['num_atoms'] = int(titl[7])
        res['space_group'] = titl[8].strip('()')
        res['enthalpy_per_atom'] = res['enthalpy'] / res['num_atoms']
    res['lattice_abc'] = [list(map(float, cell[2:5])), list(map(float, cell[5:8]))]
    # calculate lattice_cart from abc
    res['lattice_cart'] = abc2cart(res['lattice_abc'])
    if 'cell_volume' not in res:
        res['cell_volume'] = cart2volume(res['lattice_cart'])
    res['atom_types'] = []
    res['positions_frac'] = []
    res['site_occupancy'] = []
    for line_no, line in enumerate(flines):
        if 'SFAC' in line:
            i = 1
            while 'END' not in flines[line_no + i] and line_no + i < len(flines):
                cursor = flines[line_no + i].split()
                res['atom_types'].append(cursor[0])
                res['positions_frac'].append(list(map(float, cursor[2:5])))
                res['site_occupancy'].append(float(cursor[5]))
                assert len(res['positions_frac'][-1]) == 3
                i += 1
    res['positions_frac'] = wrap_frac_coords(res['positions_frac'])
    if 'num_atoms' in res:
        assert len(res['atom_types']) == res['num_atoms']
    else:
        res['num_atoms'] = len(res['atom_types'])
    # deal with implicit encapsulation
    if remark:
        if 'NTPROPS' in remark:
            res['cnt_chiral'] = [0, 0]
            res['encapsulated'] = True
            for ind, entry in enumerate(remark):
                if 'chiralN' in entry:
                    res['cnt_chiral'][0] = int(remark[ind + 1].replace(',', ''))
                if 'chiralM' in entry:
                    res['cnt_chiral'][1] = int(remark[ind + 1].replace(',', ''))
                if entry == '\'r\':':
                    res['cnt_radius'] = float(remark[ind + 1].replace(',', ''))
                if entry == '\'z\':':
                    temp_length = remark[ind + 1].replace(',', '')
                    temp_length = temp_length.replace('\n', '')
                    temp_length = temp_length.replace('}', '')
                    res['cnt_length'] = float(temp_length)

    res['stoichiometry'] = get_stoich(res['atom_types'])
    res['num_fu'] = len(res['atom_types']) / sum([elem[1] for elem in res['stoichiometry']])

    return res, True


@scraper_function
def cell2dict(seed, db=True, lattice=False, outcell=False, positions=False, **kwargs):
    """ Extract available information from .cell file; probably
    to be merged with another dict from a .param or .res file.

    Parameters:
        seed (str/list): filename or list of filenames of cell file(s)
            to scrape, with or without extension.

    Keyword arguments:
        db (bool): scrape database quality file
        lattice (bool): scrape lattice vectors
        positions (bool): scrape positions

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped
            data and True, if not, then an error string and False.

    """
    cell = dict()
    if outcell:
        lattice = True
    if seed.endswith('.cell'):
        seed = seed.replace('.cell', '')
    with open(seed + '.cell', 'r') as f:
        flines = f.readlines()
    # add cell file to source
    cell['source'] = []
    cell['source'].append(seed + '.cell')
    for line_no, line in enumerate(flines):
        if line.startswith(('#', '!')):
            continue
        if '#' or '!' in line:
            line = line.split('#')[0].split('!')[0]
        if '%block lattice_cart' in line.lower() and lattice:
            cell['lattice_cart'] = []
            i = 1
            while 'endblock' not in flines[line_no + i].lower():
                if not flines[line_no + i].strip()[0].isalpha():
                    cell['lattice_cart'].append(list(map(float, flines[line_no + i].split())))
                    assert len(cell['lattice_cart'][-1]) == 3, 'Lattice vector does not have enough elements!'
                i += 1
            assert len(cell['lattice_cart']) == 3, 'Wrong number of lattice vectors!'
            cell['cell_volume'] = cart2volume(cell['lattice_cart'])
        elif '%block species_pot' in line.lower():
            cell['species_pot'] = dict()
            i = 1
            while 'endblock' not in flines[line_no + i].lower():
                if db:
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        flines[line_no+i].split()[1].split('/')[-1]
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        cell['species_pot'][flines[line_no+i].split()[0]].replace('()', '')
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        cell['species_pot'][flines[line_no+i].split()[0]].replace('[]', '')
                else:
                    pspot_libs = ['C7', 'C8', 'C9', 'C17', 'C18', 'MS', 'HARD',
                                  'QC5', 'NCP', 'NCP18', 'NCP17', 'NCP9']
                    if flines[line_no + i].upper().split()[0] in pspot_libs:
                        cell['species_pot']['library'] = flines[line_no + i].upper().split()[0]
                    else:
                        cell['species_pot'][flines[line_no + i].split()[0]] = flines[line_no + i].split()[1]
                i += 1
        elif '%block cell_constraints' in line.lower():
            cell['cell_constraints'] = []
            for j in range(2):
                cell['cell_constraints'].append(list(map(int, flines[line_no + j + 1].split())))
        elif '%block hubbard_u' in line.lower():
            cell['hubbard_u'] = defaultdict(list)
            i = 0
            while 'endblock' not in flines[line_no + i].lower():
                line = flines[line_no + i]
                if line == 'eV' or len(line.split()) < 3:
                    i += 1
                    continue
                else:
                    atom = line.split()[0]
                    orbital = line.split()[1].replace(':', '')
                    shift = float(line.split()[-1])
                    atom = line.split()[0]
                    cell['hubbard_u'][atom] = dict()
                    cell['hubbard_u'][atom][orbital] = shift
                    i += 1
        elif '%block external_pressure' in line.lower():
            cell['external_pressure'] = []
            i = 1
            while 'endblock' not in flines[line_no + i].lower():
                if not flines[line_no + i].strip()[0].isalpha():
                    flines[line_no+i] = flines[line_no+i].replace(',', '')
                    cell['external_pressure'].append(list(map(float, flines[line_no+i].split())))
                i += 1
        # parse kpoints
        elif 'kpoints_mp_spacing' in line.lower() or 'kpoint_mp_spacing' in line.lower():
            if 'spectral_kpoints_mp_spacing' in line.lower() or 'spectral_kpoint_mp_spacing' in line.lower():
                cell['spectral_kpoints_mp_spacing'] = float(line.split()[-1])
            elif 'supercell_kpoints_mp_spacing' in line.lower() or 'supercell_kpoint_mp_spacing' in line.lower():
                cell['supercell_kpoints_mp_spacing'] = float(line.split()[-1])
            elif 'phonon_kpoints_mp_spacing' in line.lower() or 'phonon_kpoint_mp_spacing' in line.lower():
                cell['phonon_kpoint_mp_spacing'] = float(line.split()[-1])
            elif 'phonon_fine_kpoints_mp_spacing' in line.lower() or 'phonon_fine_kpoint_mp_spacing' in line.lower():
                cell['phonon_fine_kpoint_mp_spacing'] = float(line.split()[-1])
            else:
                cell['kpoints_mp_spacing'] = float(line.split()[-1])
        elif 'kpoints_mp_grid' in line.lower() or 'kpoint_mp_grid' in line.lower():
            if 'spectral_kpoints_mp_grid' in line.lower() or 'spectral_kpoint_mp_grid' in line.lower():
                cell['spectral_kpoints_mp_grid'] = list(map(int, line.split()[-3:]))
            elif 'phonon_kpoints_mp_grid' in line.lower() or 'phonon_kpoint_mp_grid' in line.lower():
                cell['phonon_kpoint_mp_grid'] = list(map(int, line.split()[-3:]))
            elif 'phonon_fine_kpoints_mp_grid' in line.lower() or 'phonon_fine_kpoint_mp_grid' in line.lower():
                cell['phonon_fine_kpoint_mp_grid'] = list(map(int, line.split()[-3:]))
            else:
                cell['kpoints_mp_grid'] = list(map(int, line.split()[-3:]))
        elif 'kpoints_mp_offset' in line.lower() or 'kpoint_mp_offset' in line.lower():
            if 'spectral_kpoints_mp_offset' in line.lower() or 'spectral_kpoint_mp_offset' in line.lower():
                cell['spectral_kpoints_mp_offset'] = list(map(float, line.split()[-3:]))
            elif 'phonon_kpoints_mp_offset' in line.lower() or 'phonon_kpoint_mp_offset' in line.lower():
                # this is a special case where phonon_kpointS_mp_offset doesn't exist
                cell['phonon_kpoint_mp_offset'] = list(map(float, line.split()[-3:]))
            elif 'phonon_fine_kpoints_mp_offset' in line.lower() or 'phonon_fine_kpoint_mp_offset' in line.lower():
                cell['phonon_fine_kpoint_mp_offset'] = list(map(float, line.split()[-3:]))
            else:
                cell['kpoints_mp_offset'] = list(map(float, line.split()[-3:]))
        elif '%block spectral_kpoints_path' in line.lower() or '%block spectral_kpoint_path' in line.lower():
            i = 1
            cell['spectral_kpoints_path'] = []
            while '%endblock' not in flines[line_no + i].lower():
                cell['spectral_kpoints_path'].append(list(map(float, flines[line_no + i].split()[:3])))
                i += 1
        elif '%block spectral_kpoints_list' in line.lower() or '%block spectral_kpoint_list' in line.lower():
            i = 1
            cell['spectral_kpoints_list'] = []
            while '%endblock' not in flines[line_no + i].lower():
                cell['spectral_kpoints_list'].append(list(map(float, flines[line_no + i].split()[:4])))
                i += 1
        elif '%block phonon_fine_kpoints_list' in line.lower() or '%block phonon_fine_kpoint_list' in line.lower():
            i = 1
            # this is a special case where phonon_fine_kpointS_list doesn't exist
            cell['phonon_fine_kpoint_list'] = []
            while '%endblock' not in flines[line_no + i].lower():
                cell['phonon_fine_kpoint_list'].append(list(map(float, flines[line_no + i].split()[:4])))
                i += 1
        elif '%block phonon_supercell_matrix' in line.lower():
            cell['phonon_supercell_matrix'] = []
            i = 1
            while 'endblock' not in flines[line_no + i].lower():
                cell['phonon_supercell_matrix'].append(list(map(int, flines[line_no + i].split())))
                assert len(cell['phonon_supercell_matrix'][-1]) == 3, 'Supercell matrix row does not have enough elements!'
                i += 1
            assert len(cell['phonon_supercell_matrix']) == 3, 'Wrong supercell matrix shape!'
        elif not db:
            if '%block positions_frac' in line.lower():
                atomic_init_spins = defaultdict(list)
                i = 1
                if positions:
                    cell['atom_types'] = []
                    cell['positions_frac'] = []
                while '%endblock positions_frac' not in flines[line_no + i].lower():
                    line = flines[line_no + i].split()
                    if positions:
                        cell['atom_types'].append(line[0])
                        cell['positions_frac'].append(list(map(float, line[1:4])))
                    if 'spin=' in flines[line_no + i].lower():
                        split_line = flines[line_no + i].split()
                        atomic_init_spins[split_line[0]] = \
                            split_line[-1].lower().replace('spin=', '')
                    i += 1
                if atomic_init_spins:
                    cell['atomic_init_spins'] = atomic_init_spins
                if positions:
                    cell['num_atoms'] = len(cell['atom_types'])
                    for ind, pos in enumerate(cell['positions_frac']):
                        for k in range(3):
                            if pos[k] > 1 or pos[k] < 0:
                                cell['positions_frac'][ind][k] %= 1
            elif 'fix_com' in line.lower():
                cell['fix_com'] = line.split()[-1]
            elif 'symmetry_generate' in line.lower():
                cell['symmetry_generate'] = True
            elif 'symmetry_tol' in line.lower():
                cell['symmetry_tol'] = float(line.split()[-1])
            elif 'snap_to_symmetry' in line.lower():
                cell['snap_to_symmetry'] = True
            elif 'quantisation_axis' in line.lower():
                cell['quantisation_axis'] = list(map(int, line.split()[1:]))
            elif 'positions_noise' in line.lower():
                cell['positions_noise'] = float(line.split()[-1])
            elif 'cell_noise' in line.lower():
                cell['cell_noise'] = float(line.split()[-1])
            elif 'kpoints_path' in line.lower() or 'kpoint_path' in line.lower():
                if 'spectral_kpoints_path_spacing' in line.lower() or 'spectral_kpoint_path_spacing' in line.lower():
                    cell['spectral_kpoints_path_spacing'] = float(line.split()[-1])
                elif 'phonon_fine_kpoints_path_spacing' in line.lower() or 'phonon_fine_kpoint_path_spacing' in line.lower():
                    cell['phonon_fine_kpoint_path_spacing'] = float(line.split()[-1])
                elif 'kpoints_path_spacing' in line.lower() or 'kpoint_path_spacing' in line.lower():
                    cell['kpoints_path_spacing'] = float(line.split()[-1])

    if 'external_pressure' not in cell or not cell['external_pressure']:
        cell['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]

    if db:
        for species in cell['species_pot']:
            if 'OTF' in cell['species_pot'][species].upper():
                pspot_seed = ''
                for directory in seed.split('/')[:-1]:
                    pspot_seed += directory + '/'
                # glob for all .usp files with format species_*OTF.usp
                pspot_seed += species + '_*OTF.usp'
                for globbed in glob.glob(pspot_seed):
                    if isfile(globbed):
                        cell['species_pot'].update(usp2dict(globbed))
                        break

    return cell, True


@scraper_function
def param2dict(seed, db=True, **kwargs):
    """ Extract available information from .param file; probably
    to be merged with other dicts from other files.

    Parameters:
        seed (str/list): param filename or list of filenames with or
            without file extension

    Keyword arguments:
        db (bool): if True, only scrape relevant info, otherwise scrape all

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped data and True,
            if not, then an error string and False.

    """
    from matador.utils.castep_params import CASTEP_PARAMS
    param = dict()
    if seed.endswith('.param'):
        seed = seed.replace('.param', '')
    with open(seed + '.param', 'r') as f:
        flines = f.readlines()
    param['source'] = []
    param['source'].append(seed + '.param')
    # exclude some useless info if importing to db
    scrub_list = ['checkpoint', 'write_bib', 'mix_history_length',
                  'fix_occupancy', 'page_wvfns', 'num_dump_cycles',
                  'backup_interval', 'geom_max_iter', 'fixed_npw',
                  'write_cell_structure', 'bs_write_eigenvalues',
                  'calculate_stress', 'opt_strategy', 'max_scf_cycles']
    false_str = ['False', 'false', '0']
    splitters = [':', '=', '\t', ' ']
    unrecognised = []
    for _, line in enumerate(flines):
        if '#' or '!' in line:
            line = line.split('#')[0].split('!')[0]
        line = line.lower()
        # skip blank lines and comments
        if line.startswith(('#', '!')) or not line.strip():
            continue
        else:
            # if scraping to db, ignore "rubbish"
            if db:
                if [rubbish for rubbish in scrub_list if rubbish in line]:
                    continue
            # read all other parameters in
            for splitter in splitters:
                if splitter in line:
                    keyword = line.split(splitter)[0].strip()
                    value = line.split(splitter)[-1].strip()
                    if keyword.lower() not in CASTEP_PARAMS:
                        unrecognised.append(keyword.lower())
                    param[keyword] = value
                    # deal with edge cases
                    if 'spin_polarised' in line:
                        param['spin_polarized'] = param['spin_polarised']
                    if 'spin_polarized' in line or 'spin_polarised' in line:
                        if [false for false in false_str if false in param['spin_polarized']]:
                            param['spin_polarized'] = False
                        else:
                            param['spin_polarized'] = True

                    if 'true' in value.lower():
                        param[keyword] = True
                    elif 'false' in value.lower():
                        param[keyword] = False

                    if 'geom_max_iter' == keyword:
                        param['geom_max_iter'] = int(param['geom_max_iter'])

                    if 'cut_off_energy' in line and 'mix_cut_off_energy' not in line:
                        temp_cut_off = (param['cut_off_energy'].lower().replace('ev', ''))
                        temp_cut_off = temp_cut_off.strip()
                        param['cut_off_energy'] = float(temp_cut_off)
                    elif 'xc_functional' in line:
                        param['xc_functional'] = param['xc_functional'].upper()
                    elif 'perc_extra_bands' in line:
                        param['perc_extra_bands'] = float(param['perc_extra_bands'])
                    elif 'geom_force_tol' in line:
                        param['geom_force_tol'] = float(param['geom_force_tol'])
                    elif 'elec_energy_tol' in line:
                        temp = (param['elec_energy_tol'].lower().replace('ev', ''))
                        temp = temp.strip()
                        param['elec_energy_tol'] = float(temp)
                    break

    if unrecognised:
        raise RuntimeError('Found several unrecognised parameters: {}'.format(unrecognised))

    return param, True


@scraper_function
def castep2dict(seed, db=True, intermediates=False, **kwargs):
    """ From seed filename, create dict of the most relevant
    information about a calculation.

    Parameters:
        seed (str/list): filename or list of filenames of castep file(s)

    Keyword arguments:
        db (bool): whether to error on missing relaxation info
        intermediates (bool): instead of a single dict containing the relaxed structure
            return a list of snapshots found in .castep file

    Returns:
        (dict/str/list, bool): if successful, a dictionary or list of dictionaries
            containing scraped data and True, if not, then an error string and False.

    """
    castep = dict()
    # read .castep, .history or .history.gz file
    if '.gz' in seed:
        with gzip.open(seed, 'r') as f:
            flines = f.readlines()
    else:
        with open(seed, 'r') as f:
            flines = f.readlines()
    # set source tag to castep file
    castep['source'] = []
    castep['source'].append(seed)
    # grab file owner
    castep['user'] = getpwuid(stat(seed).st_uid).pw_name
    if 'CollCode' in seed:
        castep['icsd'] = int(seed.split('CollCode')[-1].split('-')[0].split('.')[0].split('_')[0])
    if '-mp-' in seed:
        castep['mp-id'] = int(seed.split('-mp-')[-1].split('-')[0].split('.')[0].split('_')[0])
    # wrangle castep file for parameters in 3 passes:
    # once forwards to get number and types of atoms
    castep.update(_castep_scrape_atoms(flines, castep))
    # once backwards to get the final parameter set for the calculation
    castep.update(_castep_scrape_final_parameters(flines, castep))
    # once more forwards, from the final step, to get the final structure

    # task specific options
    if db and 'geometry' not in castep['task']:
        raise RuntimeError('CASTEP file does not contain GO calculation')

    if not db and 'thermo' in castep['task'].lower():
        castep.update(_castep_scrape_thermo_data(flines, castep))

    # only scrape snapshots/number of intermediates if requested,
    # or if not in db mode
    if intermediates or not db:
        snapshots, castep['geom_iter'] = _castep_scrape_all_snapshots(flines)
        if intermediates:
            castep['intermediates'] = snapshots

    castep.update(_castep_scrape_final_structure(flines, castep, db=db))
    castep.update(_castep_scrape_metadata(flines, castep))

    if 'positions_frac' not in castep or not castep['positions_frac']:
        raise CalculationError('Could not find positions')

    # unfortunately CASTEP does not write forces when there is only one atom
    if 'forces' not in castep and castep['num_atoms'] == 1 and 'geometry' in castep['task']:
        castep['forces'] = [[0, 0, 0]]

    # fill in any missing fields with filler
    if 'total_energy' not in castep:
        castep['total_energy'] = 'xxx'
        castep['total_energy_per_atom'] = 'xxx'
    if 'enthalpy' not in castep:
        castep['enthalpy'] = 'xxx'
        castep['enthalpy_per_atom'] = 'xxx'
    if 'pressure' not in castep:
        castep['pressure'] = 'xxx'
    if 'cell_volume' not in castep:
        castep['cell_volume'] = 'xxx'
    if 'space_group' not in castep:
        castep['space_group'] = 'xxx'

    # finally check for pseudopotential files if OTF is present in species_pot
    if db:
        for species in castep['species_pot']:
            if 'OTF' in castep['species_pot'][species].upper():
                pspot_seed = ''
                for directory in seed.split('/')[:-1]:
                    pspot_seed += directory + '/'
                # glob for all .usp files with format species_*OTF.usp
                pspot_seed += species + '_*OTF.usp'
                for globbed in glob.glob(pspot_seed):
                    if isfile(globbed):
                        castep['species_pot'].update(usp2dict(globbed))

    # check that any optimized results were saved and raise errors if not
    if not castep.get('optimised'):
        castep['optimised'] = False
        if db:
            # if importing to db, skip unconverged structure
            raise DFTError('CASTEP GO failed to converge.')

    return castep, True


@scraper_function
def bands2dict(seed, summary=False, gap=False, external_efermi=None, **kwargs):
    """ Parse a CASTEP bands file into a dictionary.

    Parameters:
        seed (str/list): filename of list of filenames to be scraped.

    Keyword arguments:
        summary (bool): print info about bandgap.
        gap (bool): re-compute bandgap info.
        external_efermi (float): override the Fermi energy with this value (eV)

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped data and True,
            if not, then an error string and False.

    """
    from matador.utils.chem_utils import HARTREE_TO_EV, BOHR_TO_ANGSTROM
    from matador.utils.cell_utils import frac2cart, real2recip
    import numpy as np

    verbosity = kwargs.get('verbosity', 0)

    bs = dict()
    with open(seed, 'r') as f:
        # read whole file into RAM, typically ~ 1 MB
        flines = f.readlines()
    header = flines[:9]
    data = flines[9:]

    seed = seed.replace('.bands', '')
    bs['source'] = [seed + '.bands']

    bs['num_kpoints'] = int(header[0].split()[-1])
    bs['num_spins'] = int(header[1].split()[-1])
    bs['num_electrons'] = float(header[2].split()[-1])
    bs['num_bands'] = int(header[3].split()[-1])
    if external_efermi is None:
        bs['fermi_energy_Ha'] = float(header[4].split()[-1])
        bs['fermi_energy'] = bs['fermi_energy_Ha'] * HARTREE_TO_EV
    else:
        bs['fermi_energy'] = external_efermi
    bs['lattice_cart'] = []
    for i in range(3):
        bs['lattice_cart'].append([BOHR_TO_ANGSTROM * float(elem) for elem in header[6 + i].split()])
    bs['kpoint_path'] = np.zeros((bs['num_kpoints'], 3))
    bs['eigenvalues_k_s'] = np.empty((bs['num_spins'], bs['num_bands'], bs['num_kpoints']))

    if verbosity > 2:
        print('Found {}'.format(bs['num_kpoints']))

    for nk in range(bs['num_kpoints']):
        kpt_ind = nk * (bs['num_spins'] * bs['num_bands'] + bs['num_spins'] + 1)
        bs['kpoint_path'][int(data[kpt_ind].split()[1]) - 1] = (
            np.asarray([float(elem) for elem in data[kpt_ind].split()[-4:-1]])
        )
        for nb in range(bs['num_bands']):
            for ns in range(bs['num_spins']):
                line_number = kpt_ind + 2 + ns + (ns * bs['num_bands']) + nb
                bs['eigenvalues_k_s'][ns][nb][int(data[kpt_ind].split()[1]) - 1] = (
                    float(data[line_number].strip()))
    bs['eigenvalues_k_s'] *= HARTREE_TO_EV
    bs['eigenvalues_k_s'] -= bs['fermi_energy']

    cart_kpts = np.asarray(frac2cart(real2recip(bs['lattice_cart']), bs['kpoint_path']))
    bs['cart_kpoints'] = cart_kpts
    kpts_diff = np.zeros((len(cart_kpts) - 1))
    kpts_diff_set = set()
    for i in range(len(cart_kpts) - 1):
        kpts_diff[i] = np.sqrt(np.sum((cart_kpts[i] - cart_kpts[i + 1])**2))
        kpts_diff_set.add(kpts_diff[i])
    bs['kpoint_path_spacing'] = np.median(kpts_diff)

    # create list containing kpoint indices of discontinuous branches through k-space
    bs['kpoint_branches'] = []
    current_branch = []
    for ind, point in enumerate(cart_kpts):
        if ind == 0:
            current_branch.append(ind)
        elif ind == len(cart_kpts) - 1:
            bs['kpoint_branches'].append(current_branch)
            continue

        if np.sqrt(np.sum((point - cart_kpts[ind + 1])**2)) < 10 * bs['kpoint_path_spacing']:
            current_branch.append(ind + 1)
        else:
            bs['kpoint_branches'].append(current_branch)
            current_branch = [ind + 1]
    assert sum([len(branch) for branch in bs['kpoint_branches']]) == bs['num_kpoints']

    if verbosity > 2:
        print('Found branch structure', [(branch[0], branch[-1]) for branch in bs['kpoint_branches']])

    if gap and bs['num_spins'] == 1:

        vbm = -1e10
        cbm = 1e10
        cbm_pos = []
        vbm_pos = []
        eps = 1e-6
        if bs['num_spins'] == 1:
            # calculate indirect gap
            for _, branch in enumerate(bs['kpoint_branches']):
                for nb in range(bs['num_bands']):
                    band = bs['eigenvalues_k_s'][0][nb][branch]
                    band_branch_min = np.min(band)
                    band_branch_max = np.max(band)
                    band_branch_argmin = np.where(band <= band_branch_min + eps)[0]
                    band_branch_argmax = np.where(band >= band_branch_max - eps)[0]
                    if band_branch_max < 0 and band_branch_max > vbm + eps:
                        vbm = band_branch_max
                        vbm_pos = [branch[max_ind] for max_ind in band_branch_argmax]
                    elif band_branch_max < 0 and band_branch_max >= vbm - eps:
                        vbm = band_branch_max
                        vbm_pos.extend([branch[val] for val in band_branch_argmax])
                    if band_branch_min > 0 and band_branch_min < cbm - eps:
                        cbm = band_branch_min
                        cbm_pos = [branch[min_ind] for min_ind in band_branch_argmin]
                    elif band_branch_min > 0 and band_branch_min <= cbm + eps:
                        cbm = band_branch_min
                        cbm_pos.extend([branch[val] for val in band_branch_argmin])
                    if band_branch_max > 0 and band_branch_min < 0:
                        vbm = 0
                        cbm = 0
                        vbm_pos = [0]
                        cbm_pos = [0]
                        break

            smallest_diff = 1e10
            for _cbm_pos in cbm_pos:
                for _vbm_pos in vbm_pos:
                    if abs(_vbm_pos - _cbm_pos) < smallest_diff:
                        tmp_cbm_pos = _cbm_pos
                        tmp_vbm_pos = _vbm_pos
                        smallest_diff = abs(_vbm_pos - _cbm_pos)
            cbm_pos = tmp_cbm_pos
            vbm_pos = tmp_vbm_pos
            bs['valence_band_min'] = vbm
            bs['conduction_band_max'] = cbm
            bs['band_gap'] = cbm - vbm
            bs['band_gap_path'] = [bs['kpoint_path'][cbm_pos], bs['kpoint_path'][vbm_pos]]
            bs['band_gap_path_inds'] = [cbm_pos, vbm_pos]
            bs['gap_momentum'] = np.sqrt(np.sum((bs['cart_kpoints'][cbm_pos] - bs['cart_kpoints'][vbm_pos])**2))

            # calculate direct gap
            direct_gaps = np.zeros((len(bs['kpoint_path'])))
            direct_cbms = np.zeros((len(bs['kpoint_path'])))
            direct_vbms = np.zeros((len(bs['kpoint_path'])))
            for ind, _ in enumerate(bs['kpoint_path']):
                direct_cbm = 1e10
                direct_vbm = -1e10
                for nb in range(bs['num_bands']):
                    band_eig = bs['eigenvalues_k_s'][0][nb][ind]
                    if band_eig < 0 and band_eig >= direct_vbm:
                        direct_vbm = band_eig
                    if band_eig > 0 and band_eig <= direct_cbm:
                        direct_cbm = band_eig
                direct_gaps[ind] = direct_cbm - direct_vbm
                direct_cbms[ind] = direct_cbm
                direct_vbms[ind] = direct_vbm
            bs['direct_gap'] = np.min(direct_gaps)
            bs['direct_conduction_band_max'] = direct_cbms[np.argmin(direct_gaps)]
            bs['direct_valence_band_min'] = direct_vbms[np.argmin(direct_gaps)]
            bs['direct_gap'] = np.min(direct_gaps)
            bs['direct_gap_path'] = 2 * [bs['kpoint_path'][np.argmin(direct_gaps)]]
            bs['direct_gap_path_inds'] = 2 * [np.argmin(direct_gaps)]

        if np.abs(bs['direct_gap'] - bs['band_gap']) < 1e-5:
            bs['valence_band_min'] = direct_vbm
            bs['conduction_band_max'] = direct_cbm
            bs['band_gap_path_inds'] = bs['direct_gap_path_inds']
            cbm_pos = bs['direct_gap_path_inds'][0]
            vbm_pos = bs['direct_gap_path_inds'][1]
            bs['band_gap_path'] = bs['direct_gap_path']
            bs['gap_momentum'] = np.sqrt(np.sum((bs['cart_kpoints'][cbm_pos] - bs['cart_kpoints'][vbm_pos])**2))
            assert bs['gap_momentum'] == 0

        if summary:
            print('Read bs for {}.'.format(seed))
            if bs['band_gap'] == 0:
                print('The structure is metallic.')
            elif bs['band_gap_path_inds'][0] == bs['band_gap_path_inds'][1]:
                print('Band gap is direct with size {:5.5f} eV'.format(bs['band_gap']), end=' ')
                print('and lies at {}'.format(bs['direct_gap_path'][0]))
            else:
                print('Band gap is indirect with size {:5.5f} eV'.format(bs['band_gap']), end=' ')
                print('between {} and {}'.format(bs['kpoint_path'][cbm_pos], bs['kpoint_path'][vbm_pos]), end=' ')
                print('corresponding to a wavenumber of {:5.5f} eV/A'.format(bs['gap_momentum']))
                print('The smallest direct gap has size {:5.5f} eV'.format(bs['direct_gap']), end=' ')
                print('and lies at {}'.format(bs['direct_gap_path'][0]))

    return bs, True


@scraper_function
def optados2dict(seed, **kwargs):
    """ Scrape optados output file (*.*.dat) or (*.pdos.*.dat)
    for DOS, projectors and projected DOS/dispersion.

    Parameters:
        seed (str/list): optados filename or list of filenames.

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped data and True,
            if not, then an error string and False.

    """
    import numpy as np
    optados = dict()
    is_pdos = False
    is_pdis = False
    is_spin_dos = False
    with open(seed, 'r') as f:
        flines = f.readlines()

    header = []
    for ind, line in enumerate(flines):
        if not line.strip().startswith('#') or 'K-point' in line:
            break
        if 'Partial' in line:
            is_pdos = True
        elif 'spin' in line:
            is_spin_dos = True
        elif 'Projected Dispersion' in line:
            is_pdis = True
        else:
            header.append(line)

    flines = flines[ind:]

    if not is_pdis:
        data = np.loadtxt(seed, comments='#')
        optados['energies'] = data[:, 0]

    if is_pdos or is_pdis:
        projectors = []
        # get pdos labels
        for ind, line in enumerate(header):
            if 'Projector:' in line:
                # skip current line and column headings
                j = 2
                elements = []
                ang_mom_channels = []
                while ind + j + 1 < len(header) and 'Projector:' not in header[ind + j + 1]:
                    elements.append(header[ind + j].split()[1])
                    ang_mom_channels.append(header[ind + j].split()[3])
                    j += 1
                projector_label = []
                if len(set(elements)) == 1:
                    projector_label.append(elements[0])
                else:
                    projector_label.append(None)

                if len(set(ang_mom_channels)) == 1:
                    projector_label.append(ang_mom_channels[0])
                else:
                    projector_label.append(None)

                projector_label = tuple(projector_label)
                projectors.append(projector_label)

        optados['num_projectors'] = len(projectors)
        print('Found projectors:', projectors)
        optados['projectors'] = projectors

    if is_pdos:
        # get pdos values
        optados['pdos'] = dict()
        optados['dos'] = np.zeros_like(data[:, 0])
        for i, projector in enumerate(projectors):
            optados['pdos'][projector] = data[:, i + 1]
            optados['dos'] += data[:, i + 1]

    elif is_spin_dos:
        optados['spin_dos'] = dict()
        optados['spin_dos']['up'] = data[:, 1]
        optados['spin_dos']['down'] = data[:, 2]
        optados['dos'] = np.abs(optados['spin_dos']['up']) + np.abs(optados['spin_dos']['down'])

    elif is_pdis:
        optados['pdis'] = []
        optados['kpoints'] = []
        optados['eigenvalues'] = []
        # get kpoints and count number of bands
        kpt_ind = -1
        for i, line in enumerate(flines):
            if 'K-point' in line:
                optados['kpoints'].append([float(val) for val in line.split()[-3:]])
                if kpt_ind == -1:
                    kpt_ind = i
                else:
                    if not optados.get('num_bands'):
                        optados['num_bands'] = i - kpt_ind - 1

        optados['num_kpoints'] = len(optados['kpoints'])

        for nk in range(optados['num_kpoints']):
            eigs = []
            pdis = []
            for nb in range(0, optados['num_bands']):
                eigs.append(float(flines[nk*(optados['num_bands']+1) + nb + 1].split()[0]))
                pdis.append([float(val) for val in flines[nk*(optados['num_bands']+1) + 1 + nb].split()[1:]])
            optados['eigenvalues'].append(eigs)
            optados['pdis'].append(pdis)

    else:
        optados['dos'] = data[:, 1]

    return optados, True


@scraper_function
def phonon2dict(seed, **kwargs):
    """ Parse a CASTEP phonon file into a dictionary.

    Parameters:
        seed (str/list): phonon filename or list of filenames.

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped data and True,
            if not, then an error string and False.

    """
    import numpy as np
    from matador.utils.cell_utils import frac2cart, real2recip

    with open(seed, 'r') as f:
        # read whole file into RAM, typically <~ 1 MB
        flines = f.readlines()

    verbosity = kwargs.get('verbosity', 0)

    ph = dict()
    seed = seed.replace('.phonon', '')
    ph['source'] = [seed + '.phonon']

    for line_no, line in enumerate(flines):
        line = line.lower()
        if 'number of ions' in line:
            ph['num_atoms'] = int(line.split()[-1])
        elif 'number of branches' in line:
            ph['num_branches'] = int(line.split()[-1])
        elif 'number of wavevectors' in line:
            ph['num_qpoints'] = int(line.split()[-1])
        elif 'frequencies in' in line:
            ph['freq_unit'] = line.split()[-1]
        elif 'unit cell vectors' in line:
            ph['lattice_cart'] = []
            for i in range(3):
                ph['lattice_cart'].append([float(elem) for elem in flines[line_no + i + 1].split()])
                assert len(ph['lattice_cart'][-1]) == 3
        elif 'fractional co-ordinates' in line:
            ph['positions_frac'] = []
            ph['atom_types'] = []
            ph['atom_masses'] = []
            i = 1
            while 'END header' not in flines[line_no + i]:
                ph['positions_frac'].append([float(elem) for elem in flines[line_no + i].split()[1:4]])
                ph['atom_types'].append(flines[line_no + i].split()[-2])
                ph['atom_masses'].append(float(flines[line_no + i].split()[-1]))
                assert len(ph['positions_frac'][-1]) == 3
                i += 1
            assert len(ph['positions_frac']) == ph['num_atoms']
            assert len(ph['atom_types']) == ph['num_atoms']
        elif 'end header' in line:
            data = flines[line_no + 1:]

    ph['phonon_kpoint_list'] = []
    ph['eigenvalues_q'] = np.zeros((1, ph['num_branches'], ph['num_qpoints']))
    line_offset = ph['num_branches'] * (ph['num_atoms'] + 1) + 3
    for qind in range(ph['num_qpoints']):
        ph['phonon_kpoint_list'].append([float(elem) for elem in data[qind * line_offset].split()[2:]])
        for i in range(1, ph['num_branches'] + 1):
            if i == 1:
                assert data[qind * line_offset + i].split()[0] == '1'
            ph['eigenvalues_q'][0][i - 1][qind] = float(data[qind * line_offset + i].split()[-1])

    ph['qpoint_path'] = np.asarray([qpt[0:3] for qpt in ph['phonon_kpoint_list']])
    ph['qpoint_weights'] = [qpt[3] for qpt in ph['phonon_kpoint_list']]
    ph['softest_mode_freq'] = np.min(ph['eigenvalues_q'])

    cart_kpts = np.asarray(frac2cart(real2recip(ph['lattice_cart']), ph['qpoint_path']))
    ph['cart_qpoints'] = cart_kpts
    kpts_diff = np.zeros((len(cart_kpts) - 1))
    kpts_diff_set = set()
    for i in range(len(cart_kpts) - 1):
        kpts_diff[i] = np.sqrt(np.sum((cart_kpts[i] - cart_kpts[i + 1])**2))
        kpts_diff_set.add(kpts_diff[i])
    ph['qpoint_path_spacing'] = np.median(kpts_diff)

    # create list containing qpoint indices of discontinuous branches through k-space
    ph['qpoint_branches'] = []
    current_branch = []
    for ind, point in enumerate(cart_kpts):
        if ind == 0:
            current_branch.append(ind)
        elif ind == len(cart_kpts) - 1:
            ph['qpoint_branches'].append(current_branch)
            continue

        if np.sqrt(np.sum((point - cart_kpts[ind + 1])**2)) < 10 * ph['qpoint_path_spacing']:
            current_branch.append(ind + 1)
        else:
            ph['qpoint_branches'].append(current_branch)
            current_branch = [ind + 1]

    assert sum([len(branch) for branch in ph['qpoint_branches']]) == ph['num_qpoints']

    if verbosity > 0:
        print('{} sucessfully scraped with {} q-points.'.format(seed, ph['num_qpoints']))

    return ph, True


@scraper_function
def usp2dict(seed, **kwargs):
    """ Extract pseudopotential string from a CASTEP
    OTF .USP file.

    Parameters:
        seed (str/list): filename of usp file, or list of filenames.

    Returns:
        dict: partial species_pot dict from usp file.

    """
    species_pot = dict()
    with open(seed, 'r') as f:
        flines = f.readlines()
        for line_no, line in enumerate(flines):
            if 'Pseudopotential Report' in line:
                i = 0
                while i + line_no < len(flines) - 3:
                    if 'Pseudopotential Report' in flines[line_no + i]:
                        i += 2
                        elem = flines[line_no + i].split(':')[1].split()[0]
                    elif 'core correction' in flines[line_no + i]:
                        i += 2
                        species_pot[elem] = flines[line_no + i].strip().split()[1]
                        # check next line for wrapped definition
                        if flines[line_no + i + 1].strip().startswith('--------'):
                            break
                        else:
                            species_pot[elem] += flines[line_no + i + 1].strip().split()[1]
                            break
                    i += 1
    species_pot[elem] = species_pot[elem].replace('"', '')
    species_pot[elem] = species_pot[elem].replace('[]', '')
    return species_pot


def _castep_scrape_thermo_data(flines, castep):
    """ Scrape the data from a CASTEP Thermodynamics claculation.

    Note:
        This only scrapes from Thermodynamics section currently,
        NOT Atomic Displacement Parameters

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Returns:
        dict: dictionary updated with scraped thermodynamics data
              in format doc['energy_type'] = {300:-0.01,200:-0.05}

    """
    for line_no, line in enumerate(flines):
        if 'Number of temperature values' in line:
            castep['num_temp_vals'] = int(line.split(':')[-1].strip())
        elif 'Initial temperature' in line:
            castep['temp_init'] = float(line.split(':')[1].strip().split(' ')[0])
        elif 'Final temperature' in line:
            castep['temp_final'] = float(line.split(':')[1].strip().split(' ')[0])
        elif 'Spacing between temperature values' in line:
            castep['temp_spacing'] = float(line.split(':')[1].strip().split(' ')[0])
        elif 'Zero-point energy' in line:
            castep['zero_point_E'] = float(line.split("=")[1].strip().split(' ')[0])
        elif 'T(K)' and 'E(eV)' in line:
            castep['thermo_temps'] = []  # temperatures calculation was done at
            castep['thermo_enthalpy'] = {}  # enthalpy E(eV)
            castep['thermo_free_energy'] = {}  # free energy F(eV)
            castep['thermo_entropy'] = {}  # entropy S(J/mol/K)
            castep['thermo_heat_cap'] = {}  # heat capacity Cv(J/mol/K)
            i = 2
            while True:
                if not flines[line_no + i + 1].strip():
                    break
                else:
                    temp_line = flines[line_no + i].split()
                    castep['thermo_temps'].append(float(temp_line[0]))
                    castep['thermo_enthalpy'][float(temp_line[0])] = float(temp_line[1])
                    castep['thermo_free_energy'][float(temp_line[0])] = float(temp_line[2])
                    castep['thermo_entropy'][float(temp_line[0])] = float(temp_line[3])
                    castep['thermo_heat_cap'][float(temp_line[0])] = float(temp_line[4])
                i += 1

    return castep


def _castep_scrape_atoms(flines, castep):
    """ Iterate forwards through flines to scrape atomic types and
    initial positions, and get a preliminary lattice.

    Parameters:
        flines (list): list of lines in file
        castep (dict): dictionary to update with scraped data

    Returns:
        dict: dictionary updated with scraped data

    """
    for line_no, line in enumerate(flines):
        if 'Real Lattice' in line:
            castep['lattice_cart'] = []
            i = 1
            while True:
                if not flines[line_no + i].strip():
                    break
                else:
                    temp_line = flines[line_no + i].split()[0:3]
                    castep['lattice_cart'].append(list(map(float, temp_line)))
                i += 1
        if 'atom types' not in castep and 'Cell Contents' in line:
            castep['atom_types'] = []
            castep['positions_frac'] = []
            i = 1
            atoms = False
            while True:
                if atoms:
                    if 'xxxxxxxxx' in flines[line_no + i]:
                        atoms = False
                        break
                    else:
                        castep['atom_types'].append(flines[line_no + i].split()[1])
                        castep['positions_frac'].append(list(map(float, (flines[line_no + i].split()[3:6]))))
                if 'x------' in flines[line_no + i]:
                    atoms = True
                i += 1
            for ind, pos in enumerate(castep['positions_frac']):
                for k in range(3):
                    if pos[k] > 1 or pos[k] < 0:
                        castep['positions_frac'][ind][k] %= 1
            castep['num_atoms'] = len(castep['atom_types'])
            castep['stoichiometry'] = get_stoich(castep['atom_types'])
            castep['num_fu'] = castep['num_atoms'] / sum([elem[1] for elem in castep['stoichiometry']])
            break
    else:
        raise CalculationError('Unable to find atoms in CASTEP file.')

    return castep


def _castep_scrape_final_parameters(flines, castep):
    """ Scrape the DFT parameters from a CASTEP file, using those listed
    last in the file (i.e. those used to make the final structure).

    Parameters:
        flines (list): list of lines in the file
        castep (dict): dictionary in which to put scraped data

    Returns:
        dict: dictionary updated with scraped data

    """

    pspot_report_dict = dict()
    for line_no, line in enumerate(reversed(flines)):
        line_no = len(flines) - 1 - line_no
        if 'task' not in castep and 'type of calculation' in line:
            castep['task'] = line.split(':')[-1].strip().replace(" ", "")
        elif 'xc_functional' not in castep and 'functional' in line:
            # convert from .castep file xc_functional to param style
            xc_string = line.split(':')[-1].strip()
            if 'Local Density Approximation' in xc_string:
                castep['xc_functional'] = 'LDA'
            elif 'Perdew Burke Ernzerhof' in xc_string:
                castep['xc_functional'] = 'PBE'
            elif 'PBE for solids' in xc_string:
                castep['xc_functional'] = 'PBESol'
            elif 'hybrid B3LYP' in xc_string:
                castep['xc_functional'] = 'B3LYP'
            elif 'hybrid HSE03' in xc_string:
                castep['xc_functional'] = 'HSE03'
            elif 'hybrid HSE06' in xc_string:
                castep['xc_functional'] = 'HSE06'
        elif 'cut_off_energy' not in castep and 'plane wave basis set' in line:
            castep['cut_off_energy'] = float(line.split(':')[-1].split()[0])
        elif 'finite_basis_corr' not in castep and 'finite basis set correction  ' in line:
            castep['finite_basis_corr'] = line.split(':')[-1].strip()
        elif 'kpoints_mp_grid' not in castep and 'MP grid size for SCF' in line:
            castep['kpoints_mp_grid'] = list(map(int, list(line.split('is')[-1].split())))
        elif 'num_kpoints' not in castep and 'Number of kpoints used' in line:
            castep['num_kpoints'] = int(line.split()[-1])
        elif 'geom_force_tol' not in castep and 'max ionic |force| tolerance' in line:
            castep['geom_force_tol'] = float(line.split()[-2])
        elif 'elec_energy_tol' not in castep and 'total energy / atom convergence tol' in line:
            castep['elec_energy_tol'] = float(line.split()[-2])
        elif 'sedc_apply' not in castep and \
                'DFT+D: Semi-empirical dispersion correction    : on' in line:
            castep['sedc_apply'] = True
            castep['sedc_scheme'] = flines[line_no + 1].split(':')[1].split()[0]
        elif 'space_group' not in castep and 'Space group of crystal' in line:
            castep['space_group'] = line.split(':')[-1].split(',')[0].strip().replace(" ", "")
        elif 'external_pressure' not in castep and 'External pressure/stress' in line:
            try:
                castep['external_pressure'] = []
                castep['external_pressure'].append(list(map(float, flines[line_no + 1].split())))
                castep['external_pressure'].append(list(map(float, flines[line_no + 2].split())))
                castep['external_pressure'].append(list(map(float, flines[line_no + 3].split())))
            except ValueError:
                castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        elif 'spin_polarized' not in castep and 'treating system as spin-polarized' in line:
            castep['spin_polarized'] = True
        elif 'hubbard_u' not in castep and 'Hubbard U values are eV' in line:
            castep['hubbard_u'] = defaultdict(list)
            i = 5
            while i < castep['num_atoms']:
                line = flines[line_no + i].strip()
                atom = line.split()[0].replace('|', '')
                shifts = list(map(float, line.split()[-5:-1]))
                for ind, shift in enumerate(shifts):
                    if shift != 0:
                        if atom not in castep['hubbard_u']:
                            castep['hubbard_u'][atom] = dict()
                        if ind == 0:
                            orbital = 's'
                        elif ind == 1:
                            orbital = 'p'
                        elif ind == 2:
                            orbital = 'd'
                        elif ind == 3:
                            orbital = 'f'
                        castep['hubbard_u'][atom][orbital] = shift
                i += 1
        elif 'Pseudopotential Report' in line:
            if 'species_pot' not in castep:
                castep['species_pot'] = dict()
            i = 0
            while i + line_no < len(flines) - 3:
                if 'Pseudopotential Report' in flines[line_no + i]:
                    i += 2
                    elem = flines[line_no + i].split(':')[1].split()[0]

                elif 'core correction' in flines[line_no + i]:
                    i += 2
                    if not pspot_report_dict.get(elem):
                        castep['species_pot'][elem] = flines[line_no + i].split('"')[1].replace('[]', '')
                        pspot_report_dict[elem] = True
                i += 1
        elif 'species_pot' not in castep and 'Files used for pseudopotentials' in line:
            if 'species_pot' not in castep:
                castep['species_pot'] = dict()
            i = 1
            while True:
                if not flines[line_no + i].strip():
                    break
                else:
                    elem = flines[line_no + i].split()[0].strip()
                    if not pspot_report_dict.get(elem):
                        castep['species_pot'][elem] = flines[line_no + i].split()[1].split('/')[-1]
                        if castep['species_pot'][elem] == 'Pseudopotential':
                            castep['species_pot'][elem] = flines[line_no + i].split()[0].strip()
                            castep['species_pot'][elem] += '_OTF.usp'
                        pspot_report_dict[elem] = False
                i += 1
    # write zero pressure if not found in file
    if 'external_pressure' not in castep:
        castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
    if 'spin_polarized' not in castep:
        castep['spin_polarized'] = False
    return castep


def _castep_scrape_final_structure(flines, castep, db=True):
    """ Scrape final structure from CASTEP file.

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Keyword arguments:
        db (bool): whether to enforce database style, e.g. geometry optimisation only

    Returns:
        dict: dictionary updated with scraped data

    """
    if 'task' in castep and castep['task'].strip() == 'geometryoptimization':
        castep['optimised'] = False
        finish_line, castep['optimised'] = _castep_find_final_structure(flines)
    elif not db:
        finish_line = 0

    final_flines = flines[finish_line + 1:]
    for line_no, line in enumerate(final_flines):
        try:
            if 'Real Lattice' in line:
                castep['lattice_cart'] = []
                i = 1
                while True:
                    if not final_flines[line_no + i].strip():
                        break
                    else:
                        temp_line = final_flines[line_no + i].split()[0:3]
                        castep['lattice_cart'].append(list(map(float, temp_line)))
                    i += 1
            elif 'Lattice parameters' in line:
                castep['lattice_abc'] = []
                i = 1
                castep['lattice_abc'].append(
                    list(map(float,
                             [final_flines[line_no+i].split('=')[1].strip().split(' ')[0],
                              final_flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                              final_flines[line_no+i+2].split('=')[1].strip().split(' ')[0]])))
                castep['lattice_abc'].append(
                    list(map(float,
                             [final_flines[line_no+i].split('=')[-1].strip(),
                              final_flines[line_no+i+1].split('=')[-1].strip(),
                              final_flines[line_no+i+2].split('=')[-1].strip()])))
            elif 'Current cell volume' in line:
                castep['cell_volume'] = float(line.split('=')[1].split()[0].strip())
            elif 'Cell Contents' in line:
                castep['positions_frac'] = []
                i = 1
                atoms = False
                while True:
                    if atoms:
                        if 'xxxxxxxxx' in final_flines[line_no + i]:
                            atoms = False
                            break
                        else:
                            temp_frac = final_flines[line_no + i].split()[3:6]
                            castep['positions_frac'].append(list(map(float, temp_frac)))
                    if 'x------' in final_flines[line_no + i]:
                        atoms = True
                    i += 1
            # don't check if final_energy exists, as this will update for each GO step
            elif 'Final energy, E' in line:
                castep['total_energy'] = float(line.split('=')[1].split()[0])
                castep['total_energy_per_atom'] = castep['total_energy'] / castep['num_atoms']
            elif 'Final free energy' in line:
                castep['free_energy'] = float(line.split('=')[1].split()[0])
                castep['free_energy_per_atom'] = castep['free_energy'] / castep['num_atoms']
            elif '0K energy' in line:
                castep['0K_energy'] = float(line.split('=')[1].split()[0])
                castep['0K_energy_per_atom'] = castep['0K_energy'] / castep['num_atoms']
            elif ' Forces **' in line:
                castep['forces'] = []
                i = 1
                max_force = 0
                forces = False
                while True:
                    if forces:
                        if '*' in final_flines[line_no + i].split()[1]:
                            forces = False
                            break
                        else:
                            force_on_atom = 0
                            castep['forces'].append([])
                            for j in range(3):
                                temp = final_flines[line_no + i].replace('(cons\'d)', '')
                                force_on_atom += float(temp.split()[3 + j])**2
                                castep['forces'][-1].append(float(temp.split()[3 + j]))
                            if force_on_atom > max_force:
                                max_force = force_on_atom
                    elif 'x' in final_flines[line_no + i]:
                        i += 1  # skip next blank line
                        forces = True
                    i += 1
                castep['max_force_on_atom'] = pow(max_force, 0.5)
            elif 'Stress Tensor' in line:
                i = 1
                while i < 20:
                    if 'Cartesian components' in final_flines[line_no + i]:
                        castep['stress'] = []
                        for j in range(3):
                            castep['stress'].append(list(map(float, (final_flines[line_no + i + j + 4].split()[2:5]))))
                    elif 'Pressure' in final_flines[line_no + i]:
                        try:
                            castep['pressure'] = float(final_flines[line_no + i].split()[-2])
                        except ValueError:
                            pass
                        break
                    i += 1
            elif 'Atomic Populations (Mulliken)' in line:
                if castep['spin_polarized']:
                    castep['mulliken_spins'] = []
                    castep['mulliken_net_spin'] = 0.0
                    castep['mulliken_abs_spin'] = 0.0
                castep['mulliken_charges'] = []
                castep['mulliken_spins'] = []
                i = 0
                while i < len(castep['atom_types']):
                    if castep['spin_polarized']:
                        castep['mulliken_charges'].append(float(final_flines[line_no + i + 4].split()[-2]))
                        castep['mulliken_spins'].append(float(final_flines[line_no + i + 4].split()[-1]))
                        castep['mulliken_net_spin'] += castep['mulliken_spins'][-1]
                        castep['mulliken_abs_spin'] += abs(castep['mulliken_spins'][-1])
                    else:
                        castep['mulliken_charges'].append(float(final_flines[line_no + i + 4].split()[-1]))
                    i += 1
            elif 'Final Enthalpy' in line:
                castep['enthalpy'] = float(line.split('=')[-1].split()[0])
                castep['enthalpy_per_atom'] = (float(line.split('=')[-1].split()[0]) / castep['num_atoms'])
            elif 'Final bulk modulus' in line:
                try:
                    castep['bulk_modulus'] = float(line.split('=')[-1].split()[0])
                except ValueError:
                    # the above will fail if bulk modulus was not printed (i.e. if it was unchanged)
                    pass

            elif 'Chemical Shielding and Electric Field Gradient Tensors'.lower() in line.lower():
                i = 5
                castep['chemical_shifts'] = []
                while True:
                    # break when the line containing just '=' is reached
                    if len(flines[line_no + i].split()) == 1:
                        break
                    castep['chemical_shifts'].append(flines[line_no + i].split()[3])
                    i += 1
                if len(castep['chemical_shifts']) != len(castep['atom_types']):
                    raise RuntimeError('Found fewer chemical shifts than atoms (or vice versa)!')

        except Exception as oops:
            msg = 'Error on line {}, contents: {}, error: {}'.format(line_no, line, oops)
            raise RuntimeError(msg)

    # calculate kpoint spacing if not found
    if 'kpoints_mp_grid' in castep and 'kpoints_mp_spacing' not in castep and 'lattice_cart' in castep:
        castep['kpoints_mp_spacing'] = calc_mp_spacing(castep['lattice_cart'], castep['kpoints_mp_grid'], prec=4)
    return castep


def _castep_scrape_metadata(flines, castep):
    """ Scrape metadata from CASTEP file.

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Returns:
        dict: dictionary updated with scraped data

    """
    # computing metadata, i.e. parallelism, time, memory, version
    for line in flines:
        if 'Release CASTEP version' in line:
            castep['castep_version'] = line.replace('|', '').split()[-1]
        try:
            if 'Run started:' in line:
                year = line.split()[5]
                month = str(strptime(line.split()[4], '%b').tm_mon)
                day = line.split()[3]
                castep['date'] = day + '-' + month + '-' + year
            elif 'Total time' in line and 'matrix elements' not in line:
                castep['total_time_hrs'] = float(line.split()[-2]) / 3600
            elif 'Peak Memory Use' in line:
                castep['peak_mem_MB'] = int(float(line.split()[-2]) / 1024)
            elif 'total storage required per process' in line:
                castep['estimated_mem_MB'] = float(line.split()[-5])

        # if any of these error, don't worry too much (at all)
        except Exception:
            pass

    return castep


def _castep_find_final_structure(flines):
    """ Search for info on final structure in .castep file.

    Parameters:
        flines (list): list of lines in file.

    Returns:
        int: line number in file where total energy of final structure is printed.

    """
    optimised = False
    finish_line = 0
    success_string = 'Geometry optimization completed successfully'
    failure_string = 'Geometry optimization failed to converge after'
    # look for final "success/failure" string in file for geometry optimisation
    for line_no, line in enumerate(reversed(flines)):
        if success_string in line:
            finish_line = len(flines) - line_no
            optimised = True
            break
        elif failure_string in line:
            finish_line = len(flines) - line_no
            optimised = False
            break

    # now wind back to get final total energies and non-symmetrised forces
    for count, line in enumerate(reversed(flines[:finish_line])):
        if 'Final energy, E' in line:
            finish_line -= count + 2
            break

    return finish_line, optimised


def _castep_scrape_all_snapshots(flines):
    """ Scrape all intermediate structures from a CASTEP file, both
    geometry optimisation snapshots, and repeated SCF calculations.

    Parameters:
        flines (list): list of lines of file.

    Returns:
        :obj:`list` of :obj:`dict`: list of dictionaries containing
            intermediate snapshots.
        int: number of completed geometry optimisation steps.

    """
    intermediates = []
    num_opt_steps = 0
    snapshot = dict()
    for line_no, line in enumerate(flines):
        try:
            # use the "Real Lattice" line as the start of a new snapshot / end of old one
            if 'Real Lattice' in line:
                # add the last snapshot only if it isn't a repeat
                if 'total_energy' in snapshot:
                    # if positions frac didn't change (and thus weren't printed, use the last value)
                    if 'positions_frac' not in snapshot:
                        snapshot['positions_frac'] = intermediates[-1]['positions_frac']
                        snapshot['atom_types'] = intermediates[-1]['atom_types']
                        snapshot['num_atoms'] = len(snapshot['positions_frac'])
                    snapshot['free_energy_per_atom'] = snapshot['free_energy'] / snapshot['num_atoms']
                    snapshot['total_energy_per_atom'] = snapshot['total_energy'] / snapshot['num_atoms']
                    snapshot['0K_energy_per_atom'] = snapshot['0K_energy'] / snapshot['num_atoms']
                    # handle single atom forces edge-case
                    if snapshot['num_atoms'] == 1:
                        snapshot['forces'] = [[0, 0, 0]]
                    if not intermediates:
                        intermediates.append(snapshot)
                    elif (snapshot['lattice_cart'] != intermediates[-1]['lattice_cart'] and
                          snapshot['total_energy'] != intermediates[-1]['total_energy']):
                        intermediates.append(snapshot)

                snapshot = dict()
                snapshot['lattice_cart'] = []
                i = 1
                while True:
                    if not flines[line_no + i].strip():
                        break
                    else:
                        temp_line = flines[line_no + i].split()[0:3]
                        snapshot['lattice_cart'].append(list(map(float, temp_line)))
                    i += 1

            elif 'Lattice parameters' in line:
                snapshot['lattice_abc'] = []
                i = 1
                snapshot['lattice_abc'].append(
                    list(map(float,
                             [flines[line_no+i].split('=')[1].strip().split(' ')[0],
                              flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                              flines[line_no+i+2].split('=')[1].strip().split(' ')[0]])))
                snapshot['lattice_abc'].append(
                    list(map(float,
                             [flines[line_no+i].split('=')[-1].strip(),
                              flines[line_no+i+1].split('=')[-1].strip(),
                              flines[line_no+i+2].split('=')[-1].strip()])))

            elif 'Current cell volume' in line:
                snapshot['cell_volume'] = float(line.split('=')[1].split()[0].strip())
            elif 'Cell Contents' in line:
                snapshot['positions_frac'] = []
                snapshot['atom_types'] = []
                i = 1
                atoms = False
                while True:
                    if atoms:
                        if 'xxxxxxxxx' in flines[line_no + i]:
                            atoms = False
                            break
                        else:
                            temp_frac = flines[line_no + i].split()[3:6]
                            snapshot['positions_frac'].append(list(map(float, temp_frac)))
                            snapshot['atom_types'].append(flines[line_no + i].split()[1])
                    if 'x------' in flines[line_no + i]:
                        atoms = True
                    i += 1
                for ind, pos in enumerate(snapshot['positions_frac']):
                    for k in range(3):
                        if pos[k] > 1 or pos[k] < 0:
                            snapshot['positions_frac'][ind][k] %= 1

                snapshot['num_atoms'] = len(snapshot['positions_frac'])
                snapshot['stoichiometry'] = get_stoich(snapshot['atom_types'])

            # don't check if final_energy exists, as this will update for each GO step
            elif 'Final energy, E' in line:
                snapshot['total_energy'] = float(line.split('=')[1].split()[0])
            elif 'Final free energy' in line:
                snapshot['free_energy'] = float(line.split('=')[1].split()[0])
            elif '0K energy' in line:
                snapshot['0K_energy'] = float(line.split('=')[1].split()[0])
            elif ' Forces **' in line:
                snapshot['forces'] = []
                i = 1
                max_force = 0
                forces = False
                while True:
                    if forces:
                        if '*' in flines[line_no + i].split()[1]:
                            forces = False
                            break
                        else:
                            force_on_atom = 0
                            snapshot['forces'].append([])
                            for j in range(3):
                                temp = flines[line_no + i].replace('(cons\'d)', '')
                                force_on_atom += float(temp.split()[3 + j])**2
                                snapshot['forces'][-1].append(float(temp.split()[3 + j]))
                            if force_on_atom > max_force:
                                max_force = force_on_atom
                    elif 'x' in flines[line_no + i]:
                        i += 1  # skip next blank line
                        forces = True
                    i += 1
                snapshot['max_force_on_atom'] = pow(max_force, 0.5)
            elif 'Stress Tensor' in line:
                i = 1
                while i < 20:
                    if 'Cartesian components' in flines[line_no + i]:
                        snapshot['stress'] = []
                        for j in range(3):
                            snapshot['stress'].append(list(map(float, (flines[line_no + i + j + 4].split()[2:5]))))
                    elif 'Pressure' in flines[line_no + i]:
                        snapshot['pressure'] = float(flines[line_no + i].split()[-2])
                        break
                    i += 1

            # use only finished iterations for counting number of complete GO steps
            elif ': finished iteration' in line and 'with enthalpy' in line:
                # don't include the "zeroth" step before anything has been moved
                if '0' not in line.split():
                    num_opt_steps += 1

        except Exception as exc:
            msg = 'Error on line {}, contents: {}, error: {}"'.format(line_no, line, exc)
            raise RuntimeError(msg)

    return intermediates, num_opt_steps
