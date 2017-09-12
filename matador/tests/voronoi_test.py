#!/usr/bin/env python
import unittest
import numpy as np
from os.path import realpath
from json import load
from matador.similarity.voronoi_similarity import get_unique_sites, collect_unique_sites
from matador.similarity.voronoi_similarity import are_sites_the_same, set_substruc_dict
from matador.similarity.voronoi_similarity import set_site_array, create_site_array


REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
plot = False


class VoronoiSimilarityTest(unittest.TestCase):
    """ Test Voronoi similarity functionality. """
    # def testVoronoiHull(self):
        # with open(REAL_PATH + 'data/voronoi_hull.json', 'r') as f:
            # cursor = load(f)
        # self.assertEqual(len(cursor), 34)
        # site_count = {'K': 0, 'P': 0}
        # for ind, doc in enumerate(cursor):
            # set_site_array(doc)

        # unique_environments = dict()
        # for elem in ['K', 'P']:
            # unique_environments[elem] = []
            # for doc in cursor:
                # get_unique_sites(doc)
                # if elem in doc['unique_substrucs']:
                    # for site in doc['unique_substrucs'][elem]:
                        # unique_environments[elem].append(site)
                # else:
                    # print(doc['stoichiometry'])
        # print(unique_environments['P'][0])
        # unique_environments, max_num_elems = create_site_array(unique_environments)
        # print(max_num_elems)
        # same = 0
        # unique = 0
        # for elem in unique_environments:
            # for i, site in enumerate(unique_environments[elem]):
                # for j in range(i, len(unique_environments[elem])):
                    # test = are_sites_the_same(unique_environments[i], unique_environments[j])
                    # print(i, j, test)


            # get_unique_sites(cursor[ind])
            # for elem in cursor[ind]['unique_sites']:
                # site_count[elem] += len(cursor[ind]['unique_sites'][elem])

        # unique_sites = collect_unique_sites(cursor)
        # self.assertIn('P', unique_sites)
        # self.assertIn('K', unique_sites)
        # for elem in unique_sites:
            # self.assertEqual(len(unique_sites[elem]), site_count[elem])

    def testOnIndividualDocs(self):
        with open(REAL_PATH + 'data/voronoi_hull.json', 'r') as f:
            cursor = load(f)
        self.assertEqual(len(cursor), 34)
        for ind, doc in enumerate(cursor):
            set_substruc_dict(cursor[ind])
        # doc = cursor[15]

        for doc in cursor:
            # print(40*'-!')
            # print(doc['stoichiometry'])

            set_site_array(doc, normalise=True)

            for elem in set(doc['atom_types']):
                self.assertTrue(elem in doc['site_array'])
                self.assertTrue(len(doc['site_array'][elem]) == sum([1 for atom in doc['atom_types'] if atom == elem]))
                self.assertTrue(isinstance(doc['site_array'][elem], list))

            for elem in doc['site_array']:
                for _elem in doc['site_array'][elem][0]:
                    shape = np.shape(doc['site_array'][elem][0][_elem])
                    for site in doc['site_array'][elem]:
                        self.assertEqual(shape, np.shape(site[_elem]))

            # check sites are properly normalised
            for elem in doc['site_array']:
                for i in range(len(doc['site_array'][elem])):
                    _sum = 0
                    for _elem in doc['site_array'][elem][0]:
                        _sum += np.sum(doc['site_array'][elem][i][_elem])
                    self.assertAlmostEqual(_sum, 1)

            rtol = 1e-3
            atol = 1e-2
            get_unique_sites(doc, rtol=rtol, atol=atol)

            # check degeneracies sum to total number of atoms
            for elem in set(doc['atom_types']):
                self.assertEqual(sum(doc['site_degeneracies'][elem]), doc['atom_types'].count(elem))

            # test stddev doesn't exceed 2*atol, and that every site is within 2*atol of mean
            for elem in doc['similar_sites']:
                for ind, site in enumerate(doc['similar_sites'][elem]):
                        for _elem in doc['unique_site_stddev'][elem][ind]:
                            try:
                                max_std = np.max(doc['unique_site_stddev'][elem][ind][_elem])
                                if isinstance(max_std, np.float64):
                                    self.assertTrue(max_std < 2*atol)
                            except(ValueError):
                                pass
                        for index in site:
                            self.assertTrue(are_sites_the_same(doc['site_array'][elem][index], doc['unique_site_array'][elem][ind], atol=2*atol, rtol=rtol))
            if plot:
                plot_doc_strucs(doc)

    def testSBCCLi(self):
        from matador.scrapers.castep_scrapers import res2dict
        rtol = 1e-3
        atol = 2e-2
        doc, s = res2dict(REAL_PATH + 'data/hull-LiP-mdm_chem_mater/Li-bcc.res')
        set_site_array(doc)
        get_unique_sites(doc, rtol=rtol, atol=atol)
        self.assertTrue(len(doc['unique_site_inds']['Li']) == 1)
        self.assertTrue(len(doc['unique_site_array']['Li'][0]['Li']) == 14)
        # print(np.shape(np.where(np.isclose(doc['unique_site_array']['Li'][0]['Li'], np.max(doc['unique_site_array']['Li'][0]['Li'])))))
        # print(np.shape(np.where(np.isclose(doc['unique_site_array']['Li'][0]['Li'], np.min(doc['unique_site_array']['Li'][0]['Li'])))))
        self.assertTrue(np.shape(np.where(np.isclose(doc['unique_site_array']['Li'][0]['Li'], np.max(doc['unique_site_array']['Li'][0]['Li']))))[-1] == 8)
        self.assertTrue(np.shape(np.where(np.isclose(doc['unique_site_array']['Li'][0]['Li'], np.min(doc['unique_site_array']['Li'][0]['Li']))))[-1] == 6)

    def testZintlLiP7(self):
        from matador.scrapers.castep_scrapers import res2dict
        rtol = 1e-3
        atol = 2e-2
        doc, s = res2dict(REAL_PATH + 'data/hull-LiP-mdm_chem_mater/LiP-ColCode23621.res')
        set_site_array(doc)
        get_unique_sites(doc, rtol=rtol, atol=atol)
        if plot:
            plot_doc_strucs(doc)
        self.assertTrue(len(doc['unique_site_inds']['Li']) == 1)

    def testLi3P_known_phases(self):
        from matador.scrapers.castep_scrapers import res2dict
        rtol = 1e-3
        atol = 2e-1
        doc, s = res2dict(REAL_PATH + 'data/hull-LiP-mdm_searches/LiP-CollCode26880.res')
        set_site_array(doc)
        get_unique_sites(doc, rtol=rtol, atol=atol)
        if plot:
            plot_doc_strucs(doc)
        doc, s = res2dict(REAL_PATH + 'data/hull-LiP-mdm_searches/LiP-CollCode81565.res')
        set_site_array(doc)
        get_unique_sites(doc, rtol=rtol, atol=atol)
        if plot:
            plot_doc_strucs(doc)

    def testAl5Y3O12Garnet(self):
        from matador.scrapers.castep_scrapers import cell2dict
        rtol = 5e-2
        atol = 1e-2
        doc, s = cell2dict(REAL_PATH + 'data/Al5Y3O12.cell', db=False, outcell=True, positions=True)
        set_site_array(doc)
        get_unique_sites(doc, atol=atol, rtol=rtol)
        # print(doc['similar_sites'])
        if plot:
            plot_doc_strucs(doc)
        self.assertTrue(len(doc['unique_site_inds']['Y']) == 1)
        self.assertTrue(len(doc['unique_site_inds']['Al']) == 2)
        self.assertTrue(len(doc['unique_site_inds']['O']) == 1)

    def testPairwiseSiteUniqueness(self):
        with open(REAL_PATH + 'data/voronoi_hull.json', 'r') as f:
            cursor = load(f)
        self.assertEqual(len(cursor), 34)
        for ind, doc in enumerate(cursor):
            set_substruc_dict(cursor[ind])
        doc = cursor[15]
        set_site_array(doc)

        site_examples = [{'P': np.array([1, 1, 1, 0.2]),
                          'K': np.array([0.3, 0.2, 0])},
                         {'P': np.array([1, 1, 1, 0.2]),
                          'K': np.array([0.3, 0.2])},
                         {'P': np.array([1, 1, 1, 0.2]),
                          'K': np.array([1e-2, 1e-2])},
                         {'P': np.array([1.1, 1, 1, 0.2]),
                          'K': np.array([1e-2, 1e-2])},
                         {'P': np.array([1, 1, 1, 0.2])},
                         {'Sn': np.array([1.1, 1, 1, 0.2]),
                          'Ba': np.array([1e-2, 1e-2])},
                         {'P': np.array([1.1, 1, 1, 0.2]),
                          'Sn': np.array([10, 10]),
                          'K': np.array([1e-2, 1e-2])},
                         {'P': np.array([1.1, 1, 1, 0.2]),
                          'Sn': np.array([0, 0]),
                          'K': np.array([1e-2, 1e-2])},
                         doc['site_array']['P'][0],
                         doc['site_array']['P'][1],
                         ]

        self.assertTrue(are_sites_the_same(site_examples[0], site_examples[0]))

        self.assertTrue(are_sites_the_same(site_examples[0], site_examples[1]))
        self.assertTrue(are_sites_the_same(site_examples[1], site_examples[0]))

        self.assertTrue(are_sites_the_same(site_examples[2], site_examples[4]))
        self.assertTrue(are_sites_the_same(site_examples[4], site_examples[2]))

        self.assertFalse(are_sites_the_same(site_examples[2], site_examples[3]))
        self.assertFalse(are_sites_the_same(site_examples[3], site_examples[2]))

        self.assertFalse(are_sites_the_same(site_examples[2], site_examples[5]))
        self.assertFalse(are_sites_the_same(site_examples[5], site_examples[2]))

        self.assertFalse(are_sites_the_same(site_examples[5], site_examples[6]))
        self.assertFalse(are_sites_the_same(site_examples[6], site_examples[5]))

        self.assertTrue(are_sites_the_same(site_examples[3], site_examples[7]))
        self.assertTrue(are_sites_the_same(site_examples[7], site_examples[3]))

        self.assertTrue(are_sites_the_same(site_examples[8], site_examples[9]))
        self.assertTrue(are_sites_the_same(site_examples[9], site_examples[8]))

        test_envs = [{'K': np.array([0.15656087, 0.10824254, 0.09639981, 0.07271435, 0.06963524, 0.00000000, 0.00000000]),
                      'P': np.array([0.23685457, 0.23661772, 0.01373757, 0.00923733, 0.00000000, 0.00000000])},
                     {'K': np.array([0.14778908, 0.10333412, 0.09553086, 0.08347127, 0.06124379, 0.00733034, 0.00118231]),
                      'P': np.array([0.23646252, 0.23457082, 0.02057224, 0.00614803, 0.00118231, 0.00118231])}]

        self.assertTrue(are_sites_the_same(test_envs[0], test_envs[1]))
        self.assertTrue(are_sites_the_same(test_envs[1], test_envs[0]))
        self.assertEqual(are_sites_the_same(test_envs[1], test_envs[0]), are_sites_the_same(test_envs[0], test_envs[1]))


def plot_doc_strucs(doc):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('whitegrid')
    colours = sns.color_palette("Dark2", 8)
    fig, ax_list = plt.subplots(len(doc['similar_sites']), len(doc['similar_sites']), figsize=(10, 5))
    for jind, elem in enumerate(doc['similar_sites']):
        print(elem)
        for ind, site in enumerate(doc['similar_sites'][elem]):
            # for cind, _elem in enumerate(doc['unique_site_array'][elem][ind]):
                # if len(doc['unique_site_array'][elem][ind][_elem]) > 0:
                    # ax.errorbar(range(len(doc['unique_site_array'][elem][ind][_elem])), doc['unique_site_array'][elem][ind][_elem],
                                # yerr=doc['unique_site_stddev'][elem][ind][_elem], ls='--', zorder=10000, lw=0, label=_elem, c=colours[cind])
            for cind, _elem in enumerate(doc['unique_site_array'][elem][ind]):
                ax = ax_list[cind, jind]
                ax.set_ylim(0, 1)
                ax.set_title('-'.join([elem, _elem]))
                if len(doc['site_array'][elem][ind][_elem]) > 1:
                    sns.tsplot(ax=ax, ci='sd', data=np.asarray([doc['site_array'][elem][index][_elem] for index in site]), color=colours[ind])
                # for index in site:
                    # if len(doc['unique_site_array'][elem][ind][_elem]) > 0:
                        # ax.plot(range(len(doc['site_array'][elem][index][_elem])), doc['site_array'][elem][index][_elem], c=colours[ind], lw=0.5, alpha=0.1)
    plt.show()


if __name__ == '__main__':
    unittest.main()
