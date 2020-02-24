from unittest import TestCase, main
from click.testing import CliRunner
from shutil import rmtree
from os.path import join
from tempfile import mkdtemp
from skbio.util import get_data_path
from skbio import DistanceMatrix
from skbio.stats.distance import DissimilarityMatrixError

from microprot.scripts.calculate_Neff import (parse_msa_file,
                                              hamming_distance_matrix,
                                              cluster_sequences,
                                              effective_family_size,
                                              _calculate_Neff)


class ProcessingTests(TestCase):
    def setUp(self):
        # temporary working directory
        self.working_dir = mkdtemp()
        self.output_fp = join(self.working_dir, 'output.txt')

        # test data files
        dir = 'test_calculate_Neff'
        self.input_a3m_fp = get_data_path(join(dir, '2phyA.a3m'))
        self.input_a3m_fp_illegal = get_data_path(join(dir, 'illegal.a3m'))
        self.input_single_a3m_fp = get_data_path(join(dir, 'single.a3m'))
        self.plain_msa_fp = get_data_path(join(dir, '2phyA.aln'))
        self.plain_msa_fp_illegal = get_data_path(join(dir, 'illegal.aln'))
        self.hamming_dm_fp = get_data_path(join(dir, '2phyA.hdm'))
        self.clusters_fp = get_data_path(join(dir, '2phyA.c80'))

    def test_parse_msa_file(self):
        obs = [str(x) for x in parse_msa_file(self.input_a3m_fp)]
        with open(self.plain_msa_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

    def test_parse_msa_file_check_illegal(self):
        obs = [str(x) for x in parse_msa_file(self.input_a3m_fp_illegal)]
        with open(self.plain_msa_fp_illegal, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

    def test_hamming_distance_matrix(self):
        msa = parse_msa_file(self.input_a3m_fp)
        obs = hamming_distance_matrix(msa)
        exp = DistanceMatrix.read(self.hamming_dm_fp)
        exp = exp.condensed_form()
        self.assertTrue((obs == exp).all())


    def test_cluster_sequences(self):
        msa = parse_msa_file(self.input_a3m_fp)
        hdm = hamming_distance_matrix(msa)
        obs = cluster_sequences(hdm, 80)
        with open(self.clusters_fp, 'r') as f:
            exp = [int(x.split('\t')[1]) for x in f.read().splitlines()]
        self.assertListEqual(obs, exp)
        obs = cluster_sequences(DistanceMatrix([[0]]), 80)
        self.assertListEqual(obs, [])

    def test_effective_family_size(self):
        msa = parse_msa_file(self.input_a3m_fp)
        hdm = hamming_distance_matrix(msa)
        clu = cluster_sequences(hdm, 80)
        obs = effective_family_size(clu, msa.shape[1])
        exp = 6.6187612134
        self.assertAlmostEqual(obs, exp)
        obs = effective_family_size(clu, msa.shape[1], Nclu=True)
        exp = 74
        self.assertEqual(obs, exp)

    def test__calculate_Neff(self):
        params = ['--infile', self.input_a3m_fp,
                  '--outfile', self.output_fp,
                  '--cutoff', 80]
        res = CliRunner().invoke(_calculate_Neff, params)
        self.assertIn('Effective family size at 80% identity: 6.619.',
                      res.output)
        with open(self.output_fp, 'r') as f:
            obs = float(f.readline().strip())
        exp = 6.6187612134
        self.assertAlmostEqual(obs, exp)
        self.assertEqual(res.exit_code, 0)
        self.assertIn('Task completed.', res.output)
        params = ['--infile', self.input_single_a3m_fp]
        res = CliRunner().invoke(_calculate_Neff, params)
        self.assertIn('Effective family size at 100% identity: 0.000.',
                      res.output)

    def tearDown(self):
        rmtree(self.working_dir)

    def test_duplicateIDs(self):
        # test that alignment is parsed correctly
        infile = get_data_path('test_calculate_Neff/GRAMNEG_T1D_899_33-87.a3m')
        msa = parse_msa_file(infile)
        self.assertEqual(msa.shape.sequence, 524)
        self.assertEqual(msa.shape.position, 55)

        # test that distance matrix has as many elements as sequences in msa
        # regardless of the fact that some sequence IDs collapse
        # Note: hdm contains elements in the upper triangle
        # of the full distance matrix
        hdm = hamming_distance_matrix(msa)
        dim = msa.shape.sequence
        exp_n = (dim ** 2 - dim) // 2
        self.assertEqual(len(hdm), exp_n)

        clu = cluster_sequences(hdm, 80)
        self.assertEqual(len(clu), msa.shape.sequence)

        Neff = effective_family_size(clu, msa.shape[1])
        self.assertAlmostEqual(Neff, 49.62, places=2)



if __name__ == '__main__':
    main()
