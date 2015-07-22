"""riboSeqR Galaxy unit tests"""
import unittest
from riboseqr import utils


class PrepareTestCase(unittest.TestCase):

    def test_process_args(self):
        """Test processing arguments. """
        # "ATG" -> c("ATG")
        rs = utils.process_args('ATG', ret_mode='charvector')
        self.assertEqual(rs, 'c("ATG")','Return string as a character vector.')

        # stop codons "TAG,TAA,TGA" -> c("TAG", "TAA", "TGA"). Also
        # replicate names, seqnames.
        rs = utils.process_args('TAG,TAA,TGA', ret_mode='charvector')
        self.assertEqual(
            rs, "c('TAG', 'TAA', 'TGA')",
            'Return comma separated strings as a character vector.')

        # "" -> None
        rs = utils.process_args('')
        self.assertIsNone(rs, 'Return empty string as None.')

        # "27,28" -> c(27, 28)
        rs = utils.process_args("27,28", ret_type='int', ret_mode='charvector')
        self.assertEqual(
            rs, "c(27, 28)", 'Return number strings as a character vector.')

        # "27,28" -> [27, 28]
        rs = utils.process_args("27,28", ret_type='int', ret_mode='list')
        self.assertEqual(rs, [27, 28], 'Return number strings as a list.')

        # "0,2" -> list(0,2)
        rs = utils.process_args("0,2", ret_type='int', ret_mode='listvector')
        self.assertEqual(
            rs, "list(0, 2)", 'Return number strings as a list vector.')

        # "50" -> 50
        rs = utils.process_args("50", ret_type='int')
        self.assertEqual(rs, 50, 'Return number string as a number.')

        # "-200" -> -200
        rs = utils.process_args("-200", ret_type='int')
        self.assertEqual(rs, -200, 'Return number string as a number.')

        # "TRUE" -> TRUE
        rs = utils.process_args("TRUE", ret_type='bool')
        self.assertEqual(rs, 'TRUE', 'Return bool string as bool.')

        # 'chlamy17,chlamy3' -> 'chlamy17,chlamy3' for ribo and rna names
        rs = utils.process_args('chlamy17,chlamy3')
        self.assertEqual(rs, 'chlamy17,chlamy3', 'Return csv string as string.')

        # 'chlamy17.idx, chlamy3.idx' -> ['chlamy17.idx', 'chlamy3.idx']
        rs = utils.process_args('chlamy17.idx, chlamy3.idx', ret_mode='list')
        self.assertEqual(rs, ['chlamy17.idx', 'chlamy3.idx'],
                         'Return files as a list.')
