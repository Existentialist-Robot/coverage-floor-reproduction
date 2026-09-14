"""Statistical invariants for the retained-data reanalysis."""
import unittest
import numpy as np
from review_reanalysis import adjust, grouped, paired_test


class ReanalysisTests(unittest.TestCase):
    def test_noise_copies_do_not_become_independent_groups(self):
        rows = [{'c':.5,'arm':'mnar','rep':i,'scen':'control',
                 'lead_fired':bool(i),'lag_fired':False} for i in range(2)]
        keys, values, _ = grouped(rows, False)
        copied_keys, copied_values, _ = grouped(rows*9, False)
        self.assertEqual(keys, copied_keys)
        np.testing.assert_array_equal(values, copied_values)
        self.assertEqual(len(keys), 2)

    def test_correction_preserves_original_order(self):
        tests = [{'p':x} for x in [.04,.001,.02,.2]]
        result = adjust(tests)
        np.testing.assert_allclose([t['holm_p'] for t in tests], [.08,.004,.06,.2])
        np.testing.assert_allclose([t['bh_p'] for t in tests], [.0533333333333,.004,.04,.2])
        self.assertEqual(result['holm_survivors'], 1)

    def test_paired_test_sign_and_null(self):
        self.assertEqual(paired_test([0]*10)['p'], 1)
        x = [-.2,.1,.4,.6]
        a,b = paired_test(x), paired_test([-v for v in x])
        self.assertAlmostEqual(a['p'],b['p'])
        self.assertAlmostEqual(a['difference'],-b['difference'])


if __name__ == '__main__':
    unittest.main()
