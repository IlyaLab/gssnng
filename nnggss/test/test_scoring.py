
import numpy as np
import pandas as pd

from nnggss.scoring import _ms_sing

DF = pd.DataFrame({
        'g1': [2, 1, 1, 0, 0],
        'g2': [1, 1, 0, 0, 0],
        'g3': [0, 0, 0, 1, 0],
        'g4': [0, 0, 0, 1, 0],
        }
    ).T



def test_mssingscore():
    gset = ['g1', 'g2']

    # note that we would usually remove the non-expressed genes.
    # but its fine for testing
    expected = {0: 0.375, 1: 0.25, 2: 0.125, 3: -0.25, 4: -0.25}
    for i in range(len(DF.columns)):
        print(_ms_sing(gset, DF[i], norm_method='standard', rankup=True)['total_score'])
        assert _ms_sing(gset, DF[i], norm_method='standard', rankup=True)['total_score'] == expected[i]

test_mssingscore()
print('test done')

