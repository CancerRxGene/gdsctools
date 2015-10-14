from gdsctools.cohens import cohens
from nose.tools import assert_almost_equal

def test_cohens():


    x = [1,2,3,1,2]
    y = [1,2,3,4,5,6,7,8]
    assert_almost_equal(cohens(x, y), 1.3378921)
    assert cohens(x, x) == 0
