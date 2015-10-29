from gdsctools.datasets import Data


def test_data():

    d = Data()
    d.description= 'toto'
    d.authors = 'toto'
    d.filename = 'here'
    print(d)

