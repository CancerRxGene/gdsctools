from gdsctools.glass import glass


def test_glass():

    x = [1,2,3,1,2]
    y = [1,2,3,4,5,6,7,8]

    g1, g2 = glass(x,y)

    assert g1 == 3.2271172452028631
    assert g2 == 1.1022703842524304
    g1, g2 = glass(x,x)

    assert g1 == g2 == 0
