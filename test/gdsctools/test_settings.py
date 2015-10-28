from gdsctools.anova import ANOVASettings


def test_settings():
    s = ANOVASettings()
    s.check()
    s.FDR_threshold = -1
    try:
        s.check()
        assert False
    except ValueError:
        assert True

