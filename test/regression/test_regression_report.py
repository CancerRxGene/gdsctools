from gdsctools import *
import pandas as pd
import pylab
import json

def test_regression_report():
    IC = gdsctools_data("IC50_v5.csv.gz")
    GF = gdsctools_data("genomic_features_v5.csv.gz")

    PREFIX = "gdsctools_regression_"
    IMAGE_DIR = "images"
    DATA_PREFIX = "data/" + PREFIX
    IMAGE_PREFIX = "images/" + PREFIX

    gd = regression.GDSCLasso(IC, GF)
    DRUGIDS = gd.drugIds[0:4]

    config = {"boxplot_n":5, "randomness":5}

    # Get best model
    inputs = []
    for drugid in DRUGIDS:
        res = gd.runCV(drugid, verbose=False, kfolds=10)
        bestmodel = gd.get_model(alpha=res.alpha)

        def _pngname(tag):
            return IMAGE_PREFIX + "%s_%s.png" % (tag, drugid)

        # Plot weights
        weights = gd.plot_weight(drugid, bestmodel)
        if len(weights):
            pylab.savefig(_pngname("weights"))
            pylab.close()

        weights = pd.DataFrame({
                "weigths": res.coefficients,
                "features":gd.feature_names})
        output = DATA_PREFIX + "weights_{}.csv".format(drugid)
        weights.to_csv(output, index=False)

        # Plot importance
        weights = gd.plot_importance(drugid, bestmodel)
        if len(weights):
            pylab.savefig(_pngname("importance"))
            pylab.close()

        # Boxplots
        boxres = gd.boxplot(drugid, model=bestmodel, n=5,
            bx_vert=False)
        if len(boxres['data']):
            pylab.savefig(_pngname("boxplot"))
            pylab.close()

        # Bayes factor
        ran = gd.check_randomness(drugid, 10, 10)
        pylab.savefig(_pngname("randomness"))
        pylab.close()
        results = {"drugid": int(drugid),
                "Rp":res.Rp,
                "alpha": res.alpha,
                "ln_alpha": res.ln_alpha,
                "ttest": ran['ttest_pval'],
                "bayes":ran['bayes_factor']}

        output = DATA_PREFIX + "results_{}.json".format(drugid)
        fh = open(output, "w")
        json.dump(results, fh)
        fh.close()
        inputs.append(output)

    # gather all results:
    data = []
    for this in inputs:
        with open(this, "r") as fh:
            data.append(json.loads(fh.read()))
    df = pd.DataFrame(data)
    df.set_index("drugid", inplace=True)
    df.to_csv(DATA_PREFIX + "results.csv")


    inputs = [DATA_PREFIX + "weights_{}.csv".format(drugid) for drugid in DRUGIDS]
    df = pd.concat(
       [pd.read_csv(this).set_index("features") for this in inputs],
            axis=1)
    df.columns = DRUGIDS
    df.to_csv(DATA_PREFIX + "weights.csv")


    from gdsctools import regression_report
    report = regression_report.RegressionReport("lasso", image_dir=IMAGE_DIR, config=config)
    report.create_html_main()
    report.create_html_drug()



