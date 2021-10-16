"""Heatmap and dendograms"""
import warnings
import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance
import easydev


__all__ = ['Linkage']


class Linkage:
    """Linkage used in other tools such as Heatmap"""

    methods = ["single", "complete", "average", "weighted", "centroid",
        "median", "ward"]
    metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching',
        'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
        'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']

    def __init__(self):
        """.. rubric:: constructor

        :param data: a dataframe or possibly a numpy matrix.

        """
        pass

    def check_metric(self, value):
        easydev.check_param_in_list(value, self.metrics)

    def check_method(self, value):
        easydev.check_param_in_list(str(value), self.methods)

    def linkage(self, df, method, metric):
        self.check_metric(metric)
        self.check_method(method)
        d = distance.pdist(df)
        D = distance.squareform(d)

        # hierarchy.ClusterWarning
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            Y = hierarchy.linkage(D, method=method, metric=metric)
            return Y
