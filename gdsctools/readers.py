# -*- python -*-
# -*- coding utf-8 -*-

#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""IO functionalities


Provides readers to read

- Matrix of IC50 data set :class:`IC50`
- Matrix of Genomic features with :class:`GenomicFeatures`
- Drug Decoder table with :class:`DrugDecode`

"""
import warnings

import pandas as pd
import pylab
import numpy as np
import easydev


__all__ = ['IC50', 'GenomicFeatures', 'Reader', 'DrugDecode']


class Reader(object):
    """Convenience base class to read CSV or TSV files (using extension)"""
    def __init__(self, data=None):
        r""".. rubric:: Constructor

        This class takes only one input parameter, however, it may be a
        filename, or a dataframe or an instance of :class:`Reader` itself. This
        means than children classes such as :class:`IC50` can also be used
        as input as long as a dataframe named :attr:`df` can be found.

        :param data: a filename in CSV or TSV format with format specified by
            child class (see e.g. :class:`IC50`), or a valid dataframe, or an
            instance of :class:`Reader`.

        The input can be a filename either in CSV (comma separated values) or
        TSV (tabular separated values). The extension will be used to interpret
        the content, so please be consistent in the naming of the file
        extensions.

        ::

            >>> from gdsctools import Reader, ic50_test
            >>> r = Reader(ic50_test.filename) # this is a CSV file
            >>> len(r.df)   # number of rows
            988
            >>> len(r)      # number of elements
            11856

        Note that :class:`Reader` is a base class and more sophisticated
        readers are available. for example, the :class:`IC50` would be
        better to read this IC50 data set.

        The data has been stored in a data frame in the :attr:`df` attribute.

        The dataframe of the object itself can be used as an input to create
        an new instance::

            >>> from gdsctools import Reader, ic50_test
            >>> r = Reader(ic50_test.filename, sep="\t")
            >>> r2 = Reader(r) # here r.df is simply copied into r2
            >>> r == r2
            True

        It is sometimes convenient to create an empty Reader that will be
        populated later on::

            >>> r = Reader()
            >>> len(r)
            0

        More advanced readers (e.g. :class:`IC50`) can also be used as input
        as long as they have a :attr:`df` attribute::

            >>> from gdsctools import Reader, ic50_test
            >>> ic = IC50(ic50_test)
            >>> r = Reader(ic)


        """
        # input data
        if data is None:
            # create an empty dataframe
            self.df = pd.DataFrame()
            self._filename = None
        elif isinstance(data, str):
            # Read a filename in TSV or CSV format
            self.read_data(data)
            self._filename = data
        elif hasattr(data, 'filename'):
            # could be a data sets from gdsctools.datasets.Data
            self.read_data(data.filename)
            self._filename = data.filename
        elif hasattr(data, 'df'):
            # an instance of a Reader (or child such as IC50, GenomicFeatures)
            self.df = data.df.copy()
            self._filename = data._filename
        elif isinstance(data, pd.DataFrame):
            # Or just a dataframe ?
            self.df = data.copy()
            self._filename = None
        else:
            raise TypeError("Input must be a filename, a IC50 instance, or " +
                            "a dataframe.")

        #: if populated, can be used to check validity of a header
        self.header = []

    def read_data(self, filename):
        # remove possible white spaces in the header's names
        if ".csv" in filename:
            try:
                # this is to cope with pandas 0.13 on ReadTheDoc
                # and newer versions
                try:
                    rawdf = pd.read_csv(filename, sep=",", comment="#")
                except:
                    rawdf = pd.read_csv(filename, sep="\t", comment="#",
                            compression='gzip')
            except:
                raise ValueError('Could not read %s' % filename)
            rawdf.rename(columns=lambda x: x.strip(), inplace=True)
        elif ".tsv" in filename or '.txt' in filename: # txt not supported
            # officialy but txt file from previous run were interepreted as tsv

            try:
                # this is to cope with pandas 0.13 on ReadTheDoc
                # and newer versions
                try:
                    rawdf = pd.read_csv(filename, sep="\t", comment="#")
                except:
                    rawdf = pd.read_csv(filename, sep="\t", comment="#",
                            compression='gzip')
            except:
                msg = 'Could not read %s' % filename
                raise ValueError(msg)
            rawdf.rename(columns=lambda x: x.strip(), inplace=True)
        else:
            raise ValueError("Only file ending in .csv or .csv.gz or .tsv"+
                " or .tsv.gz will be interpreted.")

        # let us drop columns that are unnamed and print information
        columns = [x for x in rawdf.columns if x.startswith('Unnamed')]
        if len(columns) > 0:
            print('%s  unnamed columns found and removed. ' % len(columns) +
                'Please fix your input file.')
        self.df = rawdf.drop(columns, axis=1)

        # Finally, check that names do not contain the unwanted character
        # / that was used in some old matrices.
        if len([True for x in self.df.columns if "/" in x])>0:
            print("Your input data contains unwanted / characters in " +
                    " the header. Let's remove them.")
            self.df.columns = [x.replace("/", "_") for x in self.df.columns]

    def _interpret(self):
        pass

    def _valid_header(self, df):
        for name in self.header:
            if name not in list(df.columns):
                return False
        return True

    def _read_matrix_from_r(self, name):
        """Required biokit. Will be removed"""
        print("Reading matrix %s " % (name))
        self.session.run("rnames = rownames(%s)" % name)
        self.session.run("cnames = colnames(%s)" % name)
        self.session.run("data = %s" % name)

        cnames = self.session.cnames
        rnames = self.session.rnames
        data = self.session.data
        df = pd.DataFrame(data=data.copy())
        df.columns = [x.strip() for x in cnames]
        df.index = [x.strip() for x in rnames]
        return df

    def __str__(self):
        self.df.info()
        return ""

    def __len__(self):
        return self.df.shape[0] * self.df.shape[1]

    def to_csv(self, filename, sep=",", index=False, reset_index=True):
        """Save data into a CSV file without indices"""
        #Reset the index (e.g., COSMIC ID)
        if reset_index is True:
            df = self.df.reset_index()
        else:
            df = self.df
        df.to_csv(filename, sep=sep, index=index)

    def check(self):
        """Checking the format of the matrix

        Currently, only checks that there is no duplicated column names

        """
        if len(self.df.columns.unique()) != len(self.df.columns):
            columns = list(self.df.columns)
            for this in columns:
                if columns.count(this) > 1:
                    msg ='Found identical named columns (%s)' % this
                    raise ValueError(msg)

    def __eq__(self, other):
        return all(self.df.fillna(0) == other.df.fillna(0))


class CosmicRows(object):
    """Parent class to IC50 and GenomicFeatures to handle cosmic identifiers"""

    def _get_cosmic(self):
        return list(self.df.index)
    def _set_cosmic(self, cosmics):
        for cosmic in cosmics:
            if cosmic not in self.cosmicIds:
                raise ValueError('Unknown cosmic identifier')
        self.df = self.df.ix[cosmics]
    cosmicIds = property(_get_cosmic, _set_cosmic,
            doc="return list of cosmic ids (could have duplicates)")

    def drop_cosmic(self, cosmics):
        """drop a drug or a list of cosmic ids"""
        cosmics = easydev.to_list(cosmics)
        tokeep = [x for x in self.cosmicIds if x not in cosmics]
        self.cosmicIds = tokeep


class IC50(Reader, CosmicRows):
    """Reader of IC50 data set

    This input matrix must be a comman-separated value (CSV) or
    tab-separated value file (TSV).

    The matrix must have a header and at least 2 columns. If the number of rows
    is not sufficient, analysis may not be possible.

    The header must have a column called "COSMIC_ID" or "COSMIC ID".
    This column will be used as indices (row names). All other columns will
    be considered as input data.

    The column "COSMIC_ID" contains the cosmic identifiers (cell line). The
    other columns should be filled with the IC50s corresponding to a pair
    of COSMIC identifiers and Drug. Nothing prevents you to fill the file with
    data that have other meaning (e.g. AUC).

    If at least one column starts with "Drug_", all other columns will be
    ignored. This was implemented for back compatibility.

    If no columns start with "Drug_", that prefix will be added to all column
    names.

    The order of the columns is not important.

    Here is a simple example of a valid TSV file::

        COSMIC_ID   Drug_1_IC50 Drug_20_IC50
        111111      0.5         0.8
        222222      1           2


    A test file is provided in the gdsctools package::

        from gdsctools import ic50_test

    You can read it using this class and plot information as follows:

    .. plot::
        :width: 80%
        :include-source:

        from gdsctools import IC50, ic50_test
        r = IC50(ic50_test)
        r.plot_ic50_count()

    You can get basic information using the print function::

        >>> from gdsctools import IC50, ic50_test
        >>> r = IC50(ic50_test)
        >>> print(r)
        Number of drugs: 11
        Number of cell lines: 988
        Percentage of NA 0.206569746043


    You can get the drug identifiers as follows::

        r.drugIds

    and set the drugs, which means other will be removed::

        r.drugsIds = ['Drug_1_IC50', 'Drug_1000_IC50']

    .. versionchanged:: 0.9.10
        The column **COSMIC ID** should now be **COSMIC_ID**.
        Previous name is deprecated but still accepted.

    """
    cosmic_name = 'COSMIC_ID'
    drug_prefix = 'Drug'

    def __init__(self, filename, drug_prefix=''):
        """.. rubric:: Constructor

        :param filename: input filename of IC50s. May also be an instance
            of :class:`IC50` or a valid dataframe. The data is stored as a
            dataframe in the attribute called :attr:`df`. Input file may be
            gzipped

        """
        super(IC50, self).__init__(filename)
        # interpret the raw data and check some of its contents
        self._interpret()
        self.check()

    def _interpret(self):
        # if there is at least one column that starts with Drug or drug or
        # DRUG or variant then all other columns are dropped except "COSMIC ID"
        if len(self.df) == 0:
            return

        columns = []
        drug_prefix = ''
        for col in self.df.columns:
            if col.startswith(self.drug_prefix + "_"):
                drug_prefix = self.drug_prefix
        self.drug_prefix = drug_prefix

        # Let us rename "COSMIC ID" into "COSMIC_ID" if needed
        _cols = list(self.df.columns)
        if "COSMIC ID" in _cols and self.cosmic_name not in _cols:
            warnings.warn("'COSMIC ID' column name is deprecated since " +
            "0.9.10. Please replace with 'COSMIC_ID'", DeprecationWarning)
            self.df.columns = [x.replace("COSMIC ID", "COSMIC_ID")
                    for x in self.df.columns]
        if "CL" in _cols and "COSMID_ID" not in self.df.columns:
            self.df.columns = [x.replace("CL", "COSMIC_ID")
                    for x in self.df.columns]

        # If the data has not been interpreted, COSMIC column should be
        # found in the column and set as the index
        if self.cosmic_name in self.df.columns:
            self.df.set_index(self.cosmic_name, inplace=True)
            columns = [x for x in self.df.columns
                    if x.startswith(self.drug_prefix)]
            self.df = self.df[columns]
        # If already interpreted, COSMIC name should be the index already.
        elif self.df.index.name == self.cosmic_name:
            columns = [x for x in self.df.columns
                    if x.startswith(self.drug_prefix)]
            columns = list(set(columns))
            self.df = self.df[columns]
        # Otherwise, raise an error
        else:
            raise ValueError("{0} could not be found".format(
                self.cosmic_name))

    def _get_drugs(self):
        return list(self.df.columns)
    def _set_drugs(self, drugs):
        for drug in drugs:
            if drug not in self.drugIds:
                raise ValueError('Unknown drug name')
        self.df = self.df[drugs]
    drugIds = property(_get_drugs, _set_drugs,
                       doc='list the drug identifier name or select sub set')

    def drop_drugs(self, drugs):
        """drop a drug or a list of drugs"""
        drugs = easydev.to_list(drugs)
        tokeep = [x for x in self.drugIds if x not in drugs]
        self.drugIds = tokeep

    def __contains__(self, item):
        if item in self.drugIds:
            return True
        else:
            return False

    def plot_ic50_count(self, **kargs):
        """Plots the fraction of valid/measured IC50 per drug

        :param kargs: any valid parameters accepted by pylab.plot function.
        :return: the fraction of valid/measured IC50 per drug

        """
        data = self.df.count()/len(self.df)
        pylab.clf()
        pylab.plot(data.values, **kargs)
        pylab.grid()
        pylab.xlim([0, len(self.drugIds)+1])
        pylab.xlabel('Drug index')
        pylab.ylim([0,1])
        pylab.ylabel('Percentage of valid IC50')
        return  data

    def hist(self, bins=20, **kargs):
        """Histogram of the measured IC50

        :param bins: binning of the histogram
        :param kargs: any argument accepted by pylab.hist function.
        :return: all measured IC50

        .. plot::
            :include-source:
            :width: 80%

            from gdsctools import IC50, ic50_test
            r = IC50(ic50_test)
            r.hist()

        """
        pylab.clf()
        pylab.hist(self.get_ic50(), bins=bins, **kargs)
        pylab.grid()
        pylab.xlabel('log IC50')

    def get_ic50(self):
        """Return all ic50 as a list"""
        return [x for x in self.df.values.flatten() if not np.isnan(x)]

    def __str__(self):
        txt = "Number of drugs: %s\n" % len(self.drugIds)
        txt += "Number of cell lines: %s\n" % len(self.df)
        N = len(self.drugIds) * len(self.df)
        Nna = self.df.isnull().sum().sum()
        if N != 0:
            txt += "Percentage of NA {0}\n".format(Nna / float(N))
        return txt

    def __repr__(self):
        Nc = len(self.cosmicIds)
        Nd = len(self.drugIds)
        return "IC50 object <Nd={0}, Nc={1}>".format(Nd, Nc)

    def __add__(self, other):
        print("Experimantal. combines IC50 via COSMIC IDs")
        df = pd.concat([self.df, other.df], ignore_index=True)
        df = df.drop_duplicates(cols=[self.cosmic_name])
        return df

    def copy(self):
        new = IC50(self, drug_prefix=self.drug_prefix)
        return new


class GenomicFeatures(Reader, CosmicRows):
    """Read Matrix with Genomic Features

    These are the compulsary column names required (note the spaces):

        - 'COSMIC_ID'
        - 'TISSUE_FACTOR'
        - 'MSI_FACTOR'

    If this column is found, it is removed::

        - 'SAMPLE_NAME'
        - 'Sample Name'

    and features can be also encoded with the following convention:

        - columns ending in "_mut" to encode a gene mutation (e.g., BRAF_mut)
        - columns starting with "gain_cna"
        - columns starting with "loss_cna"

    Those columns will be removed:

        - starting with `Drug_`, which are supposibly from the IC50 matrix

    ::

        >>> from gdsctools import GenomicFeatures
        >>> gf = GenomicFeatures()
        >>> print(gf)
        Genomic features distribution
        Number of unique tissues 27
        Number of unique features 677 with
        - Mutation: 270
        - CNA (gain): 116
        - CNA (loss): 291

    .. versionchanged:: 0.9.10
        The header's columns' names have changed to be more consistant.
        Previous names are deprecated but still accepted.

    """
    colnames = easydev.AttrDict()
    colnames.cosmic = 'COSMIC_ID'
    colnames.tissue = 'TISSUE_FACTOR'
    colnames.msi = 'MSI_FACTOR'
    colnames.media = 'MEDIA_FACTOR'

    def __init__(self, filename=None):
        """.. rubric:: Constructor

        If no file is provided, using the edfault file provided in the
        package that is made of 1001 cell lines times 680 features.

        """
        # first reset the filename to the shared data (if not provided)
        if filename is None:
            from gdsctools.datasets import genomic_features
            filename = genomic_features
        # used in the header so should be ser before call to super()

        super(GenomicFeatures, self).__init__(filename)

        # FIXME Remove columns related to Drug if any. Can be removed in
        # the future
        self.df = self.df[[x for x in self.df.columns
            if x.startswith('Drug_') is False]]

        for this in ['Sample Name', 'SAMPLE_NAME']:
            if this in self.df.columns:
                self.df.drop(this, axis=1, inplace=True)

        # Let us rename "COSMIC ID" into "COSMIC_ID" if needed
        for old, new in {
                    'Tissue Factor Value': 'TISSUE_FACTOR',
                    'MS-instability Factor Value': 'MSI_FACTOR',
                    'COSMIC ID': 'COSMIC_ID'}.items():
            if old in self.df.columns:
                warnings.warn("'%s' column name is deprecated " % old + 
                    " since 0.9.10. Please replace with '%s'" %  new,
                    DeprecationWarning)
                self.df.columns = [x.replace(old, new)
                        for x in self.df.columns]
        if "CL" in self.df.columns and "COSMID_ID" not in self.df.columns:
            self.df.columns = [x.replace("CL", "COSMIC_ID")
                    for x in self.df.columns]

        # There are 3 special columns to hold the factors
        self._special_names = []

        # If tissue factor is not provided, we create and fill it with dummies.
        # OTherwise, we need to change a lot in the original code in ANOVA
        if self.colnames.tissue not in self.df.columns:
            warnings.warn("column named '%s' not found" % self.colnames.tissue,
                    UserWarning)
            self.df[self.colnames.tissue] = ['unspecified'] * len(self.df)
            self._special_names.append(self.colnames.tissue)
        else:
            self._special_names.append(self.colnames.tissue)

        self.found_msi = self.colnames.msi in self.df.columns
        if self.found_msi is False:
            print("WARNING: column named '%s' not found" % self.colnames.msi)
        else:
            self._special_names.append(self.colnames.msi)

        self.found_media = self.colnames.media in self.df.columns
        if self.found_media is False:
            print("WARNING: column named '%s' not found" \
                    % self.colnames.media)
        else:
            self._special_names.append(self.colnames.media)

        # order columns and index
        self._order()

        #
        self._interpret_cosmic()

        #
        self.check()

    def _get_shift(self):
        return len(self._special_names)
    shift = property(_get_shift)

    def _interpret_cosmic(self):
        if self.colnames.cosmic in self.df.columns:
            self.df.set_index(self.colnames.cosmic, inplace=True)
        elif self.colnames.cosmic == self.df.index.name:
            pass
        else:
            error_msg = "the features input file must contains a column " +\
                " named %s" % self.colnames.cosmic
            raise ValueError(error_msg)
        self.df.sort_index(inplace=True)

    def fill_media_factor(self):
        """Given the COSMIC identifiers, fills the MEDIA_FACTOR column

        If already populated, replaced by new content.

        """
        from gdsctools import COSMICInfo
        c = COSMICInfo()
        self.df['MEDIA_FACTOR'] = [c.get(x).SCREEN_MEDIUM
                for x in self.df.index]
        self.found_media = True
        if self.colnames.media not in self._special_names:
            self._special_names.append(self.colnames.media)
        self._order()

    def _order(self):
        others = [x for x in self.df.columns if x not in self._special_names]
        self.df = self.df[self._special_names + others]

    def _get_features(self):
        return list(self.df.columns)
    def _set_features(self, features):
        for feature in features:
            if feature not in self.features:
                raise ValueError('Unknown feature name %s' % feature)
        features = [x for x in features if x.endswith('FACTOR') is False]
        features = self._special_names + features
        self.df = self.df[features]
        self._order()
    features = property(_get_features, _set_features,
                        doc="return list of features")

    def _get_tissues(self):
        return list(self.df[self.colnames.tissue])
    tissues = property(_get_tissues, doc='return list of tissues')

    def _get_unique_tissues(self):
        return list(self.df[self.colnames.tissue].unique())
    unique_tissues = property(_get_unique_tissues, doc='return set of tissues')

    def plot(self):
        """Histogram of the tissues found

        .. plot::
            :include-source:
            :width: 80%

            from gdsctools import GenomicFeatures
            gf = GenomicFeatures() # use the default file
            gf.plot()


        """
        if self.colnames.tissue not in self.df.columns:
            return
        data = pd.get_dummies(self.df[self.colnames.tissue]).sum()
        data.index = [x.replace("_", " ") for x in data.index]
        # deprecated but works for python 3.3
        try:
            data.sort_values(ascending=False)
        except:
            data.sort(ascending=False)
        pylab.figure(1)
        pylab.clf()
        labels = list(data.index)
        pylab.pie(data, labels=labels)
        pylab.figure(2)
        data.plot(kind='barh')
        pylab.grid()
        pylab.xlabel('Occurences')

        # keep the try to prevent MacOS issue
        try:pylab.tight_layout()
        except:pass
        return data

    def __str__(self):
        txt = 'Genomic features distribution\n'
        try:
            tissues = list(self.df[self.colnames.tissue].unique())
            Ntissue = len(tissues)
            txt += 'Number of unique tissues {0}'.format(Ntissue)
            if Ntissue == 1:
                 txt += ' ({0})\n'.format(tissues[0])
            elif Ntissue < 10:
                txt += '\nHere are the tissues: '
                txt += ",".join(tissues) + "\n"
            else:
                txt += '\nHere are the first 10 tissues: '
                txt += ", ".join(tissues[0:10]) + "\n"
        except:
            txt += 'No information about tissues\n'

        if self.found_msi:
            txt +=  "MSI column: yes\n"
        else:
            txt +=  "MSI column: no\n"

        if self.found_media:
            txt +=  "MEDIA column: yes\n"
        else:
            txt +=  "MEDIA column: no\n"

        # -3 since we have also the MSI, tissue, media columns
        # TODO should use shift attribute ?
        Nfeatures = len(self.features)
        txt += '\nThere are {0} unique features distributed as\n'.format(Nfeatures-self.shift)

        n_mutations = len([x for x in self.df.columns if x.endswith("_mut")])
        txt += "- Mutation: {}\n".format(n_mutations)

        n_gain = len([x for x in self.df.columns if x.startswith("gain_cna")])
        txt += "- CNA (gain): {}\n".format(n_gain)
        n_loss = len([x for x in self.df.columns if x.startswith("loss_cna")])
        txt += "- CNA (loss): {}".format(n_loss)
        return txt

    def drop_tissue_in(self, tissues):
        """Drop tissues from the list

        :param list tissues: a list of tissues to drop. If you have only
            one tissue, can be provided as a string. Since rows are removed
            some features (columns) may now be empty (all zeros). If so, those
            columns are dropped (except for the special columns (e.g, MSI).

        """
        tissues = easydev.to_list(tissues)
        mask = self.df[self.colnames.tissue].isin(tissues) == False
        self.df = self.df[mask]
        self._cleanup()

    def keep_tissue_in(self, tissues):
        """Drop tissues not in the list

        :param list tissues: a list of tissues to keep. If you have only
            one tissue, can be provided as a string. Since rows are removed
            some features (columns) may now be empty (all zeros). If so, those
            columns are dropped (except for the special columns (e.g, MSI).

        """
        tissues = easydev.to_list(tissues)
        mask = self.df[self.colnames.tissue].isin(tissues)
        self.df = self.df[mask]
        self._cleanup()

    def _cleanup(self, required_features=0):
        # FIXME: there is view/copy warning here in pandas. it should be fixed
        # or may have side-effects
        to_ignore = self._special_names
        # create a view ignoring the informative columns
        view = self.df[[x for x in self.df.columns if x not in to_ignore]]

        todrop = list(view.columns[view.sum() <= required_features])

        self.df.drop(todrop, axis=1, inplace=True)

    def __repr__(self):
        Nc = len(self.cosmicIds)
        Nf = len(self.features) - self.shift
        try:
            Nt = len(set(self.tissues))
        except:
            Nt = '?'
        return "GenomicFeatures <Nc={0}, Nf={1}, Nt={2}>".format(Nc, Nf, Nt)

    def get_TCGA(self):
        from gdsctools.cosmictools import COSMICInfo
        c = COSMICInfo()
        tcga = c.df.ix[self.df.index].TCGA
        return tcga


class PANCAN(Reader):
    """Reads RData file wit all genomic features including methylation.

    will be removed. Used to read original data in R format but
    will provide the data as CSV or TSV
    """
    def __init__(self, filename=None):
        if filename is None:
            filename = easydev.get_share_file('gdsctools', 'data',
                            'PANCAN_simple_MOBEM.rdata')
        super(PANCAN, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self._filename)
        self.df = self._read_matrix_from_r('MoBEM')


class Extra(Reader):
    def __init__(self, filename="djvIC50v17v002-nowWithRMSE.rdata"):
        super(Extra, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self._filename)

        # 3 identical matrices containing AUC, IC50 and
        self.dfAUCv17= self._read_matrix_from_r('dfAUCv17')
        self.dfIC50v17 = self._read_matrix_from_r('dfIC50v17')
        # Residual
        self.dfResv17 = self._read_matrix_from_r('dfResv17')

        # This df holds the xmid/scale parameters for each cell line
        # Can be visualised using the tools.Logistic class.
        self.dfCL= self._read_matrix_from_r('dfCL')

        # There is an extra matrix called MoBEM, which is the same as in the
        # file

    def hist_residuals(self, bins=100):
        """Plot residuals across all drugs and cell lines"""
        data = [x for x in self.dfResv17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('Residuals')
        pylab.ylabel(r'\#')

    def scatter(self):
        from biokit import scatter
        s = scatter.ScatterHist(self.dfCL)
        s.plot(kargs_histx={'color':'red', 'bins':20},
                kargs_scatter={'alpha':0.9, 's':100, 'c':'b'},
                kargs_histy={'color':'red', 'bins':20})

    def hist_ic50(self, bins=100):
        data = [x for x in self.dfIC50v17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('IC50')
        pylab.ylabel(r'\#')

    def hist_auc(self, bins=100):
        data = [x for x in self.dfAUCv17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('AUC')
        pylab.ylabel(r'\#')


class DrugDecode(Reader):
    """Reads a "drug decode" file

    The format must be 3-columns comma-separated file.
    ::

        DRUG_ID     ,DRUG_NAME   ,DRUG_TARGET
        999         ,Erlotinib   ,EGFR
        1039        ,SL 0101-1   ,"RSK, AURKB, PIM3"

    TSV file may also work out of the box. If a column name called
    'PUTATIVE_TARGET' is found, it is renamed 'DRUG_TARGET'.

    """
    def __init__(self, filename=None):
        """.. rubric:: Constructor"""
        super(DrugDecode, self).__init__(filename)
        self.header = ['DRUG_ID', 'DRUG_NAME', 'DRUG_TARGET', 'OWNED_BY',
            'WEBRELEASE']

        #self.df.drop_duplicates(inplace=True)
        try:self._interpret()
        except:pass

    def _interpret(self, filename=None):
        N = len(self.df)
        if N  == 0:
            return

        self.df.rename(columns={'PUTATIVE_TARGET': 'DRUG_TARGET'},
                inplace=True)

        if 'WEBRELEASE' not in self.df.columns:
            self.df['WEBRELEASE'] = [None] * N

        if 'OWNED_BY' not in self.df.columns:
            self.df['OWNED_BY'] = [None] * N

        for this in self.header:
            msg = "Could not read the file with expected format. "
            msg += "It should be a comma separated file."
            msg += " It should be made of 3 columns with this header"
            msg += "%s"
            if this not in self.df.columns and this != self.df.index.name:
                raise ValueError(msg % this)

        try:
            self.df.set_index('DRUG_ID', inplace=True)
        except:
            # could be done already
            pass

    def _get_names(self):
        return list(self.df.DRUG_NAME.values)
    drug_names = property(_get_names)

    def _get_target(self):
        return list(self.df.DRUG_TARGET.values)
    drug_targets = property(_get_target)

    def _get_drug_ids(self):
        return list(self.df.index)
    drugIds = property(_get_drug_ids,
            doc="return list of drug identifiers")

    def _get_row(self, drug_id, colname):
        if drug_id in self.df.index:
            return self.df.ix[drug_id][colname]
        elif str(drug_id).startswith("Drug_"):
            try:
                drug_id = int(drug_id.split("_")[1])
            except:
                print("DRUG ID %s not recognised" % drug_id)
                return None
            if drug_id in self.df.index:
                return self.df[colname].ix[drug_id]
        elif "_" in str(drug_id):
            try:
                drug_id = int(drug_id.split("_")[0])
            except:
                print("DRUG ID %s not recognised" % drug_id)
                return None
            if drug_id in self.df.index:
                return self.df[colname].ix[drug_id]
        else:
            return None

    def get_name(self, drug_id):
        return self._get_row(drug_id, 'DRUG_NAME')

    def get_target(self, drug_id):
        return self._get_row(drug_id, 'DRUG_TARGET')

    def get_public(self, drug_id):
        return self._get_row(drug_id, 'WEBRELEASE')

    def check(self):
        for x in self.drugIds:
            try:
                x += 1
            except TypeError as err:
                print("drug identifiers must be numeric values")
                raise err
        # it may happen that a drug has no target in the database ! so we
        # cannot check that for the moment:
        #if  self.df.isnull().sum().sum()>0:
        #   print(d.df.isnull().sum())
        #    raise ValueError("all values must be non-na. check tabulation")

    def get_info(self):
        # Note that there are 4 cases : Y, N, U (unknown?) and NaN
        dd = {  'N': len(self),
                'N_public': sum(self.df.WEBRELEASE == 'Y'),
                'N_prop': sum(self.df.WEBRELEASE != 'Y')}
        return dd

    def __len__(self):
        return len(self.df)

    def __str__(self):
        txt = "Number of drugs: %s\n" % len(self.df)
        return txt

    def __repr__(self):
        txt = self.__str__()
        if len(self.companies):
            txt += "Contains %s companies" % len(self.companies)
        return txt

    def _get_companies(self):
        if 'OWNED_BY' in self.df.columns:
            companies = list(self.df.OWNED_BY.dropna().unique())
        else:
            companies = []
        return sorted(companies)
    companies = property(_get_companies)

    def drug_annotations(self, df):
        """Populate the drug_name and drug_target field if possible
   
        :param df: input dataframe as given by e.g., :meth:`anova_one_drug`
        :return df: same as input but with the FDR column populated
        """
        if len(self.df) == 0:
            return df
        #      print("Nothing done. DrugDecode is empty.")

        # aliases
        if 'DRUG_ID' not in df.columns:
            raise ValueError('Expected column named DRUG_ID but not found')

        drug_names = [self.get_name(x) for x in df.DRUG_ID.values]
        drug_target = [self.get_target(x) for x in df.DRUG_ID.values]

        # this is not clean. It works but could be simpler surely.
        df['DRUG_NAME'] = drug_names
        df['DRUG_TARGET'] = drug_target
        return df

   
