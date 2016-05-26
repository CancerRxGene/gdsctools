
# fixing compatiblity python 2 and 3 related to merging or urllib and urllib2 i
# n python 3
try:     #python 3
    from urllib.request import urlopen
except:
    from urllib2  import urlopen


import pandas as pd
import easydev


__all__ = ['COSMICFetcher', 'COSMICInfo']


class COSMICFetcher(object):
    """Utility to download a flat file about cosmic identier and build a small
    dataframe with cosmic identifiers and their diseases

    The original flat file is downloaded from ftp.expasy.org/databases and
    contains records as follows::

        ID         Identifier (cell line name)     Once; starts an entry
        AC         Accession (CVCL_xxxx)           Once
        SY         Synonyms                        Optional; once
        DR         Cross-references                Optional; once or more
        RX         References identifiers          Optional: once or more
        WW         Web pages                       Optional; once or more
        CC         Comments                        Optional; once or more
        DI         Diseases                        Optional; once or more
        OX         Species of origin               Once or more
        HI         Hierarchy                       Optional; once or more
        OI         Originate from same individual  Optional; once or more
        SX         Sex (gender) of cell            Optional; once
        CA         Category                        Optional; once

    We keep only records with COSMIC cross references.
    From those records, we keep ID, AC, CA, DI (Disease) and the cosmic
    identifier.

    The resulting dataframe can then be accessed in the :attr:`df` attribute.

    ::

        >>> from gdsctools.cosmictools import COSMICFetcher
        >>> cf = COSMICFetcher() # this may take a while to download the file
        >>> cf.df.ix[0]
        ID                                       OS-A
        AC                                  CVCL_0C23
        CA                           Cancer cell line
        COSMIC_ID                             2239090
        Disease      C4917; Small cell lung carcinoma
        Name: 0, dtype: object

    """
    def __init__(self, filename=None):
        """.. rubric:: Constructor

        :param str filename: If not provided, download file
            from expasy.org and store it in :attr:`data`. Otherwise,
            if filename is provided, reads a local file. Format should be
            the same as the file downloaded from expasy

        """
        if filename is not None:
            fh = open(filename, 'r')
            self._data = fh.read()
            fh.close()
            self._scandata()
        else:
            url = 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt'
            self.url = url
            print('Downloading data. This may take a while')
            print('Consider saving the *data* attribute in a file ' +
                'for next time')
            self._data = urlopen(self.url).read()
            self._scandata()

    def _scandata(self):
        print('Parsing the data')
        self._data = self._data.split("\nID   ")[1:] # skip header
        print(len(self._data))
        self._data = [this for this in self._data if 'Cosmic' in this]
        print('Dropping records with no COSMIC cross references:')
        print('Kept %s records' % len(self._data))
        self._data2records()

    def _data2records(self):
        print("Creating records")
        self._records = {}
        for this in self._data:
            record = this.split("\n",1)
            identifier = record[0].strip()
            content = record[1].strip()
            self._records[identifier] = content

        # we want to store, AC, ID, DR if Cosmic or GDSC, OX (organism)
        # DI disease

        print("Scanning records")
        records = []
        pb = easydev.Progress(len(self._records))
        count = 0
        for ID,this in self._records.items():
            count += 1
            pb.animate(count)
            # those are to be found only once
            AC = self._scan_record_for(this, 'AC')[0] # should have only one
            OX = self._scan_record_for(this, 'OX')[0] # should have only one
            CA = self._scan_record_for(this, 'CA')[0]

            try:
                OX = OX.split("!")[1].strip()
            except:
                pass

            for line in this.split("\n"):

                # get DI. Most of the time there is only one but could have 2
                # sometimes
                DI = "__".join(self._scan_record_for(this, 'DI'))
                DI = DI.replace("NCIt;", "")
                DI = DI.strip()

                if line.startswith('DR'):
                    dummy, content = line.split(" ", 1)
                    if 'Cosmic' in content:
                        content = content.replace('Cosmic;', '').strip()
                        content = content.replace('Cosmic-CLP;', '').strip()
                        content = content.replace('CC;', '').strip()
                        records.append([ID, AC, OX, CA, int(content), DI])

        self._records_list = records
        self.df = pd.DataFrame(records, columns=['ID', 'AC', 'OX', 'CA',
                'COSMIC_ID', 'Disease'])

        # keep only homo sapiens (drop mus musculus)
        self.df = self.df[self.df.OX == 'Homo sapiens']
        del self.df['OX']
        self.df.drop_duplicates(inplace=True)
        self.df.reset_index(drop=True, inplace=True)

    def _scan_record_for(self, record, key):
        lines = [line for line in record.split("\n") if line.startswith(key)]
        content = [this.split(" ", 1)[1].strip() for this in lines]
        return content


class COSMICInfo(object):
    """Retrieve information about cell line included in GDSC1000

    This file reads a GDSCTools dataset :attr:`gdsctools.datasets.cosmic_info`.
    Its content is stored in :attr:`df`.

    In corresponds to Table S1E (List cell line samples with data 
    availability and annotations across the different omics

    The method :meth:`get` retrieves information
    contained in the dataframe :attr:`df`. Provide a known cosmic identifier
    as follows:

    .. doctest::

        >>> from gdsctools import COSMICInfo
        >>> c = COSMICInfo()
        >>> c.get(909907, 'SAMPLE_NAME')
        'ZR-75-30'

    or get all available field as follows::

        >>> c.get(909907)
        SAMPLE_NAME           ZR-75-30
        SEQ                          1
        CNA                          1
        EXP                          1
        MET                          1
        DRUG_SCR                     1
        GDSC_description_1      breast
        GDSC_description_2      breast
        Study_Abbreviation        BRCA
        MMR                      MSI-L
        SCREEN_MEDIUM                R
        GROWTH_PROPERTIES     Adherent
        Name: 909907, dtype: object

    .. note:: there are only 1000 cell lines in the :attr:`df`. Additional cell
        lines may be retrieve using :class:`COSMICFetcher`

    If a cosmic identifier is not found, the returned object has the same
    structure as above but with all fields set to False.

    .. seealso:: http://www.cancerrxgene.org/translation/CellLine
    """
    def __init__(self):
        """.. rubric:: constructor"""
        from gdsctools.datasets import cosmic_info
        #: dataframe with all information
        self.df = pd.read_csv(cosmic_info.filename, sep=',')
        self.df.set_index('COSMIC_ID', inplace=True)

    def get(self, identifier, colname=None):
        """

        :param int identifier: a cosmic identifiers. Possible values are 
            stored in :attr:`df.index` attribute
        :param colname: specific field. 

        :return: if colname is not provided, returns a time series for the 
            **identifier** with all available fields. Otherwise, returns a
            specific field.
        """
        if isinstance(identifier, str):
            identifier = int(identifier)

        if identifier not in self.df.index:
            ts = pd.Series([None]*12, index=self.df.columns, name=identifier)
        else:
            ts = self.df.ix[identifier]

        if colname is None:
            return ts.copy() # to be safe since user may change it 
        else:
            return ts[colname]

    def _get_url(self, cosmic_id):
        url = 'http://cancer.sanger.ac.uk/cell_lines/sample/overview'
        url = url + "?id={0}#overview".format(cosmic_id)
        return url

    def on_web(self, identifier):
        """Open a tab related to the COSMIC identifier (in your browser)"""
        from easydev import onweb
        url = self._get_url(identifier)
        onweb(url)
