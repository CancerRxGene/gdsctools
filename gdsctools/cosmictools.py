
# fixing compatiblity python 2 and 3 related to merging or urllib and urllib2 i
# n python 3
try:     #python 3
    from urllib.request import urlopen
except:
    from urllib2  import urlopen


import pandas as pd
import easydev


__all__ = ['COSMICFetcher', 'COSMIC']


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
            if filename is provided, reads local file. Format should be
            the same as the file downloaded from expasy

        """
        if filename is not None:
            fh = open(filename, 'r')
            self.data = fh.read()
            fh.close()
            self._scandata()
        else:
            url = 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt'
            self.url = url
            print('Downloading data. This may take a while')
            print('Consider saving the *data* attribute in a file ' +
                'for next time')
            self.data = urlopen(self.url).read()
            self._scandata()

    def _scandata(self):
        print('Parsing the data')
        self.data = self.data.split("\nID   ")[1:] # skip header
        print(len(self.data))
        self.data = [this for this in self.data if 'Cosmic' in this]
        print('Dropping records with no COSMIC cross references:')
        print('Kept %s records' % len(self.data))
        self._data2records()

    def _data2records(self):
        print("Creating records")
        self._records = {}
        for this in self.data:
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
        content = [this.split(" ",1)[1].strip() for this in lines]
        return content


class COSMIC(object):
    """A COSMIC object

    Simply hold the cosmic identifier and open browser on the relevant page.

    ::

        >>> c = COSMIC(905940)
        >>> c.on_web()


    .. seealso:: http://www.cancerrxgene.org/translation/CellLine
    """
    def __init__(self, identifier):
        self.identifier = identifier

    def _get_url(self, cosmic_id=None):
        if cosmic_id is None:
            cosmic_id = self.identifier
        url = 'http://cancer.sanger.ac.uk/cell_lines/sample/overview'
        url = url + "?id={0}".format(cosmic_id)
        return url

    def on_web(self):
        from easydev.browser import browse
        url = self._get_url(self.identifier)
        browse(url)

