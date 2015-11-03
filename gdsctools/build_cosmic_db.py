import urllib2
import pandas as pd
import easydev

class COSMICFetcher(object):
    """

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

    Skip all records without Cosmic cross references

    ::

        >>> cf = COSMICFetcher()
        >>> cf.df.ix[0]
        ID                                       OS-A
        AC                                  CVCL_0C23
        CA                           Cancer cell line
        COSMIC_ID                             2239090
        Disease      C4917; Small cell lung carcinoma
        Name: 0, dtype: object

    """
    def __init__(self):
        self.url = 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt'
        self.readdata()

    def readdata(self):
        print('Downloading data')
        self.data = urllib2.urlopen(self.url).read()
        self.data = self.data.split("\nID   ")[1:] # skip header
        print(len(self.data))
        self.data = [this for this in self.data if 'Cosmic' in this]
        print('Dropping records with no COSMIC cross references:')
        print('Kept %s records' % len(self.data))
        self._data2records()

    def _data2records(self):
        print("Creating records")
        self.records = {}
        for this in self.data:
            record = this.split("\n",1)
            identifier = record[0].strip()
            content = record[1].strip()
            self.records[identifier] = content

        # we want to store, AC, ID, DR if Cosmic or GDSC, OX (organism)
        # DI disease

        print("Scannig records")
        records = []
        pb = easydev.Progress(len(self.records))
        count = 0
        for ID,this in self.records.items():
            count += 1
            pb.animate(count)
            # those are to be found only once
            AC = self.scan_record_for(this, 'AC')[0] # should have only one
            OX = self.scan_record_for(this, 'OX')[0] # should have only one
            CA = self.scan_record_for(this, 'CA')[0]

            try:
                OX = OX.split("!")[1].strip()
            except:
                pass

            for line in this.split("\n"):

                # get DI. Most of the time there is only one but could have 2
                # sometimes
                DI = "__".join(self.scan_record_for(this, 'DI'))
                DI = DI.replace("NCIt;", "")
                DI = DI.strip()

                if line.startswith('DR'):
                    dummy, content = line.split(" ", 1)
                    if 'Cosmic' in content:
                        content = content.replace('Cosmic;', '').strip()
                        content = content.replace('Cosmic-CLP;', '').strip()
                        content = content.replace('CC;', '').strip()
                        records.append([ID, AC, OX, CA, int(content), DI])

        self.records_list = records
        self.df = pd.DataFrame(records, columns=['ID', 'AC', 'OX', 'CA',
                'COSMIC_ID', 'Disease'])

        # keep only homo sapiens (drop mus musculus)
        self.df = self.df[self.df.OX == 'Homo sapiens']
        del self.df['OX'] 
        self.df.drop_duplicates(inplace=True)
        self.df.reset_index(drop=True, inplace=True)

    def scan_record_for(self, record, key):
        lines = [line for line in record.split("\n") if line.startswith(key)]
        content = [this.split(" ",1)[1].strip() for this in lines]
        return content
