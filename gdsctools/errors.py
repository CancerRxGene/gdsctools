


# A specific error for gdsctools


class GDSCToolsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class GDSCToolsDuplicatedDrugError(GDSCToolsError):
    def __init__(self, drug_id):
        super(GDSCToolsDuplicatedDrugError, self).__init__("")
        self.value = 'Found identical named columns (%s)' % drug_id

