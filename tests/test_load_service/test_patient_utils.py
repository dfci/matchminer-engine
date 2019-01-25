import os
import unittest

from src.data_store import key_names as kn
from src.utilities.utilities import get_db
from src.utilities import settings as s
s.MONGO_URI = 'mongodb://localhost:27017'
s.MONGO_DBNAME = 'matchminer'

from src.services.load_service.patient_utils import PatientUtils

# todo enable these unit tests


class TestPatientUtils(unittest.TestCase):

    def setUp(self):
        super(TestPatientUtils, self).setUp()

        self.db = get_db(mongo_uri=s.MONGO_URI, mongo_dbname=s.MONGO_DBNAME)
        self.p = PatientUtils(db=self.db, mongo_uri=s.MONGO_URI, mongo_dbname=s.MONGO_DBNAME)

        data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../', 'data'))
        self.example_clinical_csv = os.path.join(data_path, 'example_clinical.csv')
        self.example_genomic_csv = os.path.join(data_path, 'example_genomic.csv')
        self.example_clinical_pkl = os.path.join(data_path, 'example_clinical.pkl')
        self.example_genomic_pkl = os.path.join(data_path, 'example_genomic.pkl')
        self.example_samples_bson = os.path.join(data_path, 'samples.bson')

    def tearDown(self):
        self.db.samples.drop()

    def test_load_csv(self):
        return

        self.p.load_csv(clinical=self.example_clinical_csv, genomic=self.example_genomic_csv)

        assert len(self.p.clinical_df.index) == 2, len(self.p.clinical_df.index)
        assert len(self.p.genomic_df.index) == 11

    def test_load_pkl(self):
        return

        self.p.load_pkl(clinical=self.example_clinical_pkl, genomic=self.example_genomic_pkl)

        assert len(self.p.clinical_df.index) == 2
        assert len(self.p.genomic_df.index) == 11

    def test_load_bson(self):
        return

        self.p.load_bson(samples_bson=self.example_samples_bson)
        samples = list(self.db.samples.find())
        assert samples is not None
        assert len(samples) == 2
