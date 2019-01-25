import os
import json
import yaml
import unittest

from src.utilities import settings as s
from src.utilities.utilities import get_db
from src.services.match_service.match_utils.main import main as me
s.MONGO_URI = 'mongodb://localhost:27017'
s.MONGO_DBNAME = 'matchminer'


class TestBugfix(unittest.TestCase):

    def setUp(self):
        super(TestBugfix, self).setUp()
        self.db = get_db(mongo_uri=s.MONGO_URI, mongo_dbname=s.MONGO_DBNAME)
        self.db.samples.drop()
        self.db.trial.drop()

        # insert example sample data
        test_sample_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_sample.json'))
        with open(test_sample_file) as f:
            data = json.load(f)

        self.db.samples.insert(data)

        # insert example trial data
        yml = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_trial.yml'))
        with open(yml) as f:
            t = yaml.load(f.read())
            self.db.trial.insert_one(t)

    def tearDown(self):
        self.db.samples.drop()
        self.db.trial.drop()

    def test_bugfix(self):
        setattr(self, 'mongo_uri', s.MONGO_URI)
        setattr(self, 'mongo_dbname', s.MONGO_DBNAME)
        setattr(self, 'protocol_nos', None)
        setattr(self, 'sample_ids', None)
        me(self)
