import json
import logging

from src.utilities import settings as s
from src.utilities.utilities import get_db

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s', )


class SharedUtils(object):

    def __init__(self, mongo_uri, mongo_dbname):
        self.db = get_db(mongo_uri, mongo_dbname)

    def find_trials(self, query=None):
        """
        Pulls all trials out of the database.

        :param: {null or dict}
        :return: {list of dict} Each item is a trial record
        """
        logging.info('Retrieving trials from the database')

        if query is None:
            query = {}

        proj = {
            '_id': 0,
            s.trial_protocol_no_col: 1,
            s.trial_nct_id_col: 1,
            s.trial_treatment_list_col: 1,
            s.trial_summary_col: 1
        }
        return list(self.db.trial.find(query, proj))

    def add_trial_matches(self, trial_matches_df):
        """
        Add trial matches to database

        :param trial_matches_df: {Pandas dataframe}
        :return: {null}
        """
        if trial_matches_df is None:
            return

        logging.info('Adding trial matches to database')

        if len(trial_matches_df.index) > 0:
            records = json.loads(trial_matches_df.T.to_json()).values()
            self.db.trial_match.drop()
            self.db.trial_match.insert_many(records)
