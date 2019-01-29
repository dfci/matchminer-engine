import logging

from src.utilities import settings as s
from src.utilities.utilities import get_coordinating_center
from src.services.match_service.match_utils.matchengine.matchengine import MatchEngine

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s', )


class TrialUtils:

    def __init__(self, trial, mongo_uri=None, mongo_dbname=None):
        self.trial = trial
        self.mongo_uri = mongo_uri
        self.mongo_dbname = mongo_dbname
        self.accrual_status = s.match_accrual_status_open_val

    def parse_match_trees_from_trial(self):
        """
        Iterates through all match clauses of the given trial.

        :return: {list of MatchEngine class objects} Each item is a MatchEngine object instantiated with a match tree
        """
        logging.info('Matching trial %s' % self.trial[s.trial_protocol_no_col])

        # If the trial is not open to accrual, all matches to all match trees in this trial will be marked closed
        if s.trial_summary_col in self.trial:

            trial_summary = self.trial[s.trial_summary_col]
            if s.trial_accrual_status_col in trial_summary:
                trial_status = trial_summary[s.trial_accrual_status_col]

                if isinstance(trial_status, list) and 'value' in trial_status[0]:
                    if trial_status[0]['value'].lower() != s.trial_accrual_status_open_val:
                        self.accrual_status = s.match_accrual_status_closed_val

        # match trees can be pulled from step, arm or dose levels
        match_trees = []

        # STEP #
        for step in self.trial[s.trial_treatment_list_col][s.trial_step_col]:
            if s.trial_match_tree_col in step:
                trial_info = self._get_trial_info(level=s.trial_step_col, step=step)
                matchengine = MatchEngine(match_tree=step[s.trial_match_tree_col][0],
                                          trial_info=trial_info,
                                          mongo_uri=self.mongo_uri,
                                          mongo_dbname=self.mongo_dbname)
                match_trees.append(matchengine)

            # ARM #
            for arm in step[s.trial_arm_col]:
                if s.trial_match_tree_col in arm:
                    trial_info = self._get_trial_info(level=s.trial_arm_col, step=step, arm=arm)
                    matchengine = MatchEngine(match_tree=arm[s.trial_match_tree_col][0],
                                              trial_info=trial_info,
                                              mongo_uri=self.mongo_uri,
                                              mongo_dbname=self.mongo_dbname)
                    match_trees.append(matchengine)

                # DOSE #
                for dose in arm[s.trial_dose_col]:
                    if s.trial_match_tree_col in dose:
                        trial_info = self._get_trial_info(level=s.trial_dose_col, step=step, arm=arm, dose=dose)
                        matchengine = MatchEngine(match_tree=dose[s.trial_match_tree_col][0],
                                                  trial_info=trial_info,
                                                  mongo_uri=self.mongo_uri,
                                                  mongo_dbname=self.mongo_dbname)
                        match_trees.append(matchengine)

        return match_trees

    def _get_trial_info(self, **kwargs):
        """
        Parse match tree level, step code, arm code, and dose level code for the given match tree

        :param kwargs:
            - step {dict} -- Step level information
            - arm {dict}  -- Arm level information
            - dose {dict} -- Dose level information
            - level {str}  -- (e.g. step, arm, dose)

        :return: {dict}
        """

        accrual_status = self.accrual_status
        accrual_status_dict = {True: 'open', False: 'closed'}
        accrual_suspended_key_dict = {'step': 'step', 'arm': 'arm', 'dose': 'level'}
        for key in ['step', 'arm', 'dose']:
            suspended_key = '%s_suspended' % accrual_suspended_key_dict[key]
            if key in kwargs and suspended_key in kwargs[key]:
                accrual_status = accrual_status_dict[kwargs[key][suspended_key].lower() == 'n']

        return {
            'protocol_no': self.trial[s.trial_protocol_no_col],
            'accrual_status': accrual_status,
            'level': kwargs['level'],
            'step_code': kwargs['step'][s.trial_step_code_col] if 'step' in kwargs else None,
            'arm_code': kwargs['arm'][s.trial_arm_code_col] if 'arm' in kwargs else None,
            'dose_code': kwargs['dose'][s.trial_dose_code_col] if 'dose' in kwargs else None,
            s.trial_coordinating_center_col: get_coordinating_center(self.trial)
        }
