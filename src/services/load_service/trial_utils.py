import os
import yaml
import subprocess


class TrialUtils:

    def __init__(self, db, mongo_uri, mongo_dbname):

        self.db = db
        self.mongo_uri = mongo_uri
        self.mongo_dbname = mongo_dbname
        self.load_dict = {
            'yml': self.yaml_to_mongo,
            'bson': self.bson_to_mongo,
            'json': self.json_to_mongo
        }

    def yaml_to_mongo(self, yml):
        """
        If you specify the path to a directory, all files with extension YML will be added to MongoDB.
        If you specify the path to a specific YML file, it will add that file to MongoDB.

        :param yml: Path to YML file.
        """

        # search directory for ymls
        if os.path.isdir(yml):
            for y in os.listdir(yml):
                ymlpath = os.path.join(yml, y)

                # only add files of extension ".yml"
                if ymlpath.split('.')[-1] != 'yml':
                    continue

                # convert yml to json format
                add_trial(ymlpath, self.db)
        else:
            add_trial(yml, self.db)

    def bson_to_mongo(self, bson):
        """
        If you specify the path to a directory, all files with extension BSON will be added to MongoDB.
        If you specify the path to a specific BSON file, it will add that file to MongoDB.

        :param bson: Path to BSON file.
        """
        cmd = "mongorestore --uri %s --db %s %s" % (self.mongo_uri, self.mongo_dbname, bson)
        subprocess.call(cmd.split(' '))

    def json_to_mongo(self, json):
        """
        If you specify the path to a directory, all files with extension JSON will be added to MongoDB.
        If you specify the path to a specific JSON file, it will add that file to MongoDB.

        :param json: Path to JSON file.
        """
        cmd = "mongoimport --uri %s --db %s --collection trial --file %s" % (self.mongo_uri, self.mongo_dbname, json)
        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
        if p.returncode != 0:
            cmd += ' --jsonArray'
            subprocess.call(cmd.split())


def add_trial(yml, db):
    """
    Adds file in YAML format to MongoDB

    :param yml: Path to file
    :param db: MongoDB connection
    """

    with open(yml) as f:
        t = yaml.load(f.read())
        db.trial.insert_one(t)
