## NOTE
This version of the MatchEngine is deprecated and no longer supported. Please see [MatchEngineV2](https://github.com/dfci/matchengine-V2).

# MatchEngine
The matchengine matches patient clinical and genomic information to trials. 

## Built with
* [MongoDB](https://docs.mongodb.com/) - NoSQL document database for data storage.
* [nose](http://nose.readthedocs.io/en/latest/) - Python library for unit testing.

All required python libraries can be installed by running `pip install -r requirements.txt`

### User Guide
##### Step 1: Set up MongoDB
The matchengine was initially developed using MongoDB version 3.2. For MongoDB installation instructions
for Linux, Mac OS X, and Windows please visit [their installation page](https://docs.mongodb.com/manual/administration/install-community/).

##### Step 2: Load data

###### Patient data
The matchengine expects patient data to be stored in two separate MongoDB collections:
* **clinical**: Contains clinical attributes like cancer diagnosis and age (see examples/clinical.example.bson for an example)

| MRN | SAMPLE_ID | ONCOTREE_PRIMARY_DIAGNOSIS_NAME | BIRTH_DATE | VITAL_STATUS | GENDER |
| --- | --------- | ------------------------------- | ---------- | ------------ | ------ |
| 01  | SAMPLE-01 | Breast Invasive Ductal Carcinoma| 1900-01-01 | alive        | female |

* **genomic**: Contains all genomic variants sequenced from each patient (see examples/genomic.example.csv for an example)

| SAMPLE_ID | TRUE_HUGO_SYMBOL | TRUE_PROTEIN_CHANGE | TRUE_VARIANT_CLASSIFICATION | VARIANT_CATEGORY | CNV_CALL | TRUE_TRANSCRIPT_EXON | WILDTYPE |
| --------- | ---------------- | ------------------- | --------------------------- | ---------------- | -------- | -------------------- | -------- |
| SAMPLE-01 | PIK3CA           | p.H1047R            | Missense_Mutation           | MUTATION         |          | 8                    | false    |

Clinical and genomic files can be imported to MongoDB using the matchengine in CSV, PKL, and JSON format.
MongoDB will store these collections in JSON format and is able to export the files again
in BSON, JSON, and CSV format. For more information see 
[mongodump](https://docs.mongodb.com/manual/reference/program/mongodump/) and
[mongoexport](https://docs.mongodb.com/manual/reference/program/mongoexport/)

###### Trial data
The matchengine expects trial data to also be stored in a separate MongoDB collection. Matching
information is stored in a nested structure under the root field name "treatment_list".
Trials can be imported to MongoDB using the matchengine in YML or JSON format. In YML format,
an example of the trial structure would be:
```yaml
protocol_no: 00-000
nct_id: NCT000
treatment_list:
  step:
  - arm:
    - arm_code: A
      arm_description: 'Example Arm A'
      arm_internal_id: 1
      arm_suspended: N
      dose_level: []
      match:
        - and:
          - clinical:
              oncotree_primary_diagnosis: Breast
              age_numerical: '>=18'
          - or:
            - genomic:
                hugo_symbol: PIK3CA
                variant_category: Mutation
                protein_change: p.H1047R
            - genomic:
                hugo_symbol: TP53
                variant_category: Mutation

```

There are several genomic variants that can be curated in this way. Beneath is a map
detailing how the trial field names correspond to the patient data field names:

| trial field name        | genomic field name               | example               |
| ----------------------- | -------------------------------- | --------------------- |
| hugo_symbol             | TRUE_HUGO_SYMBOL                 | ERBB2                 |
| protein_change          | TRUE_PROTEIN_CHANGE              | p.T790M               |
| wildcard_protein_change | TRUE_PROTEIN_CHANGE              | p.G719                |
| variant_classification  | TRUE_VARIANT_CLASSIFICATION 	 | In_Frame_Del          |
| variant_category 		  | VARIANT_CATEGORY 				 | Mutation              |
| exon 					  | TRUE_TRANSCRIPT_EXON 			 | 10                    |
| cnv_call 				  | CNV_CALL 						 | Heterozygous deletion |
| wildtype 				  | WILDTYPE 						 | True or False         |

| trial field name        | clinical field name              | example               |
| ----------------------- | -------------------------------- | --------------------- |
| oncotree_diagnosis      | ONCOTREE_PRIMARY_DIAGNOSIS_NAME  | Breast Invasive Ductal Carcinoma |
| age_numerical			  | BIRTH_DATE                       | 1900-01-01            | 

###### variant_classification options:
- Missense_Mutation
- In_Frame_Del
- Nonsense_Mutation
- Splice_Region
- Frame_Shift_Del
- Splice_Site
- In_Frame_Ins

###### variant_category options:
- Mutation
- Copy Number Variation
- Structural Variation
- Signature

###### cnv_call options (for '''variant_category: Copy Number Variation''' only)
- Heterozygous deletion
- Homozygous deletion
- Gain
- High level amplification

###### Our example
To import example data run:
```bash
python matchengine.py load -t examples/trial.example.yml -c examples/clinical.example.csv -g examples/genomic.example.csv --mongo-uri ${your_mongo_uri}
```

* By default, `load` inserts the data into a database named `matchminer`.
* For more information on linking your Mongo URI please see these [docs](https://docs.mongodb.com/manual/reference/connection-string/).
  For default mongo shell configurations this will likely be `mongodb://localhost:27017`
* Default trial file format is YML. To change this specify `--trial-format {yml,json,bson}`
* Default clinical file format is CSV. To change this specify `--trial-format {csv,pkl,bson}`

    
##### Step 2: Matching
Once your MongoDB is set up you can perform matching by running:
```bash
python matchengine.py match --mongo-uri ${your_mongo_uri}
```

Default output will be a csv file called "results.csv" in your current working directory.
You can specify the outpath path and filename of the results by setting the `-o` flag. <br>
***NOTE***: If using `-o`, please specify output directory **and** filename. 
You can change the file format of the output to JSON by setting the `--json` flag.

### Unit testing
The matchengine uses nose for unit testing. To run all tests from the repository's
root directory:
```bash
nosetests tests
```

## Authors
* **Zachary Zwiesler**
* **Priti Kumari**
* **James Lindsay**
