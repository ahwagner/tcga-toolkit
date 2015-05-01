#!/usr/bin/env python3
__author__ = 'Alex H Wagner'

from pathlib import Path
import os
import pickle
import datetime as dt
import requests
import json


class TCGAObject:

    def __init__(self, uuid, object_type='biospecimen'):
        self.uuid = uuid
        self.attributes = {'object_type': object_type, 'timestamp': None}
        self.queries = 0
        self.load_object()

    def load_object(self):
        root = Path(os.path.realpath(__file__)).parent
        file = (root / 'tcga_objects' / self.uuid).with_suffix('.pickle')
        self.queries = 0
        if file.exists():
            with file.open('rb') as f:
                self.attributes = pickle.load(f)
        ts = self.attributes['timestamp']
        if ts is None:
            # Load from TCGA API
            response = self.send_request()
            j = json.JSONDecoder()
            tcga_element = j.decode(response.text)['tcgaElement']
            self.attributes = self.scan_attributes(tcga_element)
            self.attributes['timestamp'] = dt.datetime.now()
            with file.open('wb') as f:
                pickle.dump(self.attributes, f)
            # Save results to pickle
        elif dt.datetime.now() - ts > dt.timedelta(days=30):
            # warn about stale data
            pass

    def send_request(self):
        print("Sending query for {0}".format(self.uuid))
        request = 'https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/json/uuid/{0}'.format(self.uuid)
        response = requests.get(request)
        self.queries += 1
        return response

    def scan_attributes(self, attributes):
        if "@href" in attributes:
            if len(attributes) > 1:
                raise AttributeError('Expected only one attribute.')
            uuid = attributes["@href"].split('/')[-1]
            return {'uuid': uuid}
        if not isinstance(attributes, dict):
            return attributes
        for attribute in attributes:
            print("Analyzing {0}...".format(attribute))
            attributes[attribute] = self.scan_attributes(attributes[attribute])
        return attributes


if __name__ == '__main__':
    biospecimen = TCGAObject('fe82fe98-0116-4657-acb0-36b2093fb639')
    print(biospecimen.queries)
    print(biospecimen.attributes)