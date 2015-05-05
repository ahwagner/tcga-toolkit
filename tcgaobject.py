#!/usr/bin/env python3
__author__ = 'Alex H Wagner'

from pathlib import Path
import os
import pickle
import datetime as dt
import requests
from requests_futures.sessions import FuturesSession
import json
import time


class Dict2(dict):

    def __getattr__(self, item):
        if item.startswith('__'):
            return super().item
        return self[item]


class TCGAObject:

    queries = list()
    timestamps = dict()
    sessions = dict()
    j = json.JSONDecoder()
    futures_session = FuturesSession(max_workers=20)

    def __init__(self, uuid, allow_update=True):
        self.attributes = Dict2(uuid=uuid, timestamp=None)
        self.load(allow_update)

    def load(self, allow_update=True):
        file = self.attrib_file()
        if not file.parent.exists():
            os.mkdir(str(file.parent))
        if file.exists():
            with file.open('rb') as f:
                self.attributes = pickle.load(f)
        if self.is_stale() and allow_update:
            # Load from TCGA API
            response = self.send_request()
            j = json.JSONDecoder()
            tcga_element = j.decode(response.text)['tcgaElement']
            self.set_attributes(tcga_element)

    def save(self):
        file = self.attrib_file()
        with file.open('wb') as f:
            pickle.dump(self.attributes, f)

    def attrib_file(self):
        root = Path(os.path.realpath(__file__)).parent
        return (root / 'tcga_objects' / self.uuid).with_suffix('.pickle')

    def send_request(self):
        print("Sending query for {0}".format(self.uuid))
        request = 'https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/json/uuid/{0}'.format(self.uuid)
        response = requests.get(request)
        return response

    def is_stale(self):
        return self.timestamp is None or dt.datetime.now() - self.timestamp > dt.timedelta(days=30)

    def set_attributes(self, attributes):
        self.attributes = self._scan_attributes(attributes)
        self.attributes['timestamp'] = dt.datetime.now()
        self.save()

    def _scan_attributes(self, attributes):
        if "@href" in attributes:
            if len(attributes) > 1:
                raise AttributeError('Expected only one attribute.')
            uuid = attributes["@href"].split('/')[-1]
            return Dict2(uuid=uuid)
        if not isinstance(attributes, dict):
            return attributes
        for attribute in attributes:
            # print("Analyzing {0}...".format(attribute))
            attributes[attribute] = self._scan_attributes(attributes[attribute])
        return Dict2(attributes)

    def __getattr__(self, item):
        val = self.attributes[item]
        if isinstance(val, dict) and 'uuid' in val:
            return TCGAObject(val['uuid'])
        else:
            return val

    def __dir__(self):
        x = super().__dir__()
        y = [z for z in x if not z.startswith('_')] + list(self.attributes)
        return y


def preload_cohort(uuids):
    stale = dict()
    linked_uuids = list()
    for uuid in uuids:
        # Check if fresh file
        o = TCGAObject(uuid, False)
        if o.is_stale():
            stale[uuid] = o
    print("Updating {0} stale UUIDs.".format(len(stale)))

    uuids = sorted(stale.keys())
    print("Quickly queueing queries.")
    for i, uuid in enumerate(uuids):
        _add_to_session(uuid, len(uuids) - i)
    print("Resubmitting erroneous responses.")
    while uuids:
        uuid = uuids.pop(0)
        try:
            response = TCGAObject.sessions[uuid].result()
            if response.status_code != 200:
                # print("Resetting {0} due to response {1}.".format(uuid, response.status_code), flush=True)
                raise ConnectionError
        except (ConnectionError, ConnectionRefusedError, requests.exceptions.ConnectionError):
            print("Retrying UUID: {0}".format(uuid), flush=True)
            uuids.append(uuid)
            _add_to_session(uuid, len(uuids))
        else:
            tcga_element = TCGAObject.j.decode(response.text)['tcgaElement']
            o = stale[uuid]
            o.set_attributes(tcga_element)
            for a in o.attributes:
                v = o.attributes[a]
                if isinstance(v, dict) and 'uuid' in v:
                    linked_uuids.append(v['uuid'])
    if linked_uuids:
        print("Examining {0} linked UUIDs.".format(len(linked_uuids)))
        preload_cohort(linked_uuids)


def _add_to_session(uuid, remaining=None):
    queries = TCGAObject.queries
    sessions = TCGAObject.sessions
    futures_session = TCGAObject.futures_session
    timestamps = TCGAObject.timestamps
    if len(queries) == 950:
        ts = queries.pop(0)
        delta = dt.datetime.now() - ts
        if delta < dt.timedelta(minutes=3):  # We have not yet passed 3 minutes since submission
            t = (dt.timedelta(minutes=3) - delta).seconds + 1
            if remaining is None:
                print("Waiting {0} seconds.".format(t))
            else:
                print("Waiting {0} seconds (server limit -- {1} objects in queue).".format(t, remaining))
            time.sleep(t)
    # make query
    s = futures_session.get('https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/json/uuid/{0}'.format(uuid))
    sessions[uuid] = s
    ts = dt.datetime.now()
    timestamps[uuid] = ts
    queries.append(ts)

if __name__ == '__main__':
    p = Path('/Users/awagner/Workspace/sclc/rnaseq/background/')
    q = p / 'unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.Level_3.1.14.0'
    gene_summary_files = [x for x in q.iterdir() if str(x).endswith('rsem.genes.normalized_results')]
    q = p / 'unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.Level_3.1.9.0'
    gene_summary_files += [x for x in q.iterdir() if str(x).endswith('rsem.genes.normalized_results')]
    test = [x.name.split('.')[2] for x in gene_summary_files]
    preload_cohort(test)
    # o = TCGAObject('0c5d09fe-e0e0-4cfa-9cd8-802b5725a3f6')
    # print(o.sample.sample.sampleType)