import json
from mrjob.job import MRJob

import corisco.quaternion

class MRWordCounter(MRJob):
    def mapper(self, _, line):
        proc_params = json.loads(line)
        for x in [1, 2]:
            yield _, (corisco.quaternion.Quat(1.0,2.0,3.0,4.0).normalize().q).tolist()
            # yield _, {'thing': x, 'file': proc_params['image_filename']}

    def reducer(self, word, results):
        for r in results:
            yield word, r

if __name__ == '__main__':
    MRWordCounter.run()
