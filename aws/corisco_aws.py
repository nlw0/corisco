import boto
import json
from mrjob.job import MRJob

import corisco.process

def get_image(filename):
    if filename[:6] != 's3n://':
        return open(filename).read()

    conn = boto.connect_s3()

    split_path = filename.split('/')
    bucket_name = split_path[2]
    remaining_path = '/'.join(split_path[3:])

    k = boto.s3.key.Key(conn.get_bucket(bucket_name))
    k.key = remaining_path
    return k.get_contents_as_string()

class MRWordCounter(MRJob):
    def mapper(self, _, input_parameters):
        input_parameters = json.loads(input_parameters)
        image_file = get_image(input_parameters['img_filename'])
        process_args = input_parameters['proc_args']
        image_id = input_parameters['img_id']

        ritL = input_parameters['ransac_iterations']
        gsL = input_parameters['grid_size']

        for rit in ritL:
            for gs in gsL:
                process_args['ransac_itr'] = rit
                process_args['grid']['size'] = gs

                # output_data, pic = corisco.process.estimate_orientation(
                #     process_args,
                #     image_file)
                while True:
                    try:
                        output_data, _ = corisco.process.estimate_orientation(
                            process_args,
                            image_file)
                        break
                    except Exception:
                        continue

                yield image_id, output_data

    def reducer(self, word, results):
        for r in results:
            yield word, r

if __name__ == '__main__':
    MRWordCounter.run()
