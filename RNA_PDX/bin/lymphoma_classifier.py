#! /usr/bin/env python

"""
Based on expected expression of a small set of genes, determine whether this
supposed tumor sample has been converted into a lymphoma.

We use two cut-offs, one for passaged tumors, the other for patient.
"""
from __future__ import print_function
import sys
import requests
import argparse
import math


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tumor', default=3, type=int,
                        help="Cut-off for tumor samples")
    parser.add_argument('-p', '--patient', default=7, type=int,
                        help="Cut-off for patient samples")
    parser.add_argument('-o', '--output', default='lymphoma_score.txt',
                        help="Output file into which the total z-score is "
                             "written")
    parser.add_argument('normalized_counts',
                        help="file containing gene name in column 0, and "
                             "normalized gene count in column 3")
    parser.add_argument('expected_expression',
                        help="file containing gene name in column 0, "
                             "expected (average) expression in column 3, "
                             "and standard deviation in column 4")
    parser.add_argument('sample_name',
                        help="first part of the fastq file name, from which "
                             "we extract the sample and model")
    return parser.parse_args()


def get_expected(fn):
    d = {}
    for line in open(fn):
        # parts[0]: gene name
        # parts[1]: expected up/down regulation
        # parts[3]: average expression
        # parts[4]: stddev
        parts = [x.strip() for x in line.split()]
        d[parts[0]] = {'updown': float(parts[1]),
                       'average': float(parts[3]),
                       'stddev': float(parts[4])}
    return d


def z_score(count, updown, average, stddev):
    arbitrary_scale_factor = 25.0
    l2 = math.log(count + 1.0, 2)
    z = updown * (l2 - average) / (arbitrary_scale_factor * stddev)
    return z


def process_counts_file(fn, expected):
    z_total = 0.0
    for line in open(fn):
        parts = [x.strip() for x in line.split()]
        if parts[0] in expected:
            e = expected[parts[0]]
            z_total += z_score(float(parts[3]),
                               e['updown'],
                               e['average'],
                               e['stddev'])
    return z_total


def get_model_sample(name):
    parts = name.split('_')
    # After naming changes were put in place in early 2017, the model is
    # always first, and the sample is always second.
    model = parts[0]
    sample = parts[1]
    return model, sample


def is_lymphoma(sample):
    """
    Try to determine whether this sample is from a model known (in ELIMS) to be
    a lymphoma. If any step fails, assume the model is not a known lymphoma.
    That is safe, because then we'll run the classifier code.
    :param sample: The sample name.
    :return: True if the model for this sample is known to be a lymphoma.
    """
    r = requests.get('http://pdx-dashboard.jax.org/elims/JSON/all',
                     {'id': sample})
    if r.status_code != 200:
        print("In is_lymphoma(): Request for {0} failed with "
              "status {1}. Assuming not lymphoma.".format(sample, r.status_code),
              file=sys.stderr)
        return False

    j = r.json()
    data = j['data']
    if len(data) > 1:
        print("Multiple sample entries returned. Assuming not lymphoma",
              file=sys.stderr)
        return False
    details = data[0]['details']
    if len(details) > 1:
        print("Multiple sets of details returned. Assuming not lymphoma",
              file=sys.stderr)
        return False
    return 'lymphoma' in details[0]['clinical_diagnosis'].lower()


def is_patient(sample):
    if sample.endswith('PT'):
        return True
    if sample[0] != 'J':
        return False

    # Last attempt: look up a J sample.
    r = requests.get('http://pdx-dashboard.jax.org/elims/JSON/all',
                     {'id': sample})
    if r.status_code != 200:
        print("In is_patient(): Request for {0} failed with status {1}. "
              "Assuming not a patient sample.".format(sample, r.status_code))
        return False
    j = r.json()
    data = j['data']
    if len(data) > 1:
        print("Multiple sample entries returned. Assuming not patient",
              file=sys.stderr)
        return False
    summary = data[0]['summary']
    return summary['passage'] == 'Patient'


def main():
    print("Starting lymphoma_classifier")
    print("Starting lymphoma_classifier", file=sys.stderr)
    args = parse_args()
    assume_statuses = False
    try:
        # We are processing data from outside organizations, which don't use
        # our naming schemes.  Handle it gracefully.  Assuming that this
        # is not a lymphoma model, and not a patient sample.
        model, sample = get_model_sample(args.sample_name)
    except:
        assume_statuses = True

    # This whole thing doesn't matter if the model is a lymphoma. Look it
    # up, based on the sample name.
    if (not assume_statuses) and is_lymphoma(sample):
        msg = "Sample is from lymphoma model. Not evaluating."
        print(msg, file=open(args.output, 'w'))
        print(msg)
        print(msg, file=sys.stderr)
        return 0

    expected = get_expected(args.expected_expression)
    z_total = process_counts_file(args.normalized_counts, expected)
    print(z_total, file=open(args.output, 'w'))
    if (not assume_statuses) and is_patient(sample):
        cutoff = args.patient
    else:
        cutoff = args.tumor
    print("Finished lymphoma_classifier")
    print("Finished lymphoma_classifier", file=sys.stderr)
    if z_total > cutoff:
        return 1
    return 0

if __name__ == '__main__':
    main()
