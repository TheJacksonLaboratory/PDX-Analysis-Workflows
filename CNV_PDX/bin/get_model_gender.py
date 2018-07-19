#! /usr/bin/env python
from __future__ import print_function
import sys
import requests
import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--details', action='store_true',
                        help="Return ID information as well as gender")
    parser.add_argument('id', help="The model ID whose gender is needed")

    return parser.parse_args()


def main():
    args = parse_args()
    # for test:
    url = 'http://JSON/gender'
    # for production:
    url = 'http://elims/JSON/gender'
    r = requests.get(url, params={'id': args.id})
    if r.status_code == 200:
        try:
            d = r.json()
        except:
            print('JSON decoding failed. Here is the returned string:',
                  file=sys.stderr)
            print(r.text, file=sys.stderr)
            print('And the request...', file=sys.stderr)
            print(r.request.__dict__, file=sys.stderr)

            sys.exit(3)

        if args.details:
            print('Query ID: {0}\tInventory Code: {1}\t'
                  'Model ID: {2}\tGender: {3}'.
                  format(d['query_id'], d['inventory_code'],
                         d['model_id'], d['gender']))
        else:
            print(d['gender'].lower())
        if d['gender'] == "NOT FOUND":
            # We're about to exit with error status. Write a reason to the log.
            print("Couldn't find model {0} in the database.".format(args.id),
                  file=sys.stderr)
            sys.exit(1)
    else:
        print('Request failed with status code:', r.status_code, file=sys.stderr)
        sys.exit(2)

if __name__ == '__main__':
    main()
