#! /usr/bin/env python
from __future__ import print_function
"""
For determining whether we have enough data to make meaningful RNA expression
estimations, we'll look at the number of human reads available after the
Xenome step.
"""
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--minimum-reads', default='1000000',
                        help="Minimum number of human reads [default: 1000000]")
    parser.add_argument('files', nargs='+',
                        help="The file[s] to test.")

    return parser.parse_args()


def process_file(fn, minimum, multiple):
    """
    Determine whether there was adequate coverage of bases in this file.
    NOTE: Returns True if the run was OK.
    :param fn:
    :param minimum:
    :param multiple:
    :return: True if meets criterion, False otherwise
    """
    line_found = False
    summary_section = False
    line = ""
    for line in open(fn):
        line = line.strip()
        if line == 'Summary':
            summary_section = True
            continue
        if summary_section and line.endswith('human'):
            line_found = True
            break
    if not line_found:
        print("Could not find coverage line in:", fn, file=sys.stderr)
        return False

    count = int(line.split()[0])
    if count < minimum:
        if multiple:
            print("{0}\t{1}".format(count, fn))
        else:
            print("Too low human read count: {0} {1}".format(
                count, fn))
        return False
    return True


def main():
    args = parse_args()

    minimum_reads = int(args.minimum_reads)
    multiple = len(args.files) > 1

    success = True
    if multiple:
        print("Human reads\tRun")

    for fn in args.files:
        success &= process_file(fn, minimum_reads, multiple)
    if not success:
        sys.exit(1)

if __name__ == '__main__':
    main()
