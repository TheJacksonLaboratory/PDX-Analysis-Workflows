#! /usr/bin/env python
from __future__ import print_function
"""
Read a Picard coverage stats file and terminate the run if coverage is below
a threshold.

A picard stats file consists of a metadata section, terminated by the line,

## METRICS CLASS        net.sf.picard.analysis.directed.HsMetrics

followed by two rows; a header and data.

we'll read both lines and then pick out the field matching the requested coverage level.

"""
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--ctp-coverage-level', default='100',
                        help="Coverage level to use for the test on CTP runs: "
                             "2, 10, 20, 30, 40, 50, or 100 [default: 100]")
    parser.add_argument('-d', '--debug', action='store_true',
                        help="Enable some debugging prints")
    parser.add_argument('-x', '--hex-coverage-level', default='20',
                        help="Coverage level to use for the test on HEX runs: "
                             "2, 10, 20, 30, 40, 50, or 100 [default: 20]")
    parser.add_argument('-p', '--ctp-percentage', type=float, default=75,
                        help="Minimum %% of bases covered at that depth. Enter"
                             "as a percentage or decimal (e.g., 75 or 0.75 "
                             "[default: 75]")
    parser.add_argument('-P', '--hex-percentage', type=float, default=75,
                        help="Minimum %% of bases covered at that depth. Enter"
                             "as a percentage or decimal (e.g., 75 or 0.75 "
                             "[default: 75]")
    parser.add_argument('files', nargs='+',
                        help="The file[s] to test.")

    return parser.parse_args()


def process_file(fn, coverage, percentage, multiple, debug):
    """
    Determine whether there was adequate coverage of bases in this file.
    NOTE: Returns True if the run was OK.
    :param fn:
    :param coverage:
    :param percentage:
    :param multiple:
    :param debug:
    :return:
    """
    metadata_ended = False
    pattern = 'PCT_TARGET_BASES_{0}X'.format(coverage)
    header = []
    data = []
    for line in open(fn):
        line = line.strip()
        if not metadata_ended:
            if line.startswith('## METRICS CLASS'):
                metadata_ended = True
            continue
        if not header:
            # process the header line
            header = line.split()
        else:
            # process the data line
            data = line.split()
            break

    try:
        idx = header.index(pattern)
    except ValueError:
        print("Could not find coverage column {0} in: {1}\n"
              "header line is\n{2}".format(pattern, fn, header),
              file=sys.stderr)
        return False

    if debug:
        print("Percentage at {0} is {1}".format(header[idx], data[idx]))
    this_percentage = float(data[idx])
    if this_percentage < percentage:
        if multiple:
            print("{0}X\t{1}\t{2}".format(coverage, this_percentage, fn))
        else:
            print("Too low coverage percentage at {0}X: {1} {2}".format(
                coverage, this_percentage, fn))
            print("Too low coverage percentage at {0}X: {1} {2}".format(
                coverage, this_percentage, fn), file=sys.stderr)
        return False
    return True


def main():
    args = parse_args()

    ctp_percentage = float(args.ctp_percentage)
    hex_percentage = float(args.hex_percentage)
    multiple = len(args.files) > 1
    debug = args.debug

    if ctp_percentage > 1.0:
        ctp_percentage /= 100.0
    if hex_percentage > 1.0:
        hex_percentage /= 100.0

    success = True
    coverage = ""
    percentage = 0.0
    for fn in args.files:
        if '_HEX' in fn:
            coverage = args.hex_coverage_level
            percentage = hex_percentage
        elif '_CTP' in fn or '_TEX' in fn:
            coverage = args.ctp_coverage_level
            percentage = ctp_percentage
        else:
            print("Couldn't determine HEX or CTP. Assuming CTP.", file=sys.stderr)

        success &= process_file(fn, coverage, percentage, multiple, debug)
    if not success:
        sys.exit(1)

if __name__ == '__main__':
    main()
