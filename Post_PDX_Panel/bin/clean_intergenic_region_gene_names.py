#! /usr/bin/env python

from __future__ import print_function

"""
A new version of snpEff started reporting the bounding genes for variants
in intergenic regions, resulting in gene "names" comprising two gene names
separated by a dash.

We don't want these.

This will be fixed in the pipeline, but we don't want to have to rerun
the pipelines just to clean this up.  So this program will go over the existing
*Annotated.tab files, and remove the gene names from any variant whose effect
is "intragenic_region".
"""
import csv
import os
from glob import glob
import argparse


def parse_args():
    """
    Parse the command line arguments
    :return: The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        "Set the gene name to '', when the variant in an intergenic_region.")
    parser.add_argument('-r', '--root',
                        help="Process all *Annotated.tab files under this "
                             "directory.")
    parser.add_argument('-f', '--file-name',
                        help="Process only this file.")
    parser.add_argument('-s', '--suffix', default='.ORIG',
                        help="append this suffix to the original file "
                             "[Default = .ORIG]")
    parser.add_argument('-d', '--delete', action='store_true',
                        help="Delete the original file.")
    args = parser.parse_args()

    if args.root and args.file_name:
        parser.error("May only specify --root or --file-name, not both.")
    if not (args.root or args.file_name):
        parser.error("Must specify exactly one of --root or --file-name")
    return args


def get_files(root):
    """
    Get the list of *Annotated.tab files.
    :param root: The base directory containing the model level directories.
    :return: a list of filepaths to the *Annotated.tab files.
    """
    files = glob(os.path.join(root, '[JT]*', '*', 'analysis', '*',
                              '*Annotated.tab'))
    return files


def get_header(fn):
    """
    Get the header columns as a list.
    :param fn: A filename from which we will use the first row to determine
        the column names.
    :return: a list of column names.
    """
    with open(fn) as f:
        line = f.readline()
        headers = line.strip().split('\t')
    return headers


def needs_processing(fn):
    """
    Check a file to see if it has the condition that we have to clean up,
    before we go to the effort to re-write the whole file.
    :param fn: The filename to check.
    :return: True if we need to process the file.
    """
    f = open(fn)
    reader = csv.DictReader(f, delimiter='\t')
    ret = False
    for row in reader:
        if row['EFF[*].EFFECT'] == 'intergenic_region' and \
                row['EFF[*].GENE'] != '':
            ret = True
            break
    f.close()
    return ret


def process_file(fn, args):
    """
    Process one file.  Make sure the headers are the same as the initial file;
    there is a possibility that different pipelines have different file formats.

    We will rename the original file with the extension .ORIG, and write a new
    file with the original name.

    :param fn: The path to the file to process.
    :param args: The parsed command line arguments
    :return:
    """
    if not needs_processing(fn):
        # Nothing to do.
        return
    headers = get_header(fn)
    orig_file = fn + args.suffix
    os.rename(fn, orig_file)
    in_f = open(orig_file)
    out_f = open(fn, 'wb')

    reader = csv.DictReader(in_f, delimiter='\t')
    writer = csv.DictWriter(out_f, headers, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in reader:
        if row['EFF[*].EFFECT'] == 'intergenic_region':
            # print(row)
            row['EFF[*].GENE'] = ''
            # print(row)
        writer.writerow(row)
    in_f.close()
    out_f.close()
    if args.delete:
        os.remove(orig_file)


def main():
    args = parse_args()
    if args.root:
        files = get_files(args.root)
        for fn in files:
            if not needs_processing(fn):
                # print('***Skipping', fn)
                continue
            print('PROCESSING', fn)
            process_file(fn, args)
    else:
        # The argument parser guarantees that we'll have either root or
        # file_name.
        process_file(args.file_name, args)


if __name__ == '__main__':
    main()
