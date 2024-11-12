#! /usr/bin/env python3

import os
import sys
import argparse
import fnmatch
import threading

import boto3
from botocore.exceptions import ClientError

# This code assumes the AWS region and access keys are available
# from the AWS config file, e.g. via the AWS_CONFIG_FILE environment
# variable.

class ProgressPercentage(object):

    def __init__(self, filename):
        self._filename = filename
        if os.path.exists (filename):
            self._size = float(os.path.getsize(filename))
        else:
            self._size = 0
        self._seen_so_far = 0
        self._lock = threading.Lock()

    def __call__(self, bytes_amount):
        with self._lock:
            self._seen_so_far += bytes_amount
            basename = os.path.basename(self._filename)
            if self._size > 0:
                percentage = (self._seen_so_far / self._size) * 100
                sys.stdout.write("\n%s  %f / %f  (%.2f%%)" % (basename, self._seen_so_far, self._size, percentage))
            else:
                sys.stdout.write("\n%s  %f" % (basename, self._seen_so_far))
            sys.stdout.flush()

class asdc_s3:

    def __init__ (self, bucket, destination_dir):
        self.bucket = bucket
        self.destination_dir = destination_dir.strip('/')
        self.s3client = boto3.client('s3')

    def upload_file (self, path, object_name=None):
        if object_name is None:
            object_name = os.path.join (self.destination_dir, os.path.basename(path))
        try:
            response = self.s3client.upload_file (path, self.bucket, object_name)
        except ClientError as e:
            print (e)
            return False
        return True

    def download_file (self, filename, local_dir=None):
        object_name = os.path.join (self.destination_dir, filename)
        if local_dir is None:
            local_filename = filename
        else:
            if not os.path.isdir(local_dir):
                print ('*** Error: nonexistent directory: {}'.format(local_dir))
                return False
            local_filename = os.path.join (local_dir, filename)
        try:
            response = self.s3client.download_file (self.bucket, object_name, local_filename)
        except ClientError as e:
            print (e)
            return False
        return True

    def list_objects (self):
        try:
            response = self.s3client.list_objects_v2 (Bucket = self.bucket, Prefix = self.destination_dir + '/')
        except ClientError as e:
            print (e)
            return None
        return response

    def list_files (self, pat):
        object_list = self.list_objects ()
        if object_list is None:
            return None
        files = []
        for obj in object_list["Contents"]:
            p = obj['Key'].split('/')
            if len(p[1]) > 0:
                files.append(p[1])
        if pat is not None:
            files = fnmatch.filter (files, pat)
        return files

    def remove_file (self, filename):
        object_name = os.path.join (self.destination_dir, filename)
        try:
            response = self.s3client.delete_object (Bucket = self.bucket, Key = object_name)
        except ClientError as e:
            print (e)
            return False
        return True

def read_files (listfile):
    with open (listfile, "r") as fp:
        filelist = fp.read()
    filelist = filelist.split ('\n')
    filelist = filter (None, filelist)
    return filelist

def main ():
    parser = argparse.ArgumentParser(description='ASDC AWS S3 bucket put/get/list/remove files')
    parser.add_argument('--bucket', metavar='BUCKET:DIR', default=None,
                        help="S3 bucket specification")
    parser.add_argument('--list', action='store_true',
                        help="List files in S3 bucket directory")
    parser.add_argument('--pattern', default=None,
                        help="fnmatch pattern to filter output of --list option")
    parser.add_argument('--put', metavar='filelist', default=None,
                        help="Upload filelist files to S3 bucket directory")
    parser.add_argument('--get', metavar='filelist', default=None,
                        help="Download filelist files from S3 bucket directory")
    parser.add_argument('--dir', default=None,
                        help="Local destination directory for --get option")
    parser.add_argument('--remove', metavar='filelist', default=None,
                        help="Remove filelist files from S3 bucket directory")

    args = parser.parse_args()

    if args.bucket is None:
        parser.print_usage(sys.stderr)
        print('*** Error: undefined S3 bucket')
        sys.exit(1)

    tok = args.bucket.split(':')
    bucket = tok[0]
    bucket_dir = tok[1]

    s3 = asdc_s3 (bucket, bucket_dir)

    if args.list:
        files = s3.list_files (args.pattern)
        for filename in files:
            print(filename)
        sys.exit(0)
    elif args.put is not None:
        filelist = read_files (args.put)
        for filename in filelist:
            if os.path.exists (filename):
                s3.upload_file (filename)
    elif args.get is not None:
        filelist = read_files (args.get)
        for filename in filelist:
            s3.download_file (filename, local_dir=args.dir)
    elif args.remove is not None:
        filelist = read_files (args.remove)
        for filename in filelist:
            s3.remove_file (filename)

if __name__ == "__main__":
    main()
