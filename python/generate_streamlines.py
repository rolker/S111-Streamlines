#!/usr/bin/env python3

import argparse
import importlib
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert vector fields to streamlines using the Jobard-Lefer algorithm", formatter_class=argparse.RawTextHelpFormatter)

    list_group = parser.add_argument_group("Time Listing")

    list_group.add_argument("--list", dest="list", action="store_const", const=True, default=False, help="Lists available times in file")
    list_group.add_argument("--list-format", dest="time_format", default=None, nargs="?", help="Specific format for time listing.\nDefault: Standard Python date to string format")

    parser.add_argument("--format", dest="file_format", default="s111_h5", help="Format of file to process.\nDefault: s111_h5", nargs="?")
    parser.add_argument("--processes", dest="process_count", type=int, default=4, help="Number of processes to use when processing whole files.\nDefault: 4", nargs="?")
    parser.add_argument("file", metavar="file", help="The file to process")
    parser.add_argument("output_directory", default=".", nargs="?", help="Output directory for GeoJSON files.\nDefault: current directory")

    generation_group=parser.add_argument_group("Generation From Subsets (Exclusive)")

    subset_group = generation_group.add_mutually_exclusive_group()
    subset_group.add_argument("--index", metavar="INDEX", type=int, dest="file_index", nargs="?", help="A specific subdataset index to process")
    subset_group.add_argument("--subdataset", metavar="NAME", dest="file_subdataset", nargs="?", help="A specific named subdataset to process")
    subset_group.add_argument("--time", metavar="TIME_STAMP", dest="file_time", nargs="?", help="A specific time to process")

    args = parser.parse_args()
    
    file_format_module = importlib.import_module(args.file_format.lower())

    if args.list:
        file_format_module.list_times(args.file, args.time_format)
    elif args.file_index != None:
        file_format_module.generate_indexed_streamlines(args.file, args.output_directory, args.file_index)
    elif args.file_subdataset:
        file_format_module.generate_subdata_streamlines(args.file, args.output_directory, args.file_subdataset)
    elif args.file_time:
        file_format_module.generate_time_streamlines(args.file, args.output_directory, args.file_time)
    else:
        file_format_module.generate_all_streamlines(args.file, args.process_count, args.output_directory)
