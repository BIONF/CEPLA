"""
Created on: 03/06/2017

@author: Dorjee Gyaltsen
@brief: calculates population coverage - standalone version
"""

import os
import sys

from argparse import ArgumentParser

from population_calculation import PopulationCoverage
from util import *


def flatten_list(lis_o_lis):
    return [val for lis in lis_o_lis for val in lis]


def calculate(args):
    population = flatten_list(args.population) if args.population else []
    population = filter(None, population)
    mhc_class = flatten_list(args.mhc_class) if args.mhc_class else []
    mhc_class = filter(None, mhc_class)

    if args.list:
        population_class_args = []
        if population:
            population_class_args.extend(population)
        if mhc_class:
            population_class_args.extend(mhc_class)

        if not population_class_args:
            get_population_list()
        else:
            get_available_allele_names(*population_class_args)
        sys.exit(0)
    else:
        _prog = ArgumentParser().prog
        if not population or not mhc_class:
            print("arguments -p/--population and -c/--mhc_class are required\n{}\n"
                  "\nfor detail usage: $ python {} --help".format(msg(), _prog))
            sys.exit(0)
        elif not args.filename:
            print("argument -f/--file: expected one argument\n{}\n"
                  "\nfor detail usage: $ python {} --help".format(msg(), _prog))
            sys.exit(0)

        pcal = PopulationCoverage()
        result, negative, every_epitope = pcal.calculate_coverage(population=population, mhc_class=mhc_class, filename=args.filename) #Modified
        """
        print_chart_table(result)

        if args.path:
            from pathlib import Path
            path = Path(args.path)
            path.mkdir(parents=True, exist_ok=True)
            generate_plot(result, args.path)

        if negative:
            print("note: data for following combinations are not available, and therefore skipped")
            for i, _neg in enumerate(negative):
                print("{} - '{}, {}'".format(i+1, _neg.get("mhc_class"), _neg.get("population")))
        """
        return every_epitope #Modified


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def msg():
    _prog = ArgumentParser().prog
    usage = """python {} [-h] -p [POPULATION] -c [MHC_CLASS] -f [FILE]""".format(_prog)

    return usage


def csv_list(s):
    return s.split(",")


def get_parser():
    """Get parser object for script calculate_population_coverage.py."""

    parser = ArgumentParser(description=__doc__, usage=msg())

    req_argument = parser.add_argument_group('required arguments')

    req_argument.add_argument("-p", "--population",
                        dest="population",
                        nargs="+",
                        type=csv_list,
                        help="select comma-separated area(s) or population(s)")

    req_argument.add_argument("-c", "--mhc_class",
                        dest="mhc_class",
                        nargs="+",
                        type=csv_list,
                        help="select one or more comma-separated mhc class option - I, II, combined")

    req_argument.add_argument("-f", "--file",
                        dest="filename",
                        type=lambda x: is_valid_file(parser, x),
                        metavar="FILE",
                        help="a file containing a list of epitopes and associated alleles (comma-separated)",)

    parser.add_argument("--list",
                        action="store_true",
                        help="list all population and ethnicity")

    parser.add_argument("--plot",
                        dest="path",
                        help="generate a plot.")

    parser.add_argument('--version', action='version', version='%(prog)s v1.0')

    return parser

""" EDITED PART AHEAD """

def calculate_dir(args):
    if not os.path.isdir(args.filename):
        calculate(args)
    else:
        all_files = [file for file in os.listdir(os.path.join(args.filename, "mhc_i")) if os.path.isfile(os.path.join(args.filename, "mhc_i", file))]
        all_files += [file for file in os.listdir(os.path.join(args.filename, "mhc_ii")) if os.path.isfile(os.path.join(args.filename, "mhc_ii", file))]
        directory = args.filename
        complete_data = {}
        for file in all_files:
            mhc_class = file[file.index("_mhc_") + 1: file.index(".txt")]
            if mhc_class == "mhc_i":
                args.mhc_class = [["I"]]
            else:
                args.mhc_class = [["II"]]
            args.filename = os.path.join(directory, mhc_class, file)
            coverage_data = calculate(args)
            seq = file[:file.index("_mhc_")]

            if seq not in complete_data:
                complete_data[seq] = [i.strip().split("\t") for i in open(os.path.join(directory, seq + ".txt"), "r").readlines()]
        
            for pep in coverage_data:
                for line in complete_data[seq]:
                    if len(line) < 10 and line[4] == mhc_class and line[7] == pep:
                        line.append(coverage_data[pep])
            
        
        for seq in complete_data:
            out_file = open(os.path.join(directory, "modified", seq + ".txt"), "a")
            for line in complete_data[seq]:
                if len(line) < 10:
                    line.append("")
                out_file.write(str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) + "\t" + str(line[6]) + "\t" + str(line[7]) + "\t" + str(line[8]) + "\t" + str(line[9]) + "\n")
                

""" EDITED PART FINISH """



if __name__ == "__main__":
    args = get_parser().parse_args()
    calculate_dir(args)
