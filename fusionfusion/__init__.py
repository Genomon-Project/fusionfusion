#! /usr/bin/env python

from .arg_parser import create_parser
from .run import fusionfusion_main 

def main():

    parser = create_parser()
    args = parser.parse_args()
    fusionfusion_main(args)

