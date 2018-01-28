#!/usr/bin/env python

import re

regRe = re.compile(r'([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')

def getCoverSize(cover_str):
    """
    get the size of covered bases from cover strings (e.g., chr1:741208-741271,chr1:745447-745482)
    """
    size = 0
    F = cover_str.split(",")
    for i in range(len(F)):
        mReg = regRe.match(F[i])
        chr = mReg.group(1)
        start = int(mReg.group(2))
        end = int(mReg.group(3))
        size = size + end - start + 1

    return size

