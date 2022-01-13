#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 07:32:11 2022

@author: buehler
"""

import sys
from pprint import pprint

import numpy as np
from loguru import logger

from uvva.field import BaseField, Field, GALEXField


def main():
    # logging
    # replacing the standard handler
    # logger.remove(0)
    # setup custom handler to make use of
    # additional fields defined by uvva
    fmt = (
        "<g>{time:YYYY-MM-DD HH:mm:ss.SSS}</> | "
        "<lvl>{level: ^8}</> | "
        "<c>{name}</>:<c>{line}</> in <c>{extra[classname]}</>:<c>{function}</> | "
        "<lvl>{message}</>"
    )
    logger.add(sys.stdout, format=fmt)
    # activate logging (is deactivated by import of uvva)
    logger.enable("uvva")

    # test BaseField
    ff = BaseField()
    print("\n******** INFO *********")
    ff.info()

    print("\n******** PRINT ********")
    print(ff)

    # test GALEX use case
    gf = GALEXField(parobs_id=42, field_name="gTest", ra=np.float64(5))
    print(gf.tt_field.info)
    pprint(gf.tt_field.meta)

    # # tt_field is the old tt_coadd
    gf = Field(6388191295067652096)
    print("GALEX field info:")


if __name__ == "__main__":
    main()
