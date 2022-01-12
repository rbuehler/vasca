#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 07:32:11 2022

@author: buehler
"""
from uvva import field


def main():
    ff = field.BaseField()
    print("\n******** INFO *********")
    ff.info()
    print("\n******** PRINT ********")
    print(ff)


if __name__ == "__main__":
    main()
