#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv

def main ():
    test("test")

def Pair(text):
    print text+" Pair"

def Single(text):
    print text+" Single"


if __name__ == '__main__':
    main()
    
if __name__ == '__test__':
    print "test"
    if argv[1] == "pe":
        Pair()
    elif argv[1] == "se":
        Single()
    
