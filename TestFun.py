#!/usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib import pyplot
from random import gauss

fig = pyplot.figure(figsize=(20, 10), dpi=100)
pyplot.title("Distribution of sonication fragment length")
pyplot.ylabel('Count')
pyplot.xlabel('Size of fragment')

list = [gauss(700, 30) for i in range (10000)]
pyplot.hist(list, bins = 100, normed=1, facecolor='green', alpha=0.5, align='mid')

# Tweak spacing to prevent clipping of ylabel
pyplot.subplots_adjust(left=0.15)

