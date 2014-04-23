from ../QualGenerator import QualGenerator
from matplotlib import pyplot

def qtester(length, quality, iteration, title):

    q = QualGenerator(length, quality)
    fig = pyplot.figure(figsize=(length/2, 10), dpi=100)
    pyplot.title(title)
    pyplot.ylabel('PHRED Quality')
    pyplot.xlabel('Position')
    position = [i+1 for i in range (length)]
    
    for i in range (iteration):
        q_list = q.random_qual_string()
        pyplot.plot(position, q_list)

    fig.savefig(title+'.png')

