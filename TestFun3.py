from QualGenerator import QualGenerator

q = QualGenerator(100,30,35,37,35,25,5,2,2,5,10,0,0.15,0.25,0.25,0.35)

with (open("test.csv", "w")) as file:
    for i in range (10000):
        s =""
        for i in q.random_qual_string(100):
            s += "{};".format(i)
        file.write(s+"\n")
