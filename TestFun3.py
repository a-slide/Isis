from QualGenerator import QualGenerator
q = QualGenerator(100, "bad")

with (open("test.csv", "w")) as file:
    for i in range (10):
        s =""
        for i in q.random_qual_string():
            s += "{};".format(i)
        file.write(s+"\n")
