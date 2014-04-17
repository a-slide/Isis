from random import gauss

def illumina_qualgen (m_start, m_mid, m_end, sd_start, sd_mid, sd_end, length):
    """"""

    start_border = int(0.2*length)
    mid_border = int(0.66*length)
    end_border = length
    mean_qual = []
    sd_qual =[]
    
    for i in range (start_border) :
        mean_qual.append    (int(m_start + (float (m_mid - m_start))/start_border*i))
        sd_qual.append      (int(sd_start + (float (sd_mid - sd_start))/start_border*i))

    for i in range (mid_border - start_border) :
        mean_qual.append(m_mid)
        sd_qual.append(sd_mid)

    for i in range(end_border - mid_border):
        mean_qual.append    (int(m_mid + (float (m_end - m_mid)/ (end_border - mid_border)*i)))
        sd_qual.append      (int(sd_mid + (float (sd_end - sd_mid)/ (end_border - mid_border)*i)))

    return (mean_qual, sd_qual)

def random_qual_string (mean_qual, sd_qual, length):
    return [valid_qual(mean_qual[i], sd_qual[i]) for i in range (length)]

def valid_qual(mean, sd):
    while True:
        qual = int(gauss(mean, sd))
        if 0 <= qual <= 40:
            return qual


mu, sd = illumina_qualgen(30, 35, 15, 2, 1, 10, 100)
h = plt.bar(range(100), random_qual_string(mu, sd, 100))
plt.show()
