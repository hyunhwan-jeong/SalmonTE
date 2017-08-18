with open("hs_origin.fa", "r") as inp:
    seqs = inp.read().strip().split(">")[1:]

with open("annotation_hs.tsv", "r") as inp:
    clades = {}
    for row in inp:
        row = row.strip().split("\t")
        clades[row[0]] = (row[1], row[2])

with open("hs.fa", "w") as oup, open("hs.clades", "w") as oup2:
    oup2.write("name\tclass\tclade\n")
    for s in seqs:
        first_line = s.split("\n")[0].split("\t")
        assert len(first_line) == 3
        if not first_line[0] in clades:
            print(first_line)
            continue
        if clades[first_line[0]][0] == "DNA transposon":
            continue
        oup.write(">{}".format(s))
        oup2.write("{},{},{}\n".format(first_line[0],
            clades[first_line[0]][0],
            clades[first_line[0]][1]))
