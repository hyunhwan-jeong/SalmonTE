with open("dm_origin.fa", "r") as inp:
    seqs = inp.read().strip().split(">")[1:]

with open("clades.tsv", "r") as inp:
    clades = {}
    for row in inp:
        row = row.strip().split("\t")
        clades[row[0]] = (row[1], row[2])

with open("dm.fa", "w") as oup, open("dm.clades", "w") as oup2:
    oup2.write("name,class,clade\n")
    for s in seqs:
        first_line = s.split("\n")[0].split("\t")
        assert len(first_line) == 3
        if first_line[-1] != "Drosophila melanogaster":
            continue
        if clades[first_line[1]][0] == "DNA transposon":
            print(first_line)
            continue
        oup.write(">{}".format(s))
        oup2.write("{},{},{}\n".format(first_line[0],
            clades[first_line[1]][0],
            clades[first_line[1]][1]))
