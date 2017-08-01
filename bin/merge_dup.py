import sys


def get_tpm(fname):
    from collections import defaultdict
    tpm = defaultdict( float )
    with open(fname, "r") as inp:
        line = inp.readline()
        for line in inp:
            line = line.strip().split()
            name = "_".join(line[0].split("_")[:-1])
            if len(name) == 0:
                name = line[0]
            tpm[name] += float(line[3])

    return tpm

def main():
    import glob
    tb = dict()
    for s in sorted(glob.glob("*_SalmonRepBase/quant.sf")):
        sid = s.split("_")[0]
        tb[sid] = get_tpm(s)

    import pandas as pd
    with open("Salmon_Repbase.out", "w") as oup:
        oup.write(pd.DataFrame(tb).to_csv(sep="\t"))


if __name__ == "__main__":
    main()
