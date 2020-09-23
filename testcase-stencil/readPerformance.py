import pandas as pd
import json
import os

def writeCsv(fileName):
    with open(fileName, encoding="utf-8") as file:
        dat_file = [l.rstrip("\n") for l in file]


    json_file = []
    dictLevel = 0
    for i,l in enumerate(dat_file):
        l = l.replace(';', '')
        l = l.replace('"', '')
        content = l.strip().split()
        # print(content)

        line = ""
        if len(content) == 2:
            #if "trigger" in content[0]:
            line = '"{}" : "{}",'.format(content[0],content[1])

        if len(content) == 1:
            #if "trigger" in content[0]:
            if not ("{" in content or "}" in content):
                line = '"{}":'.format(content[0])
            else:
                line = '{}'.format(content[0])
            if "}" in content:
                line = "},"


        endOfDict = False
        if len(json_file) >= 1:
            if "}" in line:
                endOfDict = True
        if endOfDict:
            json_file[-1] = json_file[-1][:-1]
        if line != "":
            json_file.append(line)

    # remove last comma
    json_file[-1] = json_file[-1][:-1]
    json_file = json_file[1:]
    print(json_file[-10:])

    with open("tmp.json","w") as f:
        for row in json_file:
            f.write(str(row) + '\n')


    with open("tmp.json") as f:
        prof_dict = json.load(f)
    os.remove("tmp.json")

    df = pd.DataFrame.from_dict(prof_dict)
    return df

df = writeCsv("profiling_stencilFraction.dat").T
df = df.set_index("description")
df.to_csv("profiling_stencilFraction_200.csv")
df = writeCsv("profiling_stencilOptimize.dat").T
df = df.set_index("description")
df.to_csv("profiling_stencilOptimize_200.csv")

