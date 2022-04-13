from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os
import pandas as pd

def getProfilingData():
    profData = pd.DataFrame()
    dir_name = os.path.dirname(os.path.abspath(__file__))
    print(dir_name)
    proc_dirs = [ f.path for f in os.scandir(dir_name) if "processor" in f.name]
    for proc in proc_dirs:
        print(proc)
        d = ParsedParameterFile(f"{proc}/0/uniform/profiling")
        proc_data = pd.DataFrame(d["profiling"]).T
        proc_data["processor"] = proc[2:]
        profData = profData.append(proc_data)
    return profData


df = getProfilingData()
masks = df['description'].str.contains('calcCellStencil')
print(df[masks])