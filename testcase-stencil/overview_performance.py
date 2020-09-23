import pandas as pd
import os

prof_stencilLoop = pd.DataFrame()
for f in os.listdir():
    if "profiling_stencilFraction" in f:
        name = os.path.splitext(f)[0]
        print(name)
        cellNo = name.split("_")[2]
        df = pd.read_csv(f)
        df["iterTime"] = df["totalTime"] - df["childTime"]
        df["nCells"] = int(cellNo)**3
        prof_stencilLoop = prof_stencilLoop.append(df,ignore_index=True)

prof_stencilLoop = prof_stencilLoop.dropna()
prof_stencilLoop = prof_stencilLoop[prof_stencilLoop["parentId"] == 0]
prof_stencilLoop= prof_stencilLoop.drop(['id', 'calls', 'active',  'parentId'], axis=1)
prof_stencilLoop = prof_stencilLoop.rename(columns={"childTime":"build [s]","iterTime":"iterate [s]"})
prof_stencilLoop = prof_stencilLoop[prof_stencilLoop["description"].str.contains("partial")]
prof_stencilLoop["description"] = prof_stencilLoop["description"].str.replace("partial","")
print(prof_stencilLoop)
prof_stencilLoop = prof_stencilLoop.round(4)

print(prof_stencilLoop.to_markdown(index=False))