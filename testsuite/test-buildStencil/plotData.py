import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt

def trim(val):
    val = val.split("::")[0]
    return val

profData = pd.read_csv("profDataTotal.csv")
print(profData)
profData["description"] = profData["description"].apply(trim)
sns.barplot(x="description", y="totalTime", hue="Res", data=profData)

plt.show()