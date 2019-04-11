import matplotlib.pyplot as pl
import seaborn as sns
import pandas as pd

df = pd.read_csv("./result/ica_all_qc3/ica_cluster_size.csv",header=None)
labels = df.index.values
sizes = df[0]
fig1, ax1 = pl.subplots()
ax1.pie(sizes, autopct='%1.1f%%',shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
pl.savefig("./result/ica_all_qc3/pie.png")
