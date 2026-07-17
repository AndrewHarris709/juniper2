import pandas as pd
import matplotlib.pyplot as plt
import csv

rows = []
with open('./output//hydro32_neighbours.jrep') as f:
    for row in csv.reader(f):
        d = {'index': int(row[0]), 'keyword': row[1], 'values': set([int(x) for x in row[2:]])}
        rows.append(d)

df = pd.DataFrame(rows).set_index("index")
df.loc[df['keyword'] == "true", "tree_values"] = df.loc[df['keyword'] == "tree", 'values']
df = df.loc[df['keyword'] == "true"].copy()
df['missing'] = df.apply(lambda r: len(r['values'] - r['tree_values']), axis=1)
df['extra'] = df.apply(lambda r: len(r['tree_values'] - r['values']), axis=1)
df['matching'] = df.apply(lambda r: len(r['values'] & r['tree_values']), axis=1)
df = df.sort_values('matching').reset_index(drop=True).drop(columns=["keyword", "values", "tree_values"])

fig, ax = plt.subplots()

ax.stackplot(df.index, df['matching'], df['missing'], df['extra'], labels=["Matching Neighbours", "Missing Neighbours", "Extra Tree Neighbours"])

ax.set_title("Neighbour Analysis for hydro32_00020")
ax.set_xlabel("Particle Index (sorted by increasing neighbours)")
ax.set_ylabel("# of Neighbours")
ax.legend(loc="upper left")
ax.axhline(120)

plt.show(block=True)
