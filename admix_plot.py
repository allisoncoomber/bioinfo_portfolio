#!/usr/bin/python3
  
import pandas as pd
import seaborn as sns

ped = 'names.ped'
q = 'ch4.effectome.noMM2020.plink.2.Q'

df_ped = pd.read_csv(ped, sep=' ', header=None)
df_q = pd.read_csv(q, sep=' ', header=None)


names = ["pop{}".format(i) for i in range(1, df_q.shape[1]+1)]
df_q.columns = names
df_q.insert(0, 'Sample', df_ped[0])
df_q.set_index('Sample', inplace=True)
df_q['assignment'] = df_q.idxmax(axis=1)

pal = sns.color_palette("tab10")


def sort_df_by_pops_nocat(df):
    temp_dfs = []
    for pop in sorted(df['assignment'].unique()):
        temp = df.loc[df['assignment'] == pop].sort_values(by=[pop], ascending=False)
        temp_dfs.append(temp)
    return temp_dfs
sub_dfs = sort_df_by_pops_nocat(df_q)

df_custom_sort = pd.concat([sub_dfs[1], sub_dfs[0]])

ax = df_custom_sort.plot.bar(stacked=True, 
                             figsize=(25,5), 
                             width=1,
                             color=pal, 
                             fontsize='x-small',
                             edgecolor='black', 
                             linewidth=0.5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticklabels(df_custom_sort.index, rotation=90, ha='center')
ax.legend(bbox_to_anchor=(1,1), fontsize='medium', labelspacing=0.5, frameon=False)

ax.figure.savefig('Admixture-K2-effectome.pdf', bbox_inches='tight')