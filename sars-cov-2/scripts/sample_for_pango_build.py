"""
Needs to run as local rule for internet
Input:
- pre-processed/open_pango_metadata.tsv
- https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json
Output:
- builds/pango_subsample.tsv
"""
#%%
import pandas as pd
from collections import defaultdict

#%%
N_PER_LINEAGE = 3
N_PER_CHILD = 2
#%%
df = pd.read_csv('pre-processed/open_pango_metadata.tsv', sep='\t')
#%%
# Dealias the pango column
# Calculate number of descendants per pango
# Sample specified number of lineages
# Potentially sample earliest samples of lineage

#%%
class Aliasor:
    def __init__(self):
        import pandas as pd

        aliases = pd.read_json('pre-processed/alias.json')

        self.alias_dict = {}
        for column in aliases.columns:
            if column.startswith('X'):
                self.alias_dict[column] = column
            else:
                self.alias_dict[column] = aliases[column][0]

        self.alias_dict['A'] = 'A'
        self.alias_dict['B'] = 'B'

        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self,name):
        name_split = name.split('.')
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            return letter + '.' + name_split[4]
        else:
            return letter + '.' + ".".join(name_split[4:])

    def uncompress(self,name):
        name_split = name.split('.')
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            return unaliased + '.' + name_split[1]
        else:
            return unaliased + '.' + ".".join(name_split[1:])

aliasor = Aliasor()
#%%
lineages = df['pango_designated'].unique()
#%%
lineages_unaliased = [aliasor.uncompress(lineage) for lineage in lineages]
# lineages_unaliased
# %%
# Count parents
children_count = defaultdict(int)

for lineage in lineages_unaliased:
    lineage_split = lineage.split('.')
    if len(lineage_split) > 1:
        parent = ".".join(lineage_split[0:-1])
        children_count[parent] += 1

# children_count

#%%
children_count_aliased = {aliasor.compress(k): v for k, v in children_count.items()}
# children_count_aliased
#%%
samples_per_lineage = df.groupby('pango_designated').size()
#%%

n_sample = { k: N_PER_LINEAGE for k in lineages}
for parent in children_count_aliased:
    if parent in n_sample:
        n_sample[parent] += N_PER_CHILD * children_count_aliased[parent]

for lineage in n_sample:
    n_sample[lineage] = min(n_sample[lineage], samples_per_lineage[lineage])

# sum(n_sample.values())
#%%
samples = [df.strain[df.strain.str.contains('Wuhan/Hu-1/2019')]]
for name, group in df.groupby('pango_designated'):
    samples.append(group.strain.sample(n=n_sample[name], random_state=0))

#%%
subsample = pd.concat(samples, ignore_index=True)
subsample.to_csv('builds/pango_subsample.tsv', sep='\t', index=False)
# %%