# Copyright (c) 2020 Leon Kuchenbecker <leon.kuchenbecker@uni-tuebingen.de>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Compares typing results to the ground truth

import pandas as pd

# Read the ground truth data
gtruth = pd.read_csv(snakemake.input.gtruth, '\t')
gtruth['Dataset'] = snakemake.wildcards.dataset
gtruth['Method'] = 'ground_truth'
gtruth['Options'] = ''

# Read the typing results
typing = pd.read_csv(snakemake.input.typing, '\t')

# Concatenate and pivot the data
compare = pd.concat([gtruth, typing])

# Create result column with set of called alleles per gene and drop original, individual columns
compare['Result'] = compare.apply(axis='columns', func=lambda x : {x.Allele1, x.Allele2})
compare.drop(columns=['Allele1','Allele2'], inplace = True)

# Pivot the table to have one column per method, drop constant column index level
compare = compare.pivot_table(index=['Dataset','Sample','Gene'], columns=['Method','Options'], aggfunc=lambda x : set.union(*x))
compare = compare.droplevel(0, axis='columns')

# Make the ground truth the first column
compare = compare[ [c for c in compare.columns if c[0]=='ground_truth'] + [c for c in compare.columns if c[0]!='ground_truth'] ]

def annotate(row, method):
    if pd.isna(row[method]) or len(row[method])==0:
        return 'no result'
    gt = row[('ground_truth', '')]
    if row[method] == gt:
        return 'correct'
    if pd.isna(gt) or len(gt)==0:
        return f'no ground truth: {row[method]}'
    if row[method] & gt: # Intersection
        return f'partial match: {row[method]}'
    else:
        return f'mismatch: {row[method]}'

# Create copy with annotated version of results
annotated = compare.copy()
for col in [ c for c in annotated.columns if c[0] != 'ground_truth' ]:
    annotated[col] = annotated.apply(lambda x : annotate(x, col), axis = 'columns')

# Write results
compare.to_csv(snakemake.output.raw, '\t', index = True)
annotated.to_csv(snakemake.output.annotated, '\t', index = True)
