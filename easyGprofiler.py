

import requests
import pandas as pd
import argparse
from numpy import log10
import matplotlib.pyplot as plt

""" Small program to download gProfiler term enrichment analysis"""


def args():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile',
                        help='File with the list of genes to analyse')
    parser.add_argument('organism',
                        help='''Organism name according to gProfiler list.
                        Example: mmusculus, hsapiens, ggallus,
                        dmelanogaster.''')
    parser.add_argument('--output', '-o', 
    	help="""Output file""")
    parser.add_argument('--method', '-m', default='g_SCS',
    	help="""False discovery method: 'g_SCS', 'bonferroni' and 'fdr'""")
    parser.add_argument('--threshold', '-t', default=0.05,
    	help="""p-value threshold""")
    parser.add_argument('--underrepresented', '-u', default=False, action='store_true',
    	help="""If specified, returns the under represented terms (isntead of the over represented) """)
    return parser.parse_args()


def gprofiler(namelist, organism, user_threshold=0.05, method='g_SCS',
              measure_underrepresentation=False, background=None,
              simple_out=0):
    """Run gProfiler using POST api with a json query body

    Returns a pandas DataFrame with the result

    methods = 'g_SCS', 'bonferroni' and 'fdr'
    simple_out = 3 levels of output 0|1|2
    background = TO DO!!!!!!!!
    """
    if type(namelist) is not list:
        namelist = list(namelist)
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'organism': organism,
            'query': namelist,
            'user_threshold': user_threshold,
            'measure_underrepresentation': measure_underrepresentation,
            'significance_threshold_method': method,
               },
        timeout=10
                    )
    df = pd.DataFrame(r.json()['result'])
    # failed genes
    failed = r.json()['meta']['genes_metadata']['failed']
    # succes genes
    genesdict = r.json()['meta']['genes_metadata']['query']['query_1']['mapping']
    ensdict = {v[0]: id_ for id_, v in genesdict.items()}
    enslist = r.json()['meta']['genes_metadata']['query']['query_1']['ensgs']
    geneslist = [ensdict[ens] for ens in enslist]
    # dataframe of genes and onthologies
    genesdf = pd.DataFrame(0, index=df['native'], columns=geneslist)
    for i in range(len(df)):
        intersections = df['intersections'][i]
        for j, int_ in enumerate(intersections):
            if int_:            # only not empty values are relevant
                genesdf.iloc[i, j] = 1

    # manually selecting and sorting output
    cols = ['source', 'native', 'name', 'p_value', 'description', 'query',
            'significant', 'term_size', 'intersection_size']
    extend = ['query_size', 'effective_domain_size',
              'intersections', 'parents']
    # Unkown columns, there is no direct explanation.
    wtf_cols = ['goshv', 'group_id', 'precision', 'recall',
                'source_order']
    assert simple_out in [0, 1, 2], 'Simple_out wrong value'
    if simple_out == 0:
        df = df[cols]
    elif simple_out == 1:
        df = df[cols + extend]
    elif simple_out == 2:
        df = df[cols + extend + wtf_cols]
    return df, genesdf, failed


if __name__ == '__main__':
    args = args()

    if args.output == None:
    	args.output = args.infile + '.gprofiler'

    query = []

    with open(args.infile) as inf:
    	for line in inf:
    		query.append(line.strip())

    print('Arguments')
    print(args)

    result, genes, failed = gprofiler(query, args.organism,
                                      user_threshold=args.threshold,
                                      method=args.method,
                                      measure_underrepresentation=args.underrepresented,
                                      simple_out=0)
    # save results
    result.to_csv(args.output, sep='\t')
    genes.to_csv(args.output + '.genes', sep='\t')
    with open(args.output + '.failed', 'w') as outf:
        outf.write('\n'.join(failed))


    print('Some plots')
    sources = result.source.unique()

    result = result.set_index('name')

    for source in sources:
        subdf = (-log10(result[result.source == source].p_value))
        subdf = subdf.iloc[:10]
        y = 0.2 * len(subdf)
        plt.figure(figsize=(5,y))
        # subdf.plot.barh(color=source_colors[source])
        subdf.plot.barh()
        plt.gca().invert_yaxis()
        plt.xlabel('-log10(p-value)')
        plt.title(source)
        plt.savefig(args.output + source + '.svg')


    print('[DONE] Be happy!!!')

    # # get profile
    # genes = []
    # organism = 'mmusculus'
    # profile = gprofiler(genes, organism)
    # # Total terms detected for all sources
    # print('Total detected p-val <= 0.05:')
    # print(f'{len(golargest)}')
    # # number of terms detected for each source
    # print('\nDetected by source: ')
    # print(golargest.source.value_counts())
    # # plot top 10 pval for each source
    # # source = 'GO:MF'
