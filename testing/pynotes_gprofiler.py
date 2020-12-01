import gprofiler
from gprofiler import GProfiler
GProfiler?
gp = GProfiler(return_dataframe=True)
gp.profile(organism='mmusculus', query=genes)
genes
genes = """ENSMUSG00000076488
ENSMUSG00000065231
ENSMUSG00000079120
ENSMUSG00000047222
ENSMUSG00000097494
ENSMUSG00000064419
ENSMUSG00000095668
ENSMUSG00000059606""".split()
gp.profile(organism='mmusculus', query=genes)
import requests
def mygprofiler(namelist, organism='mmusculus'):
    """Run gProfiler using POST api with a json query body
    
    Returns a pandas DataFrame with the result"""
    if type(namelist) is not list:
        namelist = list(namelist)
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'organism':organism,
            'query': namelist,
        }
        )
    df = pd.DataFrame(r.json()['result'])
    return df
myres = mygprofiler(genes)
import pandas as pd
myres = mygprofiler(genes)
res = gp.profile(organism='mmusculus', query=genes)
res
myres
res.columns
myres.columns
myres.goshv
res.pvalue
res.pval
res
res.columns
res.p_value
myres.p_value
res.significant
myres.significant
gp.profile?
myres.columns
gp.convert(organism='mmusculus', query=genes)
gp.convert(organism='mmusculus', query=genes, target_namespace='name')
gp.convert?
gp.convert(organism='mmusculus', query=genes, target_namespace='name')
gp.convert(organism='mmusculus', query=genes, )
gp.convert(organism='mmusculus', query=genes, ).namespaces
gp.convert(organism='mmusculus', query=genes, target_namespace='ENSG' )
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'UCSC',
        'query':genes,
    }
    )
x = r.json()
x
x.keys()
x['result']
df(x['result'])
pd.DataFrame(x['result'])
res
pd.DataFrame(x['result'])
pd.DataFrame(x['result']).columns
pd.DataFrame(x['result']).name
pd.DataFrame(x['result']).name
pd.DataFrame(x['result']).name
gp.convert(organism='mmusculus', query=genes, target_namespace='ENSG' )
gp.convert(organism='mmusculus', query=genes )
gp.convert(organism='mmusculus', query=genes ).columns
df.sort_values('prob1', ascending=False)
pd.DataFrame(x['result']).colymns
pd.DataFrame(x['result']).columns
pd.DataFrame(x['result']).columns
gp.convert(organism='mmusculus', query=genes ).columns
gp.convert(organism='mmusculus', query=genes )
pd.DataFrame(x['result'])
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'EMBL',
        'query':genes,
    }
    )
pd.DataFrame(x['result'])
x = r.json()
pd.DataFrame(x['result'])
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'refseq',
        'query':genes,
    }
    )
x = r.json()
pd.DataFrame(x['result'])
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'RefSeq',
        'query':genes,
    }
    )
x = r.json()
pd.DataFrame(x['result'])
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'Ensembl',
        'query':genes,
    }
    )
x = r.json()
pd.DataFrame(x['result'])
genes
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'Uniprot',
        'query':genes,
    }
    )
x = r.json()
pd.DataFrame(x['result'])
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',

        'query':genes,
    }
    )
x = r.json()
x
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'UCSC',
        'query':genes,
    }
    )
x = r.json()
df.sort_values('prob1', ascending=False)
pd.DataFrame(x['result'])
gp.convert(organism='mmusculus', query=genes )
gp.convert?
gp.convert?
gp.convert(organism='mmusculus', query=genes , target_namespace='Ensembl')
gp.convert(organism='mmusculus', query=genes , target_namespace='RefSeq')
gp.convert(organism='mmusculus', query=genes , target_namespace='UCSC')
gp.convert(organism='mmusculus', query=genes , target_namespace='Uniprot')
gp.convert(organism='mmusculus', query=genes , target_namespace='HUGO')
gp.convert(organism='mmusculus', query=genes , target_namespace='IPI')
gp.convert(organism='mmusculus', query=genes , target_namespace='ENSG')
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'mmusculus',
        'target':'ENSG',
        'query':genes,
    }
    )
x = r.json()
pd.DataFrame(x['result'])
gp.snpense(query=genes)
gp.snpense(query=['rs11734132', 'rs7961894', 'rs4305276', 'rs17396340'])
gp.orth(organism='mmusculus', query=genes, target='chicken')
gp.orth(organism='mmusculus', query=genes, target='human')
gp.orth(organism='mmusculus', query=genes, target='msapiens')
gp.orth(organism='mmusculus', query=genes, target='hsapiens')
gp.orth(organism='mmusculus', query=genes, target='hsapiens')
ls
%history -f pynotes_gprofiler.py
