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

import helpers

def get_samples(dataset):
    """Report samples per dataset"""
    return list({ record['SampleName'] for record in helpers.read_sample_sheet(dataset) })

####################################################################################################
# INPUT HELPER FOR RULES THAT WORK ON RAW READS
####################################################################################################

def input_reads(wildcards):
    dataset = wildcards.dataset
    sample  = wildcards.sample
    trim    = wildcards.trim == 'trim'

    records = [ record for record in helpers.read_sample_sheet(dataset) if record['SampleName'] == sample ]

    fastqs_r1 = [ record['FileNameR1'] for record in records ]
    fastqs_r2 = [ record['FileNameR2'] for record in records if record['FileNameR2'] != '']

    if len(fastqs_r1) != 1:
        print(fastqs_r1)
        raise RuntimeError('Exactly one R1 reads file has to be specified for a sample')
    if len(fastqs_r2) not in [0,1] :
        raise RuntimeError('One or no R2 reads file has to be specified for a sample')

    reads = fastqs_r1 + fastqs_r2

    if trim:
        reads = [ os.path.join('trim', re.sub('fastq.gz$', 'trimmed.fastq.gz', fastq)) for fastq in reads ]
    else:
        reads = [ os.path.join('fastq', fastq) for fastq in reads ]

    return reads


####################################################################################################
# XHLA
####################################################################################################

rule xhla_conversion:
    """Converts xHLA output files to the MultiHLA standard format"""
    input:
        'typing/xhla/{dataset}_{sample}_{trim}_{ref}.xhla.json',
    output:
        'typing/xhla/{dataset}_{sample}_{trim}_{ref}.xhla.multihla',
    params:
        opts = lambda wildcards : f'version=0 trim={wildcards.trim} ref={wildcards.ref}',

        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    run:
        # Parse JSON output of xHLA
        import json
        with open(input[0], 'r') as infile:
            alleles = json.load(infile)['hla']['alleles']

        # Write multihla output
        import collections
        with open(output[0], 'w') as outfile:
            print('Dataset\tSample\tMethod\tOptions\tGene\tAllele1\tAllele2', file=outfile)
            data = collections.defaultdict(lambda : [])
            for allele in alleles:
                gene, allele = allele.split('*')
                data[gene].append(allele)
            for gene, alleles in data.items():
                if len(alleles) == 1:
                    s = allele[0] + '\tNA'
                elif len(alleles) == 2:
                    s = '\t'.join(alleles)
                else:
                    raise RuntimeError(f'More than two alleles reported for gene {gene}.')
                print(f'{wildcards.dataset}\t{wildcards.sample}\txHLA\t{params["opts"]}\t{gene}\t{s}', file = outfile)

rule hla_xhla_dataset:
    input:
        lambda wildcards : [
            f'typing/xhla/{wildcards.dataset}_{sample}_{trim}_{ref}.xhla.multihla'
            for sample in get_samples(wildcards.dataset)
            for ref in [ 'hg38.noalt' ] # We only map against hg38 w/o alt contigs
            for trim in [ 'trim' ]      # We only work with adapter trimmed reads
            ]
    output:
        'typing/xhla/{dataset}.xhla.ds.multihla'
    params:
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    script:
        '../scripts/concat_tables.py'

####################################################################################################
# HLA-LA
####################################################################################################

rule hla_la_conversion:
    """Converts HLA VBSeq output files to the MultiHLA standard format"""
    input:
        'typing/hla_la/{dataset}_{sample}_{trim}_{ref}.hla_la.txt',
    output:
        'typing/hla_la/{dataset}_{sample}_{trim}_{ref}.hla_la.multihla',
    params:
        opts = lambda wildcards : f'version=1.0.1 trim={wildcards.trim} ref={wildcards.ref}',

        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    run:
        import csv
        from collections import defaultdict
        ddict = defaultdict(lambda : set())
        with open(input[0], 'r') as infile:
            reader = csv.DictReader(infile, delimiter = '\t')
            for rec in reader:
                # extract first two fields and drop gene name
                allele = ':'.join(rec['Allele'].split(':')[:2]).split('*')[1]
                # record allele
                ddict[rec['Locus']].add(allele)
        with open(output[0], 'w') as outfile:
            print('Dataset\tSample\tMethod\tOptions\tGene\tAllele1\tAllele2', file = outfile)
            for gene, alleles in ddict.items():
                alleles = list(alleles)
                if len(alleles)==2:
                    print(f'{wildcards.dataset}\t{wildcards.sample}\tHLA-LA\t{params.opts}\t{gene}\t{alleles[0]}\t{alleles[1]}', file = outfile)
                elif len(alleles)==1:
                    print(f'{wildcards.dataset}\t{wildcards.sample}\tHLA-LA\t{params.opts}\t{gene}\t{alleles[0]}\t{alleles[0]}', file = outfile)
                else:
                    raise RuntimeError('hla_la_conversion expects either one or two alleles to be reported per gene')

rule hla_la_dataset:
    input:
        lambda wildcards : [
            f'typing/hla_la/{wildcards.dataset}_{sample}_{trim}_{ref}.hla_la.multihla'
            for sample in get_samples(wildcards.dataset)
            for ref in [ 'hg38.noalt' ] # We only map against hg19 w/o alt contigs
            for trim in [ 'trim' ]      # We only work with adapter trimmed reads
            ]
    output:
        'typing/hla_la/{dataset}.hla_la.ds.multihla'
    params:
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    script:
        '../scripts/concat_tables.py'

####################################################################################################
# HLA-VBSEQ
####################################################################################################

rule hla_vbseq_conversion:
    """Converts HLA VBSeq output files to the MultiHLA standard format"""
    input:
        'typing/vbseq/{dataset}_{sample}_{trim}_{ref}.v{version}.vbseq.tsv',
    output:
        'typing/vbseq/{dataset}_{sample}_{trim}_{ref}.v{version}.vbseq.multihla',
    params:
        opts = lambda wildcards : f'version={wildcards.version} trim={wildcards.trim} ref={wildcards.ref}',

        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    shell:
        r"""
        echo -e 'Dataset\tSample\tMethod\tOptions\tGene\tAllele1\tAllele2' >{output:q}
        sed 's/\t[^*]\+[*]/\t/g' < {input:q} | egrep '^(A|B|C|DP|DQ|DRB)' | sed 's/^/{wildcards.dataset}\t{wildcards.sample}\tHLA-VBSeq\t{params[opts]}\t/' >>{output:q}
        """

rule hla_vbseq_dataset:
    input:
        lambda wildcards : [
            f'typing/vbseq/{wildcards.dataset}_{sample}_{trim}_{ref}.v{ver}.vbseq.multihla'
            for sample in get_samples(wildcards.dataset)
            for ref in [ 'hg19.noalt' ] # We only map against hg19 w/o alt contigs
            for trim in [ 'trim' ]      # We only work with adapter trimmed reads
            for ver in ['1', '2']       # We perfom the analysis with VBSeq version 1 and 2
            ]
    output:
        'typing/vbseq/{dataset}.vbseq.ds.multihla'
    params:
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    script:
        '../scripts/concat_tables.py'

####################################################################################################
# OPTITYPE
####################################################################################################

rule optitype_conversion:
    """Converts OptiType output files to the MultiHLA standard format"""
    input:
        'typing/optitype/{dataset}_{sample}_{trim}_{ref}_result.tsv',
    output:
        'typing/optitype/{dataset}_{sample}_{trim}_{ref}.optitype.multihla',
    params:
        opts = lambda wildcards : f'version=1.3.4 trim={wildcards.trim} filt={wildcards.ref}',
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    run:
        import csv
        from collections import defaultdict
        ddict = defaultdict(lambda : set())
        with open(input[0], 'r') as infile:
            reader = csv.DictReader(infile, delimiter = '\t')
            for rec in reader:
                if ddict:
                    raise RuntimeError('optitype_conversion requires the result file to contain exactly one record')
                ddict['A'] = { rec['A1'].split('*')[1], rec['A2'].split('*')[1] }
                ddict['B'] = { rec['B1'].split('*')[1], rec['B2'].split('*')[1] }
                ddict['C'] = { rec['C1'].split('*')[1], rec['C2'].split('*')[1] }
        with open(output[0], 'w') as outfile:
            print('Dataset\tSample\tMethod\tOptions\tGene\tAllele1\tAllele2', file = outfile)
            for gene, alleles in ddict.items():
                alleles = list(alleles)
                print(f'{wildcards.dataset}\t{wildcards.sample}\tHLA-LA\t{params.opts}\t{gene}\t{alleles[0]}\t{alleles[1]}', file = outfile)

rule optitype_dataset:
    input:
        lambda wildcards : [
            f'typing/optitype/{wildcards.dataset}_{sample}_{trim}_{filt}.optitype.multihla'
            for sample in get_samples(wildcards.dataset)
            for filt in [ 'nofilt' ] # No prior read filtering against hg19 or hg38
            for trim in [ 'trim' ]  # We only work with adapter trimmed reads
            ]
    output:
        'typing/optitype/{dataset}.optitype.ds.multihla'
    params:
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    script:
        '../scripts/concat_tables.py'


####################################################################################################
# GLOBAL COLLECTOR
####################################################################################################

rule collect:
    input:
        'typing/vbseq/{dataset}.vbseq.ds.multihla',
        'typing/xhla/{dataset}.xhla.ds.multihla',
        'typing/hla_la/{dataset}.hla_la.ds.multihla',
        'typing/optitype/{dataset}.optitype.ds.multihla',
    output:
        'typing/{dataset}.all.multihla'
    params:
        # Parameters for cluster execution
        cluster_mem = '1G',
        cluster_rt = '0:15:00',
    script:
        '../scripts/concat_tables.py'
