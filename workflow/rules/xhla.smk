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

import os
import re
import shlex
from csv import DictReader

def read_sample_sheet(dataset):
    try:
        with open(os.path.join('datasets', dataset + '.tsv'), 'r') as sample_sheet:
            reader = DictReader(sample_sheet, delimiter = '\t')
            return [ record for record in reader  ]
    except FileNotFoundError as e:
        raise RuntimeError(f'Could not find a sample sheet for the dataset \'{dataset}\'\n{e}.')

def bwamem_input(wildcards):
    dataset = wildcards.dataset
    sample  = wildcards.sample

    records = [ record for record in read_sample_sheet(dataset) if record['SampleName'] == sample ]

    fastqs = [ record['FileNameR1'] for record in records ] + [ record['FileNameR2'] for record in records ]
    fastqs = { 'fastq_' + str(idx) : os.path.join('trim', re.sub('fastq.gz$', 'trimmed.fastq.gz', fastq)) for idx, fastq in zip(range(len(fastqs)), fastqs) }

    return fastqs

def bwamem_params(input):
    # We extract all input keys that start with fastq_
    fastqs = [ input[key] for key in input.keys() if key.startswith('fastq_') ]

    # The list of FASTQ files is expected to contain the R1 fastq files first
    # and then the R2 fastq files in the corresponding order. Hence we split
    # the list in two even halfs. We return two strings, one for all R1 files
    # and one for all R2 files, each providing a space separated list of shell
    # escaped paths.
    if not len(fastqs) or len(fastqs) % 2 != 0:
        raise RuntimeError('None or uneven number of FASTQ files supplied')
    r1_fastqs = ' '.join([ shlex.quote(path) for path in fastqs[:len(fastqs)//2] ])
    r2_fastqs = ' '.join([ shlex.quote(path) for path in fastqs[len(fastqs)//2:] ])
    return r1_fastqs, r2_fastqs

rule xhla_bwamem:
    # Read mapping with BWQ mem default params as instructed by xhla authors
    conda:
        "../envs/bwa.yaml"
    input:
        unpack(bwamem_input)
    output:
        bam = 'typing/xhla/{dataset}_{sample}.bam',
        bai = 'typing/xhla/{dataset}_{sample}.bam.bai',
    log:
        bwa = 'typing/xhla/{dataset}_{sample}.bwa.log',
        sort = 'typing/xhla/{dataset}_{sample}.sort.log',
        index = 'typing/xhla/{dataset}_{sample}.index.log',
    params:
        r1_fastqs = lambda wildcards, input : bwamem_params(input)[0],
        r2_fastqs = lambda wildcards, input : bwamem_params(input)[1],

        # Parameters for cluster execution
        cluster_mem = '16G',
        cluster_rt = '4:00:00',
    threads:
        6
    shell:
        """
        # Read mapping using BWA
        bwa mem -t {threads} \
                ./ref/hg38.main.fa \
                <(unpigz -c {params[r1_fastqs]}) \
                <(unpigz -c {params[r2_fastqs]}) \
                2>{log.bwa:q} | samtools sort -T /share/scratch/kuchenb.tmp/samtools.$HOSTNAME.$$ -o {output.bam} &>{log.sort:q}

        # BAM Indexing
        samtools index -@ {threads} {output.bam:q} &>{log.index:q}
        """

rule xhla_typing:
    container:
        'docker://humanlongevity/hla'
    input:
        bam = '{base}.bam',
        bai = '{base}.bam.bai',
    output:
        json = '{base}.xhla',
        tdir = temp(directory('{base}.xhla.workdir')),
        bindir = temp(directory('{base}.xhla.bindir')),
    log:
        '{base}.xhla.log',
    params:
        # Parameters for cluster execution
        # Some times xHLA just leaks infinite amounts of memory and will be
        # killed by the OOM killer or grid engine. Oh well...
        cluster_mem = '200G',
        cluster_rt = '2:00:00',
    threads:
        6
    shell:
        """
        mkdir -p {output.bindir:q}
        mkdir -p {output.tdir:q}

        # This is a hack to circumvent that xHLA does not allow for any control
        # over the number of threads it consumes - it will always consume
        # everything the machine reports. Two components seem to use
        # multithreading, the invocation of diamond and the lpsolver in R. Here
        # we fix the diamond issue by integrating a wrapper script called
        # 'diamond' in the PATH that sets the max threads paramter of diamond
        # accordingly. Currently there is no fix for the R code, which makes
        # this rule vulnerable to cluster engines that enforce the SMP limit.

        echo '#!/bin/bash' > {output.bindir:q}/diamond
        echo 'shift' >> {output.bindir:q}/diamond
        echo '/usr/bin/diamond blastx -p {threads} "$@"' >> {output.bindir:q}/diamond
        chmod +x {output.bindir:q}/diamond
        export PATH="{output.bindir}:$PATH"

        export LC_ALL=C
        ecode=0
        if python /opt/bin/run.py --sample_id sample.$HOSTNAME.$$ --input_bam_path {input.bam:q} --output_path {output.tdir:q} &>{log:q}
        then
            mv -v {output.tdir:q}/report-sample.$HOSTNAME.$$-hla.json {output.json:q} &>>{log:q}
        else
            ecode=1
        fi
        rm -rfv hla-sample.$HOSTNAME.$$ &>>{log:q}
        exit $ecode
        """

def collector_dataset(wildcards, prefix, suffix):
    samples = { record['SampleName'] for record in read_sample_sheet(wildcards.dataset) }
    return [ f'{prefix}{wildcards.dataset}_{sample}{suffix}' for sample in samples ]

rule xhla_collector_bam:
    input:
        lambda wildcards : collector_dataset(wildcards, 'typing/xhla/', '.bam')
    output:
        'typing/xhla/{dataset}.mapping'
    shell:
        """
        touch {output:q}
        """

rule xhla_collector_xhla:
    input:
        lambda wildcards : collector_dataset(wildcards, 'typing/xhla/', '.xhla')
    output:
        'typing/xhla/{dataset}.tsv'
    shell:
        """
        touch {output:q}
        """
