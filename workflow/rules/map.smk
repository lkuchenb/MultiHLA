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

import shlex
import helpers
import os

def bwa_mem_input(wildcards):
    dataset = wildcards.dataset
    sample  = wildcards.sample
    trim    = wildcards.trim == 'trim'
    ref     = wildcards.ref

    records = [ record for record in helpers.read_sample_sheet(dataset) if record['SampleName'] == sample ]

    fastqs = [ record['FileNameR1'] for record in records ] + [ record['FileNameR2'] for record in records ]

    if trim:
        inputs = { 'fastq_' + str(idx) : os.path.join('trim', re.sub('fastq.gz$', 'trimmed.fastq.gz', fastq)) for idx, fastq in zip(range(len(fastqs)), fastqs) }
    else:
        inputs = { 'fastq_' + str(idx) : os.path.join('fastq', fastq) for idx, fastq in zip(range(len(fastqs)), fastqs) }

    # Add the dependecy to the specified reference
    inputs['ref'] = os.path.join('ref', ref + '.fa')

    return inputs

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

rule bwa_mem:
    # Read mapping with BWQ mem default params as instructed by xhla authors
    conda:
        "../envs/bwa.yaml"
    input:
        unpack(bwa_mem_input)
    output:
        bam = 'map/{dataset}_{sample}_{trim}_{ref}.bam',
        bai = 'map/{dataset}_{sample}_{trim}_{ref}.bam.bai',
    log:
        bwa = 'map/{dataset}_{sample}_{trim}_{ref}.bwa.log',
        sort = 'map/{dataset}_{sample}_{trim}_{ref}.sort.log',
        index = 'map/{dataset}_{sample}_{trim}_{ref}.index.log',
    wildcard_constraints:
        trim = 'trim|orig'
    params:
        r1_fastqs = lambda wildcards, input : bwamem_params(input)[0],
        r2_fastqs = lambda wildcards, input : bwamem_params(input)[1],

        # Parameters for cluster execution
        cluster_mem = '16G',
        cluster_rt = '8:00:00',
    threads:
        6
    shell:
        """
        # Read mapping using BWA and BAM sorting using samtools
        bwa mem -t {threads} \
                {input.ref:q} \
                <(unpigz -c {params[r1_fastqs]}) \
                <(unpigz -c {params[r2_fastqs]}) \
                2>{log.bwa:q} | samtools sort -T /share/scratch/kuchenb.tmp/samtools.$HOSTNAME.$$ -o {output.bam} &>{log.sort:q}

        # BAM indexing using samtools
        samtools index -@ {threads} {output.bam:q} &>{log.index:q}
        """

def bwa_mem_dataset_collector_input(wildcards):
    return [
        os.path.join('map', f'{wildcards.dataset}_{sample}_{wildcards.trim}_{wildcards.ref}.bam')
        for sample in { record['SampleName'] for record in helpers.read_sample_sheet(wildcards.dataset) }
        ]

rule bwa_mem_dataset_collector:
    input:
        bwa_mem_dataset_collector_input
    output:
        'map/{dataset}_{trim}_{ref}.map.all'
    shell:
        """
        touch {output:q}
        """
