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

rule bwa_index:
    input:
        'ref/{genome}'
    output:
        'ref/{genome}.0123',
        'ref/{genome}.amb',
        'ref/{genome}.ann',
        'ref/{genome}.bwt.2bit.64',
        'ref/{genome}.bwt.8bit.32',
        'ref/{genome}.pac',
    log:
        'ref/{genome}.bwaindex.log'
    params:
        prefix='ref/{genome}',
        algorithm='bwtsw'
    resources:
        mem_mb = '100G',
        time   = '1:00:00',
    wrapper:
        '0.65.0/bio/bwa-mem2/index'

def bwa_mem_input(wildcards):
    return {
            'reads' : input_reads(wildcards),
            'ref'   : os.path.join('ref', wildcards.ref + '.fa'),
            'pac'   : os.path.join('ref', wildcards.ref + '.fa.pac'),
        }

rule samtools_index:
    conda:
        '../envs/samtools.yaml'
    input:
        'map/{sample}.bam'
    output:
        'map/{sample}.bam.bai'
    resources:
        mem_mb = '16G',
        time   = '1:00:00',
    threads:
        1
    shell:
        """
        samtools index {input:q} {output:q}
        """

rule bwa_mem_2:
    input:
        unpack(bwa_mem_input)
    output:
        bam = 'map/{dataset}_{sample}_{trim}_{ref}.bam',
    log:
        bwa = 'map/{dataset}_{sample}_{trim}_{ref}.log',
    wildcard_constraints:
        trim = 'trim|orig'
    params:
        # Wrapper parameters
        index=lambda wildcards, input : input.ref,
        extra="",
        sort="samtools",         # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",           # Extra args for samtools/picard.
    resources:
        mem_mb = '70G',
        time   = '8:00:00',
    threads:
        6
    wrapper:
        "0.65.0/bio/bwa-mem2/mem"

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
    resources:
        mem_mb = '1GB',
        time   = '0:15:00',
    shell:
        """
        touch {output:q}
        """
