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

def star_input(wildcards):
    reads = input_reads(wildcards)
    return {
            'fq1'     : reads[0],
            'fq2'     : reads[1],
            'ref'     : os.path.join('ref', wildcards.ref + '.fa'),
            'ref_idx' : os.path.join('ref', wildcards.ref + '.star'),
        }

rule star_index:
    input:
        fasta = "ref/{base}.fa"
    output:
        directory("ref/{base}.star"),
    resources:
        mem_mb = '40G',
        time   = '2:00:00',
    threads:
        8
    params:
        extra = ""
    log:
        "ref/{base}.star.log"
    wrapper:
        "0.67.0/bio/star/index"

rule star_map_pe:
    input:
        unpack(star_input)
    output:
        bam = 'map/star/{dataset}_{sample}_{trim}_{ref}/Aligned.sortedByCoord.out.bam',
    log:
        'map/star/{dataset}_{sample}_{trim}_{ref}.log'
    resources:
        mem_mb = '40G',
        time   = '2:00:00',
    params:
        # path to STAR reference genome index
        index=lambda wildcards, input : input.ref_idx,
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate"
    threads: 8
    wrapper:
        "0.67.0/bio/star/align"

