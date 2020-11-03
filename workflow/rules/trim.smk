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

rule trim:
    # Executes Trim Galore! for adapter and low-quality base trimming
    conda:
        "../envs/trim.yaml"
    input:
        r1 = 'fastq/{prefix}_R1{suffix}.fastq.gz',
        r2 = 'fastq/{prefix}_R2{suffix}.fastq.gz',
    params:
        r1_tmp = 'trim/{prefix}_R1{suffix}_val_1.fq.gz',
        r2_tmp = 'trim/{prefix}_R2{suffix}_val_2.fq.gz',
        r1fastqc_tmp = 'trim/{prefix}_R1{suffix}_val_1_fastqc.html',
        r2fastqc_tmp = 'trim/{prefix}_R2{suffix}_val_2_fastqc.html',
    resources:
        mem_mb = '8G',
        time   = '1:00:00',
    output:
        r1 = 'trim/{prefix}_R1{suffix}.trimmed.fastq.gz',
        r2 = 'trim/{prefix}_R2{suffix}.trimmed.fastq.gz',
        r1fastqc = 'trim/{prefix}_R1{suffix}.trimmed.fastqc.html',
        r2fastqc = 'trim/{prefix}_R2{suffix}.trimmed.fastqc.html',
    wildcard_constraints:
        suffix=".*"
    log:
        'trim/{prefix}{suffix}.trim.log',
    threads:
        6
    shell:
        """
        trim_galore --output_dir trim --fastqc -j {threads} --paired {input.r1:q} {input.r2:q} &>{log:q}
        mv -f {params.r1_tmp:q} {output.r1:q}
        mv -f {params.r2_tmp:q} {output.r2:q}
        mv -f {params.r1fastqc_tmp:q} {output.r1fastqc:q}
        mv -f {params.r2fastqc_tmp:q} {output.r2fastqc:q}
        """
