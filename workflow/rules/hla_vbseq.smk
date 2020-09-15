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

def check_bin(fname):
    path = os.path.join('./bin/HLA-VBSeq', fname)
    if not os.path.isfile(path):
        raise RuntimeError(f'Missing HLA VBSeq program \'{path}\'. Make sure all required HLA VBSeq files have been installed to bin/HLA-VBSeq')
    return path

def get_reference(wildcards, ftype):
    infix = '_v2' if wildcards.version == '2' else ''

    if ftype=='fasta':
        return f'hla_all{infix}.fasta'
    if ftype=='txt':
        return f'Allelelist{infix}.txt'
    raise RuntimeError(f'get_reference() unknown ftype "{ftype}"')

rule hla_vbseq_typing:
    input:
        bam = 'map/{base}.bam',
        bai = 'map/{base}.bam.bai',
    output:
        result_txt      = 'typing/vbseq/{base}.v{version}.vbseq.txt',
        result_tsv      = 'typing/vbseq/{base}.v{version}.vbseq.tsv',
        idx_dir         = temp(directory('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.bam.idx')),
        locus_ids       = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.loc_ids'),
        partial_sam     = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.partial.sam'),
        partial_fq_1    = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.partial_1.fastq'),
        partial_fq_2    = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.partial_2.fastq'),
        unmapped_fq_1   = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.unmapped_1.fastq'),
        unmapped_fq_2   = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.unmapped_2.fastq'),
        cat_fq_1        = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.cat_1.fastq'),
        cat_fq_2        = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.cat_2.fastq'),
        bam_unmapped    = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.unmapped.bam'),
        sam_final       = temp('typing/vbseq/{base}.v{version}.vbseq.workdir/{base}.v{version}.vbseq.final.sam'),
    log:
        'typing/vbseq/{base}.v{version}.vbseq.log',
    conda:
        '../envs/hla_vbseq.yaml'
    params:
        hla_vbseq_main_jar               = lambda wildcards : check_bin('HLAVBSeq.jar'),
        hla_vbseq_bam_name_index_jar     = lambda wildcards : check_bin('bamNameIndex.jar'),
        hla_vbseq_call_hla_digits        = lambda wildcards : check_bin('call_hla_digits.py'),

        hla_vbseq_ref_fa                 = lambda wildcards : check_bin(get_reference(wildcards, 'fasta')),
        hla_vbseq_ref_txt                = lambda wildcards : check_bin(get_reference(wildcards, 'txt')),
    resources:
        # Some times xHLA just leaks infinite amounts of memory and will be
        # killed by the OOM killer or grid engine. Oh well...
        mem_mb = '25G',
        time   = '6:00:00',
    threads:
        6
    shell:
        """
        # Extract a list of read name that were aligned to HLA loci (HLA-A, B, C, DM, DO, DP, DQ, DR, E, F, G, H, J, K, L, P, V, MIC, and TAP)
        samtools view -@ {threads} {input.bam:q} chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863 \
            chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825 \
            chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664 \
            chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613 \
            chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106 \
            chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009 \
            chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202 \
            chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588 \
            2> {log:q} | awk '{{print $1}}' | LC_ALL=C sort | uniq > {output.locus_ids:q}

        echo "'samtools view' finished extracting chromosome 6 reads." >>{log:q}

        mkdir -p {output.idx_dir:q}

        java -jar -Xmx32g -Xms32g '{params[hla_vbseq_bam_name_index_jar]}' index {input.bam:q} --indexFile {output.idx_dir:q}/index &>>{log:q}
        echo "bamNameIndex.jar [1/2] finished." >>{log:q}

        java -jar '{params[hla_vbseq_bam_name_index_jar]}' search {input.bam:q} --indexFile {output.idx_dir:q}/index --name {output.locus_ids:q} --output {output.partial_sam:q} &>>{log:q}
        echo "bamNameIndex.jar [2/2] finished." >>{log:q}

        picard SamToFastq I={output.partial_sam:q} F={output.partial_fq_1:q} F2={output.partial_fq_2:q} &>>{log:q}
        echo "picard SamToFastq [1/2] finished." >>{log:q}

        samtools view -@ {threads} -bh -f 12 {input.bam:q} >{output.bam_unmapped:q} 2>>{log:q}
        echo "'samtools view' finished extracting unmapped reads." >>{log:q}

        picard SamToFastq I={output.bam_unmapped:q} F={output.unmapped_fq_1:q} F2={output.unmapped_fq_2:q} &>>{log:q}
        echo "picard SamToFastq [2/2] finished." >>{log:q}

        cat {output.partial_fq_1:q} {output.unmapped_fq_1:q} >{output.cat_fq_1} 2>>{log:q}
        cat {output.partial_fq_2:q} {output.unmapped_fq_2:q} >{output.cat_fq_2} 2>>{log:q}
        echo "Generating input FASTQ files finished." >>{log:q}

        bwa index '{params[hla_vbseq_ref_fa]}' &>>{log:q}
        bwa mem -t {threads} -P -L 10000 -a '{params[hla_vbseq_ref_fa]}' {output.cat_fq_1} {output.cat_fq_2} >{output.sam_final:q} 2>>{log:q}
        echo "Mapping against HLA reference finished." >>{log:q}

        # Run VBSeq
        java -jar '{params[hla_vbseq_main_jar]}'  '{params[hla_vbseq_ref_fa]}' {output.sam_final:q} {output.result_txt:q} --alpha_zero 0.01 --is_paired &>>{log:q}
        echo "VBSeq finished." >>{log:q}

        # Call alleles
        python '{params[hla_vbseq_call_hla_digits]}' -v {output.result_txt:q} -a '{params[hla_vbseq_ref_txt]}'  -r 140 -d 4 --ispaired >{output.result_tsv:q} 2>>{log:q}
        echo "Calling HLA alleles finished." >>{log:q}

        echo "" >>{log:q}
        echo "Rule 'hla_vbseq_typing' finished." >>{log:q}
        """
