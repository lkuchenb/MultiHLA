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

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule arcashla_software_installation:
    conda:
        '../envs/arcasHLA.yaml'
    input:
        HTTP.remote('https://github.com/RabadanLab/arcasHLA/archive/v{version}.tar.gz', keep_local = False)
    output:
        outdir = directory("typing/arcasHLA/arcasHLA.v{version}.ref{refversion}"),
        outbin = "typing/arcasHLA/arcasHLA.v{version}.ref{refversion}/arcasHLA"
    log:
        "typing/arcasHLA/arcasHLA.v{version}.ref{refversion}.log"
    resources:
        mem_mb = '5G',
        time   = '1:00:00',
    shell:
        """
        base="$PWD"
        log="$base"/{log:q}

        # DOWNLOAD AN UNPACK ARCASHLA
        mkdir -p {output.outdir:q} &>"$log"
        cd {output.outdir:q} &>>"$log"
        tar --strip-components=1 -xf "$base"/{input:q} &>>"$log"

        # OBTAIN REFERENCE
        ./arcasHLA reference --version {wildcards.refversion:q} &>>"$log"
        """

rule arcashla_extract:
    conda:
        '../envs/arcasHLA.yaml'
    input:
        arcas_bin = "typing/arcasHLA/arcasHLA.v{version}.ref{refversion}/arcasHLA",
        bam = "map/star/{base}/Aligned.sortedByCoord.out.bam",
        bai = "map/star/{base}/Aligned.sortedByCoord.out.bam.bai",
    output:
        r1 = "typing/arcasHLA/{base}_R1.arcasHLA.v{version}.ref{refversion}.extract.fastq.gz",
        r2 = "typing/arcasHLA/{base}_R2.arcasHLA.v{version}.ref{refversion}.extract.fastq.gz",
        temp = temp(directory("typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.extract"))
    log:
        main = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.extract.main.log",
        cli  = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.extract.cli.log",
    resources:
        mem_mb = '5G',
        time   = '1:00:00',
    threads:
        8
    shell:
        """
        rm -rf {output.temp} &>>{log.cli:q}
        mkdir -p {output.temp:q} &>{log.cli:q}
        {input.arcas_bin:q} extract -t {threads} --paired -o {output.temp:q} --log {log.main:q} {input.bam:q} &>>{log.cli:q}
        mv {output.temp:q}/*1.fq.gz {output.r1:q} &>>{log.cli:q}
        mv {output.temp:q}/*2.fq.gz {output.r2:q} &>>{log.cli:q}
        """

rule arcashla_genotype:
    conda:
        '../envs/arcasHLA.yaml'
    input:
        arcas_bin = "typing/arcasHLA/arcasHLA.v{version}.ref{refversion}/arcasHLA",
        r1 = "typing/arcasHLA/{base}_R1.arcasHLA.v{version}.ref{refversion}.extract.fastq.gz",
        r2 = "typing/arcasHLA/{base}_R2.arcasHLA.v{version}.ref{refversion}.extract.fastq.gz",
    output:
        genes    = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.genes.json",
        genotype = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.genotype.json",
        align    = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.alignment.p",
        temp     = temp(directory("typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.genotype")),
    log:
        main = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.extract.main.log",
        cli  = "typing/arcasHLA/{base}.arcasHLA.v{version}.ref{refversion}.extract.cli.log",
    resources:
        mem_mb = '5G',
        time   = '1:00:00',
    threads:
        8
    shell:
        """
        rm -vrf {output.temp} &>>{log.cli:q}
        mkdir -vp {output.temp:q} &>{log.cli:q}
        {input.arcas_bin:q} genotype -t {threads} -o {output.temp:q} --log {log.main:q} {input.r1:q} {input.r2:q} &>>{log.cli:q}
        mv -v {output.temp:q}/*genes.json {output.genes:q} &>>{log.cli}
        mv -v {output.temp:q}/*genotype.json {output.genotype:q} &>>{log.cli}
        mv -v {output.temp:q}/*alignment.p {output.align:q} &>>{log.cli}
        """
