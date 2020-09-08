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
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule hla_la_graph_compilation:
    conda:
        '../envs/hla_la.yaml'
    input:
        HTTP.remote('http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz', keep_local = True)
    output:
        'typing/hla_la/hla_la.graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH'
    log:
        'typing/hla_la/hla_la.graphs.log',
    params:
        checksum='525a8aa0c7f357bf29fe2c75ef1d477d',
        # Parameters for cluster execution
        cluster_mem = '40G',
        cluster_rt = '5:00:00',
    threads:
        1
    shell:
        """
        logpath="$PWD/{log}"
        origpath="$PWD"
        date >"$logpath"

        # CHECK DOWNLOADED FILE
        md5sum < {input:q} 2>>"$logpath" | fgrep -q {params.checksum} || echo 'Checksum failure for PRG_MHC_GRCh38_withIMGT.tar.gz' | tee -a "$logpath" >&2

        # UNPACK DOWNLOADED FILE
        cd 'typing/hla_la/hla_la.graphs' &>>"$logpath"
        tar xf "$origpath"/{input:q} &>>"$logpath"

        # PREPARE THE INFERENCE GRAPH
        "$CONDA_PREFIX"/opt/hla-la/bin/HLA-LA --action prepareGraph --PRG_graph_dir "$PWD"/PRG_MHC_GRCh38_withIMGT &>>"$logpath"
        """

rule hla_la_typing:
    conda:
        '../envs/hla_la.yaml'
    input:
        bam            = 'map/{base}.bam',
        bai            = 'map/{base}.bam.bai',
        graphs         = 'typing/hla_la/hla_la.graphs/',
        graphs_preproc = 'typing/hla_la/hla_la.graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH',
    output:
        main = 'typing/hla_la/{base}.hla_la.txt',
        outdir = directory('typing/hla_la/{base}.hla_la.out'),
    params:
        tempdir = lambda wildcards : f'typing/hla_la/{wildcards.base}.hla_la.tempdir'
    log:
       'typing/hla_la/{base}.hla_la.log',
    shell:
        """
        LOGPATH="$PWD/{log}"
        date > "$LOGPATH"
        env >> "$LOGPATH"
        # Check if a proper conda environment is available
        if [ -z "${{CONDA_PREFIX+_}}" ] || [ ! -d "$CONDA_PREFIX"/opt/hla-la ]; then
            echo "--use-conda is required to run the HLA-LA rule." | tee -a "$LOGPATH" >&2
            exit 1
        fi

        # Create and prepare the temp directory
        rm -rf {params.tempdir:q}
        mkdir -p {params.tempdir:q}/src {params.tempdir:q}/bin >>"$LOGPATH"

        # PREPARE DIRECTORY SRC
        cd {params.tempdir:q}/src
        ln -s "$CONDA_PREFIX"/opt/hla-la/src/* . &>>"$LOGPATH"
        mv HLA-LA.pl HLA-LA.pl~
        cp "$CONDA_PREFIX"/opt/hla-la/src/HLA-LA.pl . &>>"$LOGPATH"
        cd ..

        # PREPARE DIRECTORY BIN
        ln -s "$CONDA_PREFIX"/opt/hla-la/bin/HLA-LA bin/ &>>"$LOGPATH"

        # PREPARE DIRECTORY GRAPHS
        ln -s ../../../{input.graphs:q} graphs &>>"$LOGPATH"

        # RUN MAIN PROGRAM
        ./src/HLA-LA.pl --sampleID sample --BAM ../../../{input.bam:q} --workingDir "$PWD" --graph PRG_MHC_GRCh38_withIMGT &>>"$LOGPATH"

        # MOVE RESULTS INTO PLACE
        cd ../../../
        rm -f {output.outdir:q} &>>"$LOGPATH"
        mv {params.tempdir:q}/sample {output.outdir:q} &>>"$LOGPATH"
        ln -s ../../{output.outdir:q}/hla/R1_bestguess_G.txt {output.main:q} &>>"$LOGPATH"

        # REMOVE TEMPORARY DIRECTORY
        rm -rf {params.tempdir:q} &>>"$LOGPATH"
        """
