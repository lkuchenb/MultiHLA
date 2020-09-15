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

rule xhla_typing:
    container:
        'docker://humanlongevity/hla'
    input:
        bam = 'map/{base}.bam',
        bai = 'map/{base}.bam.bai',
    output:
        json = 'typing/xhla/{base}.xhla.json',
        tdir = temp(directory('typing/xhla/{base}.xhla.workdir')),
        bindir = temp(directory('typing/xhla/{base}.xhla.bindir')),
    log:
        'typing/xhla/{base}.xhla.log',
    resources:
        # Some times xHLA just leaks infinite amounts of memory and will be
        # killed by the OOM killer or grid engine. Oh well...
        mem_mb = '200G',
        time   = '2:00:00',
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
