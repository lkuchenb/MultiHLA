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

rule optitype:
    input:
        unpack(lambda wildcards : { 'reads' : input_reads(wildcards) })
    output:
        multiext("typing/optitype/{dataset}_{sample}_{trim,trim|orig}_nofilt", "_coverage_plot.pdf", "_result.tsv")
    log:
        "typing/optitype/{dataset}_{sample}_{trim,trim|orig}_nofilt.log"
    params:
        # Type of sequencing data. Can be 'dna' or 'rna'. Default is 'dna'.
        sequencing_type="dna",
        # optiype config file, optional
        config="",
        # additional parameters
        extra=lambda wildcards : "--prefix " + '_'.join([wildcards.dataset, wildcards.sample, wildcards.trim, 'nofilt'])
    resources:
        mem_mb = '30G',
        time   = '1:00:00',
    wrapper:
        "0.66.0/bio/optitype"
