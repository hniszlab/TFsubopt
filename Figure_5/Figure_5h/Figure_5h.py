#!/usr/bin/env python

import coolbox
from coolbox.api import *
%config InlineBackend.figure_formats = ['svg']

from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

data_dir = "/Volumes/hnisz_seq_data-1/ChIPseq/ChIP_CEBPA/PAPER/Julian_ChIP_CEBPA/ToMerge/ToUpload/"
IS15_24_I = f"{data_dir}/IS15_24_ChIP_input_sorted.bw"
IS15_24_ChIP = f"{data_dir}/IS15_24_ChIP_merged.bw"

IS15_48_I = f"{data_dir}/IS15_48_ChIP_input_sorted.bw"
IS15_48_ChIP = f"{data_dir}/IS15_48_ChIP_merged.bw"

WT_24_I = f"{data_dir}/WT_24_ChIP_input_sorted.bw"
WT_24_ChIP = f"{data_dir}/WT_24_ChIP_merged.bw"

WT_48_I = f"{data_dir}/WT_48_ChIP_input_sorted.bw"
WT_48_ChIP = f"{data_dir}/WT_48_ChIP_merged.bw"

WT_48_P = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/ToMerge/WT_48h_merged.bed"
IS15_24_P = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/ToMerge/IS15_24h_merged.bed"
IS15_48_P = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/ToMerge/IS15_48h_merged.bed"
WT_24_P = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/ToMerge/WT_24h_merged.bed"

HG38_bed = f"/Users/magalhae/Desktop/HumanDatabase/HG38_proteincoding_UCSC.sorted.bed"
HG38_gtf = f"/Users/magalhae/Desktop/HumanDatabase/gencode.v41.annotation.sorted.gtf"

IS24_fimo = f"/Users/magalhae/Desktop/fimo.gff"


with MinValue(0), MaxValue(20), Color("#ef4448"):
    frame1 = XAxis() +\
             BigWig(WT_24_ChIP) +\
             Title("Wild type 24h") +\
             BigWig(WT_48_ChIP) +\
             Title("Wild type 48h")
             

with MinValue(0), MaxValue(20), Color("#673ab7"):
    frame2 = BigWig(IS15_24_ChIP) +\
             Title("AroPERFECT IS15 24h") +\
             BigWig(IS15_48_ChIP) +\
             Title("AroPERFECT IS15 48h")

             
frame = frame1 + frame2 + GTF(HG38_gtf) + BED(HG38_bed)  + TrackHeight(2)

#frame.plot("chr1:161,460,000-161,560,000")

frame.plot("chr2:33,592,728-33,608,048")