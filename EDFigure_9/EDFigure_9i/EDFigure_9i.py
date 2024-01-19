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



CEBP_promoters = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/SE/REs_BtoiMac/promoters_sorted_cut.bed"
CEBP_MacroSE = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/SE/REs_BtoiMac/iMac_SE_cut.bed"
CEBP_Enh = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/SE/REs_BtoiMac/enhancers_sorted_cut.bed"
CEBP_BcellSE = f"/Users/magalhae/Desktop/Julian_ChIP_CEBPA/SE/REs_BtoiMac/BLaER_SE_cut.bed"

With MinValue(0), MaxValue(35), Color("#ef4448"):
    frame1 = XAxis() +\
             BigWig(WT_24_I) +\
             Title("Wild type Input") +\
             BigWig(WT_24_ChIP) +\
             Title("Wild type 24h") +\
             BigWig(WT_48_ChIP) +\
             Title("Wild type 48h")


with MinValue(0), MaxValue(35), Color("#673ab7"):
    frame2 = BigWig(IS15_24_I) +\
             Title("AroPERFECT IS15 Input") +\
             BigWig(IS15_24_ChIP) +\
             Title("AroPERFECT IS15 24h") +\
             BigWig(IS15_48_ChIP) +\
             Title("AroPERFECT IS15 48h")

with Color("#000000"):
    frame3 = BED(CEBP_promoters) +\
             Title("CEBP_promoters") +\
             BED(CEBP_Enh) +\
             Title("CEBP_Enh") +\
             BED(CEBP_MacroSE) +\
             Title("CEBP_MacroSE") +\
             BED(CEBP_BcellSE) +\
             Title("CEBP_BcellSE")

frame = frame1 + frame2 + GTF(HG38_gtf) + BED(HG38_bed) + frame3  + TrackHeight(2)

frame.plot("chr19:42,359,074-42,676,713")

frame.plot("chr1:161,460,000-161,660,000")
