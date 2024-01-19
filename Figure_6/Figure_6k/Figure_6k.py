#!/usr/bin/env python

import coolbox
from coolbox.api import *
%config InlineBackend.figure_formats = ['svg']

from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

data_dir = "Tracks"

AroLITE_12h_f = f"{data_dir}/ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_forward.bw"
AroLITE_12h_r = f"{data_dir}/ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_reverse.bw"
AroLITE_24h_f = f"{data_dir}/ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_forward.bw"
AroLITE_24h_r = f"{data_dir}/ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_reverse.bw"
AroPERFECT_12_r = f"{data_dir}/ZIP13K2_NGN2_AroPERFECT_12_TCReads_cov_CPM_reverse.bw"
AroPERFECT_12_f = f"{data_dir}/ZIP13K2_NGN2_AroPERFECT_12h_TCReads_cov_CPM_forward.bw"
AroPERFECT_24_f = f"{data_dir}/ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_forward.bw"
AroPERFECT_24_r = f"{data_dir}/ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_reverse.bw"
WT_12h_r = f"{data_dir}/ZIP13K2_NGN2_WT_12h_r1_TCReads_cov_CPM_reverse.bw"
WT_12h_f = f"{data_dir}/ZIP13K2_NGN2_WT_12h_TCReads_cov_CPM_forward.bw"
WT_24h_f = f"{data_dir}/ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_forward.bw"
WT_24h_r = f"{data_dir}/ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_reverse.bw"
Parental_f = f"{data_dir}/ZIP13K2_TCReads_cov_CPM_forward.bw"
Parental_r = f"{data_dir}/ZIP13K2_TCReads_cov_CPM_reverse.bw"

HG38_bed = f"HG38_proteincoding_UCSC.sorted.bed"
HG38_gtf = f"gencode.v41.annotation.sorted.gtf"
```


AroLITE_12h = f"{data_dir}/bw/ZIP13K2_NGN2_24h_AroLITE_merged.bw"
AroPERFECT_24h = f"{data_dir}/bw/ZIP13K2_NGN2_24h_AroPERFECT_merged.bw"
WT_24h = f"{data_dir}/bw/ZIP13K2_NGN2_24h_WT_merged.bw"
AroLITE_48h = f"{data_dir}/bw/ZIP13K2_NGN2_48h_AroLITE_merged.bw"
AroPERFECT_48h = f"{data_dir}/bw/ZIP13K2_NGN2_48h_AroPERFECT_merged.bw"
WT_48h = f"{data_dir}/bw/ZIP13K2_NGN2_48h_WT_merged.bw"

with MinValue(0), MaxValue(25), Color("#ef4548"):
    frame1 = XAxis() +\
             Title("Wild type 24h") +\
             BigWig(WT_24h) +\
             Title("Wild type 24h") +\
             BigWig(WT_48h) +\
             Title("Wild type 48h")


with MinValue(0), MaxValue(25), Color("#664a9e"):
    frame2 = BigWig(AroPERFECT_24h) +\
             Title("AroPERFECT 24h") +\
             BigWig(AroPERFECT_48h) +\
             Title("AroPERFECT 48h")

with MinValue(0), MaxValue(25), Color("#1d59a8"):
    frame3 = BigWig(AroLITE_12h) +\
             Title("AroLITE_12h") +\
             BigWig(AroLITE_48h) +\
             Title("AroLITE_48h")



frame = frame1 + frame2 + frame3 + GTF(HG38_gtf) + BED(HG38_bed)  + TrackHeight(2)

frame.plot("chr2:30,100,000-30,390,000")
