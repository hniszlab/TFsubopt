#!/usr/bin/env python

import coolbox
from coolbox.api import *
%config InlineBackend.figure_formats = ['svg']

from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

data_dir = "Tracks"

ZIP13K2 = f"{data_dir}/ZIP13K2_WT_merged.bw"

AroLITE_24 = f"{data_dir}/ZIP13K2_NGN2_24h_AroLITE_merged.bw"
AroPERFECT_24 = f"{data_dir}/ZIP13K2_NGN2_24h_AroPERFECT_merged.bw"
NGN2_24 = f"{data_dir}/ZIP13K2_NGN2_24h_WT_merged.bw"

AroLITE_48 = f"{data_dir}/ZIP13K2_NGN2_48h_AroLITE_merged.bw"
AroPERFECT_48 = f"{data_dir}/ZIP13K2_NGN2_48h_AroPERFECT_merged.bw"
NGN2_48 = f"{data_dir}/ZIP13K2_NGN2_24h_WT_merged.bw"
NGN2_482 = f"{data_dir}/ZIP13K2_NGN2_48h_WT_merged.bw"

HG38_bed = f"HG38_proteincoding_UCSC.sorted.bed"
HG38_gtf = f"gencode.v41.annotation.sorted.gtf"

with MinValue(0), MaxValue(25), Color("#ef4448"):
    frame1 = XAxis() +\
             BigWig(NGN2_24) +\
             Title("ZIP13K2 Parental")

with MinValue(0), MaxValue(25), Color("#ef4448"):
    frame2 = BigWig(NGN2_24) +\
             Title("NGN2 24h") +\
             BigWig(NGN2_48) +\
             Title("NGN2 48h")


with MinValue(0), MaxValue(25), Color("#ef4448"):
    frame3 = BigWig(AroLITE_24) +\
             Title("AroLITE 24h") +\
             BigWig(AroLITE_48) +\
             Title("AroLITE 48h")


with MinValue(0), MaxValue(25), Color("#673ab7"):
    frame4 = BigWig(AroPERFECT_24) +\
             Title("AroPERFECT 24h") +\
             BigWig(AroPERFECT_48) +\
             Title("AroPERFECT 48h")


frame = frame1 + frame2 + frame3 + frame4 + GTF(HG38_gtf) + BED(HG38_bed)  + TrackHeight(1)

#frame.plot("chr19:41770000-41830761")

frame.plot("chr17:28,276,513-28,371,372")
