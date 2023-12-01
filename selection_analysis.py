#!/usr/bin/python

import dendropy
from dendropy.calculate import popgenstat
import pandas as pd
import glob
import os
import dnds
import sys 
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

########################################################################

## Restructure the fasta files so they are one file for each locus with
## all samples included.

########################################################################

path='path'
r_files= glob.glob(os.path.join(path, "*_SETA_stub.fasta"))
effector_files = glob.glob(os.path.join(path, "*_SETA_pinf.fasta"))


rli = []
eli = []


for filename in r_files:
    seqs = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        seqs.append((str(seq_record.id)+'\t'+str(seq_record.seq)))
    name = os.path.basename(filename)
    sample = name.split("_",1)[0]
    df = pd.DataFrame([sub.split("\t") for sub in seqs])
    df["Sample"] = sample
    rli.append(df)



for filename in effector_files:
    seqs=[]
    for seq_record in SeqIO.parse(filename, "fasta"):
        seqs.append((str(seq_record.id)+'\t'+str(seq_record.seq)))
    name = os.path.basename(filename)
    sample = name.split("_",1)[0]
    df = pd.DataFrame([sub.split("\t") for sub in seqs])
    df["Sample"] = sample
    eli.append(df)

rframe = pd.concat(rli, axis=0, ignore_index=True)
eframe = pd.concat(eli, axis=0, ignore_index=True)

headers=["Location","Sequence","Sample"]

rframe.columns=headers
eframe.columns=headers

rframe = rframe.drop(rframe[rframe.Sample.isin(["K10", "K48", "K52", "US0186686", "US0186932", "US0186856","US0186860", "US0186841","US0186843"])].index)
eframe = eframe.drop(eframe[eframe.Sample.isin(["K10", "K48", "K52", "US0186856","US0186860"])].index)

print(rframe)
print(eframe)

rgroups=rframe.groupby("Location")
egroups=eframe.groupby("Location")



for name, group in rgroups:
    group = group.reset_index()
    with open("{}_rgenesbysamples.fasta".format(name), 'w') as f:
        for i in range(len(group)):
            f.write(">" + str(group["Sample"][i]) + "\n")
            f.write(str(group["Sequence"][i]) + "\n")

for name, group in egroups:
    group = group.reset_index()
    with open("{}_egenesbysamples.fasta".format(name), 'w') as f:
        for i in range(len(group)):
            f.write(">" + str(group["Sample"][i]) + "\n")
            f.write(str(group["Sequence"][i]) + "\n")

#######################################################################

## Load in the resturctured files for use with DendroPy.

#######################################################################



r_bytaxa= glob.glob(os.path.join(path, "*_rgenesbysamples.fasta"))
e_bytaxa = glob.glob(os.path.join(path, "*_egenesbysamples.fasta"))


locus_list=[]
avg_num_pairwise_diff_list=[]
num_seg_sites_list=[]
nucleotide_diversity_list=[]
wattersons_theta_list=[]
tajimas_d_list=[]


for filename in r_bytaxa:
    locus_fasta = dendropy.DnaCharacterMatrix.get(file=open(filename), schema="fasta")
    name = os.path.basename(filename)
    locus = name.split("_",2)[1]
    avg_num_pairwise_diff = dendropy.calculate.popgenstat.average_number_of_pairwise_differences(locus_fasta)
    num_seg_sites = dendropy.calculate.popgenstat.num_segregating_sites(locus_fasta)
    nucleotide_diversity = dendropy.calculate.popgenstat.nucleotide_diversity(locus_fasta)
    wattersons_theta = dendropy.calculate.popgenstat.wattersons_theta(locus_fasta)
    tajimas_d = dendropy.calculate.popgenstat.tajimas_d(locus_fasta)

    locus_list.append(locus)
    avg_num_pairwise_diff_list.append(avg_num_pairwise_diff)
    num_seg_sites_list.append(num_seg_sites)
    nucleotide_diversity_list.append(nucleotide_diversity)
    wattersons_theta_list.append(wattersons_theta)
    tajimas_d_list.append(tajimas_d)

rselection_analysis = pd.DataFrame({"Locus": locus_list, "Average Number of Pairwise Differences": avg_num_pairwise_diff_list, "Nucleotide Diversity": nucleotide_diversity_list, "Number of Segregating Sites": num_seg_sites_list, "Tajima's D": tajimas_d_list, "Watterson's Theta": wattersons_theta_list})

rselection_analysis.to_csv("04212023_SETA_rgene_selection_analysis.csv")


locus_list=[]
avg_num_pairwise_diff_list=[]
num_seg_sites_list=[]
nucleotide_diversity_list=[]
wattersons_theta_list=[]
tajimas_d_list=[]


for filename in e_bytaxa:
    locus_fasta = dendropy.DnaCharacterMatrix.get(file=open(filename), schema="fasta")
    name = os.path.basename(filename)
    locus = name.split("_",1)[0]
    try:
        avg_num_pairwise_diff = dendropy.calculate.popgenstat.average_number_of_pairwise_differences(locus_fasta)
    except:
        avg_num_pairwise_diff = "ERROR"
    try:
        num_seg_sites = dendropy.calculate.popgenstat.num_segregating_sites(locus_fasta)
    except:
        num_seg_sites = "ERROR"
    try:
        nucleotide_diversity = dendropy.calculate.popgenstat.nucleotide_diversity(locus_fasta)
    except:
        nucleotide_diversity = "ERROR"
    try:
        wattersons_theta = dendropy.calculate.popgenstat.wattersons_theta(locus_fasta)
    except:
        wattersons_theta = "ERROR"
    try:
        tajimas_d = dendropy.calculate.popgenstat.tajimas_d(locus_fasta)
    except:
        wattersons_theta = "ERROR"

    locus_list.append(locus)
    avg_num_pairwise_diff_list.append(avg_num_pairwise_diff)
    num_seg_sites_list.append(num_seg_sites)
    nucleotide_diversity_list.append(nucleotide_diversity)
    wattersons_theta_list.append(wattersons_theta)
    tajimas_d_list.append(tajimas_d)

eselection_analysis = pd.DataFrame({"Locus": locus_list, "Average Number of Pairwise Differences": avg_num_pairwise_diff_list, "Nucleotide Diversity": nucleotide_diversity_list, "Number of Segregating Sites": num_seg_sites_list, "Tajima's D": tajimas_d_list, "Watterson's Theta": wattersons_theta_list})

eselection_analysis.to_csv("04212023_SETA_effector_selection_analysis.csv")
    














