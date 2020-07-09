#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cBioPortal.py

Usage: ./cBioPortal.py --gene 'HUGO_gene_name'
"""

# Michael Kesling, July 1, 2020

##########
# Import packages
##########


import pandas as pd
import argparse
from bravado.client import SwaggerClient


# very useful for trouble-shooting:
# https://www.cbioportal.org/api/swagger-ui.html#/Mutations/getMutationsInMolecularProfileBySampleListIdUsingGET


##########
# Set up argparse and activate the swagger client
##########
parser = argparse.ArgumentParser(description="Query all somatic variants for an inputted \
                        HUGO gene name in the cBioPortal website. \
                            Returns each variant allele and the number of \
                                distinct samples it's been found in.")
parser.add_argument('--gene', type=str)
args = parser.parse_args()
gene = args.gene.upper()


cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False, \
                                        "validate_swagger_spec":False})

    
#########
# for QA:
# gene = 'KRAS'
# for a in dir(cbioportal):
#     cbioportal.__setattr__(a.replace(' ', '_').lower(), cbioportal.__getattr__(a))
#########


##########
# convert HUGO gene name to entrezID:
##########
try:
    gene_tmp = cbioportal.Genes.getGeneUsingGET(geneId=gene).response().result
except:
    print("\nGene Name Entered is not a valid HUGO name\n")
    
entrezID = gene_tmp['entrezGeneId']



###########3
# Gather 'MUTATION_EXTEDED Molecular Profiles:
###########
MPs = cbioportal.Molecular_Profiles.getAllMolecularProfilesUsingGET().result()

MP_mutations = {}
for i in range(len(MPs)):
    if MPs[i]['molecularAlterationType'] == 'MUTATION_EXTENDED':
        MP_mutations[i] = MPs[i]
#print(len(MP_mutations))  
# down to 286 molecular profiles from 1186


##########
# for QA:
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
##########


##########
# Query gene-specific mutations one Molecular Profile at a time.
# Store individual mutations in mutDF dataframe.
##########
mutDF = pd.DataFrame(columns=['chr','start','end','gene', \
                              'refAllele','varAllele', \
                              'sampleID', 'molecAltType'])
    

printerCounter = 0
for idx in MP_mutations.keys():
    mutProfileID = MPs[idx]['studyId'] + "_mutations"
    smplListID = MPs[idx]['studyId'] + '_all'
    muts = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        molecularProfileId = mutProfileID,
        sampleListId = smplListID,
        entrezGeneId=entrezID,
        projection='DETAILED'          # needed for gene-level info
        ).result()


    # grab desired fields and place in dataframe
    for jdx in range(len(muts)):
        gene = muts[jdx]['gene']['hugoGeneSymbol']                             # UNCOMMENT
        chrm = muts[jdx]['chr']
        start = muts[jdx]['startPosition']
        end = muts[jdx]['endPosition']
        refAllele = muts[jdx]['referenceAllele']
        varAllele = muts[jdx]['variantAllele']
        sampleID = muts[jdx]['sampleId']
        lst =[[chrm, start, end, gene, refAllele,varAllele, sampleID,\
                MP_mutations[idx]['molecularAlterationType']]]
        tmpDF = pd.DataFrame(lst, columns=['chr', 'start', 'end','gene', 'refAllele', \
                              'varAllele', 'sampleID','molecAltType'])
        mutDF = pd.concat([mutDF, tmpDF], axis=0)
    printerCounter += 1
    print("Molecular Profile #:", printerCounter, "of", len(MP_mutations), \
          ". Number of possibly redundant", gene, "mutations in this Mol. Profile= ",len(muts), \
          "     MP name =", mutProfileID)
    



##########
# Group the mutation data frame by exact coordinates and 
# refAllele and varAllele values, then gather unique samples,
# and convert to output dataframe:
##########
grouped = mutDF.groupby(['chr', 'start', 'end', \
                          'refAllele', 'varAllele', 'gene'])['sampleID']

outputDF = pd.DataFrame(columns=['chr','start','end','gene', \
                              'refAllele','varAllele','numUniqSamples'])
for name, group in grouped:
    if name[1] < 1:           # remove if not a coordinate
        continue
    lst = [list(name) + list([len(group.unique())])]
    tmpDF = pd.DataFrame(lst, columns=['chr', 'start', 'end','gene', \
                                        'refAllele', 'varAllele', 'numUniqSamples'])

    outputDF = pd.concat([outputDF, tmpDF], axis=0)
    

########## 
# output file:
##########
outputFile = gene + "-mutationTable.tsv"
outputDF.to_csv(outputFile, sep='\t', mode='w', \
                index=False)


    

    
################################################
#
#
# QA -- not used in the running of the completed script
#
#
###############################################

# check to see how many molecularAlterationTypes
# there are in cBioPortal:
# mp = set()
# for i in range(len(MPs)):
#     mp.add(MPs[i]['molecularAlterationType'])
# print(mp)


# ########## 
# Check how many Molecular Profiles the various  non-MUTATION_EXTENDED classes encompass:
# MP_mutations = {}
# for i in range(len(MPs)):
#     if MPs[i]['molecularAlterationType'] == 'STRUCTURAL_VARIANT': 
#             MP_mutations[i] = MPs[i]
# print(len(MP_mutations))
# 286 for MUTATION_EXTENDED, and it has SNVs.
#   1 for MUTATION_UNCALLED, and it has SNVs.
# 103 for STRUCTURAL_VARIANT, and it has SNVs.
# 398 for MRNA_EXPRESSION, and it has SNVs.
#   3 FOR GENERIC_ASSAY, and it has SNVs.
#  65 for METHYLATION, and it has SNVs.
# 260 for COPY_NUMBER_ALTERATION, and it has SNVs.
#  80 for PROTEIN_LEVEL, and it has SNVs.
# (showing that these have SNVs is not shown here.)



##########
# Confirm number of samples in Molecular Profile = MUTATION_EXTENDED
# are equal to the number can be compared with one on cBioPortal website
#########

# keep track of all samples -- not unique across MPs.
# totalSamples = []

# printerCounter = 0
# for idx in MP_mutations.keys():
#     samples = set()
#     mutProfileID = MPs[idx]['studyId'] + "_mutations"
#     smplListID = MPs[idx]['studyId'] + '_all'
#     muts = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
#         molecularProfileId = mutProfileID,
#         sampleListId = smplListID,
#         ).result()

#     # grab desired fields and add to set:
#     for jdx in range(len(muts)):
#         sampleID = muts[jdx]['sampleId']
#         # refAllele = muts[jdx]['referenceAllele']
#         # print(sampleID, refAllele)
#         samples.add(sampleID)    # want non-redundant samples within each MP
#         # I should count if there are any mutations
        

#     # now convert samples in samples set into totalSamples list
#     totalSamples = totalSamples + list(samples)
    
#     # print out status to stderr
#     printerCounter += 1
#     print("Molecular Profile #:", printerCounter, "of", len(MP_mutations), \
#           "Number of unique samples in this MP:", mutProfileID, "is", len(samples), \
#           ". QA index:", idx)

# len(totalSamples)       # total samples
# # 78,812
# len(set(totalSamples))  # unique samples
# # 52357 samples




########################
# checking more ideas
########################

# studies = cbioportal.Studies.getAllStudiesUsingGET().result()
# print("The total number of samples in all studies is: {}".format(sum([x.allSampleCount for x in studies])))
# # 91,340 total samples

# # total genes:
# allGenes = cbioportal.Genes.getAllGenesUsingGET().result()
# len(allGenes)
######################################


