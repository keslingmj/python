# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import requests
import json
# import numpy as np
from bravado.client import SwaggerClient
# from pprint import pprint

# very useful for trouble-shooting:
# https://www.cbioportal.org/api/swagger-ui.html#/Mutations/getMutationsInMolecularProfileBySampleListIdUsingGET



cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False, \
                                        "validate_swagger_spec":False})
            # rm validation requests is essential


for a in dir(cbioportal):
    cbioportal.__setattr__(a.replace(' ', '_').lower(), cbioportal.__getattr__(a))



#gene = 'BRCA1'
gene = 'KRAS'

gene_tmp = cbioportal.Genes.getGeneUsingGET(geneId=gene).response().result
# might not need this.  Can use gene name directly in Molecular Profiles
entrezID = gene_tmp['entrezGeneId']

# we next want to connect the entrezID to the mutation_event table
#cbioportal.Mutations.fetchMutationsInMultipleMolecularProfilesUsingPOST()
# dir(cbioportal)
# dir(cbioportal.Mutations)




###########3
# Molecular Profiles:
###########
MPs = cbioportal.Molecular_Profiles.getAllMolecularProfilesUsingGET().result()
# mp = set()
# for i in range(len(MPs)):
#     mp.add(MPs[i]['molecularAlterationType'])
# print(mp)

################
# Need to verify for other classes!!!!!!!   DONE, they all work.
# want to also verify the total number of studies I'm querying and 
# compare that to the website numbers
################
# for right now, I'm keeping Molecular Profiles whose 'molecularAlterationType'
# is either 'MUTATION_UNCALLED', 'MUTATION_EXTENDED', 'STRUCTURAL_VARIANT'
# MP_mutations = {}
# for i in range(len(MPs)):
#     if MPs[i]['molecularAlterationType'] == 'MUTATION_EXTENDED' or \
#        MPs[i]['molecularAlterationType'] == 'MUTATION_UNCALLED' or \
#        MPs[i]['molecularAlterationType'] == 'STRUCTURAL_VARIANT':
#            MP_mutations[i] = MPs[i]
# print(len(MP_mutations))  # down to 390 molecular profiles

# ############################### Try other classes
# MP_mutations = {}
# for i in range(len(MPs)):
#     if MPs[i]['molecularAlterationType'] == 'PROTEIN_LEVEL': #or \
#        # MPs[i]['molecularAlterationType'] == 'MUTATION_UNCALLED' or \
#        # MPs[i]['molecularAlterationType'] == 'STRUCTURAL_VARIANT':
#            MP_mutations[i] = MPs[i]
# print(len(MP_mutations))
# 398 for MRNA_EXPRESSION, and works for it!
#   3 FOR GENERIC_ASSAY, and it works!
#  65 for METHYLATION, and it works!
# 260 for COPY_NUMBER_ALTERATION, and it works!
#  80 for PROTEIN_LEVEL, and it works!


MP_mutations = MPs





# now I'm going to look into these and see if they're the kind I
# want

# mutProfileID = MP_mutations[0]['studyId'] + '_mutations'
# smplListID = MP_mutations[0]['studyId'] + '_all'
# muts = cbioportal.Mutations.fetchMutationsInMolecularProfileUsingPOST(
#     molecularProfileId = mutProfileID
#     ).result()

####### This works, but I need the sampleListId: #########
# Right now, I'm guessing at sampleListId.  Better to look up
#####################3
# could write an outer loop to go through all MP_mutations
######################

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

mutDF = pd.DataFrame(columns=['chr','start','end','gene', \
                              'refAllele','varAllele', \
                              'sampleID', 'molecAltType'])
    
# keep track of all unique samples
samples = set()


# response_spec = 'default'               # don't know if this is proper way or not
# response = requests.get(url)
# if response.status_code != 200: #could also check == requests.codes.ok
#    continue

for idx in range(200):  #(len(MPs)): #MP_mutations.keys():
    #mutProfileID = MPs[idx]['studyId'] + '_mutations'  # 'mutations not always correct
    # this should be looked up via:
    # MPs[idx]['molecularProfileId']
    mutProfileID = MPs[idx]['studyId'] + "_mutations"
    smplListID = MPs[idx]['studyId'] + '_all'
    muts = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        # molecularProfileId = mutProfileID,
        molecularProfileId = mutProfileID,
        sampleListId = smplListID,
        entrezGeneId=entrezID,
        # responses='default',
        projection='DETAILED'
        ).result()


    # projection='DETAILED' is essential to have gene-level info
    # grab desired fields and place in dataframe
    for jdx in range(len(muts)):
        #sampleID = muts[jdx]['sampleId']
        #samples.add(sampleID)
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
    print("Molecular Profile #: ", idx, ". Number of ", gene, " mutations in this Mol. Profile= ",len(muts))


##########################################################


# next, remove all cases where either start=='' or end==''                           # UNCOMMENT
#mutDF['refAllele']

###########
# Now we'll groups the mutDF by the first 5 fields
grouped = mutDF.groupby(['chr', 'start', 'end', \
                          'refAllele', 'varAllele', 'gene'])['sampleID']


##########
# Group the mutation data frame, gather unique samples,
# then convert to output dataframe:
##########
# this is working!  Just need to change it back to a dataframe!
outputDF = pd.DataFrame(columns=['chr','start','end','gene', \
                              'refAllele','varAllele','numUniqSamples'])
for name, group in grouped:
    lst = [list(name) + list([len(group.unique())])]
    tmpDF = pd.DataFrame(lst, columns=['chr', 'start', 'end','gene', \
                                        'refAllele', 'varAllele', 'numUniqSamples'])

    outputDF = pd.concat([outputDF, tmpDF], axis=0)
    

                                                                                   ####################3


# next, remove all cases where 
# len(grouped.groups)   # 48 groups for 419 examples
# groupIDs = list(grouped.groups.keys())
# groupIDs[0]
# grouped.groups[groupIDs[0]]

# singleBool = (mutDF['start']==25398284) & (mutDF['end']==25398284) &\
#     (mutDF['refAllele']=='C') & (mutDF['varAllele']=='A')
# # QA:
# mutDF.loc[singleBool,]['sampleID'].sort_values()
# # the command I want is unique().  count() tells me how many total events
# # even if some events are replicated.  ¿What is freq()?



# checking more ideas

studies = cbioportal.Studies.getAllStudiesUsingGET().result()
print("The total number of samples in all studies is: {}".format(sum([x.allSampleCount for x in studies])))
# 91,340 total samples

# total genes:
allGenes = cbioportal.Genes.getAllGenesUsingGET().result()
len(allGenes)


# goodies = {}
# for j in range(len(muts)):
#     if muts[j]['entrezGeneId'] == entrezID:
#         goodies[j] = muts[j]

#     sampleListId = smplListID

# but I want the mutation_event table and not the mutation or the structural_variant
# table!!!

# MP_mutations = {}
# for i in range(len(MPs)):
#     if MPs[i]['molecularAlterationType'] == 'MUTATION_UNCALLED':
#         MP_mutations[i] = MPs[i]
# there is only a single study that is 'MUTATION_UNCALLED
# its ID is 1094 and the study_id is 'glioma_msk_2018'











          
        

# dir(cbioportal)

# dir(cbioportal.Genes)
# cbioportal.Genes.getGeneUsingGET(
#     geneId='672',
#     #geneticEntityId='hg19'
#     ).response().result

# cbioportal.Genes.getGeneUsingGET(geneId="BRCA1").response().result
# # this works, but I get ValidationError: 'geneticEntityId' is a required property




# # dir(cbioportal.Mutations.fetchMutationCountsByPositionUsingPOST)
# # muts = cbioportal.Mutations.fetchMutationCountsByPositionUsingPOST

# allMolecularProfiles = cbioportal.Molecular_Profiles.getAllMolecularProfilesUsingGET().result()
# # type(allMolecularProfiles[0]$molecularProfileId)
# type(allMolecularProfiles)


# cbioportal.Molecular_Profiles.fetchMolecularProfilesUsingPOST()
# # SwaggerMappingError: molecularProfileFilter is a required parameter

# gene_tmp can be used to query mutations directlly
# mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
#     molecularProfileId='msk_impact_2017_mutations',
#     sampleListId='msk_impact_2017_all',
#     projection='DETAILED'
# ).result()
# mutations[0]

# pprint(dir(mutations[0])['_Model__dict']) 
# mutations[-1]


# cbioportal.Genes.getGeneUsingGET.__repr__()


# cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
#     molecularProfileId='msk_impact_2017_mutations',
#     sampleListId='msk_impact_2017_all',
#     projection='DETAILED',
    
#     ).result()

# cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET.__getattribute__()

# from inspect import signature
# cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET.__getattr__()
