#!/usr/bin/python
# this script takes

import re, urllib2
from bs4 import BeautifulSoup

################################################################################
# FUNCTIONS

## A.
def locus_finder(ref_id):
    #grab the GI of the top blast hit and use it to look up the gene ID
    GI = re.match( r'gi\|(\d+)', ref_id, re.M|re.I)
    GI = GI.group(1)
    req = urllib2.Request('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id={}' .format(GI))
    response = urllib2.urlopen(req)
    elink = response.read()
    
    #use geneID to fetch the gene record for the top blast hit (putative homologue for our query sequence)
    soup = BeautifulSoup(elink, 'html.parser')
    ID = soup.elinkresult.linksetdb.link.id
    geneID = re.match( r'\<id\>(\d+)', str(ID), re.M|re.I)
    geneID = geneID.group(1)
    req = urllib2.Request('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&rettype=xml' .format(geneID))
    response = urllib2.urlopen(req)
    efetch = response.read()
    soup = BeautifulSoup(efetch, 'html.parser')
    
    #use the gene record to print the location information for our sequence
    soup_ref_locus = soup.find('gene-ref_locus')
    ref_locus = re.search( r'\>(\w+)\<', str(soup_ref_locus), re.M|re.I)
    ref_locus_name = ref_locus.group(1)
    soup_ref_maploc = soup.find('gene-ref_maploc')
    ref_maploc = re.search( r'\>(\S+)\<', str(soup_ref_maploc), re.M|re.I)
    ref_maploc = ref_maploc.group(1)
    
    return ref_locus_name, ref_maploc


