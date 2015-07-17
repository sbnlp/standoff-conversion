import pickle
import os
import sys
import urllib, urllib2
from xml.dom import minidom

#get the Uniprot identifier from Text query
def map_entity(entity_string, organism = '9606', base_url = "http://www.uniprot.org/uniprot/?", special_characters = ["alpha","beta","gamma","kappa"]):
    
    query = { 'query' : entity_string, 'organism' : organism}
    for character in special_characters:
        if character in query['query']:
            query['query'] = ((' ' + character + ' ').join( query['query'].split( character))).strip()
    params={
        'query': "\"" + query['query'] + "\"" + " AND organism:\"" + query['organism'] + "\" AND reviewed:yes",
        'format':'tab',
        'sort':'score',
        }
    url = base_url + urllib.urlencode(params)
    data = urllib.urlopen(url).read()
    return data
    
#From a Uniprot ID get the other details of the protein
def get_details( uniprot_id, base_url = "http://www.uniprot.org/uniprot/"):
    data = { "ensembl_ids" : [], "pdb_ids" : [], "uniprot_ids" : [], "refseq_ids" : [], "alternate_names" : [], "gene_name_synonyms" : [], 'gene_id' : None, "official_name" : None, "gene_name" : None}
    
    url = base_url + uniprot_id + ".xml"
    received_data = minidom.parseString( urllib.urlopen( url).read())
    
    data['official_name'] = received_data.getElementsByTagName('entry')[0].getElementsByTagName('name')[0].childNodes[0].data
    
    for i in received_data.getElementsByTagName('entry')[0].getElementsByTagName('accession'):
        data["uniprot_ids"].append(i.childNodes[0].data)
        
    data["alternate_names"].append(received_data.getElementsByTagName('entry')[0].getElementsByTagName('protein')[0].getElementsByTagName('recommendedName')[0].getElementsByTagName('fullName')[0].childNodes[0].data)
    for i in received_data.getElementsByTagName('entry')[0].getElementsByTagName('protein')[0].getElementsByTagName('alternativeName'):
        data["alternate_names"].append(i.getElementsByTagName('fullName')[0].childNodes[0].data)
    
    data['gene_name'] = received_data.getElementsByTagName('entry')[0].getElementsByTagName('gene')[0].getElementsByTagName('name')[0].childNodes[0].data
    
    for i in received_data.getElementsByTagName('entry')[0].getElementsByTagName('gene')[0].getElementsByTagName('name'):
        data['gene_name_synonyms'].append(i.childNodes[0].data)
    
    for i in received_data.getElementsByTagName('entry')[0].getElementsByTagName('dbReference'):
        if i.getAttribute('type') == 'Ensembl':
            data['ensembl_ids'].append(i.getAttribute('id'))
            for j in i.getElementsByTagName('property'):
                data['ensembl_ids'].append(j.getAttribute('value'))
        if i.getAttribute('type') == 'PDB':
            data['pdb_ids'].append(i.getAttribute('id'))
        if i.getAttribute('type') == 'RefSeq':
            data['refseq_ids'].append(i.getAttribute('id'))
        if i.getAttribute('type') == 'HGNC' or i.getAttribute('type')=='MGI':
            data['gene_id'] = i.getAttribute('id')

    return data

def get_uniprot_data( entity_string, uniprot_pickle = None):
    data = None
    if uniprot_pickle is not None:
        # use the file
        try:
            with open( uniprot_pickle, 'rb') as f:
                data = pickle.load( f)
            if data.get( entity_string):
                return data[entity_string]
        except:
            pass

    #query online
    uniprot_id = ""
    try:
        uniprot_id = map_entity( entity_string)
    except:
        pass;
    details = {}
    if len( uniprot_id)>0:
        try:
            details = get_details( [x.split('\t') for x in uniprot_id.split('\n')][1][0])
        except:
            pass
    
    if data is not None:
        data[entity_string] = details
        with open( uniprot_pickle, 'rb') as f:
                pickle.dump( data, f)
        
    return details
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Retrieve uniprot information.')
    parser.add_argument('string', help="the string used for looking up information")
    parser.add_argument('--use-uniprot-pickle',
                    action="store_true",
                    dest="uniprot_pickle",
                    default=None,
                    help="whether to add uniprot information")
    cmd_args = parser.parse_args()
    

    print "Handling " + cmd_args.string
    uniprot_data = get_uniprot_data( cmd_args.string)
        
    data = {}
    if cmd_args.uniprot_pickle != None and not os.path.isfile( "uniprot.pickle"):
        print "Creating new data"
    elif cmd_args.uniprot_pickle != None: 
        print "Loading existing data"
        with open( cmd_args.uniprot_pickle, 'rb') as f:
            data = pickle.load( f)
    data[cmd_args.string] = uniprot_data
    
    print(uniprot_data)
    
    if cmd_args.uniprot_pickle != None: 
        print "Saving data"
        with open( cmd_args.uniprot_pickle, "wb") as f:
            pickle.dump( data, f)
