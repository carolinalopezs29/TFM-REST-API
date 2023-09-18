import sys
import pandas as pd
from pymongo import MongoClient

import json


rute = sys.argv[1]
bbddip = sys.argv[2]
client = MongoClient(bbddip, 27017)
db = client["CLINVAR"]
collection_hg19 = db['hg19']
collection_hg38 = db['hg38']

def obtener_bd_segura():
    host = "marie"
    puerto = "27017"
    usuario = "homerJ"
    palabra_secreta = "zlevnz24"
    base_de_datos = "phenopackets"
    cliente = MongoClient("mongodb://{}:{}@{}:{}".format(usuario, palabra_secreta, host, puerto))
    return cliente[base_de_datos]

db_segura = obtener_bd_segura()

def queryVEP(db, uplodad_variation) :    
    myquery = [{'$match':{"Uploaded_variation":uplodad_variation}},
        {'$group': {"_id": {"id_variant":"$Uploaded_variation", "gene_ensembl":"$Gene", "rsid": "$Existing_variation"}, 
        "consequences":{'$push':{"gene_transcript":"$Feature", "consequence":"$Consequence"}}}},
        {'$project':{"_id": 0, "id_variant":'$_id.id_variant', "gene_ensembl":'$_id.gene_ensembl', "rsid":'$_id.rsid', "consequences": 1}}]    
    mycol = db["vep"]    
    mydoc = mycol.aggregate(myquery)    
    return mydoc


#Cargamos las 10.000 primeras filas de la base de datos variant summary.
data = pd.read_csv(rute, sep="\t", chunksize=1000000, low_memory=False)

for variants in data:

    #Echamos un vistazo a los datos del dataframe.
    variants.head()

    #Comprobamos que nos quedamos con las 10.000 primeras filas y todas las columnas, 34.
    variants.shape

    #Imprimimos el nombre de las columnas para llevar una idea de los datos que contiene el dataframe.
    variants.columns

    #Vamos a modificar los nombres de las columnas para que no contengan carácteres especiales.

    variants.rename(columns={'#AlleleID':'alleleid' , 'Type':'type', 'Name':'name', 'GeneID':'geneid', 'GeneSymbol':'genesymbol',
                            'HGNC_ID':'hgncid',
        'ClinicalSignificance':'clinicalsignificance', 'ClinSigSimple':'clinsigsimple', 'LastEvaluated':'lastevaluated',
                            'RS# (dbSNP)':'dbsnp',
        'nsv/esv (dbVar)':'dbvar', 'RCVaccession':'rcvaccession', 'PhenotypeIDS':'phenotypeids', 'PhenotypeList':'phenotypelist',
        'Origin':'origin', 'OriginSimple':'originsimple', 'Assembly':'assembly', 'ChromosomeAccession':'chromosomeaccession',
        'Chromosome':'chromosome', 'Start':'start', 'Stop':'stop', 'ReferenceAllele':'referenceallele', 'AlternateAllele':'alternateallele',
        'Cytogenetic':'cytogenetic', 'ReviewStatus':'reviewstatus', 'NumberSubmitters':'numbersubmitters', 'Guidelines':'guidelines',
        'TestedInGTR':'testedingtr', 'OtherIDs':'otherids', 'SubmitterCategories':'submittercategories', 'VariationID':'variationid',
        'PositionVCF':'positionvcf', 'ReferenceAlleleVCF':'referenceallelevcf', 'AlternateAlleleVCF':'alternateallelevcf'},
                inplace=True)

    #### Tenemos que separar las variantes con la versión 37 del genoma humano de las que tienen la versión 38.

    #Filtar por el valor de la columna assembly = a GRCh37
    variants37=variants.query('assembly == "GRCh37"')
    #Filtrar por el valor de la columna = a GRCh38

    variants38=variants.query('assembly == "GRCh38"')

    #Una vez separadas las variantes en las versiones del genoma humano, vamos a crear una lista de diccionarios
    #con las variantes de cada una de las versiones.
    variants37dic = variants37.to_dict(orient='records')
    variants38dic = variants38.to_dict(orient='records')

    ##### CREANDO JSON PARA LA VERSION 37 DEL GENOMA HUMANO #####

    #Vamos a meter en listas los campos que lo requieren: 
    # -rcvaccession -> separamos por |
    # -phenotypeids -> separamos por ,
    # -phenotypelist -> separamos por |
    # -origin -> separamos por ;
    # -otherids -> separamos por ,

    for variant in variants37dic:
        #RCVACCESSION
        rcvaccession = variant['rcvaccession']
        list_rcvaccession = rcvaccession.split('|')
        variant['rcvaccession'] = list_rcvaccession
        #PHENOTYPEIDS
        phenotypeids = variant['phenotypeids']
        list_phenotypeids = phenotypeids.split(',')
        list2_phenotypeids = [i for phenotypeid in list_phenotypeids for i in phenotypeid.split("|")]
        variant['phenotypeids'] = list2_phenotypeids
        #PHENOTYPELIST
        phenotypelist = variant['phenotypelist']
        list_phenotypelist = phenotypelist.split('|')
        variant['phenotypelist'] = list_phenotypelist
        #ORIGIN
        origin = variant['origin']
        list_origin = origin.split(';')
        variant['origin'] = list_origin
        #OTHERIDS
        otherids = variant['otherids']
        list_otherids = otherids.split(',')
        variant['otherids'] = list_otherids
        #Cargar información de VEP
        mydoc = queryVEP(db_segura, str(variant['variationid']))
        for docvep in mydoc :                        
            variant.update(docvep)            
            break        
        #print(variant)


        #En las siguientes lineas de código iteramos por las variantes y eliminamos las columnas lastevaluated, 
    # dbvar y guidelines si no contienen información.

    #CUIDADO! si guidelines contiene información debemos meterla como una lista.

    #Además eliminamos la columna assembly de todas las variantes ya que no es necesaria tras separarlas por la 
    #versión del genoma.

    for i in variants37dic:
        if i['lastevaluated'] == "-":
            del i['lastevaluated']

    for i in variants37dic:
        if i['dbvar'] == "-":
            del i['dbvar']

    for i in variants37dic:
        if i['guidelines'] == "-":
            del i['guidelines']
        else:
            guidelines = i['guidelines']
            list_guidelines = guidelines.split(',')
            i['guidelines'] = list_guidelines

    for i in variants37dic:
        del i['assembly']

    #Guardamos un en un archivo json el array de diccionarios.
    with open("variants37.json", "w") as file:
        for i in variants37dic:
            if i['testedingtr'] == "N":
                i['testedingtr'] = False
            else:
                i['testedingtr'] = True
        json.dump(variants37dic, file, indent=3, default=bool)
    
    
    collection_hg19.insert_many(variants37dic)

    ##### CREANDO JSON PARA LA VERSION 38 DEL GENOMA HUMANO #####

    #Seguimos los mismos pasos que para la versión 38.

    for variant in variants38dic:
        #RCVACCESSION
        rcvaccession = variant['rcvaccession']
        list_rcvaccession = rcvaccession.split('|')
        variant['rcvaccession'] = list_rcvaccession
        #PHENOTYPEIDS
        phenotypeids = variant['phenotypeids']
        list_phenotypeids = phenotypeids.split(',')
        list2_phenotypeids = [i for phenotypeid in list_phenotypeids for i in phenotypeid.split("|")]
        variant['phenotypeids'] = list2_phenotypeids
        #PHENOTYPELIST
        phenotypelist = variant['phenotypelist']
        list_phenotypelist = phenotypelist.split('|')
        variant['phenotypelist'] = list_phenotypelist
        #ORIGIN
        origin = variant['origin']
        list_origin = origin.split(';')
        variant['origin'] = list_origin
        #OTHERIDS
        otherids = variant['otherids']
        list_otherids = otherids.split(',')
        variant['otherids'] = list_otherids
        mydoc = queryVEP(db_segura, str(variant['variationid']))
        for docvep in mydoc :                        
            variant.update(docvep)            
            break        

    for i in variants38dic:
        if i['lastevaluated'] == "-":
            del i['lastevaluated']

    for i in variants38dic:
        if i['dbvar'] == "-":
            del i['dbvar']

    for i in variants38dic:
        if i['guidelines'] == "-":
            del i['guidelines']
        else:
            guidelines = i['guidelines']
            list_guidelines = guidelines.split(',')
            i['guidelines'] = list_guidelines

    for i in variants38dic:
        del i['assembly']

    #Guardamos un en un archivo json el array de diccionarios.
    with open("variants38.json", "w") as file:
        for i in variants38dic:
            if i['testedingtr'] == "N":
                i['testedingtr'] = False
            else:
                i['testedingtr'] = True
        json.dump(variants38dic, file, indent=3, default=bool)
    
    
    collection_hg38.insert_many(variants38dic)  
