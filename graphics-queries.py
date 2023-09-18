import pymongo
import pandas as pd
from gprofiler import GProfiler
import csv
from pycirclize import Circos
from pycirclize.parser import Matrix
import matplotlib.pyplot as plt


client = pymongo.MongoClient("localhost:27017")

db = client["CLINVAR"]
collection_hg19 = db['hg19']
collection_hg38 = db['hg38']

def busqueda_genes(omim, clinical): 

       cursor = collection_hg19.aggregate([{'$match': {"phenotypeids":{'$in':[omim]}, "clinicalsignificance":clinical}},
       {'$group':{"_id":"$genesymbol"}}
       ])
       
       gene_list = []
       for document in cursor:
              for values in document.values():
                     gene_list.append(values)
              
       
       return gene_list

### CÁNCER DE PULMÓN
 #211980

resultado = busqueda_genes("OMIM:211980","Pathogenic")
print(resultado)


gp = GProfiler(return_dataframe=True)
result = gp.profile(organism='hsapiens',
            query=resultado, no_evidences=False)

reac = result[result['source'] == 'REAC'] 
dfplot = pd.DataFrame(columns=["from", "to", "value"])


for index in range(0, len(reac["native"]) - 1) :
      for index2 in range(index + 1, len(reac["native"])) :
           lista = [reac["native"].iloc[index], reac["native"].iloc[index2]]
           intersections1 = reac["intersections"].iloc[index]
           intersections2 = reac["intersections"].iloc[index2]
           listresult = set(intersections1).intersection(intersections2)
           if (len(listresult) > 0) :
              lista.append(len(listresult))
              dfplot.loc[len(dfplot)] = lista

matrix = Matrix.parse_fromto_table(dfplot)

circos = Circos.initialize_from_matrix(
    matrix,
    space=3,
    cmap="tab10",
    label_kws=dict(size=8, r=110, orientation="vertical"),
)

fig = circos.plotfig()
fig.savefig("/tmp/cancer_pulmon.png")

### HEMOFILIA ###

# 306900

resultado2 = busqueda_genes("OMIM:306900","Pathogenic")
print(resultado2)


gp2 = GProfiler(return_dataframe=True)
result2 = gp2.profile(organism='hsapiens',
            query=resultado2, no_evidences=False)

reac2 = result2[result2['source'] == 'REAC'] 
dfplot2 = pd.DataFrame(columns=["from", "to", "value"])


for index in range(0, len(reac2["native"]) - 1) :
      for index2 in range(index + 1, len(reac2["native"])) :
           lista = [reac2["native"].iloc[index], reac2["native"].iloc[index2]]
           intersections1 = reac2["intersections"].iloc[index]
           intersections2 = reac2["intersections"].iloc[index2]
           listresult = set(intersections1).intersection(intersections2)
           if (len(listresult) > 0) :
              lista.append(len(listresult))
              dfplot2.loc[len(dfplot2)] = lista

matrix2 = Matrix.parse_fromto_table(dfplot2)

circos2 = Circos.initialize_from_matrix(
    matrix2,
    space=3,
    cmap="tab10",
    label_kws=dict(size=8, r=110, orientation="vertical"),
)

fig2 = circos2.plotfig()
fig2.savefig("/tmp/hemofilia.png")

### BUSCAMOS VARIANTES PARA HACER EL GRAFICO DE BARRAS

def busqueda_var(chromosome,start,stop, review_status):
    
    chromosome2 = str(chromosome)
    
    lista_review =  ["practice guideline", "reviewed by expert panel", 
                                          "criteria provided, multiple submitters, no conflicts",
                                          "criteria provided, conflicting interpretations",
                                          "criteria provided, single submitter"
                                         ]
    
    #Este condicional nos permite que si nuestro cliente quiere las review desde una estrella, tenga las de una estrella
    #hasta 4 estrellas y viceversa.
    
    #Mapeo a partir de las estrellas que quiere mi cliente, su equivalencia en el reviewstatus de Clinvar
    if review_status == 1:
        review_status = ["criteria provided, single submitter",  "criteria provided, conflicting interpretations"]
    elif review_status == 2:
        review_status = "criteria provided, multiple submitters, no conflicts"
    elif review_status == 3: 
        review_status = "reviewed by expert panel"
    else:
        review_status = "practice guideline"
    
    #Tras mapearlo indico los resultados que quiero que me devuelva.

    if review_status == ["criteria provided, single submitter",  "criteria provided, conflicting interpretations"] :
        result = lista_review
    elif review_status == "criteria provided, multiple submitters, no conflicts" :
        result = lista_review[:-2]
    elif review_status == "reviewed by expert panel" :
        result = lista_review[:-3]
    else:
        result = lista_review[:-4]

    #Consulta a MongoDB
    
    myquery2 ={"chromosome": chromosome2,'start': {'$gte': start}, 'stop': {'$lte': stop},'reviewstatus' : { '$in': result }}

    #Búsqueda de la query en la colección indicada.
    
    mydoc2 = collection_hg38.find(myquery2)
    
    #Output de la función con la información necesaria.
    
    lista_variantes= []

    for x2 in mydoc2:
        lista_variantes.append(x2)
    dataframe_final = pd.DataFrame(lista_variantes)
    return dataframe_final

#GEN CROMOSOMA 17

#GEN TTN EL DEL CROMOSOMA 2
resultado = busqueda_var(2,178527298,178527498,1)
#GEN TTN EL DEL CROMOSOMA 2 desde 2
resultado = busqueda_var(2,178527298,178527498,2)

#GEN CROMOSOMA 17
resultado = busqueda_var(17,47284362,47284662,1)
#GEN CROMOSOMA 17 desde 3
resultado = busqueda_var(17,47284362,47284662,3)

#def grafico_barras(frame, ruta):
    #Cuento las veces que se repite cada significado clínico. 
column = resultado[["clinicalsignificance"]] 
lista_pat = column.values.tolist()

fig, ax = plt.subplots()

categories = ["Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic", "Others"]
colors = ["green", "lightgreen", "lightblue", "#FF9999", "red", "grey"]
values = [0, 0, 0, 0, 0, 0]

for row in lista_pat :    
    categoria = row[0]
    if categoria == categories[0]:
            values[0] = values[0] + 1
    elif categoria == categories[1]:
        values[1] = values[1] + 1 
    elif categoria == categories[2]:
        values[2] = values[2] + 1 
    elif categoria == categories[3]:
        values[3] = values[3] + 1       
    elif categoria == categories[4]:
        values[4] = values[4] + 1 
    else :
        values[5] = values[5] + 1 
            


ax.bar(categories, values, label=categories, color=colors)

ax.set_ylabel('Número de variantes')
ax.set_title('Variantes patogénicas / benignas cercanas a las de interés')

# Personalizar la leyenda
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels, title='Color según el significado clínico', bbox_to_anchor=(1, 0.5), loc='center left')
plt.setp(legend.get_title(), multialignment='center')

    # Rotar las etiquetas del eje x
ax.set_xticklabels(categories, rotation='vertical')

    # Ajustar el espacio entre los subplots para evitar que se superpongan los ejes y la leyenda
plt.subplots_adjust(right=0.85)

    plt.show()

plt.savefig('/tmp/gencr17-desde3.png')

grafico_barras(resultado_var_todas, '/tmp/ttn-todas.png')
grafico_barras(resultado_var_desde3, '/tmp/ttn-desde3.png')
