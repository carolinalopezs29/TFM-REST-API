# TFM-REST-API
Web services for variant prioritization and efficient exploitation of clinical information registered in ClinVar.

The present study seeks to develop a REST Application Programming Interface (REST API) to exploit the Clinvar database, which includes information on genetic variants and enrich its information to know the effect of variants in proteins or genes.

The methodology followed in this study is divided into two parts: (1) to have automatic mechanisms in place to dump and integrate the information from ClinVar and Ensembl Variant Effect Predictor (VEP) and (2) to put in place an efficient REST API to exploit these results. 

In the first part we have used as Python technology and its respective libraries to create a program that integrates ClinVar and VEP and allows to load them into the and VEP and allows to load them into the MongoDB database.

In the second part, Python has also been used to create the REST API, together with bioinformatics tools such as Gprofiler and databases such as REACTOME.

To validate the methodology proposed in the paper, two REST API use cases have been proposed. On one hand, to apply some of the ACMG recommendations for the prioritization of variants, we have consulted the variants close in DNA to the ITGB3 (related to congenital platelet disorders) and TTN (related to heart disease) genes based on their evidence and clinical significance. On the other hand, we have used Clinvar to query implicated genes (based on pathogenic variants identified in these) of two diseases of interest such as lung cancer and hemophilia. In addition, another REST API provided by third parties has been integrated to obtain metabolic pathways in which the obtained genes are involved.
