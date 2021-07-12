## **What is Done**

1. Data-extraction for the Concrete-usecase
2. Creation of a Graph ready for upload to some triplestore

## **tools needet to run the script**

* Python 3.9 *atleast that is what i am running*
* Owlready2
* xldr
  * Version: 1.2.0 `pip install xlrd==1.2.0`
* csv
  * a library includet in python 3.9

## **tools needet for ontologie editing & creation**

* Protègè
* Python 3.9
* Owlready2
* Mobi `https://mobi.matolab.org`

## **How to create the graph**

1. Run Creategraph.py `py CreateGraph.py`
2. Thats it

## **What is the Output of CreateGraph.py**

1. The graph containing all Experiments from the Concrete-usecase
    * as .rdf `Concrete_data.rdf`
    * as .owl `Concrete_data.owl`
2. The Ontologie used
    * as .owl `Concrete.owl`

All outputs are createt in the "/GraphCreation" directory.      
I Am thinking about moving the `Concrete.owl` file into the "/Ontologies" directory

## **what is inside the Ontologies folder**

All ontologies needet.    
In .ttl and .owl format

* CommonCore
* MSEO
* PeriodicTable
* WCTmid (some extentions to the MSEO)
* ConcreteStressTest (classes specifically relevant to the concrete example )
