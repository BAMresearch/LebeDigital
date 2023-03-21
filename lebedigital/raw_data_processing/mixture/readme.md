There are two scripts for extracting mixture metadata:

# The first script: mixture_metadata_extraction.py

This script was the first one and is designed for the raw data we have ("..\usecases\MinimumWorkingExample\Data\Mischungen").
This raw data contains only one type of cement and up to two additions. We don't use the resulting metadata-yaml-file at the 
moment for any ontology. 

# The second script: mixdesign_metadata_extraction.py

This script is made to make the existing raw data fit the MixDesign ontology, so that we can produce some
data for the Minimum Working Example. The MixDesign ontology is refering to different types of cement, as 
stated in this [paper](https://www.sciencedirect.com/science/article/pii/S0008884608000884).

Main differences to the raw data that we have: 
1. There are placeholders for Cem I and Cem II, tho our raw data only has one type of cement.
2. No placeholder(s) for addition ('Zusatzstoff'). 
3. Aggregate is always of type "Okrilla". 

The following information is attempted to be extracted from raw data:
- Raw data file name ("$$RawDataFile_Value$$"^^xsd:string)
- Water cement ratio ("$$WaterCementRatio_Value$$"^^xsd:decimal), calculated in the script
- Mixing Date ("$$MixingDate_Value$$"^^xsd:dateTimeStamp)
- Location of the lab ("$$Lab_Value$$"^^xsd:string)
- Cement 1: Quantity and density ("$$CEMIQtyInMix_Value$$"^^xsd:decimal, "$$CEMIDensity_Value$$"^^xsd:decimal)
- Water: Quantity and density ("$$MixingWaterQtyInMix_Value$$"^^xsd:decimal, "$$WaterDensity_Value$$"^^xsd:decimal)
- Aggregates: Quantity and density ("$$OkrillaQtyInMix_Value$$"^^xsd:decimal, "$$OkrillaDensity_Value$$"^^xsd:decimal)
- Admixture: Quantity and density ("$$PlasticizerQtyInMix_Value$$"^^xsd:decimal, "$$PlasticizerDensity_Value$$"^^xsd:decimal)
- AirContent: Quantity and density ("$$AirContent_Value$$"^^xsd:decimal, "$$AirDensity_Value$$"^^xsd:decimal)

Please note:
- The following information is not existent in the raw data: 
Cement 2: Quantity and density ("$$CEMIIQtyInMix_Value$$"^^xsd:decimal, "$$CEMIIDensity_Value$$"^^xsd:decimal)
- For calculating the water-cement-ratio, Cem I is used.
- 'OkrillaDensity' ('Zuschlag (gesamt)') seems to have no value given in most raw data files. Instead the volume is given.
- 'AirContent' seems to have no value given for quantity and density in most raw files. Only the volume is given.