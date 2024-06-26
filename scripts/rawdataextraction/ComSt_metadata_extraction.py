import re
import json
import uuid
import datetime


# the function read each line and return metadata as key and value
def get_metadata_in_one_line(line):
    s = re.sub('\t+', '\t', line)
    s = s.replace('\n', '\t')
    result = s.split('\t')[:-1]
    return result


# function to convert german formatting to english
def replace_comma(string):
    string = string.replace(',', '.')
    return string


def extract_metadata_ComSt(blob,mix_json,processed_data):
    """Returns two dictionaries: one with extracted emodule-metadata and one with
    extracted specimen metadata.

    Parameters
    ----------
    blob : dat BLOB file from db
    mix_json : Json from the Mixture used (binary format)
    processed_data : dataFrame from ComSt_generate_processed_data

    Returns
    -------
    [metadata_ComSt, metadata_specimen_ComSt] : Array containing metadata of compressive strength and specimen

    """

    # create empty dictionary for metadata
    metadata_ComSt = {}
    metadata_specimen_ComSt = {}

    # Write the blob data to a .dat file
    specimen_file = 'specimen.dat'
    with open(specimen_file, 'wb') as file:
        file.write(blob)


    # read specimen file
    with open(specimen_file, 'r', encoding="utf8", errors='ignore') as data:
        # read data file line by line
        lines = data.readlines()

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0] + 1, emptyLineIndex[1] + 4):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))

        ###########  D A T A   A B O U T    E X P E R I M E N T  #######

        # humanreadable ID = name of the data file
        metadata_ComSt['humanreadableID'] = serviceInformation[3][1]

        # ID of this experiment
        ComStID = str(uuid.uuid4())
        metadata_ComSt['ID'] = ComStID

        # get experiment date and time in protege format YYYY-MM-DDTHH:mm:SS
        date = serviceInformation[11][4]  # datetime.datetime.strptime(,'%d.%m.%y')
        date_only = datetime.datetime.strptime(date.split(" ")[0], '%d.%m.%Y')
        date_protegeformat = date_only.strftime('%Y-%m-%d') + "T" + date.split(" ")[1]
        metadata_ComSt['ExperimentDate'] = str(date_protegeformat)

        # operator name - This data has no placeholder yet.
        # metadata_ComSt['tester_name'] = serviceInformation[2][1]

        # remarks - This data has no placeholder yet.
        # metadata_ComSt['remark'] = serviceInformation[4][1]

        # set experiment lab location to BAM
        metadata_ComSt['Lab'] = 'BAM'

        # set Compression and Transducer Column
        metadata_ComSt['CompressionColumn'] = [4]
        metadata_ComSt['CompressionForce_Unit'] = "kN"
        metadata_ComSt['TransducerColumn'] = [5]
        metadata_ComSt['Extensometer_Unit'] = "mm"  # Transducer messen eine Verschiebung.


        # name of specimen (humanreadable)
        metadata_specimen_ComSt['humanreadableID'] = serviceInformation[3][1]
        # set size of specimen
        metadata_specimen_ComSt['SpecimenDiameter'] = float(replace_comma(serviceInformation[5][1]))  # diameter
        metadata_specimen_ComSt['SpecimenDiameter_Unit'] = 'mm'
        metadata_specimen_ComSt['SpecimenHeight'] = float(replace_comma(serviceInformation[6][1]))  # Height
        metadata_specimen_ComSt['SpecimenHeight_Unit'] = 'mm'

        # set the length
        metadata_specimen_ComSt['SpecimenLength'] = float(replace_comma(serviceInformation[7][1]))  # Length
        metadata_specimen_ComSt['SpecimenLength_Unit'] = 'mm'

        # weight
        metadata_specimen_ComSt['SpecimenMass'] = float(replace_comma(serviceInformation[8][1]))
        metadata_specimen_ComSt['SpecimenMass_Unit'] = 'g'

        

        # ID of this specimen
        #specimenID = str(uuid.uuid4())
        metadata_ComSt['specimenID'] = metadata_specimen_ComSt['ID'] = ComStID
        # save Mixdesign ID to specimen metadata
        try:
            mixdesign = json.loads(mix_json.decode('utf-8'))
            mixtureID = mixdesign['ID']
            metadata_specimen_ComSt['MixtureID'] = mixtureID
            # Extract mixing date
            mixing_date_str = mixdesign['MixingDate']
            mixing_date = datetime.datetime.strptime(mixing_date_str, '%Y-%m-%dT%H:%M:%S')
        except AttributeError:
            raise Exception("No mixdesign json-file found! Can't import the ID and save it to the output!")

        # Calculate specimen age
        if 'ExperimentDate' in metadata_ComSt and mixing_date:
            # Convert dates to midnight
            experiment_date = datetime.datetime.strptime(metadata_ComSt['ExperimentDate'], '%Y-%m-%dT%H:%M:%S').replace(
                hour=0, minute=0, second=0)
            mixing_date = mixing_date.replace(hour=0, minute=0, second=0)
            # Calculate age
            specimen_age = (experiment_date - mixing_date).days
            metadata_ComSt['SpecimenAge'] = specimen_age
            metadata_ComSt['SpecimenAge_Unit'] = 'day'

        # set shape
        if 'SpecimenHeight' in metadata_specimen_ComSt:
            metadata_specimen_ComSt['SpecimenShape'] = 'Cube'

        # set paths
        #metadata_ComSt['ProcessedFile'] = os.path.join('../../../usecases/MinimumWorkingExample/Druckfestigkeit/processeddata')  # path to csv file with values extracted by ComSt_generate_processed_data.py
        metadata_ComSt['RawDataFile'] = "Download"

        try:
            #my_data = pd.read_csv(metadata_ComSt['ProcessedFile'])
            # print(my_data)
            min_force = processed_data['Force [kN]'].min()
            #diameter = 100.0
            #height = 100.3
            area = metadata_specimen_ComSt['SpecimenDiameter'] ** 2
            normalValue = (min_force * -1)
            CompressiveStrength = (normalValue / area) * 1000
            metadata_ComSt['CompressiveStrength'] = CompressiveStrength
        except:

            raise Exception("No processed_file found!")

        metadata_ComSt['CompressiveStrength_Unit'] = "GPa"

    return [metadata_ComSt, metadata_specimen_ComSt]

