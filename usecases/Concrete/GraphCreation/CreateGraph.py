from owlready2 import *
import xlrd  ### pip install xlrd==1.2.0
import csv



Compression = 'usecases/Concrete/Data/Druckfestigkeit'
Tensile = 'usecases/Concrete/Data/E-modul'
Mixture = 'usecases/Concrete/Data/Mischung'


def get_weight(data):
    with open('{}'.format(data), newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        next(spamreader)
        next(spamreader)
        for row in spamreader:
            print(row)
            if "Masse" in row[0]:
                return(row[-1])


def ConcreteCreation(Mischung, Probekoerper):

    book = xlrd.open_workbook('../Data/{}.xlsx'.format('Versuche'), encoding_override = "utf-8")
    sheet = book.sheet_by_index(0)

    wolf = xlrd.open_workbook('../Data/Mischungen/2014_12_10 Wolf.xls', encoding_override = "utf-8")
    WolfRezept = wolf.sheet_by_index(1)

    scharfelanke = xlrd.open_workbook('../Data/Mischungen/2014_03_11 Reprofilierungsplatten_Scharfe_Lanke.xls', encoding_override = "utf-8")
    scharfelankeRezept = scharfelanke.sheet_by_index(1)

    Werner = xlrd.open_workbook('../Data/Mischungen/2014_08_04 Rezepturen_auf 85 Liter_Werner_Losert.xlsx', encoding_override = "utf-8")
    WernerM1Rezept = Werner.sheet_by_index(2)
    WernerM2Rezept = Werner.sheet_by_index(3)
    WernerM3Rezept = Werner.sheet_by_index(4)
    WernerM4Rezept = Werner.sheet_by_index(5)
    WernerM5Rezept = Werner.sheet_by_index(6)
    WernerM6Rezept = Werner.sheet_by_index(7)


    M1 = ['BA-Losert MI E-Modul 28d v. 04.08.14 Probe 6', 'BA-Losert MI E-Modul 28d v. 04.08.14 Probe 5', 'BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4']
    M2 = ['BA-Losert MII E-Modul 28d v. 04.08.14 Probe 6', 'BA-Losert MII E-Modul 28d v. 04.08.14 Probe 4', 'BA-Losert MII E-Modul 28d v. 04.08.14 Probe 5-']
    M3 = ['Werner 7.0 M III E-Modul 28d v. 06.08.14 Probe 5', 'Werner 7.0 M III E-Modul 28d v. 06.08.14 Probe 4', 'Werner 7.0 M III E-Modul 28d v. 06.08.14 Probe 6']
    M4 = ['Werner 7.0 M IV E-Modul 28d v. 06.08.14 Probe 5', 'Werner 7.0 M IV E-Modul 28d v. 06.08.14 Probe 4', 'Werner 7.0 M IV E-Modul 28d v. 06.08.14 Probe 6']
    M5 = ['BA Los M V-4', 'BA Los M V-5', 'BA Los M V-6']
    M6 = ['BA Los M VI-6', 'BA Los M VI-4', 'BA Los M VI-5-']

    if Mischung == "2014_12_10 Wolf":
        worksheet = WolfRezept
    if Mischung == "2014_03_11 Reprofilierungsplatten_Scharfe_Lanke":
        worksheet = scharfelankeRezept
    if Probekoerper in M1:
        worksheet = WernerM1Rezept
    if Probekoerper in M2:
        worksheet = WernerM2Rezept
    if Probekoerper in M3:
        worksheet = WernerM3Rezept
    if Probekoerper in M4:
        worksheet = WernerM4Rezept
    if Probekoerper in M5:
        worksheet = WernerM5Rezept
    if Probekoerper in M6:
        worksheet = WernerM6Rezept

    ConcreteCreation = CST.ActOfCreation()

    CementInput = CST.Cement()

    CementWeight = CCO.Weight()
    CementWeightMeasurement = CCO.MeasurementInformationContentEntity()
    CementWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(21,2).value])
    CementWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
    CementWeightMeasurement.RO_0010001.append(CementWeightValue)
    CementWeight.is_measured_by.append(CementWeightMeasurement)
    CementInput.RO_0000086.append(CementWeight)

    CementVolume = CST.Volume()
    CementVolumeMeasurement = CCO.MeasurementInformationContentEntity()
    CementVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(21,6).value])
    CementVolumeValue.uses_measurement_unit.append(CCO.CubicDecimeterMeasurementUnit)
    CementVolumeMeasurement.RO_0010001.append(CementVolumeValue)
    CementVolume.is_measured_by.append(CementVolumeMeasurement)
    CementInput.RO_0000086.append(CementVolume)


    WaterInput = CST.Water()

    WaterWeight = CCO.Weight()
    WaterWeightMeasurement = CCO.MeasurementInformationContentEntity()
    WaterWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(22,2).value])
    WaterWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
    WaterWeightMeasurement.RO_0010001.append(WaterWeightValue)
    WaterWeight.is_measured_by.append(WaterWeightMeasurement)
    WaterInput.RO_0000086.append(WaterWeight)

    WaterVolume = CST.Volume()
    WaterVolumeMeasurement = CCO.MeasurementInformationContentEntity()
    WaterVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(22,6).value])
    WaterVolumeValue.uses_measurement_unit.append(CCO.CubicDecimeterMeasurementUnit)
    WaterVolumeMeasurement.RO_0010001.append(WaterVolumeValue)
    WaterVolume.is_measured_by.append(WaterVolumeMeasurement)
    WaterInput.RO_0000086.append(WaterVolume)


    AirInput = CST.Air()
    AirVolume = CST.Volume()
    AirVolumeMeasurement = CCO.MeasurementInformationContentEntity()
    AirVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(24,6).value])
    AirVolumeValue.uses_measurement_unit.append(CCO.CubicDecimeterMeasurementUnit)
    AirVolumeMeasurement.RO_0010001.append(AirVolumeValue)
    AirVolume.is_measured_by.append(AirVolumeMeasurement)
    AirInput.RO_0000086.append(AirVolume)

    if Mischung == "2014_03_11 Reprofilierungsplatten_Scharfe_Lanke":
        WaterAbsorbed = CST.AbsorbedAmount()
        WaterAbsorbedMeasurement = CCO.MeasurementInformationContentEntity()
        WaterAbsorbedValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(23,2).value])
        WaterAbsorbedValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        WaterAbsorbedMeasurement.RO_0010001.append(WaterAbsorbedValue)
        WaterAbsorbed.is_measured_by.append(WaterAbsorbedMeasurement)
        WaterInput.RO_0000086.append(WaterAbsorbed)



    if Mischung == "2014_03_11 Reprofilierungsplatten_Scharfe_Lanke":
        ZuschlagInput = CST.Grain()

        ZuschlagInputWeight = CCO.Weight()
        ZuschlagInputWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,2).value])
        ZuschlagInputWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagInputWeightMeasurement.RO_0010001.append(ZuschlagInputWeightValue)
        ZuschlagInputWeight.is_measured_by.append(ZuschlagInputWeightMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputWeight)

        ZuschlagInputVolume = CST.Volume()
        ZuschlagInputVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,6).value])
        ZuschlagInputVolumeValue.uses_measurement_unit.append(CCO.CubicDecimeterMeasurementUnit)
        ZuschlagInputVolumeMeasurement.RO_0010001.append(ZuschlagInputVolumeValue)
        ZuschlagInputVolume.is_measured_by.append(ZuschlagInputVolumeMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputVolume)


        #----------------------------------------------------------------
        ZuschlagOne = CST.Grain()

        ZuschlagOneWeight = CCO.Weight()
        ZuschlagOneWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,2).value])
        ZuschlagOneWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagOneWeightMeasurement.RO_0010001.append(ZuschlagOneWeightValue)
        ZuschlagOneWeight.is_measured_by.append(ZuschlagOneWeightMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneWeight)

        ZuschlagOneVolume = BWMD.BWMD_00276()
        ZuschlagOneVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,7).value / 100])
        ZuschlagOneVolumeMeasurement.RO_0010001.append(ZuschlagOneVolumeValue)
        ZuschlagOneVolume.is_measured_by.append(ZuschlagOneVolumeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneVolume)

        ZuschlagOneMinSize = CST.MinimumDiameter()
        ZuschlagOneMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,0).value.split('/')[0].replace(',','.')])
        ZuschlagOneMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagOneMinSizeMeasurement.RO_0010001.append(ZuschlagOneMinSizeValue)
        ZuschlagOneMinSize.is_measured_by.append(ZuschlagOneMinSizeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneMinSize)

        ZuschlagOneMaxSize = CST.MaximumDiameter()
        ZuschlagOneMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,0).value.split('/')[1].replace(',','.')])
        ZuschlagOneMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagOneMaxSizeMeasurement.RO_0010001.append(ZuschlagOneMaxSizeValue)
        ZuschlagOneMaxSize.is_measured_by.append(ZuschlagOneMaxSizeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagOne)


        #------------------------------------------
        ZuschlagTwo = CST.Grain()

        ZuschlagTwoWeight = CCO.Weight()
        ZuschlagTwoWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,2).value])
        ZuschlagTwoWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagTwoWeightMeasurement.RO_0010001.append(ZuschlagTwoWeightValue)
        ZuschlagTwoWeight.is_measured_by.append(ZuschlagTwoWeightMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoWeight)

        ZuschlagTwoVolume = BWMD.BWMD_00276()
        ZuschlagTwoVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,7).value / 100])
        ZuschlagTwoVolumeMeasurement.RO_0010001.append(ZuschlagTwoVolumeValue)
        ZuschlagTwoVolume.is_measured_by.append(ZuschlagTwoVolumeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoVolume)

        ZuschlagTwoMinSize = CST.MinimumDiameter()
        ZuschlagTwoMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,0).value.split('/')[0].replace(',','.')])
        ZuschlagTwoMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagTwoMinSizeMeasurement.RO_0010001.append(ZuschlagTwoMinSizeValue)
        ZuschlagTwoMinSize.is_measured_by.append(ZuschlagTwoMinSizeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoMinSize)

        ZuschlagTwoMaxSize = CST.MaximumDiameter()
        ZuschlagTwoMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,0).value.split('/')[1].replace(',','.')])
        ZuschlagTwoMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagTwoMaxSizeMeasurement.RO_0010001.append(ZuschlagTwoMaxSizeValue)
        ZuschlagTwoMaxSize.is_measured_by.append(ZuschlagTwoMaxSizeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagTwo)
        #-----------------------------------------------
        ZuschlagThree = CST.Grain()

        ZuschlagThreeWeight = CCO.Weight()
        ZuschlagThreeWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,2).value])
        ZuschlagThreeWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagThreeWeightMeasurement.RO_0010001.append(ZuschlagThreeWeightValue)
        ZuschlagThreeWeight.is_measured_by.append(ZuschlagThreeWeightMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeWeight)

        ZuschlagThreeVolume = BWMD.BWMD_00276()
        ZuschlagThreeVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,7).value / 100])
        ZuschlagThreeVolumeMeasurement.RO_0010001.append(ZuschlagThreeVolumeValue)
        ZuschlagThreeVolume.is_measured_by.append(ZuschlagThreeVolumeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeVolume)

        ZuschlagThreeMinSize = CST.MinimumDiameter()
        ZuschlagThreeMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,0).value.split('/')[0].replace(',','.')])
        ZuschlagThreeMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagThreeMinSizeMeasurement.RO_0010001.append(ZuschlagThreeMinSizeValue)
        ZuschlagThreeMinSize.is_measured_by.append(ZuschlagThreeMinSizeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeMinSize)

        ZuschlagThreeMaxSize = CST.MaximumDiameter()
        ZuschlagThreeMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,0).value.split('/')[1].replace(',','.')])
        ZuschlagThreeMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagThreeMaxSizeMeasurement.RO_0010001.append(ZuschlagThreeMaxSizeValue)
        ZuschlagThreeMaxSize.is_measured_by.append(ZuschlagThreeMaxSizeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagThree)

        #--------------------------------
        ZuschlagFour = CST.Grain()

        ZuschlagFourWeight = CCO.Weight()
        ZuschlagFourWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,2).value])
        ZuschlagFourWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagFourWeightMeasurement.RO_0010001.append(ZuschlagFourWeightValue)
        ZuschlagFourWeight.is_measured_by.append(ZuschlagFourWeightMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourWeight)

        ZuschlagFourVolume = BWMD.BWMD_00276()
        ZuschlagFourVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,7).value / 100])
        ZuschlagFourVolumeMeasurement.RO_0010001.append(ZuschlagFourVolumeValue)
        ZuschlagFourVolume.is_measured_by.append(ZuschlagFourVolumeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourVolume)

        ZuschlagFourMinSize = CST.MinimumDiameter()
        ZuschlagFourMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,0).value.split('/')[0].replace(',','.')])
        ZuschlagFourMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFourMinSizeMeasurement.RO_0010001.append(ZuschlagFourMinSizeValue)
        ZuschlagFourMinSize.is_measured_by.append(ZuschlagFourMinSizeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourMinSize)

        ZuschlagFourMaxSize = CST.MaximumDiameter()
        ZuschlagFourMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,0).value.split('/')[1].replace(',','.')])
        ZuschlagFourMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFourMaxSizeMeasurement.RO_0010001.append(ZuschlagFourMaxSizeValue)
        ZuschlagFourMaxSize.is_measured_by.append(ZuschlagFourMaxSizeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagFour)

        #---------------------------------------------
        ZuschlagFive = CST.Grain()

        ZuschlagFiveWeight = CCO.Weight()
        ZuschlagFiveWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(36,2).value])
        ZuschlagFiveWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagFiveWeightMeasurement.RO_0010001.append(ZuschlagFiveWeightValue)
        ZuschlagFiveWeight.is_measured_by.append(ZuschlagFiveWeightMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveWeight)

        ZuschlagFiveVolume = BWMD.BWMD_00276()
        ZuschlagFiveVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(36,7).value / 100])
        ZuschlagFiveVolumeMeasurement.RO_0010001.append(ZuschlagFiveVolumeValue)
        ZuschlagFiveVolume.is_measured_by.append(ZuschlagFiveVolumeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveVolume)

        ZuschlagFiveMinSize = CST.MinimumDiameter()
        ZuschlagFiveMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(36,0).value.split('/')[0].replace(',','.')])
        ZuschlagFiveMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFiveMinSizeMeasurement.RO_0010001.append(ZuschlagFiveMinSizeValue)
        ZuschlagFiveMinSize.is_measured_by.append(ZuschlagFiveMinSizeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveMinSize)

        ZuschlagFiveMaxSize = CST.MaximumDiameter()
        ZuschlagFiveMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(36,0).value.split('/')[1].replace(',','.')])
        ZuschlagFiveMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFiveMaxSizeMeasurement.RO_0010001.append(ZuschlagFiveMaxSizeValue)
        ZuschlagFiveMaxSize.is_measured_by.append(ZuschlagFiveMaxSizeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagFive)

        #-------------------

        ZuschlagInputK_Number = CST.KNumber()
        ZuschlagInputK_NumberMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputK_NumberValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(9,8).value])
        ZuschlagInputK_NumberMeasurement.RO_0010001.append(ZuschlagInputK_NumberValue)
        ZuschlagInputK_Number.is_measured_by.append(ZuschlagInputK_NumberMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputK_Number)

        #-------------------------------------------------

    if Mischung == "2014_12_10 Wolf" or Mischung == "2014_08_04 Rezepturen_auf 85 Liter_Werner_Losert":

        ZuschlagInput = CST.Grain()

        ZuschlagInputWeight = CCO.Weight()
        ZuschlagInputWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(28,2).value])
        ZuschlagInputWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagInputWeightMeasurement.RO_0010001.append(ZuschlagInputWeightValue)
        ZuschlagInputWeight.is_measured_by.append(ZuschlagInputWeightMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputWeight)

        ZuschlagInputVolume = CST.Volume()
        ZuschlagInputVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(28,6).value])
        ZuschlagInputVolumeValue.uses_measurement_unit.append(CCO.CubicDecimeterMeasurementUnit)
        ZuschlagInputVolumeMeasurement.RO_0010001.append(ZuschlagInputVolumeValue)
        ZuschlagInputVolume.is_measured_by.append(ZuschlagInputVolumeMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputVolume)


        #----------------------------------------------------------------
        ZuschlagOne = CST.Grain()

        ZuschlagOneWeight = CCO.Weight()
        ZuschlagOneWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(29,2).value])
        ZuschlagOneWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagOneWeightMeasurement.RO_0010001.append(ZuschlagOneWeightValue)
        ZuschlagOneWeight.is_measured_by.append(ZuschlagOneWeightMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneWeight)

        ZuschlagOneVolume = BWMD.BWMD_00276()
        ZuschlagOneVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(29,7).value / 100])
        ZuschlagOneVolumeMeasurement.RO_0010001.append(ZuschlagOneVolumeValue)
        ZuschlagOneVolume.is_measured_by.append(ZuschlagOneVolumeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneVolume)

        ZuschlagOneMinSize = CST.MinimumDiameter()
        ZuschlagOneMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(29,0).value.split('/')[0].replace(',','.')])
        ZuschlagOneMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagOneMinSizeMeasurement.RO_0010001.append(ZuschlagOneMinSizeValue)
        ZuschlagOneMinSize.is_measured_by.append(ZuschlagOneMinSizeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneMinSize)

        ZuschlagOneMaxSize = CST.MaximumDiameter()
        ZuschlagOneMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagOneMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(29,0).value.split('/')[1].replace(',','.')])
        ZuschlagOneMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagOneMaxSizeMeasurement.RO_0010001.append(ZuschlagOneMaxSizeValue)
        ZuschlagOneMaxSize.is_measured_by.append(ZuschlagOneMaxSizeMeasurement)
        ZuschlagOne.RO_0000086.append(ZuschlagOneMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagOne)


        #------------------------------------------
        ZuschlagTwo = CST.Grain()

        ZuschlagTwoWeight = CCO.Weight()
        ZuschlagTwoWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(30,2).value])
        ZuschlagTwoWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagTwoWeightMeasurement.RO_0010001.append(ZuschlagTwoWeightValue)
        ZuschlagTwoWeight.is_measured_by.append(ZuschlagTwoWeightMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoWeight)

        ZuschlagTwoVolume = BWMD.BWMD_00276()
        ZuschlagTwoVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(30,7).value / 100])
        ZuschlagTwoVolumeMeasurement.RO_0010001.append(ZuschlagTwoVolumeValue)
        ZuschlagTwoVolume.is_measured_by.append(ZuschlagTwoVolumeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoVolume)

        ZuschlagTwoMinSize = CST.MinimumDiameter()
        ZuschlagTwoMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(30,0).value.split('/')[0].replace(',','.')])
        ZuschlagTwoMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagTwoMinSizeMeasurement.RO_0010001.append(ZuschlagTwoMinSizeValue)
        ZuschlagTwoMinSize.is_measured_by.append(ZuschlagTwoMinSizeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoMinSize)

        ZuschlagTwoMaxSize = CST.MaximumDiameter()
        ZuschlagTwoMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagTwoMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(30,0).value.split('/')[1].replace(',','.')])
        ZuschlagTwoMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagTwoMaxSizeMeasurement.RO_0010001.append(ZuschlagTwoMaxSizeValue)
        ZuschlagTwoMaxSize.is_measured_by.append(ZuschlagTwoMaxSizeMeasurement)
        ZuschlagTwo.RO_0000086.append(ZuschlagTwoMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagTwo)
        #-----------------------------------------------
        ZuschlagThree = CST.Grain()

        ZuschlagThreeWeight = CCO.Weight()
        ZuschlagThreeWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,2).value])
        ZuschlagThreeWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagThreeWeightMeasurement.RO_0010001.append(ZuschlagThreeWeightValue)
        ZuschlagThreeWeight.is_measured_by.append(ZuschlagThreeWeightMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeWeight)

        ZuschlagThreeVolume = BWMD.BWMD_00276()
        ZuschlagThreeVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,7).value / 100])
        ZuschlagThreeVolumeMeasurement.RO_0010001.append(ZuschlagThreeVolumeValue)
        ZuschlagThreeVolume.is_measured_by.append(ZuschlagThreeVolumeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeVolume)

        ZuschlagThreeMinSize = CST.MinimumDiameter()
        ZuschlagThreeMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,0).value.split('/')[0].replace(',','.')])
        ZuschlagThreeMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagThreeMinSizeMeasurement.RO_0010001.append(ZuschlagThreeMinSizeValue)
        ZuschlagThreeMinSize.is_measured_by.append(ZuschlagThreeMinSizeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeMinSize)

        ZuschlagThreeMaxSize = CST.MaximumDiameter()
        ZuschlagThreeMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagThreeMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(31,0).value.split('/')[1].replace(',','.')])
        ZuschlagThreeMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagThreeMaxSizeMeasurement.RO_0010001.append(ZuschlagThreeMaxSizeValue)
        ZuschlagThreeMaxSize.is_measured_by.append(ZuschlagThreeMaxSizeMeasurement)
        ZuschlagThree.RO_0000086.append(ZuschlagThreeMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagThree)

        #--------------------------------
        ZuschlagFour = CST.Grain()

        ZuschlagFourWeight = CCO.Weight()
        ZuschlagFourWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,2).value])
        ZuschlagFourWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagFourWeightMeasurement.RO_0010001.append(ZuschlagFourWeightValue)
        ZuschlagFourWeight.is_measured_by.append(ZuschlagFourWeightMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourWeight)

        ZuschlagFourVolume = BWMD.BWMD_00276()
        ZuschlagFourVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,7).value / 100])
        ZuschlagFourVolumeMeasurement.RO_0010001.append(ZuschlagFourVolumeValue)
        ZuschlagFourVolume.is_measured_by.append(ZuschlagFourVolumeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourVolume)

        ZuschlagFourMinSize = CST.MinimumDiameter()
        ZuschlagFourMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,0).value.split('/')[0].replace(',','.')])
        ZuschlagFourMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFourMinSizeMeasurement.RO_0010001.append(ZuschlagFourMinSizeValue)
        ZuschlagFourMinSize.is_measured_by.append(ZuschlagFourMinSizeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourMinSize)

        ZuschlagFourMaxSize = CST.MaximumDiameter()
        ZuschlagFourMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFourMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(32,0).value.split('/')[1].replace(',','.')])
        ZuschlagFourMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFourMaxSizeMeasurement.RO_0010001.append(ZuschlagFourMaxSizeValue)
        ZuschlagFourMaxSize.is_measured_by.append(ZuschlagFourMaxSizeMeasurement)
        ZuschlagFour.RO_0000086.append(ZuschlagFourMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagFour)

        #---------------------------------------------
        ZuschlagFive = CST.Grain()

        ZuschlagFiveWeight = CCO.Weight()
        ZuschlagFiveWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,2).value])
        ZuschlagFiveWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagFiveWeightMeasurement.RO_0010001.append(ZuschlagFiveWeightValue)
        ZuschlagFiveWeight.is_measured_by.append(ZuschlagFiveWeightMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveWeight)

        ZuschlagFiveVolume = BWMD.BWMD_00276()
        ZuschlagFiveVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,7).value / 100])
        ZuschlagFiveVolumeMeasurement.RO_0010001.append(ZuschlagFiveVolumeValue)
        ZuschlagFiveVolume.is_measured_by.append(ZuschlagFiveVolumeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveVolume)

        ZuschlagFiveMinSize = CST.MinimumDiameter()
        ZuschlagFiveMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,0).value.split('/')[0].replace(',','.')])
        ZuschlagFiveMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFiveMinSizeMeasurement.RO_0010001.append(ZuschlagFiveMinSizeValue)
        ZuschlagFiveMinSize.is_measured_by.append(ZuschlagFiveMinSizeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveMinSize)

        ZuschlagFiveMaxSize = CST.MaximumDiameter()
        ZuschlagFiveMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagFiveMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(33,0).value.split('/')[1].replace(',','.')])
        ZuschlagFiveMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagFiveMaxSizeMeasurement.RO_0010001.append(ZuschlagFiveMaxSizeValue)
        ZuschlagFiveMaxSize.is_measured_by.append(ZuschlagFiveMaxSizeMeasurement)
        ZuschlagFive.RO_0000086.append(ZuschlagFiveMaxSize)


        ZuschlagInput.is_made_of.append(ZuschlagFive)
        #------------------------------------------------------------

        ZuschlagSix = CST.Grain()

        ZuschlagSixWeight = CCO.Weight()
        ZuschlagSixWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSixWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,2).value])
        ZuschlagSixWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagSixWeightMeasurement.RO_0010001.append(ZuschlagSixWeightValue)
        ZuschlagSixWeight.is_measured_by.append(ZuschlagSixWeightMeasurement)
        ZuschlagSix.RO_0000086.append(ZuschlagSixWeight)

        ZuschlagSixVolume = BWMD.BWMD_00276()
        ZuschlagSixVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSixVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,7).value / 100])
        ZuschlagSixVolumeMeasurement.RO_0010001.append(ZuschlagSixVolumeValue)
        ZuschlagSixVolume.is_measured_by.append(ZuschlagSixVolumeMeasurement)
        ZuschlagSix.RO_0000086.append(ZuschlagSixVolume)

        ZuschlagSixMinSize = CST.MinimumDiameter()
        ZuschlagSixMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSixMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,0).value.split('/')[0].replace(',','.')])
        ZuschlagSixMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagSixMinSizeMeasurement.RO_0010001.append(ZuschlagSixMinSizeValue)
        ZuschlagSixMinSize.is_measured_by.append(ZuschlagSixMinSizeMeasurement)
        ZuschlagSix.RO_0000086.append(ZuschlagSixMinSize)

        ZuschlagSixMaxSize = CST.MaximumDiameter()
        ZuschlagSixMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSixMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(34,0).value.split('/')[1].replace(',','.')])
        ZuschlagSixMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagSixMaxSizeMeasurement.RO_0010001.append(ZuschlagSixMaxSizeValue)
        ZuschlagSixMaxSize.is_measured_by.append(ZuschlagSixMaxSizeMeasurement)
        ZuschlagSix.RO_0000086.append(ZuschlagSixMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagSix)
        #---------------------------------------------------------------------

        ZuschlagSeven = CST.Grain()

        ZuschlagSevenWeight = CCO.Weight()
        ZuschlagSevenWeightMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSevenWeightValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,2).value])
        ZuschlagSevenWeightValue.uses_measurement_unit.append(CCO.KilogramMeasurementUnit)
        ZuschlagSevenWeightMeasurement.RO_0010001.append(ZuschlagSevenWeightValue)
        ZuschlagSevenWeight.is_measured_by.append(ZuschlagSevenWeightMeasurement)
        ZuschlagSeven.RO_0000086.append(ZuschlagSevenWeight)

        ZuschlagSevenVolume = BWMD.BWMD_00276()
        ZuschlagSevenVolumeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSevenVolumeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,7).value / 100])
        ZuschlagSevenVolumeMeasurement.RO_0010001.append(ZuschlagSevenVolumeValue)
        ZuschlagSevenVolume.is_measured_by.append(ZuschlagSevenVolumeMeasurement)
        ZuschlagSeven.RO_0000086.append(ZuschlagSevenVolume)

        ZuschlagSevenMinSize = CST.MinimumDiameter()
        ZuschlagSevenMinSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSevenMinSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,0).value.split('/')[0].replace(',','.')])
        ZuschlagSevenMinSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagSevenMinSizeMeasurement.RO_0010001.append(ZuschlagSevenMinSizeValue)
        ZuschlagSevenMinSize.is_measured_by.append(ZuschlagSevenMinSizeMeasurement)
        ZuschlagSeven.RO_0000086.append(ZuschlagSevenMinSize)

        ZuschlagSevenMaxSize = CST.MaximumDiameter()
        ZuschlagSevenMaxSizeMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagSevenMaxSizeValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(35,0).value.split('/')[1].replace(',','.')])
        ZuschlagSevenMaxSizeValue.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
        ZuschlagSevenMaxSizeMeasurement.RO_0010001.append(ZuschlagSevenMaxSizeValue)
        ZuschlagSevenMaxSize.is_measured_by.append(ZuschlagSevenMaxSizeMeasurement)
        ZuschlagSeven.RO_0000086.append(ZuschlagSevenMaxSize)

        ZuschlagInput.is_made_of.append(ZuschlagSeven)

        #------------------------------------------------------------------

        ZuschlagInputK_Number = CST.KNumber()
        ZuschlagInputK_NumberMeasurement = CCO.MeasurementInformationContentEntity()
        ZuschlagInputK_NumberValue = BWMD.BWMD_00342(has_decimal_value = [worksheet.cell(9,8).value])
        ZuschlagInputK_NumberMeasurement.RO_0010001.append(ZuschlagInputK_NumberValue)
        ZuschlagInputK_Number.is_measured_by.append(ZuschlagInputK_NumberMeasurement)
        ZuschlagInput.RO_0000086.append(ZuschlagInputK_Number)

    ConcreteCreation.has_input = [WaterInput,CementInput,ZuschlagInput,AirInput]

    return(ConcreteCreation)





def include_Experiment(Experiment):


    test_DATFile = BWMD.BWMD_00032()



    if Experiment[7].value == 1:
        Experiment_individual = CST.ConcreteTensileTest()
        test_DATFile.has_URI_value = [Tensile + '/{}/specimen.dat'.format(Experiment[2].value)]
    if Experiment[8].value == 1:
        Experiment_individual = CST.ConcreteCompressionTest()
        test_DATFile.has_URI_value = [Compression + '/{}/specimen.dat'.format(Experiment[2].value)]


    Stamp = WCTmid.Stamp()
    Stamp.RO_0000086 = []


    Timepoint_pre = WCTmid.PointInTime( has_text_value = ["Pre"])
    Timepoint_post = WCTmid.PointInTime( has_text_value = ["Post"])


    test_Sample = BWMD.BWMD_00048()
    test_Sample.RO_0000086 = []
    test_Sample.RO_0000086.append(CCO.Cylindrical())


    test_Sample_hight = CCO.Height()
    height_Measurement = CCO.RatioMeasurementInformationContentEntity()
    height_Value = BWMD.BWMD_00342(has_decimal_value = [Experiment[5].value])
    height_Value.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
    height_Measurement.RO_0010001.append(height_Value)
    test_Sample_hight.is_measured_by.append(height_Measurement)
    test_Sample_hight.is_output_of.append(CCO.StasisOfQuality( occurs_on = [Timepoint_pre]))
    test_Sample.RO_0000086.append(test_Sample_hight)


    test_Sample_diameter = CCO.Diameter()
    diameter_Measurement = CCO.RatioMeasurementInformationContentEntity()
    diameter_Value = BWMD.BWMD_00342(has_decimal_value = [Experiment[4].value])
    diameter_Value.uses_measurement_unit.append(CCO.MillimeterMeasurementUnit)
    diameter_Measurement.RO_0010001.append(diameter_Value)
    test_Sample_diameter.is_measured_by.append(diameter_Measurement)
    test_Sample_diameter.is_output_of.append(CCO.StasisOfQuality( occurs_on = [Timepoint_pre]))
    test_Sample.RO_0000086.append(test_Sample_diameter)


    Test_sample_weight = CCO.Weight()
    weight_Measurement = CCO.RatioMeasurementInformationContentEntity()
    weight_Value = BWMD.BWMD_00342(has_decimal_value = [get_weight("../.."+test_DATFile.has_URI_value[0][17:])])
    weight_Value.uses_measurement_unit.append(CCO.GramMeasurementUnit)
    weight_Measurement.RO_0010001.append(weight_Value)
    Test_sample_weight.is_measured_by.append(weight_Measurement)
    Test_sample_weight.is_output_of.append(CCO.StasisOfQuality( occurs_on = [Timepoint_pre]))
    test_Sample.RO_0000086.append(Test_sample_weight)

    test_Concrete = CST.Concrete()

    if Experiment[0].value != '':
        Concrete_data = CCO.QualitySpecification()
        Concrete_data_file = BWMD.BWMD_00032()
        Concrete_data_file.has_URI_value = [Mixture + '/{}'.format(Experiment[0].value)]
        Concrete_data.RO_0010001 = [Concrete_data_file]
        test_Concrete.prescribed_by = [Concrete_data]
        print(Experiment[0], Experiment[2])
        test_Concrete.is_output_of = [ConcreteCreation(Experiment[0].value, Experiment[2].value)]

    test_Sample.is_made_of.append(test_Concrete)

    test_Output = BWMD.BWMD_00067()
    test_Output.RO_0010001 = [test_DATFile]


    test_Sample.is_affected_by.append(Experiment_individual)
    test_Output.is_output_of.append(Experiment_individual)
    Stamp.is_object_of.append(Experiment_individual)                          #maybe wanna change these to functional properties --> the objects do not have to be in a list

    Experiment_individual.RO_0000057 = [test_Output, Stamp, test_Sample]




#file_errors_location = '../Data/{}.xlsx'.format('Versuche')
#df = pd.read_excel(file_errors_location)
#print(df)

#print('../Data/{}.xlsx'.format('Versuche'))


onto_path.append(".")
My_world = World()

ccoonto = My_world.get_ontology("Ontologies/MergedAllCoreOntology.owl").load()
onto = My_world.get_ontology("Ontologies/MSEO_mid.owl").load()
PTonto = My_world.get_ontology("Ontologies/PeriodicTable.owl").load()
WCTmidonto = My_world.get_ontology("Ontologies/WCTmid.owl").load()
CSTonto = My_world.get_ontology("Ontologies/ConcreteStressTestOntologie.owl").load()

BWMD = My_world.get_namespace("https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#")
WCTmid = My_world.get_namespace("https://mobi.com/ontologies/6/2021/WCTmid#")
WCT = My_world.get_namespace("https://mobi.com/ontologies/6/202https://mobi.com/ontologies/6/2021/WoodCompressionTest#1/WoodCompressionTest#")
CCO = My_world.get_namespace("http://www.ontologyrepository.com/CommonCoreOntologies/")
CST = My_world.get_namespace("https://mobi.com/ontologies/7/2021/ConcreteStressTestOntologie#")
onto.imported_ontologies.append(ccoonto)
onto.imported_ontologies.append(PTonto)
onto.imported_ontologies.append(WCTmidonto)
onto.imported_ontologies.append(CSTonto)

onto.save('Concrete.owl')

book = xlrd.open_workbook('../Data/{}.xlsx'.format('Versuche'), encoding_override = "utf-8")
sheet = book.sheet_by_index(0)
for i in range(2, sheet.nrows):
    Experiment = []
    for j in range(10):
        Experiment.append(sheet.cell(i,j))
    with onto:
        include_Experiment(Experiment)
onto.save('Concrete_data.owl')
onto.save('Concrete_data.rdf', format = 'rdfxml')
