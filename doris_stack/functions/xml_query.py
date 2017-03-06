# This function performs an xml query on a provided xml file.

import sys
import collections


def xml_query(input_xml):

    try:
        import xml.etree.cElementTree as etree
    except:
        try:
            from lxml import etree
        except:
            #import xml.etree.ElementTree as etree
            print 'Failed to load lxml.etree or xml.etree.cElementTree'
            sys.exit(1)

    inTree = etree.parse(input_xml)

    queryList = collections.OrderedDict([
        ('Volume_file'                                  , 'dummy'),
        ('Volume_ID'                                    , './/adsHeader/missionDataTakeId'),
        ('Volume_identifier'                            , 'dummy'),
        ('Volume_set_identifier'                        , 'dummy'),
        ('Number of records in ref. file'               , 'dummy'),
        ('SAR_PROCESSOR'                                , 'update_1'),
        ('SWATH'                                        , './/adsHeader/swath'),
        ('PASS'                                         , './/generalAnnotation/productInformation/pass'),
        ('IMAGE_MODE'                                   , './/adsHeader/mode'),
        ('polarisation'                                 , './/adsHeader/polarisation'),
        ('Product type specifier'                       , './/adsHeader/missionId'),
        ('Logical volume generating facility'           , 'dummy'),
        ('Location and date/time of product creation'   , 'dummy'),
        ('Number_of_lines_Swath'                        , './/imageAnnotation/imageInformation/numberOfLines'),
        ('number_of_pixels_Swath'                       , './/imageAnnotation/imageInformation/numberOfSamples'),
        ('rangePixelSpacing'                            , './/imageAnnotation/imageInformation/rangePixelSpacing'),
        ('azimuthPixelSpacing'                          , './/imageAnnotation/imageInformation/azimuthPixelSpacing'),
        ('total_Burst'                                  , 'update_1'),
        ('Burst_number_index'                           , 'update_2'),
        ('RADAR_FREQUENCY (HZ)'                         , './/generalAnnotation/productInformation/radarFrequency'),
        ('Scene identification'                         , 'update_1'),
        ('Scene location'                               , 'update_1'),
        ('Sensor platform mission identifer'            , './/adsHeader/missionId'),
        ('Scene_center_heading'                         , './/generalAnnotation/productInformation/platformHeading'),
        ('Scene_centre_latitude'                        , 'update_2'),
        ('Scene_centre_longitude'                       , 'update_2'),
        ('Radar_wavelength (m)'                         , 'update_1'),
        ('Azimuth_steering_rate (deg/s)'                , './/generalAnnotation/productInformation/azimuthSteeringRate'),
        ('Pulse_Repetition_Frequency_raw_data(TOPSAR)'  , './/generalAnnotation/downlinkInformationList/downlinkInformation/prf'),
        ('First_pixel_azimuth_time (UTC)'               , 'update_2'),
        ('Pulse_Repetition_Frequency (computed, Hz)'    , './/imageAnnotation/imageInformation/azimuthFrequency'),
        ('Azimuth_time_interval (s)'                    , './/imageAnnotation/imageInformation/azimuthTimeInterval'),
        ('Total_azimuth_band_width (Hz)'                , './/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/azimuthProcessing/totalBandwidth'),
        ('Weighting_azimuth'                           , './/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/azimuthProcessing/windowType'),
        ('Range_time_to_first_pixel (2way) (ms)'        , 'update_1'),
        ('Range_sampling_rate (computed, MHz)'          , 'update_1'),
        ('Total_range_band_width (MHz)'                 , 'update_1'),
        ('Weighting_range'                              , './/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/windowType'),
        ('DC_reference_azimuth_time'                    , 'update_2'),
        ('DC_reference_range_time'                      , 'update_2'),
        ('Xtrack_f_DC_constant (Hz, early edge)'        , 'update_2'),
        ('Xtrack_f_DC_linear (Hz/s, early edge)'        , 'update_2'),
        ('Xtrack_f_DC_quadratic (Hz/s/s, early edge)'   , 'update_2'),
        ('FM_reference_azimuth_time'                    , 'update_2'),
        ('FM_reference_range_time'                      , 'update_2'),
        ('FM_polynomial_constant_coeff (Hz, early edge)', 'update_2'),
        ('FM_polynomial_linear_coeff (Hz/s, early edge)', 'update_2'),
        ('FM_polynomial_quadratic_coeff (Hz/s/s, early edge)', 'update_2'),
        ('Datafile'                                     , 'update_2'),
        ('Dataformat'                                   , 'update_2'),
        ('Number_of_lines_original'                     , 'update_2'),
        ('Number_of_pixels_original'                    , 'update_2')
    ])

    queryList_aux = collections.OrderedDict([
        ('Swath_startTime'                      , './/adsHeader/startTime'),
        ('Swath_stopTime'                       , './/adsHeader/stopTime'),
        ('imageLines'                           , './/swathTiming/linesPerBurst'),
        ('imagePixels'                          , './/swathTiming/samplesPerBurst'),
        ('firstValidSample'                     , './/swathTiming/burstList/burst/firstValidSample'),
        ('lastValidSample'                      , './/swathTiming/burstList/burst/lastValidSample'),
        ('productSpec'                          , './/generalHeader/referenceDocument'),
        ('productVolDate'                       , './/setup//IOCSAuxProductGenerationTimeUTC'),
        ('productDate'                          , './/generalHeader/generationTime'),
        ('productFacility'                      , './/productInfo/generationInfo/level1ProcessingFacility'),
        ('scenePol'                             , './/adsHeader/polarisation'),
        ('sceneMode'                            , './/adsHeader/mode'),
        ('sceneCenLine_number'                  , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/line'),
        ('sceneCenPixel_number'                 , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/pixel'),
        ('sceneCenLat'                          , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/latitude'),
        ('sceneCenLon'                          , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/longitude'),
        ('height'                               , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/height'),
        ('azimuthTime'                          , './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/azimuthTime'),
        ('sceneRecords'                         , './/imageDataInfo/imageRaster/numberOfRows'),
        ('orbitABS'                             , './/adsHeader/absoluteOrbitNumber'),
        ('orbitTime'                            , './/generalAnnotation/orbitList/orbit/time'),
        ('orbitX'                               , './/generalAnnotation/orbitList/orbit/position/x'),
        ('orbitY'                               , './/generalAnnotation/orbitList/orbit/position/y'),
        ('orbitZ'                               , './/generalAnnotation/orbitList/orbit/position/z'),
        ('rangeRSR'                             , './/generalAnnotation/productInformation/rangeSamplingRate'),
        ('rangeBW'                              , './/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/processingBandwidth'),
        ('rangeTimePix'                         , './/imageAnnotation/imageInformation/slantRangeTime'),
        ('azimuthTimeStart'                     , './/swathTiming/burstList/burst/azimuthTime'),
        ('heading'                              , './/generalAnnotation/productInformation/platformHeading'),
        ('doppler_azimuth_Time'                 , './/dopplerCentroid/dcEstimateList/dcEstimate/azimuthTime'),
        ('doppler_range_Time'                   , './/dopplerCentroid/dcEstimateList/dcEstimate/t0'),
        ('dopplerCoeff'                         , './/dopplerCentroid/dcEstimateList/dcEstimate/dataDcPolynomial'),
        ('azimuthFmRate_reference_Azimuth_time' , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthTime'),
        ('azimuthFmRate_reference_Range_time'   , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/t0'),
        ('azimuthFmRate_c0'                     , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/c0'),
        ('azimuthFmRate_c1'                     , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/c1'),
        ('azimuthFmRate_c2'                     , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/c2'),
        ('azimuthFmRatePolynomial'              , './/generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthFmRatePolynomial')

    ])

    # Now find the variables of queryList and queryList_aux in the xml data.

    for key in queryList.keys():
        try:
            vars()[key]
        except KeyError or NameError:
            vars()[key] = []

        if queryList[key][0:3] == './/':
            for nodes in inTree.findall(queryList[key]):
                vars()[key].append(nodes.text)
            dat = vars()[key]
            if isinstance(dat,list):
                dat = dat[0]
            if isinstance(dat,(int,float)):
                dat = str(dat)
            queryList[key] = dat

    for key in queryList_aux.keys():
        try:
            vars()[key]
        except KeyError or NameError:
            vars()[key] = []

        for nodes in inTree.findall(queryList_aux[key]):
            vars()[key].append(nodes.text)
        queryList_aux[key] = vars()[key]

    # Finally do the first update
    queryList['SAR_PROCESSOR'] = 'Sentinel-' + queryList['Sensor platform mission identifer'][-2:]
    queryList['total_Burst'] = str(len(queryList_aux['azimuthTimeStart']))
    queryList['Scene identification'] = 'Orbit: '+ queryList_aux['orbitABS'][0]
    queryList['Scene location'] = 'lat: ' + queryList_aux['sceneCenLat'][0] + ' lon:' + queryList_aux['sceneCenLon'][0]
    queryList['Radar_wavelength (m)'] = "{:.9f}".format(299792458.0/float(queryList['RADAR_FREQUENCY (HZ)']))
    queryList['Range_time_to_first_pixel (2way) (ms)'] = "{:.15f}".format(float(queryList_aux['rangeTimePix'][0])*1000)
    queryList['Range_sampling_rate (computed, MHz)'] = "{:.9f}".format(float(queryList_aux['rangeRSR'][0])/1000000)
    queryList['Total_range_band_width (MHz)'] = "{:.9f}".format(float(queryList_aux['rangeBW'][0])/1000000)

    deldict = ['orbitABS','rangeBW']
    for d in deldict:
        queryList_aux.pop(d)

    queryList['aux'] = queryList_aux

    return queryList
