# File: river_transect.py
# Created by: Ron Dalumpines, turugban@yahoo.com
# Date created: December 8, 2012
# Last revised: June 9, 2014
# Requirements: ArcGIS 10.1 (Advanced License), Python 2.7.
# Description: Extract distance and elevation information via automated river transect (renamed ValleyMorphTool).
# Note: Tested in ArcGIS 10.1; some geoprocessing functions need to be changed for earlier versions from 9.3.
# Known issue: Unlike the standalone (IDLE) version, this version randomly suffers from "ERROR 010213: Error in
# reading raster", which is not fully addressed at the moment. Temporary solution would be to rerun the tool
# from the WID (watershed ID) where the error occurred.

# Imports.
import os
import sys
import arcgisscripting
import math
import traceback
import shutil

# Create the geoprocessor object.
gp = arcgisscripting.create()

# Set the necessary product code.
gp.SetProduct("ArcInfo")

# Check out necessary licenses.
gp.CheckOutExtension("Spatial")
gp.CheckOutExtension("Network")

# Overwrite output.
gp.overwriteoutput = True

def point2line(inPts, outLine, IDField, sortField):
    # Input point FC: inPts
    # Output polyline FC: outLine
    # LineID Field: IDField
    # Sort Field
    if sortField == "#":
        sortField = ""

    if sortField == "":
        cursorSort = IDField
    else:
        cursorSort = IDField + ";" + sortField

    createLinesFromPoints(inPts, outLine, IDField, cursorSort) 

def createLinesFromPoints(inPts, outLine, IDField, cursorSort):
    try:
        # Assign empty values to cursor and row objects
        iCur, sRow, sCur, feat = None, None, None, None

        shapeName = gp.Describe(inPts).ShapeFieldName

        # Create the output feature class
        #
        outPath, outFC = os.path.split(outLine)
        gp.CreateFeatureClass(outPath, outFC, "Polyline", inPts, "", "", inPts)

        # Open an insert cursor for the new feature class
        #
        iCur = gp.InsertCursor(outLine)
        sCur = gp.SearchCursor(inPts, "", None, cursorSort, cursorSort)

        sRow = sCur.Next()

        # Create an array and point object needed to create features
        #
        lineArray = gp.CreateObject("Array")
        pt = gp.CreateObject("Point")

        # Initialize a variable for keeping track of a feature's ID.
        #
        ID = -1
        while sRow:
            pt = sRow.GetValue(shapeName).GetPart(0)            
            currentValue = sRow.GetValue(IDField)

            if ID == -1:
                ID = currentValue

            if ID <> currentValue:
                if lineArray.count > 1:
                    feat = iCur.NewRow()
                    if ID: #in case the value is None/Null
                        feat.SetValue(IDField, ID)
                    feat.SetValue(shapeName, lineArray)
                    iCur.InsertRow(feat)
                else:
                    #gp.addwarning("Not enough points to create a line for %s: %s" % (IDField, str(ID)))
                    raise Exception("Not enough points to create a line for %s: %s" % (IDField, int(ID)))
                lineArray.RemoveAll()

            lineArray.Add(pt)
            ID = currentValue
            sRow = sCur.Next()

        # Add the last feature
        #
        if lineArray.count > 1:
            feat = iCur.NewRow()
            if ID: #in case the value is None/Null
                feat.SetValue(IDField, currentValue)
            feat.SetValue(shapeName, lineArray)
            iCur.InsertRow(feat)
        else:
            #gp.addwarning("Not enough points to create a line for %s: %s" % (IDField, str(ID)))
            raise Exception("Not enough points to create a line for %s: %s" % (IDField, int(ID)))
        lineArray.RemoveAll()

    except Exception, err:
        #gp.addwarning("Converting points to line not successful. %s" % err.message)
        raise
        
    finally:
        if iCur:
            del iCur
        if sRow:
            del sRow
        if sCur:
            del sCur
        if feat:
            del feat

def get_featureclass_name(workspace, filename):
    """
    Returns a feature class name appropriate to a specified workspace.
    """
    if gp.Describe(workspace).DataType in ["Folder", "Workspace"]:
        if isinstance(filename, str):
            pass
        else:
            raise Exception("File name must be in string format.")
    else:
        raise Exception("Workspace must be a full path to a folder or workspace.")
    try:
        if gp.Describe(workspace).WorkspaceType == "FileSystem":
            outputName = "%s.shp" % filename
        else:
            outputName = filename
        featureClassName = os.path.join(workspace, outputName)
    except Exception, err:
        gp.AddError("Error in get_featureclass_name - %s" % err.message)
    else:
        return featureClassName

def get_minvalue(in_raster):
    return gp.getrasterproperties_management(in_raster, "MINIMUM")

def get_maxvalue(in_raster):
    return gp.getrasterproperties_management(in_raster, "MAXIMUM")

def list_watershedid(watershedbound, watershedfid, interval_profile):
    """
    Inputs:
    @watershedfid = watershed ID fieldname
    """
    widlist = []
    rows = gp.searchcursor(watershedbound)
    row = rows.next()
    while row:
        wid_value = row.getvalue(watershedfid)
        # Check watershed ID if integer.
        try:
            if isinstance(wid_value, int):
                pass
            elif isinstance(wid_value, float):
                if wid_value == int(wid_value):
                    wid_value = int(wid_value)
                else:
                    raise IOError("watershed ID must be an integer, not a float")
            elif isinstance(wid_value, str):
                try:
                    ewid_value = eval(wid_value)
                except:
                    raise IOError("watershed ID must be an integer, not a string")
                else:
                    if isinstance(ewid_value, float):
                        if ewid_value == int(ewid_value):
                            wid_value = int(ewid_value)
                        else:
                            raise IOError("watershed ID must be an integer, not a float")
                    elif isinstance(ewid_value, int):
                        wid_value = ewid_value
                    else:
                        raise IOError("watershed ID must be an integer, not unknown type")
            else:
                raise
        except Exception, err:
            gp.adderror(err.message)
            sys.exit()
        # Add to list if no duplicates.
        if wid_value not in widlist:
            if interval_profile > row.shape.length:
                gp.addwarning("skip watershed ID %i - too small for river profile interval" % wid_value)
            else:
                widlist.append(wid_value)
        else:
            raise IOError("some watershed IDs have duplicates")
        # Next row.
        row = rows.next()
    del rows, row
    # Return list.
    try:
        if widlist:
            return widlist
        else:
            raise IOError("no valid watersheds to process")
    except Exception, err:
        gp.adderror(err.message)
        sys.exit

def get_watershedid(watershedpolygon, watershedfid):
    """
    Inputs:
    @watershedfid = watershed ID fieldname
    """
    rows = gp.searchcursor(watershedpolygon)
    row = rows.next()
    while row:
        wid_value = row.getvalue(watershedfid)
        if isinstance(wid_value, int):
            pass
        elif isinstance(wid_value, float):
            wid_value = int(wid_value)
        elif isinstance(wid_value, str):
            wid_value = int(eval(wid_value))
        else:
            raise IOError("watershed ID must be an integer")
        break
    del rows, row
    try:
        return wid_value
    except:
        raise IOError("watershed is empty")

def remove_duplicate_measure(fc):
    """
    Ensures only one selection for each interval measure.
    """
    n = 0
    rows = gp.updatecursor(fc, "", "", "", "MEASURE;Shape_Length")
    row = rows.next()
    while row:
        m = row.getvalue("MEASURE")
        if n == m:
            row.setvalue("DIST_AB", 0)
            rows.updaterow(row)
        n = m
        row = rows.next()
    del rows, row

def update_fieldvalue(in_fc, fieldname, fieldvalue):
    #fieldname = gp.validatefieldname(fieldname, gp.workspace)
    rows = gp.updatecursor(in_fc)
    row = rows.next()
    while row:
        row.setvalue(fieldname, fieldvalue)
        rows.updaterow(row)
        row = rows.next()
    del rows, row

def get_valueraster(in_raster, type):
    """
    Inputs:
    @type = ["min"|"max"] value to extract
    """
    # Extract min or max value, and select expression.
    out_raster = gp.createuniquename("vgrid", gp.workspace)
    if type == "min":
        cvalue = get_minvalue(in_raster)
        whereclause = '"VALUE" > %f' % cvalue
        #gp.con_sa(in_raster, in_raster, out_raster, "#", ('"VALUE" <= %f' % cvalue))  # problematic if 32-bit raster
    elif type == "max":
        cvalue = get_maxvalue(in_raster)
        whereclause = '"VALUE" < %f' % cvalue
    else:
        raise IOError("type must be either min or max")
    # Generate raster based on min or max value.
    try:
        gp.setnull_sa(in_raster, in_raster, out_raster, whereclause)
    except Exception, err:
        raise IOError(err.message)
    return out_raster

def get_lowestpoint(in_raster):
    # Determine the lowest elevation.
    out_pointfc = get_featureclass_name(gp.workspace, "mylowpoint")
    out_raster = get_valueraster(in_raster, "min")
    gp.rastertopoint_conversion(out_raster, out_pointfc, "VALUE")
    # Select only one lowest point.
    gp.makefeaturelayer_management(out_pointfc, "mylowpoint_lyr")
    gp.selectlayerbyattribute_management("mylowpoint_lyr", "NEW_SELECTION", '"pointid" > 1')
    gp.deleterows_management("mylowpoint_lyr")
    # Delete intermediate data.
    gp.delete_management(gp.describe(out_raster).catalogpath)
    return out_pointfc

def get_farthestpoint(in_raster):
    # Set lowest elevation as reference point.
    out_pointfc = get_featureclass_name(gp.workspace, "myfarpoint")
    cellsize = gp.getrasterproperties_management(in_raster, "CELLSIZEX")
    minvalueraster = get_valueraster(in_raster, "min")
    # Create distance raster from reference point.
    eucdistraster = get_featureclass_name(gp.workspace, "myeucdistraster")
    gp.eucdistance_sa(minvalueraster, eucdistraster, "#", cellsize, "#")
    # Get farthest cell then convert to point.
    farthestraster = get_valueraster(eucdistraster, "max")
    gp.rastertopoint_conversion(farthestraster, out_pointfc, "VALUE")
    # Delete intermediate data.
    for pointfc in [minvalueraster, eucdistraster, farthestraster]:
        try:
            gp.delete_management(gp.describe(pointfc).catalogpath)
        except:
            gp.delete_management(pointfc)
    return out_pointfc

def merge(in_fcs, out_merge="out_merge"):
    """
    Returns the merge output of all input feature classes.
    
    Inputs:
    @in_fcs = list of input feature classes
    @out_merge = output feature class
    """
    out_merge = get_featureclass_name(gp.workspace, out_merge)
    fieldmappings = gp.createobject("FieldMappings")
    vt = gp.createobject("ValueTable")  # table to hold inputs to merge
    for input in in_fcs:
        fieldmappings.addtable(input)
        vt.addrow(input)
    gp.merge(vt, out_merge, fieldmappings)  # run the merge tool
    return out_merge

def get_coordpriority(in_raster, lowpoint):
    rows = gp.searchcursor(lowpoint)
    row = rows.next()
    while row:
        xcoord = row.shape.getpart().x
        ycoord = row.shape.getpart().y
        break
    del rows, row
    
    top = gp.getrasterproperties(in_raster, "TOP")
    bottom = gp.getrasterproperties(in_raster, "BOTTOM")
    right = gp.getrasterproperties(in_raster, "RIGHT")
    left = gp.getrasterproperties(in_raster, "LEFT")
    if abs(top - ycoord) < abs(bottom - ycoord):
        key1 = "UPPER"
    else:
        key1 = "LOWER"
    if abs(left - xcoord) < abs(right - xcoord):
        key2 = "LEFT"
    else:
        key2 = "RIGHT"
    return "%s_%s" % (key1, key2)
    
def get_longaxis(in_raster, wid, lowpoint, farpoint):
    line_field="WID"
    sort_field="ENDPOINT"
    
    x = 1
    for pointfc in [lowpoint, farpoint]:
        gp.addfield_management(pointfc, line_field, "LONG", "", "", "", "WATERSHED ID")
        gp.addfield_management(pointfc, sort_field, "SHORT")
        gp.calculatefield_management(pointfc, line_field, wid, "PYTHON", "#")
        gp.calculatefield_management(pointfc, sort_field, x, "PYTHON", "#")
        x += 1
    
    endpoints = merge([lowpoint, farpoint], "myendpoints")
    out_fc = get_featureclass_name(gp.workspace, "mylongaxis")
    try:
        point2line(endpoints, out_fc, line_field, sort_field)  # ArcGIS v9.3.1 or earlier
    except:
        gp.pointstoline_management(endpoints, out_fc, line_field, sort_field)  # ArcGIS 10.0 or later
    gp.calculatefield_management(out_fc, line_field, wid, "PYTHON", "#")
    gp.delete_management(gp.describe(endpoints).catalogpath)
    return out_fc

def get_axisroute(in_raster, wid, lowpoint, longaxis):
    line_field = "WID"
    axisroute = get_featureclass_name(gp.workspace, "myaxisroute")
    coordpriority = get_coordpriority(in_raster, lowpoint)
    gp.createroutes_lr(longaxis, line_field, axisroute, "", "", "", coordpriority)
    return axisroute

def create_sidefield(in_raster, lowpoint, axismidlines, sidefield="SIDE", sidealias="SIDE OF MIDLINE"):
    # Determine side values based on bearings.
    coordpriority = get_coordpriority(in_raster, lowpoint)
    sidekey = coordpriority.split("_")[1].upper()
    if sidekey == "RIGHT":
        fieldname = "LBEARING"
    elif sidekey == "LEFT":
        fieldname = "RBEARING"
    else:
        raise IOError("sidekey not in the expected list")
    # Create field name and alias then update field values.
    whereclause = "\"%s\" IS NOT NULL" % fieldname
    fclayer = "axismidlineslayer"
    gp.makefeaturelayer_management(axismidlines, fclayer)
    gp.addfield_management(fclayer, sidefield, "TEXT", "", "", "", sidealias)
    gp.selectlayerbyattribute_management(fclayer, "NEW_SELECTION", whereclause)
    gp.calculatefield_management(fclayer, sidefield, "\"RIGHT\"", "PYTHON", "#")
    gp.selectlayerbyattribute_management(fclayer, "SWITCH_SELECTION")
    gp.calculatefield_management(fclayer, sidefield, "\"LEFT\"", "PYTHON", "#")
    gp.selectlayerbyattribute_management(fclayer, "CLEAR_SELECTION")

def get_axislength(axisroute):
    rows = gp.searchcursor(axisroute)
    row = rows.next()
    length = 0
    while row:
        length = length + row.shape.length
        row = rows.next()
    del rows, row
    try:
        return length
    except Exception, err:
        raise IOError(err.message)

def get_minbound(in_raster):
    top = gp.getrasterproperties(in_raster, "TOP")
    bottom = gp.getrasterproperties(in_raster, "BOTTOM")
    right = gp.getrasterproperties(in_raster, "RIGHT")
    left = gp.getrasterproperties(in_raster, "LEFT")
    return int(min([abs(top - bottom), abs(right - left)]))

def get_maxbound(in_raster):
    top = gp.getrasterproperties(in_raster, "TOP")
    bottom = gp.getrasterproperties(in_raster, "BOTTOM")
    right = gp.getrasterproperties(in_raster, "RIGHT")
    left = gp.getrasterproperties(in_raster, "LEFT")
    return int(max([abs(top - bottom), abs(right - left)]))

def create_intervaltable(out_name, length, interval, wid, maxbound, count_interval=None):
    gp.createtable_management(gp.workspace, out_name)
    gp.addfield_management(out_name, "WID", "LONG", "", "", "", "WATERSHED ID")
    gp.addfield_management(out_name, "MEASURE", "DOUBLE", "", "", "", "INTERVAL DISTANCE")
    gp.addfield_management(out_name, "DISTANCE", "DOUBLE")
    gp.addfield_management(out_name, "RBEARING", "DOUBLE")
    gp.addfield_management(out_name, "LBEARING", "DOUBLE")
    
    tcount = int(math.ceil(length/interval))
    if count_interval and (count_interval <= tcount - 1):
        intervals = count_interval + 1
    elif count_interval and (count_interval > tcount - 1):
        raise IOError("requested number of midline intervals exceeds the length of longaxis")
    elif tcount < 2:
        raise IOError("midline interval exceeds the length of longaxis")  # watershed too narrow for user-specified interval
    else:
        intervals = tcount
    
    rows = gp.insertcursor(out_name)
    for x in xrange(1, intervals):
        row = rows.newrow()
        row.setvalue("WID", wid)
        row.setvalue("MEASURE", interval*x)
        row.setvalue("DISTANCE", maxbound)
        rows.insertrow(row)
    del rows, row
    return out_name

def create_csectionpoints(in_raster, lowpoint, riverline, wid, interval, count_interval, outputname):
    intervalpoints = create_intervalpoints(in_raster, lowpoint, riverline, wid, interval, count_interval)
    csectionpoints = get_csectionbearings(intervalpoints, riverline, interval, outputname)
    return csectionpoints

def create_axismidpoints(in_raster, lowpoint, longaxis, wid, interval, count_interval):
    intervalpoints = create_intervalpoints(in_raster, lowpoint, longaxis, wid, interval, count_interval)
    append_axisbearings(intervalpoints, lowpoint, interval)
    return intervalpoints

def create_intervalpoints(in_raster, lowpoint, longaxis, wid, interval, count_interval):
    # Convert longaxis line to route.
    axisroute = get_axisroute(in_raster, wid, lowpoint, longaxis)
    # Create axismidpoints via route events.
    length = get_axislength(axisroute)
    maxbound = get_maxbound(in_raster)
    intervaltable = create_intervaltable("mytable", length, interval, wid, maxbound, count_interval)
    props = "WID POINT MEASURE"
    out_layer = "mymidpointslayer"
    gp.makerouteeventlayer_lr(axisroute, "WID", intervaltable, props, out_layer)
    out_fc = get_featureclass_name(gp.workspace, "myaxismidpoints")
    gp.copyfeatures_management(out_layer, out_fc)
    gp.addxy_management(out_fc)
    # Delete intermediate data.
    gp.delete_management(gp.describe(axisroute).catalogpath)
    gp.delete_management(gp.describe(intervaltable).catalogpath)
    gp.delete_management(out_layer)
    return out_fc

def calculate_bearings(intervalpoints, anglefield="NEAR_ANGLE"):
    rows = gp.updatecursor(intervalpoints)
    row = rows.next()
    while row:
        angle = float(row.getvalue(anglefield))
        forward = (90 - angle) + 90
        backward = (90 - angle) - 90
        row.setvalue("RBEARING", backward)  # right bearing
        row.setvalue("LBEARING", forward)  # left bearing
        rows.updaterow(row)
        row = rows.next()
    del rows, row

##def get_azimuth(angle):
##    if angle <= 0: angle += 360
##    return angle
##
##def convert_nearangle(angle):
##    if angle <= 0:
##        angle += 360
##    elif angle >= 360:
##        angle -= 360
##    return (-1*angle) + 270

##def calculate_bearings(intervalpoints, anglefield="NEAR_ANGLE"):
##    rows = gp.updatecursor(intervalpoints)
##    row = rows.next()
##    while row:
##        angle = float(row.getvalue(anglefield))
##        bearing1 = (90 - angle) + 90
##        bearing2 = (90 - angle) - 90
##        absbear = abs(bearing1)
##        if (absbear > 90) and (absbear <= 270):
##            row.setvalue("RBEARING", bearing1)  # right bearing
##            row.setvalue("LBEARING", bearing2)  # left bearing
##        else:
##            row.setvalue("RBEARING", bearing2)
##            row.setvalue("LBEARING", bearing1)
##        rows.updaterow(row)
##        row = rows.next()
##    del rows, row 

def append_axisbearings(intervalpoints, lowpoint, interval):
    # Get side angle orientation.
    maxmeasure = (gp.getcount_management(intervalpoints) + 1)*interval   
    radius = "%i" % maxmeasure
    gp.near_analysis(intervalpoints, lowpoint, radius, "LOCATION", "ANGLE")  # north = 90 degrees, east = 0 degree

    gp.addmessage("near analysis for intervalpoints calculated ...")
    
    # Calculate bearings.
    calculate_bearings(intervalpoints)

def get_neardistance(baseline, targetline):
    radius = get_axislength(targetline)
    gp.near_analysis(baseline, targetline, radius, "LOCATION", "ANGLE")
    
    rows = gp.searchcursor(baseline, "", "", "", "NEAR_DIST D")
    row = rows.next()
    neardist = None
    while row:
        neardist = math.ceil(row.getvalue("NEAR_DIST"))
        break
    del rows, row
    if neardist is None: raise IOError("baseline is an empty feature, failed to get neardistance")
    return neardist

def get_nearpoint(riverline, targetpoint, outputpoint):
    radius = get_axislength(riverline)
    fvertices = gp.createuniquename("riverline_vertices")
    gp.featureverticestopoints_management(riverline, fvertices, "DANGLE")
    gp.near_analysis(fvertices, targetpoint, radius, "LOCATION", "ANGLE")
    
##    rows = gp.searchcursor(fvertices)
##    row = rows.next()
##    oidfieldname = gp.describe(fvertices).oidfieldname
##    dictdistance = {}
##    while row:
##        oid = row.getvalue(oidfieldname)
##        distance = row.getvalue("NEAR_DIST")
##        dictdistance[oid] = distance
##        row = rows.next()
##    del rows, row
##
##    if not dictdistance: raise IOError("riverline is an empty feature, failed to generate vertices")
##    mindist = min(dictdistance.values())
##    minoid = [x for x, y in dictdistance.items() if y == mindist][0]
    rows = gp.searchcursor(fvertices, "", "", "", "NEAR_DIST")
    row = rows.next()
    oidfieldname = gp.describe(fvertices).oidfieldname
    minoid = None
    while row:
        minoid = row.getvalue(oidfieldname)
        break
    del rows, row
    if minoid is None: raise IOError("riverline is an empty feature, failed to generate vertices")
    
    gp.makefeaturelayer_management(fvertices, "riverline_verticeslyr")
    whereclause = '"%s" = %i' % (oidfieldname, minoid)
    gp.selectlayerbyattribute_management("riverline_verticeslyr", "NEW_SELECTION", whereclause)
    gp.copyfeatures_management("riverline_verticeslyr", outputpoint)
    gp.delete_management(gp.describe(fvertices).catalogpath)
    return outputpoint

def get_streamheadpoint(in_raster, farpoint, wid):
    riverline, rivergrid = extract_mainstream(in_raster, wid, order="LOW_HIGH")
    streamheadpoint = get_nearpoint(riverline, farpoint, "streamheadpoint")
    return streamheadpoint

def get_geometricnetwork_trace(geonetworkpath, rivername, flags, tracetask, outputname):
    """
    Inputs:
    @tracetask = ["FIND_PATH", "FIND_CONNECTED"]
    """
    outputgeonetwork = "geometricrivertrace_net"
    gp.tracegeometricnetwork_management(geonetworkpath, outputgeonetwork, flags, tracetask)
    geonetworkpath = os.path.join(outputgeonetwork, rivername)
    gp.copyfeatures_management(geonetworkpath, outputname)
    gp.delete_management(outputgeonetwork)
    return outputname

def get_longest_riverline(riverline, sinkpoint, sourcepoint, interval_profile):
    
    # Create geometric network from riverline.
    featuredataset, riveredge = get_featuredataset_pair(riverline)
    geonetworkname = "geometricriver_net"
    rivername = os.path.basename(riveredge)
    infeatureclass = "%s SIMPLE_EDGE NO" % rivername
    gp.creategeometricnetwork_management(featuredataset, geonetworkname, infeatureclass)
    geonetworkpath = os.path.join(featuredataset, geonetworkname)

    gp.addmessage("geometric network for streams created ...")
    
    # Create flags using sinkpoint and sourcepoint.
    flag1 = get_nearpoint(riveredge, sinkpoint, "sinkpoint")  # makes flag1 connected to network
    riverline_connected = get_geometricnetwork_trace(geonetworkpath, rivername, flag1, "FIND_CONNECTED", "riverline_connected")
    
    if interval_profile > get_axislength(riverline_connected): #????
        raise IOError("riverline_connected too short for the river profile interval")
    
    flag2 = get_nearpoint(riverline_connected, sourcepoint, "sourcepoint")  # makes flag2 connected to network
    pointstomerge = '"%s;%s"' % (flag1, flag2)
    flags = "flags"
    gp.merge_management(pointstomerge, flags)

    gp.addmessage("flags consisting of source and sink points created ...")
    
    # Trace the path between flags and generate longest riverline.
    longestpath = get_geometricnetwork_trace(geonetworkpath, rivername, flags, "FIND_PATH", "longestpath")
    gp.dissolve_management(longestpath, "longest_riverline", "WID")

    # Check if longest river exceeds river profile interval.
    if interval_profile > get_axislength("longest_riverline"):
        raise IOError("main (longest) river too short for the river profile interval")

    gp.addmessage("longest_riverline created ...")
    
##    for fc in [flag1, flag2, flags, featuredataset, longestpath]:
##        gp.delete_management(gp.describe(fc).catalogpath)
    return "longest_riverline"

def get_featuredataset_pair(riverline):
    gp.createfeaturedataset_management(gp.workspace, "myfeaturedataset", gp.describe(riverline).spatialreference)
    featuredataset = os.path.join(gp.workspace, "myfeaturedataset")
    riveredge = os.path.join(featuredataset, "myriveredge")
    gp.copyfeatures_management(riverline, riveredge)
    return featuredataset, riveredge

def get_csectionbearings(intervalpoints, riverline, interval, outputname, sequencefield="SEQUENCE", anglefield="ANGLE"):
    # Get bearings from riverline segments.
    rsegments = get_segments_with_bearings(riverline, interval, sequencefield, anglefield)

    fieldstoretain = ["MEASURE", "ELEVATION", "WID", "SIDE", "DIST_AB", "NEAR_ANGLE", "ANGLE", "POINT_X", "POINT_Y", "DISTANCE", "RBEARING", "LBEARING"]  # NEAR_ANGLE???
    fieldmappings = get_customfieldmappings(intervalpoints, rsegments, fieldstoretain)
    gp.spatialjoin_analysis(intervalpoints, rsegments, outputname, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "#", "#")
    # Calculate bearings and add to output.
    calculate_bearings(outputname, anglefield)
    # Delete intermediate data.
    #gp.delete_management(gp.describe(rsegments).catalogpath)
    return outputname

def get_segments_with_bearings(inputline, interval, sequencefield="SEQUENCE", anglefield="ANGLE"):
    rsegments, rvertices = get_sequenced_features(inputline, sequencefield)
    add_bearings(rvertices, interval, sequencefield)
    gp.joinfield_management(rsegments, sequencefield, rvertices, sequencefield, anglefield)
    gp.delete_management(gp.describe(rvertices).catalogpath)
    return rsegments

def add_bearings(rvertices, interval, sequencefield="SEQUENCE", anglefield="ANGLE"):
    maxmeasure = (gp.getcount_management(rvertices) + 1)*interval   
    radius = "%i" % maxmeasure
    gp.makefeaturelayer_management(rvertices, "pointlyr")
    gp.makefeaturelayer_management(rvertices, "targetlyr")
    gp.addfield_management("pointlyr", anglefield, "DOUBLE", "", "", "", "ANGLE TO NEXT POINT")

    gp.setprogressorlabel("calculating bearings for each point ...")

    sids = gp.getcount_management("pointlyr")
    for x in xrange(1, sids):
        whereclause = "\"%s\" = %i" % (sequencefield, x)
        gp.selectlayerbyattribute_management("pointlyr", "NEW_SELECTION", whereclause)

        whereclause = "\"%s\" = %i" % (sequencefield, x + 1)
        gp.selectlayerbyattribute_management("targetlyr", "NEW_SELECTION", whereclause)
        
        gp.near_analysis("pointlyr", "targetlyr", radius, "LOCATION", "ANGLE")
        gp.calculatefield_management("pointlyr", anglefield, "!NEAR_ANGLE!", "PYTHON", "#")
    gp.delete_management("pointlyr")
    gp.delete_management("targetlyr")

def get_sequenced_features(inputline, sequencefield="SEQUENCE"):
    gp.dissolve_management(inputline, "riverdis")
    gp.featureverticestopoints_management("riverdis", "rvertices", "ALL")
    sinkpoint = get_nearpoint("riverdis", "rvertices", "sinkpoint")
    gp.splitline_management("riverdis", "rsegments")
    
    gp.makefeaturelayer_management("rvertices", "pointlyr")
    gp.makefeaturelayer_management("rsegments", "linelyr")
    gp.addfield_management("pointlyr", sequencefield, "SHORT")
    gp.addfield_management("linelyr", sequencefield, "SHORT")
    sids = gp.getcount_management("linelyr")
    gp.selectlayerbylocation_management("pointlyr", "INTERSECT", sinkpoint, "#", "NEW_SELECTION")
    gp.selectlayerbylocation_management("linelyr", "INTERSECT", sinkpoint, "#", "NEW_SELECTION")
    sid = 1
    gp.calculatefield_management("pointlyr", sequencefield, sid, "PYTHON", "#")
    gp.calculatefield_management("linelyr", sequencefield, sid, "PYTHON", "#")
    gp.selectlayerbylocation_management("pointlyr", "INTERSECT", "linelyr", "#", "NEW_SELECTION")
    whereclause = "\"%s\" IS NULL" % sequencefield
    gp.selectlayerbyattribute_management("pointlyr", "SUBSET_SELECTION", whereclause)
    gp.calculatefield_management("pointlyr", sequencefield, sid + 1, "PYTHON", "#")

    gp.setprogressorlabel("sorting features by sequence ...")
    
    for x in xrange(sid + 1, sids + 1):
        whereclause = "\"%s\" IS NOT NULL" % sequencefield
        gp.selectlayerbyattribute_management("linelyr", "NEW_SELECTION", whereclause)
        
        gp.selectlayerbylocation_management("linelyr", "INTERSECT", "linelyr", "#", "NEW_SELECTION")
        whereclause = "\"%s\" IS NULL" % sequencefield
        gp.selectlayerbyattribute_management("linelyr", "SUBSET_SELECTION", whereclause)
        gp.calculatefield_management("linelyr", sequencefield, x, "PYTHON", "#")

        gp.selectlayerbylocation_management("pointlyr", "INTERSECT", "linelyr", "#", "NEW_SELECTION")
        whereclause = "\"%s\" IS NULL" % sequencefield
        gp.selectlayerbyattribute_management("pointlyr", "SUBSET_SELECTION", whereclause)
        gp.calculatefield_management("pointlyr", sequencefield, x + 1, "PYTHON", "#")

    gp.delete_management("pointlyr")
    gp.delete_management("linelyr")
    a, b = "riverdis", "sinkpoint"
    for letter in [a, b]:
        gp.delete_management(gp.describe(letter).catalogpath)
    return "rsegments", "rvertices"

def get_bearingtoline(axismidpoints, bearing_field):
    """
    Inputs:
    @bearing_field = ["RBEARING"|"LBEARING"]
    """
    out_line = get_featureclass_name(gp.workspace, "my%s" % bearing_field.lower())
    gp.bearingdistancetoline_management(axismidpoints, out_line, "POINT_X", "POINT_Y", "DISTANCE",
                                        "METERS", bearing_field, "", "GEODESIC", "MEASURE")  # north = 0 degree, east = 90 degrees
    return out_line

def get_axismidlines(in_raster, axismidpoints, watershedpolygon, wid, longaxis, lowpoint, crosstype, outputname, sidefield="SIDE", sidealias="SIDE OF MIDLINE"):
    # Create lines to right and left of interval points.
    rightmidlines = get_bearingtoline(axismidpoints, "RBEARING")
    leftmidlines = get_bearingtoline(axismidpoints, "LBEARING")
    out_merge = merge([rightmidlines, leftmidlines], "mergemidlines")

    gp.addmessage("bearing lines created and merged ...")
    
    out_clip = get_featureclass_name(gp.workspace, outputname)
    gp.clip_analysis(out_merge, watershedpolygon, out_clip)
    # Assign watershed ID.
    gp.addfield_management(out_clip, "WID", "LONG", "", "", "", "WATERSHED ID")
    gp.calculatefield_management(out_clip, "WID", wid, "PYTHON", "#")

    gp.addmessage("watershed ID added ...")
    
    # Remove dangling or broken lines.
    axismidlines = remove_danglinglines(out_clip, longaxis)

    gp.addmessage("dangling lines removed ...")
    
    # Update sidefield if right or left.
    create_sidefield(in_raster, lowpoint, axismidlines, sidefield, sidealias)

    gp.addmessage("sidefield created and updated ...")

    # Assign crossline type.
    gp.addfield_management(axismidlines, "PERTO", "TEXT", "", "", "", "PERPENDICULAR TO")
    gp.addmessage("crosstype = %s" % crosstype)
    crosstypevalue = "\"%s\"" % crosstype
    gp.calculatefield_management(axismidlines, "PERTO", crosstypevalue, "PYTHON", "#")

    gp.addmessage("crosstype created and updated ...")
    
##    # Delete intermediate data.
##    for infc in [rightmidlines, leftmidlines, out_merge]:
##        gp.delete_management(gp.describe(infc).catalogpath)
        
    return axismidlines

def get_flowdirection(in_raster):
    gp.flowdirection_sa(in_raster, "flowdir", "NORMAL", "#")
    return "flowdir"

def get_flowaccumulation(flowdir):
    gp.flowaccumulation_sa(flowdir, "flowacc", "#", "FLOAT")
    return "flowacc"

def get_streamorder(flowacc, flowdir):
    gp.streamorder_sa(flowacc, flowdir, "streamorder", "STRAHLER")
    return "streamorder"

def get_hiflowaccpoint(in_raster):
    try:
        # Create flow accumulation based on DEM.
        flowdir = get_flowdirection(in_raster)
        flowacc = get_flowaccumulation(flowdir)
        # Get point with highest accumulation.
        hiflowacc = get_valueraster(flowacc, "max")
        gp.rastertopoint_conversion(hiflowacc, "hiflowaccpoint", "VALUE")
        # Delete intermediate data.
        for inras in [flowdir, flowacc, hiflowacc]:
            gp.delete_management(gp.describe(inras).catalogpath)
    except Exception, err:
        raise IOError(err.message)
    else:
        return "hiflowaccpoint"

def extract_mainstream(in_raster, wid, order="HIGH_ONLY"):
    """
    Returns main stream (in vector and raster formats) from DEM based on max stream order.

    Inputs:
    @order = ["HIGH_ONLY", "LOW_HIGH"]
    """
    # Define raster names.
    rivergrid = "rivergrid"  # raster format
    riverline = "riverline"  # vector format
    # Create stream order based on DEM.
    flowdir = get_flowdirection(in_raster)
    flowacc = get_flowaccumulation(flowdir)

    try:
        streamorder = get_streamorder(flowacc, flowdir)
    except:
        raise IOError("unknown error in reading raster")
    
    # Create stream (river) based on stream order.
    if order == "HIGH_ONLY":
        maxorder = get_maxvalue(streamorder)
        saquery = '"VALUE" < %i' % maxorder
    elif order == "LOW_HIGH":
        minorder = get_minvalue(streamorder)
        saquery = '"VALUE" < %i' % (minorder + 1)
    elif order == "ALL":
        minorder = get_minvalue(streamorder)
        saquery = '"VALUE" < %i' % minorder
    else:
        raise IOError('stream order must be specified either "HIGH_ONLY", "LOW_HIGH", or "ALL"')
    gp.setnull_sa(streamorder, streamorder, rivergrid, saquery)
    # Check if non-empty output.
    if gp.getcount_management(rivergrid) == 0: raise IOError("rivergrid is empty, uncheck option for edge effects then try again")
    # Convert stream raster to line.
    gp.streamtofeature_sa(rivergrid, flowdir, riverline, "SIMPLIFY")
    # Append watershed ID field.
    gp.addfield_management(riverline, "WID", "LONG", "", "", "", "WATERSHED ID")
    gp.calculatefield_management(riverline, "WID", wid, "PYTHON", "#")
    # Delete intermediate data.
    for inras in [flowdir, flowacc, streamorder]:
        gp.delete_management(gp.describe(inras).catalogpath)
    return riverline, rivergrid
    
def extract_valleyfloor(in_raster, riverline, wid, pthreshold=12.5):
    """
    Returns the valley floor areas (vector format) extracted from DEM.
    Inputs:
    @in_raster = DEM
    @rivergrid = river in raster format
    @pthreshold = default value is 12.5 percent  # 1 out of 8 cells (rule-of-thumb for valley floor)
    """
    # Identify flat areas based on percent slope value.
    percentslope = "percentslope"
    flatslope = "flatslope"
    gp.slope_sa(in_raster, percentslope, "PERCENT_RISE", "1")
    saquery = '"VALUE" > %f' % pthreshold
    gp.setnull_sa(percentslope, percentslope, flatslope, saquery)
    # Combine stream and flat areas.
    vafloor1 = "vafloor1"
    vafloor2 = "vafloor2"
    vafloor12 = "vafloor12"
    riverlinegrid = "riverlinegrid"
    valleygrid = "valleygrid"
    valleyfloor = "valleyfloor"
    gp.isnull_sa(flatslope, vafloor1)
    
    gp.extent = gp.describe(in_raster).extent
    gp.polylinetoraster_conversion(riverline, "OBJECTID", riverlinegrid, "MAXIMUM_LENGTH", "NONE", in_raster)
    if gp.describe(riverlinegrid).spatialreference.name == "Unknown":
        gp.defineprojection_management(riverlinegrid, gp.describe(in_raster).spatialreference)
    gp.extent = None
    
    gp.isnull_sa(riverlinegrid, vafloor2)
    gp.times_sa(vafloor1, vafloor2, vafloor12)
    gp.setnull_sa(vafloor12, "1", valleygrid, '"VALUE" = 1')
    # Convert areas to polygon.
    gp.rastertopolygon_conversion(valleygrid, valleyfloor, "SIMPLIFY", "Value")
    # Delete intermediate data.
    for inras in [percentslope, flatslope, vafloor1, vafloor2, vafloor12, riverlinegrid, valleygrid]:
        gp.delete_management(gp.describe(inras).catalogpath)
    # Adjust valleyfloor minimum width to one cell size.
    valleyfloor = update_valleyfloor(in_raster, valleyfloor, riverline)

    gp.addfield_management(valleyfloor, "WID", "LONG", "", "", "", "WATERSHED ID")
    gp.calculatefield_management(valleyfloor, "WID", wid, "PYTHON", "#")
    
    return valleyfloor

def update_valleyfloor(in_raster, valleyfloor, riverline):
    gp.makefeaturelayer_management(valleyfloor, "valleyfloor_lyr")
    gp.selectlayerbylocation_management("valleyfloor_lyr", "INTERSECT", riverline, "#", "NEW_SELECTION")
    gp.copyfeatures_management("valleyfloor_lyr", "valleyfloor_red")
    cellsize = float(gp.getrasterproperties_management(in_raster, "CELLSIZEX"))
    buffdist = "%f" % (cellsize/2)
    gp.buffer_analysis(riverline, "riverline_buf", buffdist, "FULL", "ROUND", "NONE", "#")
    gp.update_analysis("valleyfloor_red", "riverline_buf", "valleyfloor_upd", "NO_BORDERS", "#")
    gp.dissolve_management("valleyfloor_upd", "valleyfloor_dis")
    gp.copyfeatures_management("valleyfloor_dis", valleyfloor)
    a, b, c, d = "valleyfloor_red", "riverline_buf", "valleyfloor_upd", "valleyfloor_dis"
    for fc in [a, b, c, d]:
        gp.delete_management(gp.describe(fc).catalogpath)
    return valleyfloor

def dissolve_axismidlines(axismidlines, outputname):
    gp.dissolve_management(axismidlines, outputname, "MEASURE", "#", "SINGLE_PART", "DISSOLVE_LINES")
    return outputname

def get_fieldtype(in_table, fieldname):
    fieldtype = None
    f = gp.listfields(in_table)
    g = f.next()
    while g:
        if g.name.lower() == fieldname.lower():
            fieldtype = g.type
            break
        g = f.next()
    del f, g
    if fieldtype is None:
        raise IOError("fieldname '%s' cannot be found" % fieldname)
    return fieldtype

def get_watershedmidpoints(axisdissolve, longaxis, lowpoint, farpoint):
    # Convert dissolved axismidlines to points (axismidpoints).
##    gp.makefeaturelayer_management(axisdissolve, "axisdissolvelayer")
##    gp.selectlayerbylocation_management("axisdissolvelayer", "INTERSECT", longaxis, "#", "NEW_SELECTION")
##    gp.featuretopoint_management("axisdissolvelayer", "axismidpoints", "INSIDE")
##    gp.delete_management("axisdissolvelayer")

    gp.featuretopoint_management(axisdissolve, "axismidpoints", "INSIDE")
    linecount = gp.getcount_management("axismidpoints")  # used to assign ID number to farpoint
    # Merge axismidpoints with first (lowpoint) and last (farpoint) points.
    pointstomerge = '"%s;%s;%s"' % (lowpoint, "axismidpoints", farpoint)
    gp.merge_management(pointstomerge, "watershedmidpoints")
    # Sort all points in order.
    gp.addfield_management("watershedmidpoints", "SORT_ID", "SHORT")
    rows = gp.updatecursor("watershedmidpoints", "", "", "", "MEASURE")
    row = rows.next()
    n = 1
    while row:
        endpoint = row.getvalue("ENDPOINT")
        if endpoint == 1:
            rvalue = n
        elif endpoint == 2:
            rvalue = linecount + 2
        else:
            n += 1
            rvalue = n
        row.setvalue("SORT_ID", rvalue)
        rows.updaterow(row)
        row = rows.next()
    del rows, row
    # Delete intermediate data.
    axismidpoints = "axismidpoints"
    gp.delete_management(gp.describe(axismidpoints).catalogpath)
    return "watershedmidpoints"

def get_watershedmidline(watershedmidpoints, watershedpolygon, wid):
    # Convert points to line.
    try:
        point2line(watershedmidpoints, "watershedmidline", "", "SORT_ID")  # ArcGIS v9.3.1 or earlier
    except:
        gp.pointstoline_management(watershedmidpoints, "watershedmidline", "", "SORT_ID")  # ArcGIS 10.0 or later
    # Check if midline not outside watershed.
    gp.makefeaturelayer("watershedmidline", "watershedmidlinelyr")
    gp.selectlayerbylocation_management("watershedmidlinelyr", "COMPLETELY_WITHIN", watershedpolygon, "#", "NEW_SELECTION")
    if gp.getcount_management("watershedmidlinelyr") == 0:
        gp.addwarning("part of midline is outside of watershed, consider reducing the longaxis interval then try again ...")
    gp.selectlayerbyattribute("watershedmidlinelyr", "CLEAR_SELECTION")
    gp.delete_management("watershedmidlinelyr")
    # Add watershed ID field.
    gp.addfield_management("watershedmidline", "WID", "LONG", "", "", "", "WATERSHED ID")
    gp.calculatefield_management("watershedmidline", "WID", wid, "PYTHON", "#")
    return "watershedmidline"

def check_nodatavalue(elevationpoints, fieldname="RASTERVALU"):
    rows = gp.searchcursor(elevationpoints)
    row = rows.next()
    while row:
        rasterv = row.getvalue(fieldname)
        if rasterv <= -9999:
            gp.addwarning("some elevation points have nodata values")
            #raise IOError("some elevation points have nodata values")
        break
    row = rows.next()
    del rows, row

def retain_fieldnames(in_table, fieldstoretain):
    if isinstance(fieldstoretain, str):
        fieldstodelete = None
    else:
        fieldlist = [fs.lower() for fs in fieldstoretain]
        removelist = []
        f = gp.listfields(in_table)
        g = f.next()
        while g:
            if g.name.lower() not in fieldlist and g.required.upper() == "FALSE":
                removelist.append(g.name)
            g = f.next()
        del f, g
        if float(gp.getinstallinfo()["Version"]) >= 10.0:
            fieldstodelete = removelist
        else:
            fieldstodelete = '"%s"' % (";".join(removelist))
    if fieldstodelete: gp.deletefield_management(in_table, fieldstodelete)

def get_tlooklines(in_raster, axismidpoints, axisdissolve, axismidlines, valleyfloor, riverline, outputname):
    # Get endpoints from dissolved axismidlines.
    gp.featureverticestopoints_management(axisdissolve,"axisendpoints","BOTH_ENDS")
    
    # Exract elevation values for each endpoint.
    origmask = gp.mask  # get original mask
    gp.mask = in_raster  # set new mask (larger than origmask to avoid edge effects)
    gp.extractvaluestopoints_sa("axisendpoints", in_raster, "wendpoints")
    check_nodatavalue("wendpoints", fieldname="RASTERVALU")  # check for nodata values
    fieldtype = get_fieldtype("wendpoints", "RASTERVALU")
    gp.addfield_management("wendpoints", "ELEVATION", fieldtype, "", "", "", "ENDPOINT ELEVATION")
    gp.calculatefield_management("wendpoints", "ELEVATION", "!RASTERVALU!", "PYTHON", "#")
    gp.mask = origmask  # return to original mask
    
    # Split axisdissolve by axismidpoints.
##    if float(gp.getinstallinfo()["Version"]) >= 10.0:
##        ntsinputs = [axisdissolve, riverline]
##    else:
##        ntsinputs = '"%s;%s"' % (axisdissolve, riverline)
##    gp.intersect_analysis(ntsinputs, "rivercrosslines_mts_points", "ALL", "#", "POINT")
    gp.splitlineatpoint_management(axisdissolve, axismidpoints, "axisdissolve_upd", "#")
    
    # Transfer attribute info from axismidlines to wendpoints.

    fieldmappings = "#"
    
##    retainfield = ["MEASURE", "WID", "SIDE_RIVER", "ELEVATION"]  # fieldnames to appear in wendpoints
##    fieldmappings = get_customfieldmappings("wendpoints", axismidlines, retainfield)
    gp.spatialjoin_analysis("wendpoints", axismidlines, "wendpoints_upd", "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "#", "#")
    # Transfer attribute info from wendpoints to tlooklines.
##    retainfield = ["MEASURE", "WID", "SIDE_RIVER", "ELEVATION", "Shape_Length"]  # fieldnames to appear in tlooklines or riverprofile
##    fieldmappings = get_customfieldmappings("axisdissolve_upd", "wendpoints_upd", retainfield)
    gp.spatialjoin_analysis("axisdissolve_upd", "wendpoints_upd", outputname, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "#", "#")
    gp.deletefield_management(outputname, "Join_Count;TARGET_FID")
    # Delete intermediate data.
    del fieldmappings
    a, b, c, d = "axisendpoints", "wendpoints", "axisdissolve_upd", "wendpoints_upd"
    for fc in [a, b, c, d]:
        gp.delete_management(gp.describe(fc).catalogpath)
    
    # Measure valleyfloor and update riverprofile.
    measure_valleyfloorwidth(axisdissolve, axismidpoints, valleyfloor, outputname)
    
    return outputname

def get_whereclause(fieldname, fieldvalue, operator="="):
    if isinstance(fieldvalue, str):
        whereclause = "\"%s\" %s '%s'" % (fieldname, operator, fieldvalue)
    else:
        whereclause = "\"%s\" %s %s" % (fieldname, operator, fieldvalue)
    return whereclause

def get_fieldvalue_expression(fieldvalue):
    if isinstance(fieldvalue, str):
        fieldparam = '"%s"' % fieldvalue
    else:
        fieldparam = fieldvalue
    return fieldparam
    
def update_field(fc, referencefield, targetfield, valuelist=None):
    """
    Inputs:
    @valuelist = [(referencevalue, targetvalue1), (referencevalue2, targetvalue2), ...]
    """
    if valuelist is None: raise IOError("valuelist must not be empty")
    fclayer = gp.createuniquename("pointlyr")
    gp.makefeaturelayer_management(fc, fclayer)

    for x, y in valuelist:
        whereclause = get_whereclause(referencefield, x, operator="=")
        gp.selectlayerbyattribute_management(fclayer, "NEW_SELECTION", whereclause)
        fieldvalue = get_fieldvalue_expression(y)
        gp.calculatefield_management(fclayer, targetfield, fieldvalue, "PYTHON", "#")
    gp.selectlayerbyattribute(fclayer, "CLEAR_SELECTION")
    gp.delete_management(fclayer)

def get_elevation_values(in_raster, pointfc, outfc):
    # Exract elevation values for each endpoint.
    origmask = gp.mask  # get original mask
    gp.mask = in_raster  # set new mask (larger than origmask to avoid edge effects)
    gp.extractvaluestopoints_sa(pointfc, in_raster, outfc)
    check_nodatavalue(outfc, fieldname="RASTERVALU")  # check for nodata values
    fieldtype = get_fieldtype(outfc, "RASTERVALU")
    gp.addfield_management(outfc, "ELEVATION", fieldtype, "", "", "", "ELEVATION")
    gp.calculatefield_management(outfc, "ELEVATION", "!RASTERVALU!", "PYTHON", "#")
    gp.mask = origmask  # return to original mask
    return outfc

def get_elevationpoints(in_raster, crosslines, crosslinedissolve, valleyfloor, csectionpoints, sidefield="SIDE_RIVER"):
    # Generate endpoints from crosslinedissolve -> elevendpoints.
    gp.featureverticestopoints_management(crosslinedissolve,"diss_endpoints","BOTH_ENDS")
    # Classify LOCATION = [WSR, WSL] for elevendpoints.
    gp.spatialjoin_analysis("diss_endpoints", crosslines, "elevendpoints", "JOIN_ONE_TO_ONE", "KEEP_COMMON", "#", "INTERSECT", "#", "#")
    gp.addfield_management("elevendpoints", "LOCATION", "TEXT", "", "", "", "POINT LOCATION")
    update_field("elevendpoints", sidefield, "LOCATION", valuelist=[("RIGHT","WSR"), ("LEFT","WSL")])
    # Delete intermediate data.
    gp.delete_management("diss_endpoints")

    gp.addmessage("elevendpoints created ...")
    
    # Clip crosslinedissolve by valleyfloor.
    gp.clip_analysis(crosslinedissolve, valleyfloor, "clipcross")
    #gp.intersect_analysis([crosslinedissolve, valleyfloor], "clipcross", "ALL", "#", "LINE")

    # Select by location clipcross that intersect with csectionpoints.
    gp.multiparttosinglepart_management("clipcross", "clipcross_sgp")
    gp.unsplitline_management("clipcross_sgp", "clipcross_tsm", "MEASURE", "#")
    
    gp.makefeaturelayer_management("clipcross_tsm", "clipcrosslyr")
    gp.selectlayerbylocation_management("clipcrosslyr", "INTERSECT", csectionpoints, "#", "NEW_SELECTION")
    
    # Generate endpoints from valleycross.
    gp.featureverticestopoints_management("clipcrosslyr","vfendpoints","BOTH_ENDS")
    # Classify LOCATION = [VFR, VFL] for elevendpoints.
    gp.spatialjoin_analysis("vfendpoints", crosslines, "elevflrpoints", "JOIN_ONE_TO_ONE", "KEEP_COMMON", "#", "INTERSECT", "#", "#")
    gp.addfield_management("elevflrpoints", "LOCATION", "TEXT", "", "", "", "POINT LOCATION")
    update_field("elevflrpoints", sidefield, "LOCATION", valuelist=[("RIGHT","VFR"), ("LEFT","VFL")])
    # Delete intermediate data.
    a, b, c = "clipcross", "clipcross_sgp", "vfendpoints"
    for letter in [a, b, c]:
        gp.delete_management(gp.describe(letter).catalogpath)
    
    gp.addmessage("elevflrpoints created ...")
    
    # Generate midpoints from valleycross.
    gp.featureverticestopoints_management("clipcrosslyr","vfmidpoints","MID")
    gp.spatialjoin_analysis("vfmidpoints", crosslines, "elevmidpoints", "JOIN_ONE_TO_ONE", "KEEP_COMMON", "#", "INTERSECT", "#", "#")
    # Classify LOCATION = VFM for elevmidpoints.
    gp.addfield_management("elevmidpoints", "LOCATION", "TEXT", "", "", "", "POINT LOCATION")
    fieldvalue = get_fieldvalue_expression("VFM")
    gp.calculatefield_management("elevmidpoints", "LOCATION", fieldvalue, "PYTHON", "#")
    # Delete intermediate data.
    a, b = "vfmidpoints", "clipcross_tsm"
    for letter in [a, b]:
        gp.delete_management(gp.describe(letter).catalogpath)
    
    gp.addmessage("elevmidpoints created ...")
    
    # Merge elevendpoints, elevflrpoints, elevmidpoints.
    merge_elevpoints = merge(["elevendpoints", "elevflrpoints", "elevmidpoints"], "merge_elevpoints")
    # Extract elevation from in_raster.
    elevationpoints = get_elevation_values(in_raster, merge_elevpoints, "elevationpoints")
    # Delete intermediate data.
    a, b, c = "elevendpoints", "elevflrpoints", "elevmidpoints"
    for letter in [a, b, c]:
        gp.delete_management(gp.describe(letter).catalogpath)

    gp.addmessage("extracted elevation values to points ...")

    fieldstoretain = ["MEASURE", "WID", "SIDE_RIVER", "PERTO", "LOCATION", "ELEVATION"]
    retain_fieldnames(elevationpoints, fieldstoretain)
    
    return elevationpoints

def get_customfieldmappings(infc1, infc2, fieldstoretain=None):
    # Combine fields from two inputs.
    fieldmappings = gp.createobject("fieldmappings")
    fieldmappings.addtable(infc1)
    fieldmappings.addtable(infc2)
    # Remove fields not in fieldstoretain.
    if fieldstoretain:
        f = fieldmappings.fields
        g = f.next()
        while g:
            if g.name not in fieldstoretain:
                k = fieldmappings.findfieldmapindex(g.name)
                fieldmappings.removefieldmap(k)
            g = f.next()
        del f, g
    return fieldmappings

def create_vfsegment_field(tlooklines):
    gp.addfield_management(tlooklines, "VFSEGMENT", "TEXT", "", "", "", "VALLEY FLOOR WIDTH SEGMENT")
    rows = gp.updatecursor(tlooklines)
    row = rows.next()
    while row:
        side = row.getvalue("SIDE_RIVER")
        distab = row.getvalue("WID_VF")
        #gp.addmessage("side: %s, distab: %i" % (side, distab))
        if side == "RIGHT" and distab == 1:
            row.setvalue("VFSEGMENT", "RVSEG_Y")  # VALLEY FLOOR WIDTH SEGMENT, RIGHT SIDE OF RIVER
        elif side == "LEFT" and distab == 1:
            row.setvalue("VFSEGMENT", "LVSEG_Y")  # VALLEY FLOOR WIDTH SEGMENT, LEFT SIDE OF RIVER
        elif side == "RIGHT" and distab == 0:
            row.setvalue("VFSEGMENT", "RVSEG_N")  # NOT A VALLEY FLOOR WIDTH SEGMENT, RIGHT SIDE OF RIVER
        elif side == "LEFT" and distab == 0:
            row.setvalue("VFSEGMENT", "LVSEG_N")  # NOT A VALLEY FLOOR WIDTH SEGMENT, LEFT SIDE OF RIVER
        else:
            row.setvalue("VFSEGMENT", "ERROR")
        rows.updaterow(row)
        row = rows.next()
    del rows, row

def create_distclass_field(tlooklines):
    gp.addfield_management(tlooklines, "DISTCLASS", "TEXT", "", "", "", "DISTANCE CLASS")
    rows = gp.updatecursor(tlooklines)
    row = rows.next()
    while row:
        side = row.getvalue("SIDE")
        distab = row.getvalue("DIST_AB")
        #gp.addmessage("side: %s, distab: %i" % (side, distab))
        if side == "RIGHT" and distab == 1:
            row.setvalue("DISTCLASS", "RDIST_AB")  # DISTANCE FROM MIDLINE TO RIVER, RIGHT SIDE OF MIDLINE
        elif side == "LEFT" and distab == 1:
            row.setvalue("DISTCLASS", "LDIST_AB")  # DISTANCE FROM MIDLINE TO RIVER, LEFT SIDE OF MIDLINE
        elif side == "RIGHT" and distab == 0:
            row.setvalue("DISTCLASS", "RDIST_BC")  # DISTANCE FROM RIVER TO EDGE OF WATERSHED, RIGHT SIDE OF MIDLINE
        elif side == "LEFT" and distab == 0:
            row.setvalue("DISTCLASS", "LDIST_BC")  # DISTANCE FROM RIVER TO EDGE OF WATERSHED, LEFT SIDE OF MIDLINE
        else:
            row.setvalue("DISTCLASS", "ERROR")
        rows.updaterow(row)
        row = rows.next()
    del rows, row

def create_distclass_field2(tlooklines):
    # Add field.
    gp.addfield_management(tlooklines, "DISTCLASS", "TEXT", "", "", "", "DISTANCE CLASS")
    # Update field based on the following codes:
    #   RDIST_AB = DISTANCE FROM MIDLINE TO RIVER, RIGHT SIDE OF MIDLINE
    #   LDIST_AB = DISTANCE FROM MIDLINE TO RIVER, LEFT SIDE OF MIDLINE
    #   RDIST_BC = DISTANCE FROM RIVER TO EDGE OF WATERSHED, RIGHT SIDE OF MIDLINE
    #   LDIST_BC = DISTANCE FROM RIVER TO EDGE OF WATERSHED, LEFT SIDE OF MIDLINE
    codeblock = '''def getclass(x, y):
        if x == "RIGHT" and y == 1:
            return "RDIST_AB"
        elif x == "LEFT" and y == 1:
            return "LDIST_AB"
        elif x == "RIGHT" and y == 0:
            return "RDIST_BC"
        elif x == "LEFT" and y == 0:
            return "LDIST_BC"
        else:
            return "ERROR"'''
    expression = "getclass( !SIDE!, !DIST_AB!)"
    gp.calculatefield_management(tlooklines, "DISTCLASS", expression, "PYTHON", codeblock)
    
def measure_midlinetoriver(tlooklines, riverline, watershedmidpoints):
    # Add field if segment connects from midline to riverline.
    gp.addfield_management(tlooklines, "DIST_AB", "SHORT", "", "", "", "MIDLINE TO RIVER")
    gp.calculatefield_management(tlooklines, "DIST_AB", "0", "PYTHON", "#")
    # Split riverprofile by riverline.
    if float(gp.getinstallinfo()["Version"]) >= 10.0:
        ntsinputs = [tlooklines, riverline]
    else:
        ntsinputs = '"%s;%s;%s"' % (tlooklines, riverline)
    gp.intersect_analysis(ntsinputs, "tlooklines_nts", "ALL", "#", "POINT")
    
    gp.splitlineatpoint_management(tlooklines, "tlooklines_nts", "riverprofile_vup", "#")
    
    # Select and update riverprofile that connects to midline.
    gp.makefeaturelayer_management("riverprofile_vup", "tlooklines_ftlk")
    gp.selectlayerbylocation_management("tlooklines_ftlk", "INTERSECT", riverline, "#", "NEW_SELECTION")
    gp.selectlayerbylocation_management("tlooklines_ftlk", "INTERSECT", watershedmidpoints, "#", "SUBSET_SELECTION")
    gp.calculatefield_management("tlooklines_ftlk", "DIST_AB", "1", "PYTHON", "#")
    remove_duplicate_measure("tlooklines_ftlk")  # remove duplicate measure with high midlinetoriver distance

    gp.addmessage("measured distance from midline to river ...")
    
    fieldstoretain = ["MEASURE", "WID", "SIDE", "PERTO", "DIST_AB"]
    #fieldstoretain = "#"
    retain_fieldnames("tlooklines_ftlk", fieldstoretain)
    
    # Add distance class field needed for summary table.
    gp.selectlayerbyattribute_management("tlooklines_ftlk", "CLEAR_SELECTION")
    create_distclass_field("tlooklines_ftlk")

    gp.addmessage("distclass field added and updated ...")
    
    # Delete intermediate data.
    gp.copyfeatures_management("riverprofile_vup", tlooklines)
    a, b = "tlooklines_nts", "riverprofile_vup"
    for fc in [a, b]:
        gp.delete_management(gp.describe(fc).catalogpath)

def get_pivot_table(inputtable, outtable, pivotfield, valuefield, inputfields="MEASURE;WID;PERTO"):
    gp.pivottable_management(inputtable, inputfields, pivotfield, valuefield, outtable)
    return outtable

def measure_valleyfloorwidth(axisdissolve, axismidpoints, valleyfloor, riverprofile):
    # Add field if segment is within valley floor.
    gp.addfield_management(riverprofile, "WID_VF", "SHORT", "", "", "", "VALLEY FLOOR")
    gp.calculatefield_management(riverprofile, "WID_VF", "0", "PYTHON", "#")
    # Split riverprofile by valley floor.
    if float(gp.getinstallinfo()["Version"]) >= 10.0:
        ntsinputs = [axisdissolve, valleyfloor]
    else:
        ntsinputs = '"%s;%s"' % (axisdissolve, valleyfloor)
    gp.intersect_analysis(ntsinputs, "tlooklines_nts", "ALL", "#", "LINE")
    gp.multiparttosinglepart_management("tlooklines_nts", "tlooklines_sgp")
    gp.unsplitline_management("tlooklines_sgp", "tlooklines_tsm", "MEASURE", "#")
    
    gp.makefeaturelayer_management("tlooklines_tsm", "tlooklinestsm_lyr")
    gp.selectlayerbylocation_management("tlooklinestsm_lyr", "INTERSECT", axismidpoints, "#", "NEW_SELECTION")
    
    gp.featureverticestopoints_management("tlooklinestsm_lyr", "valleyfloor_endpoints", "BOTH_ENDS")
    gp.splitlineatpoint_management(riverprofile, "valleyfloor_endpoints", "riverprofile_vup", "#")
    # Select and update riverprofile within valley floor.
    gp.makefeaturelayer_management("riverprofile_vup", "riverprofilevup_lyr")
    gp.selectlayerbylocation_management("riverprofilevup_lyr", "INTERSECT", axismidpoints, "#", "NEW_SELECTION")
    gp.calculatefield_management("riverprofilevup_lyr", "WID_VF", "1", "PYTHON", "#")
    gp.copyfeatures_management("riverprofile_vup", riverprofile)

    fieldstoretain = ["MEASURE", "WID", "SIDE_RIVER", "PERTO", "WID_VF"]
    retain_fieldnames(riverprofile, fieldstoretain)

    # Add valley floor segment classification for summary table.
    create_vfsegment_field(riverprofile)
    
    # Delete intermediate data.
    a, b, c, d, e = "tlooklines_nts", "tlooklines_sgp", "tlooklines_tsm", "valleyfloor_endpoints", "riverprofile_vup"
    for fc in [a, b, c, d, e]:
        gp.delete_management(gp.describe(fc).catalogpath)
    
def get_reducedinputs(in_raster, watershedpolygon, adjust_edge, reduceby=2):
    """
    Inputs:
    @adjust_edge = boolean [True, False]
    @reduceby = number of cells, to move back from watershed boundary
    """
    gp.mask = watershedpolygon
    if adjust_edge:
        # Create an inner buffer to compensate for edge effects.
        gp.featuretoline_management(watershedpolygon, "wshedoutline", "#", "ATTRIBUTES")
        cellsize = gp.getrasterproperties_management(in_raster, "CELLSIZEX")
        buffdist = "%f" % (cellsize*reduceby)
        gp.buffer_analysis("wshedoutline", "wshedbuff", buffdist, "FULL", "ROUND", "NONE", "#")
        # Extract adjusted versions of DEM and watershed area.
        gp.erase_analysis(watershedpolygon, "wshedbuff", "wshedpolygon_adj", "#")
        gp.extractbymask_sa(in_raster, "wshedpolygon_adj", "wsheddem_adj")
        # Delete intermediate data.
        x, y = "wshedoutline", "wshedbuff"
        gp.delete_management(gp.describe(x).catalogpath)
        gp.delete_management(gp.describe(y).catalogpath)
    else:
        gp.copyfeatures_management(watershedpolygon, "wshedpolygon_adj")
        gp.copyraster_management(in_raster, "wsheddem_adj")
    return "wsheddem_adj", "wshedpolygon_adj"

def remove_danglinglines(rivercrosslines, riverline):
    gp.multiparttosinglepart_management(rivercrosslines, "rivercrosslines_mts")
    gp.makefeaturelayer_management("rivercrosslines_mts", "rivercrosslines_mtsk")

    neardist = get_neardistance(rivercrosslines, riverline)
    gp.selectlayerbylocation_management("rivercrosslines_mtsk", "WITHIN_A_DISTANCE", riverline, neardist, "NEW_SELECTION")
    
    #gp.selectlayerbylocation_management("rivercrosslines_mtsk", "INTERSECT", riverline, "#", "NEW_SELECTION")
    
    gp.selectlayerbyattribute_management("rivercrosslines_mtsk", "SWITCH_SELECTION", "#")
    gp.deleterows_management("rivercrosslines_mtsk")
    
    gp.copyfeatures_management("rivercrosslines_mts", rivercrosslines)
    gp.delete_management("rivercrosslines_mts")
    return rivercrosslines

def get_gdbfiles_list(data_type, wildcard=None):
    fcslist = []
    if data_type == "featureclasses":
        fcs = gp.listfeatureclasses(wildcard)
    elif data_type == "tables":
        fcs = gp.listtables(wildcard)
    elif data_type == "datasets":
        fcs = gp.listdatasets(wildcard)
    else:
        raise IOError("cannot generate list for unknown data type")
    fc = fcs.next()
    while fc:
        fcslist.append(fc)
        fc = fcs.next()
    del fcs, fc
    return fcslist

def get_fcs_list(workspace, wildcard=None):
    if workspace:
        origws = gp.workspace  # get current workspace
        gp.workspace = workspace
    fcs_list = []
    features = get_gdbfiles_list("featureclasses", wildcard)
    datasets = get_gdbfiles_list("datasets", wildcard)
    tables = get_gdbfiles_list("tables", wildcard)
    fcs_list.extend(features)
    fcs_list.extend(datasets)
    fcs_list.extend(tables)
    if fcs_list:
        catalogpathlist = []
        for fc in fcs_list:
            try:
                fpath = gp.describe(fc).catalogpath
            except Exception, err:
                gp.addwarning("cannot add path to list - %s" % err.message)
            else:
                catalogpathlist.append(fpath)
        fcs_list = catalogpathlist
    if workspace:
        gp.workspace = origws  # return to current workspace
    return fcs_list

def delete_gisdata(list_data_to_retain=[], workspace=None):
    """
    Deletes only feature classes (including shapefiles), tables, and raster datasets.
    """
    if workspace:
        origws = gp.workspace  # get current workspace
        gp.workspace = workspace
    if not isinstance(list_data_to_retain, list):
        raise IOError("input must be a list of files that cannot be removed")
    elif list_data_to_retain:
        orig_fcs = []
        for retain in list_data_to_retain:
            try:
                retainpath = gp.describe(retain).catalogpath
            except:
                pass
            else:
                if retainpath:
                    orig_fcs.append(retainpath)
                else:
                    orig_fcs.append(retain)
                    raise ValueError("no catalogpath returned")
        list_data_to_retain = orig_fcs
    else:
        pass
        #gp.addwarning("all feature classes, tables, and raster datasets will be removed from %s" % gp.workspace)
    final_fcs = get_fcs_list(workspace, wildcard=None)
    for fc in final_fcs:
        try:
            if fc not in list_data_to_retain:
                gp.delete_management(fc)
        except:
            gp.addwarning("cannot delete %s" % fc)
    if workspace:
        gp.workspace = origws  # return to current workspace

def get_tempgdb(name):
    #origws = gp.workspace
    #gp.workspace = None
    #folder = os.path.dirname(gp.workspace)
    folder = gp.getsystemenvironment("TEMP")
    basename = "%s.gdb" % name
    tempgdb = gp.createuniquename(basename, folder)
    gdbname = os.path.basename(tempgdb)
    gp.createfilegdb_management(folder, gdbname)
    #gp.workspace = origws
    return tempgdb

####################
def get_riverprofile(in_raster, localdem, localwatershed, wid, interval, interval_profile, count_intervalpr, adjust_edge, count_interval):
    
    # Reduce watershed area to account for edge effects.
    reduced_raster, watershedpolygon = get_reducedinputs(localdem, localwatershed, adjust_edge, reduceby=2)
    
    # Set mask for raster analysis.
    gp.mask = watershedpolygon
    
    # Get lowpoint and farpoint.
    #lowpoint = get_lowestpoint(reduced_raster)
    lowpoint = get_hiflowaccpoint(reduced_raster)
    gp.setprogressorlabel("lowpoint created")
    farpoint = get_farthestpoint(reduced_raster)
    gp.setprogressorlabel("farpoint created")
    
    # Get main (longest) river.
    riverline, rivergrid = extract_mainstream(reduced_raster, wid, order="LOW_HIGH")
    gp.setprogressorlabel("mainstream extracted")
    longest_riverline = get_longest_riverline(riverline, lowpoint, farpoint, interval_profile)
    valleyfloor = extract_valleyfloor(reduced_raster, longest_riverline, wid, pthreshold=12.5)
    gp.setprogressorlabel("valleyfloor extracted")

    # Create watershedmidline.
    longaxis = get_longaxis(reduced_raster, wid, lowpoint, farpoint)
    gp.setprogressorlabel("longaxis created")
    axismidpoints = create_axismidpoints(reduced_raster, lowpoint, longaxis, wid, interval, count_interval)
    gp.setprogressorlabel("axismidpoints created")
    axismidlines = get_axismidlines(reduced_raster, axismidpoints, watershedpolygon, wid, longaxis, lowpoint, "LONGAXIS", "myaxismidlines")
    gp.setprogressorlabel("axismidlines created")
    
    axisdissolve = dissolve_axismidlines(axismidlines, "axisdissolve")
    gp.setprogressorlabel("axisdissolve created")
    watershedmidpoints = get_watershedmidpoints(axisdissolve, longaxis, lowpoint, farpoint)
    gp.setprogressorlabel("watershedmidpoints created")
    watershedmidline = get_watershedmidline(watershedmidpoints, watershedpolygon, wid)
    gp.setprogressorlabel("watershedmidline created")

    # Create midlineprofile.
    midcrosspoints = create_csectionpoints(reduced_raster, lowpoint, watershedmidline, wid, interval_profile, count_intervalpr, "midcrosspoints")
    midcrosslines = get_axismidlines(reduced_raster, midcrosspoints, watershedpolygon, wid, watershedmidline, lowpoint, "MIDLINE", "midcrosslines")
    measure_midlinetoriver(midcrosslines, longest_riverline, midcrosspoints)
    # Create summary table for midcrosslines.
    midcrosslinestable = get_pivot_table(midcrosslines, "midcrosslines_sumtable", "DISTCLASS", "Shape_Length")
    gp.addmessage("midlinetable created ...")

    
    # Create riverprofile.
    csectionpoints = create_csectionpoints(reduced_raster, lowpoint, longest_riverline, wid, interval_profile, count_intervalpr, "csectionpoints")
    gp.addmessage("csectionpoints created ...")
    gp.setprogressorlabel("csectionpoints created")
    
    rivercrosslines = get_axismidlines(reduced_raster, csectionpoints, watershedpolygon, wid, longest_riverline, lowpoint, "RIVER", "rivercrosslines", sidefield="SIDE_RIVER", sidealias="SIDE OF RIVER")
    gp.addmessage("rivercrosslines created ...")
    gp.setprogressorlabel("rivercrosslines created")
    crosslinedissolve = dissolve_axismidlines(rivercrosslines, "crosslinedissolve")
    gp.addmessage("crosslinedissolve created ...")
    gp.setprogressorlabel("crosslinedissolve created")

    # Extract elevation points.
    elevationpoints = get_elevationpoints(in_raster, rivercrosslines, crosslinedissolve, valleyfloor, csectionpoints, sidefield="SIDE_RIVER")
    gp.setprogressorlabel("elevationpoints created")
    # Create summary table for elevationpoints.
    elevtable = get_pivot_table(elevationpoints, "elevationpoints_sumtable", "LOCATION", "ELEVATION")
    gp.addmessage("elevationpoints_sumtable created ...")

    measure_valleyfloorwidth(crosslinedissolve, csectionpoints, valleyfloor, rivercrosslines)
    gp.addmessage("valley floor width measured ...")
    # Create summary table for valley floor width.
    vfwtable = get_pivot_table(rivercrosslines, "vfw_sumtable", "VFSEGMENT", "Shape_Length")
    gp.addmessage("vfw_sumtable created ...")
    
##    riverprofile = get_tlooklines(in_raster, csectionpoints, crosslinedissolve, rivercrosslines, valleyfloor, longest_riverline, "riverprofile")
##    gp.setprogressorlabel("riverprofile created")

    # Delete intermediate data.
##    fctodelete = [reduced_raster, lowpoint, farpoint, longaxis, axismidpoints, axismidlines, axisdissolve,
##                  watershedmidpoints, rivergrid, csectionpoints, rivercrosslines, crosslinedissolve]
##    for fc in fctodelete:
##        gp.delete_management(gp.describe(fc).catalogpath)
    
    # Return output feature class.
    return rivercrosslines, midcrosslines, watershedmidline, longest_riverline, valleyfloor, elevationpoints, midcrosslinestable, elevtable, vfwtable

####################   
def main(workspace, in_raster, watershedbound, watershedfid, widlist, interval, interval_profile, count_intervalpr, adjust_edge, count_interval):

    # Set workspaces.
    origws = gp.workspace  # get current workspace
    origscratchws = gp.scratchworkspace  # get current scratchworkspace
    tempgdb = get_tempgdb("rivertransect")
    gp.workspace = tempgdb
    gp.scratchworkspace = tempgdb
    
##    gp.addmessage("gp.workspace = %s" % gp.workspace)
##    gp.addmessage("gp.scratchworkspace = %s" % gp.scratchworkspace)
    
    # Generate river profile for all listed watersheds.
    gp.setprogressor("step", "generating river profiles...", 0, len(widlist), 1)

    #######    
    start = True
    gp.makefeaturelayer_management(watershedbound, "watershedbound_lyr")
    for wid in widlist:

##        if wid != 6:
##            continue

        gp.addmessage("processing watershed ID %i ..." % wid)
        # Extract local watershed.
        whereclause = '"%s" = %i' % (watershedfid, wid)
        gp.selectlayerbyattribute_management("watershedbound_lyr", "NEW_SELECTION", whereclause)
        gp.copyfeatures_management("watershedbound_lyr", "local_watershed")
        gp.mask = None
        gp.extractbymask_sa(in_raster, "local_watershed", "local_dem")

        # Generate river profile for extracted watershed.
        try:
            rivercrosslines, midcrosslines, watershedmidline, longest_riverline, valleyfloor, elevationpoints, midcrosslinestable, elevtable, vfwtable = get_riverprofile(in_raster, "local_dem", "local_watershed", wid, interval, interval_profile, count_intervalpr, adjust_edge, count_interval)
        except IOError, err:
            gp.addwarning("%s - skip watershed ID %i" % (err.message, wid))
            
            # Delete intermediate data.
            gp.setprogressorlabel("removing intermediate data ...")
            delete_gisdata(workspace=tempgdb)
            #gp.delete_management(gp.workspace)
            gp.setprogressorposition()
            continue
        except:
            raise

        # Save output data.
        if start:
            gp.copyfeatures_management(rivercrosslines, os.path.join(workspace, "rt_rivercrosslines"))
            gp.copyfeatures_management(midcrosslines, os.path.join(workspace, "rt_midcrosslines"))
            gp.copyfeatures_management(watershedmidline, os.path.join(workspace, "rt_watershedmidline"))
            gp.copyfeatures_management(longest_riverline, os.path.join(workspace, "rt_longest_riverline"))
            gp.copyfeatures_management(valleyfloor, os.path.join(workspace, "rt_valleyfloor"))
            gp.copyfeatures_management(elevationpoints, os.path.join(workspace, "rt_elevationpoints"))
            gp.copyrows_management(midcrosslinestable, os.path.join(workspace, "rt_midcrosslinestable"))
            gp.copyrows_management(elevtable, os.path.join(workspace, "rt_elevationpointstable"))
            gp.copyrows_management(vfwtable, os.path.join(workspace, "rt_valleyfloorwidthtable"))
            start = False
        else:
            gp.append_management(rivercrosslines, os.path.join(workspace, "rt_rivercrosslines"), "NO_TEST", "#", "#")
            gp.append_management(midcrosslines, os.path.join(workspace, "rt_midcrosslines"), "NO_TEST", "#", "#")
            gp.append_management(watershedmidline, os.path.join(workspace, "rt_watershedmidline"), "NO_TEST", "#", "#")
            gp.append_management(longest_riverline, os.path.join(workspace, "rt_longest_riverline"), "NO_TEST", "#", "#")
            gp.append_management(valleyfloor, os.path.join(workspace, "rt_valleyfloor"), "NO_TEST", "#", "#")
            gp.append_management(elevationpoints, os.path.join(workspace, "rt_elevationpoints"), "NO_TEST", "#", "#")
            gp.append_management(midcrosslinestable, os.path.join(workspace, "rt_midcrosslinestable"), "NO_TEST", "#", "#")
            gp.append_management(elevtable, os.path.join(workspace, "rt_elevationpointstable"), "NO_TEST", "#", "#")
            gp.append_management(vfwtable, os.path.join(workspace, "rt_valleyfloorwidthtable"), "NO_TEST", "#", "#")
        
        gp.setprogressorposition()

    # Restore original workspaces.
    gp.workspace = origws  # restore current workspace
    gp.scratchworkspace = origscratchws  # restore current scratchworkspace

    # Delete intermediate data.
    gp.setprogressorlabel("removing intermediate data ...")
    delete_gisdata(workspace=tempgdb)
    #gp.delete_management(gp.workspace)
    
    gp.addmessage("process successfully completed")
    
####################
if __name__ == '__main__':
    try:
        # Set input parameters.
        workspace = gp.getparameterastext(0)  # workspace (file geodatabase)
        in_raster = gp.getparameterastext(1)  # watershed DEM (raster layer)
        watershedbound = gp.getparameterastext(2)  # watershed boundaries (feature layer)
        watershedfid = gp.getparameterastext(3)  # watershed ID (fieldname)
        interval = gp.getparameterastext(4)  # midline intervals in meters (float)
        interval_profile = gp.getparameterastext(5)  # river profile intervals in meters (float)
        count_intervalpr = gp.getparameterastext(6)  # number of river profile intervals (integer)
        adjust_edge = gp.getparameterastext(7)  # adjustment for edge effects (boolean)
        
        # Convert inputs to proper formats.
        try:
            if not interval or float(interval) == 0:
                interval = 500  # default value in meters
                gp.addwarning("midline interval set to default value: %s meters" % interval)
            else:
                interval = float(interval)
        except ValueError:
            gp.adderror("midline interval must be a number")
            sys.exit()
        
        try:
            if not interval_profile or float(interval_profile) == 0:
                interval_profile = 500  # default value in meters
                gp.addwarning("river profile interval set to default value: %s meters" % interval_profile)
            else:
                interval_profile = float(interval_profile)
        except ValueError:
            gp.adderror("river profile interval must be a number")
            sys.exit()
        
        if not count_intervalpr or int(count_intervalpr) == 0:
            count_intervalpr = None
        else:
            count_intervalpr = int(count_intervalpr)
        
        if adjust_edge.upper() == "TRUE":
            adjust_edge = True
            gp.addwarning("reduced inputs by 2 cells to account for edge effects")
        else:
            adjust_edge = False
        
        count_interval = None  # number of midline intervals
        
        # List watershed IDs if valid.
        widlist = list_watershedid(watershedbound, watershedfid, interval_profile)
        
        # Generate riverprofile based on inputs.
        main(workspace, in_raster, watershedbound, watershedfid, widlist, interval, interval_profile, count_intervalpr, adjust_edge, count_interval)
    except:
        gp.addmessage(traceback.format_exc())
        del gp
        sys.exit()


