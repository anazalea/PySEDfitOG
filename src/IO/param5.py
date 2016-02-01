from copy import deepcopy
import re, os

#FIXME: add in other optional parameters



def SetParams(pfile, paramList, args=None):
    """Read in parameters from a file and overwrite default params. 
    Parameter file must have lines of the form 'keyword value'.
    Anything after Hashtags are ignored. Whitespace is ignored.
    Parameters can also be input directly from command line with command line
    over-riding text file over-riding default values.

    Inputs:
        pfile: String of file path to parameter text file
        paramList: List of Param objects. Says what keys to us and allows 
            for setting default values and error checking.

    Outputs:
        Dictionary containing parameter values.
    """
    paramList = deepcopy(paramList)
    allowedKeys = []
    usedKeys = []
    keyList = []
    valueList = []
    for x in paramList:
        allowedKeys.append(x.key)
 
    with open(pfile, "r") as f:
        for line in f:
            line = line.split('#')[0].strip() #remove comments & whitespace
            if len(line) > 0: #ignore blank lines
                stuff = line.split(None, 1) #split into keyword, value
                if len(stuff) != 2:
                    raise ValueError("Line in param file is missing either "
                                     "value or keyword: %s" % line)
                key, value = stuff
                keyList.append(key.lower())
                valueList.append(value)

    if args is not None and len(args) > 0:
        args = map(str.lower, args)
        clKeys = []
        clValues = []
        for x in args:
            if x[0] == '-' and x[1] not in '1234567890':
                clKeys.append(x)
        clValues = re.split("-[a-z_]+\s+", ' '.join(args))[1:]
        clValues = map(str.strip, clValues)
        if len(clKeys) != len(clValues):
            raise ValueError("Mismatched key/value pairs read in from command "
                             "line. Correct useage is -KEYWORD1 VALUE1 "
                             "-KEYWORD2 VALUE2 etc...")
        for i in xrange(len(clKeys)):
            value = clValues[i]
            key = clKeys[i].lower() 
            if key[0] == '-':
                key = key[1:]
            keyList.append(key)
            valueList.append(value)

    keyList, valueList = CompatabilityHandler(keyList, valueList)

    for i in xrange(len(keyList)):
        key = keyList[i]
        value = valueList[i]
        if key not in allowedKeys: 
            raise KeyError("%s is not a valid keyword" % key)
        for x in paramList:
            if x.key == key:
                x.setValue(value)
                break
        if key in usedKeys and not x.multipleKey:
            print("WARNING: Keyword %s used twice. Using value of %s" % 
                  (key, value))
        else:
            usedKeys.append(key)
  
    for x in paramList:

        #check that all mandatory keys are present
        if x.mandatory and x.key not in usedKeys:
            print usedKeys
            raise KeyError("Mandatory keyword %s is not present." % x.key)

        #check that associated keys are present
        if len(x.assocKeys) > 0:
            for key in x.assocKeys:
                if key not in usedKeys:
                    raise KeyError("When using %s, %s key and value must "
                                   "also be declared." % (x.key, key))

        #check that appropriate keys have the same length lists
        if len(x.sameLengthKeys) > 0 and hasattr(x, "value"):
            for key in x.sameLengthKeys:
                for y in paramList:
                    if key == y.key:
                       if len(x.value) != len(y.value): 
                           raise ValueError("%s and %s must have the same "
                                            "number of elements." % 
                                            (x.key, y.key))

        #set values of optional keywords to default if not present
        if not hasattr(x, "value"):
            x.value = x.defaultValue

    paramDict = {}
    for x in paramList:
        paramDict[x.key] = x.value

    return paramDict

#############################################################################

def CompatabilityHandler(keyList, valueList):
    """Make sed_fit.pl and make_bbsed.pl param files compatible with python
    versions.

    Inputs:
        keyList -- list of keywords all in lower case
        valueList -- list of values as strings

    Outputs:
        updated keyList and valueList
    """

    deprecatedKeys = {"model_bbsed_file" : "model_file",
                      "model_mags" : "model_flux_columns",
                      "restrict_param" : "restrict_model",
                      "bestfit_spectra_yn" : "output_bestfit_spectra",
                      "chisq_matrix_yn" : "output_chisq_matrix",
                      "fit_errorbars" : "errorbar_method",
                      "save_mc_results_yn" : "output_mc_results",
                      "data_mags" : "data_flux_columns",
                      "data_uncertainties" : "data_error_columns"}

    obsoleteKeys = ["model_dir", "verbose", "bestfit_spectra_mkfile", 
                    "chisq_matrix_dir", "data_dir", "data_upperlim_yn",
                    "data_upperlim_flag_cols", "data_upperlim_flag_vals",
                    "data_upperlim_flag_oper", "data_upperlim_lim_cols",
                    "data_upperlim_nsigma", "fitmethod_upperlim"]

    removeDex = []
    for i in xrange(len(keyList)):
        if keyList[i] in deprecatedKeys:
            print("WARNING: The use of %s is deprecated. Please use %s "
                  "instead." % (keyList[i], deprecatedKeys[keyList[i]]))
            keyList[i] = deprecatedKeys[keyList[i]]
        if keyList[i] in obsoleteKeys:
            print("WARNING: The keyword %s is obsolete and no longer has any "
                  "functionality.")
            removeDex.append(i)

    for i in xrange(len(keyList)):
        if keyList[i] == "model_dir":
            if keyList[i][-1] != "/":
                keyList[i] += "/"
            for j in xrange(len(keyList)):
                if keyList[j] == "model_file":
                    valueList[j] = valueList[i] + valueList[j]

        if keyList[i] == "data_dir":
            if keyList[i][-1] != "/":
                keyList[i] += "/"
            for j in xrange(len(keyList)):
                if keyList[j] == "data_file":
                    valueList[j] = valueList[i] + valueList[j]

        if keyList[i] == "model_param":
            valueList[i] = valueList[i].replace("fluxscale","true")

        if keyList[i] == "errorbar_method":
            valueList[i] = valueList[i].replace("fraction", "montecarlo")
            temp = valueList[i].split()
            valueList[i] = temp[0]
            if len(temp) > 1:
                keyList.append("errorbar_range")
                valueList.append(temp[1])
            if len(temp) > 2:
                keyList.append("montecarlo_iters")
                valueList.append(temp[2])

    for index in sorted(removeDex, reverse=True):
        del keyList[index]
        del valueList[index]

    for i in range(len(keyList)):
        print keyList[i], valueList[i]

    return keyList, valueList

#########################################################################

class Param(object):
    """Contains all information about a parameter that will be used to 
    control the flow of the main program. Paramaters are initiated with just a
    keyword (and other options for format and error handling), but no value.
    The value must be explicitly set using the setValue() method.

    Inputs:
        key -- The key word name. Must be all lower case containing only 
            letters and underscores.
        defaultValue -- Optional default value
        dataType -- Data type the value should be. Can be list of types if 
            value will be a list. Default is str
        isList -- Boolean flag stating if the value can be a list
        multipleKey -- Boolean flag allowing for a key to be input multiple 
            times. Each value will be saved in a list, rather than 
            overwriting the previous value.
        mandatory -- Boolean flag stating whether it is necessary for the user
            to input a value for this keyword. i.e. the program cannot run
            without it and there is no default value.
        maxSize -- If the value will be a list, this is the max size allowed.
            Can only be used if isList is True.
        allowedValues -- Optional list of possible values for this key. Any 
            other values will raise an exception.
        canBeNeg -- Boolean flag stating whether a value can be less than zero
        canHaveSpace -- Boolean flag stating whether a value can contain a 
            whitespace character.
        assocKeys -- List of associated keywords that must be present if this 
            key is used.
        sameLengthKeys -- List of keywords that must have the same number of 
            elements as this keyword. Can only be used if isList is True.
        flattenList -- For the rare case when a multipleKey is a list, this 
            will ensure that the value is one list, rather than a list of 
            lists. Can only be used if isList is True.
    """
    def __init__(self, key, defaultValue=None, dataType=str,  isList=False, 
                 multipleKey=False, mandatory=False, maxSize=None, 
                 allowedValues=None, canBeNeg=True, canHaveSpace=True,
                 assocKeys=[], sameLengthKeys=[], flattenList=False,outputfmt='%s '):
        self.key = key
        self.defaultValue = defaultValue
        self.dataType = dataType
        self.isList = isList
        self.multipleKey = multipleKey
        self.mandatory = mandatory
        self.maxSize = maxSize
        self.allowedValues = allowedValues
        self.canBeNeg = canBeNeg
        self.canHaveSpace = canHaveSpace
        self.assocKeys = assocKeys
        self.sameLengthKeys = sameLengthKeys
        self.flattenList = flattenList
        self.outputFmt = outputfmt

        if self.isList:
            if not isinstance(self.dataType, list):
                self.dataType = [self.dataType]
            if not isinstance(self.defaultValue, list):
                self.defaultValue = [self.defaultValue]
            if not isinstance(self.canBeNeg, list):
                self.canBeNeg = [self.canBeNeg]
            if not isinstance(self.allowedValues, list):
                self.allowedValues = [self.allowedValues]

        if not self.isList and self.maxSize is not None:
            raise ValueError("maxSize can only be set if isList is True")

        if not self.isList and len(self.sameLengthKeys) > 0:
            raise ValueError("sameLengthKeys can only be applied if isList "
                             "is set to True.")

        if self.flattenList and (not self.isList or not self.multipleKey):
            raise ValueError("flattenList can only be applied if isList "
                             "and multipleKey are set to True.")

    def setValue(self, value):
        """Set the value associated with a keyword. Raise an exception if 
        the value is incompatible with the options chosen during initiation.

        Inputs:
            value -- The value to be set. Must be a string. If the value should
                be a list, then the elements should be separated by whitespace.
        """
        if not self.canHaveSpace and len(value.split()) > 1:
            raise ValueError("Please remove space in value of %s" % self.key)
        if self.isList:
            values = value.split()
        else:
            values = [value]
            self.defaultValue = [self.defaultValue]
            self.dataType = [self.dataType]
            self.canBeNeg = [self.canBeNeg]
            self.allowedValues = [self.allowedValues]

        #check that max list size is not violated
        if self.maxSize is not None:
            if len(values) > self.maxSize:
                raise ValueError("Too many values in %s. Expected at most %i. "
                                 "Recieved %i: %s" % 
                                 (self.key, self.maxSize, len(values), value))

        passedValues = []
        for i in xrange(len(values)):
            
            #check if value is allowed
            if i < len(self.allowedValues):
                j = i 
            else:
                j = -1
            if self.allowedValues[j] is not None:
                if values[i].lower() not in self.allowedValues[j]:
                    raise ValueError("Value for %s must be one of %s. Instead "
                                     "value is %s." % 
                                     (self.key, self.allowedValues[j], 
                                      values[i]))
                else:
                    values[i] = values[i].lower()

            #attempt to convert string to appropriate datatype
            if i < len(self.defaultValue):
                j = i 
            else:
                j = -1
            if self.dataType[j].__name__ == 'bool':
                self.dataType[j] = boolean
            try:
                values[i] = self.dataType[j](values[i])
            except ValueError:
                raise ValueError("Something went wrong with key %s with value "
                                 "%s. Unable to convert %s to format %s." %
                                 (self.key, value, values[i], 
                                  self.dataType[j].__name__))

            #check if value can be negative
            if i < len(self.canBeNeg):
                j = i 
            else:
                j = -1
            if not self.canBeNeg[j] and values[i] < 0:
                raise ValueError("Value of %s cannot be negative." % self.key)

            #passed all tests
            passedValues.append(values[i])

        #if fewer values given than expected, fill in rest with defaults
        if len(passedValues) < len(self.defaultValue):
            for x in self.defaultValue[len(passedValues):]:
                passedValues.append(x)

        if self.isList:
            value = passedValues
        else:
            value = passedValues[0]

        if self.multipleKey:
            
            if not hasattr(self, "value"):
                self.value = [value]
            else:
                self.value.append(value)
            
            if self.flattenList and isinstance(self.value[0], list):
                temp = []
                for sublist in self.value:
                    for item in sublist:
                        temp.append(item)
                self.value = temp
        else:
            self.value = value

###################################################################

def boolean(x):
    if x.lower() in ["y", "yes", "t", "true"]:
        return True
    elif x.lower() in ["n", "no", "f", "false"]:
        return False
    else:
        string = ("%s is an invalid string to turn into a boolean value. "
                  "Must be one of: y, yes, t, true, n, no, f, false.")
        raise ValueError(string)

def formatConversion(x):
    try:
        test = 2
        x % test
    except:
        string = "Unsupported format conversion string %s." % x
        raise ValueError(string)
    return x

#########################################################################

fitsedParams = [Param("fitting_method", defaultValue="brutefluxspace", 
                      allowedValues=["brutedaisychain","brutecolorspace",
                      "brutefluxspace","brutefiterrorbars"]),
                Param("model_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param("data_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param("model_param", isList=True, multipleKey=True,
                      dataType=[int,str,formatConversion,bool],
                      defaultValue=[0,'model_param','%.4f',0]),
                Param("output_overwrite", defaultValue=False, dataType=bool),
                Param('dchi2',dataType=float),
                Param("data_param", isList=True, multipleKey=True,
                      defaultValue=[ "data_param",0, "%.4f"],
                      dataType=[ int, str, formatConversion],
                      canBeNeg=False, maxSize=3),
                Param("mag_softening", defaultValue=0.03, dataType=float,
                      canBeNeg=False),
                Param("model_file", mandatory=True, canHaveSpace=False),
                Param("output_file", mandatory=True, canHaveSpace=False),
                Param("output_dir", defaultValue=os.getcwd()),
                Param("data_file", mandatory=True, canHaveSpace=False),
                Param("model_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False, 
                      sameLengthKeys=["data_flux_columns"]),
                Param("data_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False, 
                      sameLengthKeys=["model_flux_columns"]),
                Param("data_error_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False,
                      sameLengthKeys=["data_flux_columns"]),
                Param("mcits",dataType=int,defaultValue=1),
                Param("oldschoolmc",dataType=bool,defaultValue=False),
                Param("restrict_model", isList=True, multipleKey=True,
                      defaultValue=[0, "range", 0., 1.], maxSize=4,
                      dataType=[int, str, float, float],
                      canBeNeg=[False, True, True, True],
                      allowedValues=[None, ["range", "closecol", "closeval"], 
                                     None, None]),
                Param("output_bestfit_spectra", dataType=bool),
                Param("output_chisq_matrix", dataType=bool),
                Param("output_mc_results", dataType=bool),
                Param("restrict_data", isList=True, multipleKey=True,
                      defaultValue=[0, "range", 0., 1.], maxSize=4,
                      dataType=[int, str, float, float],
                      canBeNeg=[False, True, True, True],
                      allowedValues=[None, ["range", "closecol", "closeval"], 
                                     None, None]),
                Param("data_mag_offsets", isList=True, dataType=float,
                      sameLengthKeys=["data_flux_columns"]),
                Param("data_wavelengths", isList=True, dataType=float,
                      sameLengthKeys=["data_flux_columns"]),
                Param("yourname",defaultValue="playa")
                      ]
#------------------------------------------------------------------------------ 
makesedParams = [
                Param("rffmt", mandatory=True,  
                      dataType=str, allowedValues=['galaxev']),
                Param("dotsed", dataType=str),
                Param("dotfourcolor", dataType=str),
                Param("filter_dir",dataType=str),
                Param("filter_names",dataType=[str,str],  
                       multipleKey=True,mandatory=True,isList=True),
                Param("redshifts", isList=True, dataType=[str,float,float,float],
                      allowedValues=[['range','values'],None,None,None],
                      ),
                Param("igm_law", dataType=str, allowedValues=['madau','inoue']),
                Param("igm_opacities",dataType=[str,float,float,float],isList=True,
                    allowedValues=[['range','values'],None,None,None]),
                Param("dust_law",dataType=str, allowedValues=['calzetti2000',\
                'calzetti1997','lmc','smc','mw','dor30']),
                Param("ebvs",isList=True,dataType=[str,float,float,float],
                      allowedValues=[['range','values'],None,None,None]),
                Param('returnmags',dataType=bool,defaultValue=False),
                Param('output_file',dataType=str),
                Param('cosmology',dataType=[str,float,float,float],allowedValues=\
                    [['lcdm','wmap'],None,None,None],isList=True),
                Param('models',dataType=[str,int,int],\
                      isList=True,allowedValues=[['all','range','values'],None,None])
        
                ]
                
#------------------------------------------------------------------------------ 
                
spaceCheckParams = [
                Param("model_file", mandatory=True, canHaveSpace=False),
                Param("data_file", mandatory=True, canHaveSpace=False),
                Param("filter_names",isList=True,dataType=str),
                Param("model_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False, 
                      sameLengthKeys=["data_flux_columns"]),
                Param("data_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False, 
                      sameLengthKeys=["model_flux_columns"]),
                Param("model_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param("data_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param('nbins',mandatory=False,defaultValue=None,dataType=int),
                Param('data_plot_filename',dataType=str,defaultValue=None),
                Param('model_plot_filename',dataType=str,defaultValue=None),
                Param('outsider_flag_file',dataType=str,defaultValue=None)
                ]                
                
                
#------------------------------------------------------------------------------               
def SetMakeSedParams(pfile,args=None):
    params = SetParams(pfile, makesedParams, args)
    return(params)
    
def SetFitSedParams(pfile,args=None):
    params = SetParams(pfile, fitsedParams, args)
    print('...')
    return(params)
    
def SetSpaceCheckParams(pfile,args=None):
    params = SetParams(pfile, spaceCheckParams, args)
    return(params)

