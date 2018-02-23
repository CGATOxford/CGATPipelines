"""Parameters.py - Parameter handling for ruffus pipelines
==========================================================

Reference
---------

"""

import re
import collections
import os

import configparser

import sys

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

from CGATPipelines.Pipeline.Utils import getCallerLocals, isTest

# sort out script paths

# root directory of CGAT Code collection
CGATSCRIPTS_ROOT_DIR = os.path.dirname(
    os.path.dirname(E.__file__))

# CGAT Code collection scripts
CGATSCRIPTS_SCRIPTS_DIR = os.path.join(CGATSCRIPTS_ROOT_DIR, "CGAT", "scripts")

# root directory of CGAT Pipelines
CGATPIPELINES_ROOT_DIR = os.path.dirname(os.path.dirname(
    os.path.dirname(__file__)))

# CGAT Pipeline scripts
CGATPIPELINES_SCRIPTS_DIR = os.path.join(CGATPIPELINES_ROOT_DIR,
                                         "scripts")
# Directory of CGAT pipelines
CGATPIPELINES_PIPELINE_DIR = os.path.join(CGATPIPELINES_ROOT_DIR,
                                          "CGATPipelines")
# CGAT Pipeline R scripts
CGATPIPELINES_R_DIR = os.path.join(CGATPIPELINES_ROOT_DIR, "R")

# Location of conda installation folder
CGATPIPELINES_CONDA_DIR = os.path.abspath(CGATPIPELINES_ROOT_DIR + "/../conda-install")

# if Pipeline.py is called from an installed version, scripts are
# located in the "bin" directory.
if not os.path.exists(CGATSCRIPTS_SCRIPTS_DIR):
    SCRIPTS_DIR = os.path.join(sys.exec_prefix, "bin")

if not os.path.exists(CGATPIPELINES_SCRIPTS_DIR):
    PIPELINE_SCRIPTS_DIR = os.path.join(sys.exec_prefix, "bin")

# Global variable for configuration file data
CONFIG = configparser.SafeConfigParser()


class TriggeredDefaultFactory:
    with_default = False

    def __call__(self):
        if TriggeredDefaultFactory.with_default:
            return str()
        else:
            raise KeyError("missing parameter accessed")

# Global variable for parameter interpolation in commands
# This is a dictionary that can be switched between defaultdict
# and normal dict behaviour.
PARAMS = collections.defaultdict(TriggeredDefaultFactory())

# patch - if --help or -h in command line arguments,
# switch to a default dict to avoid missing paramater
# failures
if isTest() or "--help" in sys.argv or "-h" in sys.argv:
    TriggeredDefaultFactory.with_default = True

# A list of hard-coded parameters within the CGAT environment
# These can be overwritten by command line options and
# configuration files
HARDCODED_PARAMS = {
    'scriptsdir': CGATSCRIPTS_SCRIPTS_DIR,
    'toolsdir': CGATSCRIPTS_SCRIPTS_DIR,
    'pipeline_scriptsdir': CGATPIPELINES_SCRIPTS_DIR,
    'pipelinedir': CGATPIPELINES_PIPELINE_DIR,
    'pipeline_rdir': CGATPIPELINES_R_DIR,
    'pipelines_conda_dir': CGATPIPELINES_CONDA_DIR,
    # script to perform map/reduce like computation.
    'cmd-farm': """python %(pipeline_scriptsdir)s/farm.py
                --method=drmaa
                --bashrc=%(pipeline_scriptsdir)s/bashrc.cgat
                --cluster-options=%(cluster_options)s
                --cluster-queue=%(cluster_queue)s
                --cluster-num-jobs=%(cluster_num_jobs)i
                --cluster-priority=%(cluster_priority)i
                --cluster-queue-manager=%(cluster_queue_manager)s
                --cluster-memory-resource=%(cluster_memory_resource)s
                --cluster-memory-default=%(cluster_memory_default)s
    """,
    # command to get tab-separated output from database
    'cmd-sql': """sqlite3 -header -csv -separator $'\\t' """,
    # DEPRECATED: options to use for csv2db upload
    "csv2db_options": "--backend=sqlite --retry --map=gene_id:str "
    "--map=contig:str --map=transcript_id:str",
    # database backend
    'database_backend': "sqlite",
    # database host
    'database_host': "",
    # name of database
    'database_name': "csvdb",
    # database connection options
    'database_username': "cgat",
    # database password - if required
    'database_password': "",
    # database port - if required
    'database_port': 3306,
    # wrapper around non-CGAT scripts
    'cmd-run': """%(pipeline_scriptsdir)s/run.py""",
    # legacy directory used for temporary local files
    #     Use of this var can be problematic (issue #174)
    #     - it may be depreciated.
    'tmpdir': os.environ.get("TMPDIR", '/scratch'),
    # directory used for temporary local tempory files on compute nodes
    # *** passed directly to the shell      ***
    # *** may not exist on login/head nodes ***
    # default matches 'tmpdir' only for backwards compatibility
    # typically a shell environment var is expected, e.g.
    # 'local_tmpdir': '$SCRATCH_DIR',
    'local_tmpdir': os.environ.get("TMPDIR", '/scratch'),
    # directory used for temporary files shared across machines
    'shared_tmpdir': os.environ.get("SHARED_TMPDIR", "/ifs/scratch"),
    # queue manager (supported: sge, slurm, torque, pbspro)
    'cluster_queue_manager': 'sge',
    # cluster queue to use
    'cluster_queue': 'all.q',
    # priority of jobs in cluster queue
    'cluster_priority': -10,
    # number of jobs to submit to cluster queue
    'cluster_num_jobs': 100,
    # name of consumable resource to use for requesting memory
    'cluster_memory_resource': "mem_free",
    # amount of memory set by default for each job
    'cluster_memory_default': "2G",
    # general cluster options
    'cluster_options': "",
    # parallel environment to use for multi-threaded jobs
    'cluster_parallel_environment': 'dedicated',
    # ruffus job limits for databases
    'jobs_limit_db': 10,
    # ruffus job limits for R
    'jobs_limit_R': 1,
}

# After all configuration files have been read, some
# parameters need to be interpolated with other parameters
# The list is below:
INTERPOLATE_PARAMS = ('cmd-farm', 'cmd-run')


def configToDictionary(config):
    """convert the contents of a :py:class:`ConfigParser.ConfigParser`
    object to a dictionary

    This method works by iterating over all configuration values in a
    :py:class:`ConfigParser.ConfigParser` object and inserting values
    into a dictionary. Section names are prefixed using and underscore.
    Thus::

        [sample]
        name=12

    is entered as ``sample_name=12`` into the dictionary. The sections
    ``general`` and ``DEFAULT`` are treated specially in that both
    the prefixed and the unprefixed values are inserted: ::

       [general]
       genome=hg19

    will be added as ``general_genome=hg19`` and ``genome=hg19``.

    Numbers will be automatically recognized as such and converted into
    integers or floats.

    Returns
    -------
    config : dict
        A dictionary of configuration values

    """
    p = {}
    for section in config.sections():
        for key, value in config.items(section):
            try:
                v = IOTools.str2val(value)
            except TypeError:
                E.error("error converting key %s, value %s" % (key, value))
                E.error("Possible multiple concurrent attempts to "
                        "read configuration")
                raise

            p["%s_%s" % (section, key)] = v
            if section in ("general", "DEFAULT"):
                p["%s" % (key)] = v

    for key, value in config.defaults().items():
        p["%s" % (key)] = IOTools.str2val(value)

    return p


def inputValidation(PARAMS, pipeline_script=""):
    '''Inspects the PARAMS dictionary looking for problematic input values.

    So far we just check that:

        * all required 3rd party tools are on the PATH

        * input parameters are not empty

        * input parameters do not contain the "?" character (used as a
          placeholder in different pipelines)

        * if the input is a file, check whether it exists and
          is readable
    '''

    E.info('''=== Input Validation starts ===''')
    E.info(''' Checking 3rd party dependencies ''')

    ### check 3rd party dependencies ###
    if len(pipeline_script) > 0:
        # this import requires the PYTHONPATH in the following order
        # PYTHONPATH=<src>/CGATPipelines:<src>/cgat
        import scripts.cgat_check_deps as cd
        deps, check_path_failures = cd.checkDepedencies(pipeline_script)
        # print info about dependencies
        if len(deps) == 0:
            print('\nNo dependencies found.\n')
        else:
        # print dictionary ordered by value
            for k in sorted(deps, key=deps.get, reverse=True):
                print('\nProgram: {0!s} used {1} time(s)'.format(k, deps[k]))

            n_failures = len(check_path_failures)
            if n_failures == 0:
                print('\nCongratulations! All required programs are available on your PATH\n')
            else:
                print('\nThe following programs are not on your PATH')
                for p in check_path_failures:
                    print('\n{0!s}'.format(p))
                print

    ## check PARAMS
    num_missing = 0
    num_questions = 0

    E.info(''' Checking pipeline configuration ''')

    for key, value in sorted(PARAMS.items()):

        key   = str(key)
        value = str(value)

        # check for missing values
        if value == "":
            print('\n"{}" is empty, is that expected?'.format(key))
            num_missing += 1

        # check for a question mark in the dictironary (indicates
        # that there is a missing input parameter)
        if "?" in value:
            print('\n"{}" is not defined (?), is that expected?'.format(key))
            num_questions += 1

        # validate input files listed in PARAMS
        if (value.startswith("/") \
           or value.endswith(".gz") \
           or value.endswith(".gtf")) \
           and "," not in value:

            if os.access(value, os.R_OK):
                pass
            else:
                print('\n"{}": "{}" is not readable'.format(key, value))

    if num_missing == 0 and num_questions == 0:
        while True:
            confirmation = input('''
            ##########################################################

            Your input data seems all correct, congratulations!

            Do you want to continue running the pipeline? (y/n)

            ##########################################################

            ''')
            if confirmation.lower() == "y":
                E.info('=== Input Validation ends ===')
                break
            elif confirmation.lower() == "n":
                E.info('=== Input Validation ends ===')
                E.info('Pipeline aborted.')
                sys.exit(0)
    else:
        while True:
            start_pipeline = input('''
            ###########################################################

            Please check the WARNING messages and if you are
            happy then enter "y" to continue or "n" to abort running
            the pipeline.

            ###########################################################
            ''')
            if start_pipeline.lower() == "y":
                E.info('=== Input Validation ends ===')
                break
            if start_pipeline.lower() == "n":
                E.info('=== Input Validation ends ===')
                E.info('Pipeline aborted.')
                sys.exit(0)


def getParameters(filenames=["pipeline.ini", ],
                  defaults=None,
                  site_ini=True,
                  user_ini=True,
                  default_ini=True,
                  only_import=None):
    '''read a config file and return as a dictionary.

    Sections and keys are combined with an underscore. If a key
    without section does not exist, it will be added plain.

    For example::

       [general]
       input=input1.file

       [special]
       input=input2.file

    will be entered as { 'general_input' : "input1.file",
    'input: "input1.file", 'special_input' : "input2.file" }

    This function also updates the module-wide parameter map.

    The section [DEFAULT] is equivalent to [general].

    The order of initialization is as follows:

    1. hard-coded defaults
    2. pipeline specific default file in the CGAT code installation
    3. :file:`.cgat` in the users home directory
    4. files supplied by the user in the order given

    If the same configuration value appears in multiple
    files, later configuration files will overwrite the
    settings form earlier files.

    Path names are expanded to the absolute pathname to avoid
    ambiguity with relative path names. Path names are updated
    for parameters that end in the suffix "dir" and start with
    a "." such as "." or "../data".

    Arguments
    ---------
    filenames : list
       List of filenames of the configuration files to read.
    defaults : dict
       Dictionary with default values. These will be overwrite
       any hard-coded parameters, but will be overwritten by user
       specified parameters in the configuration files.
    default_ini : bool
       If set, the default initialization file will be read from
       'CGATPipelines/configuration/pipeline.ini'
    user_ini : bool
       If set, configuration files will also be read from a
       file called :file:`.cgat` in the user`s home directory.
    only_import : bool
       If set to a boolean, the parameter dictionary will be a
       defaultcollection. This is useful for pipelines that are
       imported (for example for documentation generation) but not
       executed as there might not be an appropriate .ini file
       available. If `only_import` is None, it will be set to the
       default, which is to raise an exception unless the calling
       script is imported or the option ``--is-test`` has been passed
       at the command line.

    Returns
    -------
    config : dict
       Dictionary with configuration values.
    '''

    global CONFIG
    global PARAMS
    old_id = id(PARAMS)

    caller_locals = getCallerLocals()

    # check if this is only for import
    if only_import is None:
        only_import = isTest() or \
            "__name__" not in caller_locals or \
            caller_locals["__name__"] != "__main__"

    # important: only update the PARAMS variable as
    # it is referenced in other modules. Thus the type
    # needs to be fixed at import. Raise error where this
    # is not the case.
    # Note: Parameter sharing in the Pipeline module needs
    # to be reorganized.
    if only_import:
        # turn on default dictionary
        TriggeredDefaultFactory.with_default = True

    # Clear up ini files on the list that do not exist.
    # Please note the use of list(filenames) to create
    # a clone to iterate over as we remove items from
    # the original list (to avoid unexpected results)
    for fn in list(filenames):
        if not os.path.exists(fn):
            filenames.remove(fn)

    if site_ini:
        # read configuration from /etc/cgat/pipeline.ini
        fn = "/etc/cgat/pipeline.ini"
        if os.path.exists(fn):
            filenames.insert(0, fn)

    if default_ini:
        # The link between CGATPipelines and Pipeline.py
        # needs to severed at one point.
        # 1. config files into CGAT module directory?
        # 2. Pipeline.py into CGATPipelines module directory?
        filenames.insert(0,
                         os.path.join(CGATPIPELINES_PIPELINE_DIR,
                                      'configuration',
                                      'pipeline.ini'))

    if user_ini:
        # read configuration from a users home directory
        fn = os.path.join(os.path.expanduser("~"),
                          ".cgat")
        if os.path.exists(fn):
            index = filenames.index('pipeline.ini')
            filenames.insert(index,fn)

    # IMS: Several legacy scripts call this with a string as input
    # rather than a list. Check for this and correct

    if isinstance(filenames, str):
        filenames = [filenames]

    PARAMS['pipeline_ini'] = filenames

    try:
        CONFIG.read(filenames)
        p = configToDictionary(CONFIG)
    except configparser.InterpolationSyntaxError as ex:
        # Do not log, as called before logging module is initialized -
        # this will mess up loging configuration in Control.py and Experiment.py
        # E.debug(
        #     "InterpolationSyntaxError when reading configuration file, "
        #     "likely due to use of '%'. "
        #     "Please quote '%' if ini interpolation is required. "
        #     "Orginal error: {}".format(str(ex)))
        CONFIG = configparser.RawConfigParser()
        CONFIG.read(filenames)
        p = configToDictionary(CONFIG)
    
    # update with hard-coded PARAMS
    PARAMS.update(HARDCODED_PARAMS)

    if defaults:
        PARAMS.update(defaults)
    PARAMS.update(p)

    # interpolate some params with other parameters
    for param in INTERPOLATE_PARAMS:
        try:
            PARAMS[param] = PARAMS[param] % PARAMS
        except TypeError as msg:
            raise TypeError('could not interpolate %s: %s' %
                            (PARAMS[param], msg))

    # expand pathnames
    for param, value in list(PARAMS.items()):
        if param.endswith("dir"):
            if value.startswith("."):
                PARAMS[param] = os.path.abspath(value)

    # make sure that the dictionary reference has not changed
    assert id(PARAMS) == old_id

    return PARAMS


def loadParameters(filenames):
    '''load parameters from one or more files.

    Parameters are processed in the same way as :func:`getParameters`,
    but the global parameter dictionary is not updated.

    Arguments
    ---------
    filenames : list
       List of filenames of the configuration files to read.

    Returns
    -------
    config : dict
       A configuration dictionary.

    '''
    try:
        config = configparser.SafeConfigParser()
        config.read(filenames)
        p = configToDictionary(config)
    except configparser.InterpolationSyntaxError as ex:
        E.warn(
            "InterpolationSyntaxError when reading configuration file, "
            "likely due to use of '%'. "
            "Please quote '%' if ini interpolation is required. "
            "Orginal error: {}".format(str(ex)))
        config = configparser.RawConfigParser()
        config.read(filenames)
        p = configToDictionary(config)

    return p


def matchParameter(param):
    '''find an exact match or prefix-match in the global
    configuration dictionary param.

    Arguments
    ---------
    param : string
        Parameter to search for.

    Returns
    -------
    name : string
        The full parameter name.

    Raises
    ------
    KeyError if param can't be matched.

    '''
    if param in PARAMS:
        return param

    for key in list(PARAMS.keys()):
        if "%" in key:
            rx = re.compile(re.sub("%", ".*", key))
            if rx.search(param):
                return key

    raise KeyError("parameter '%s' can not be matched in dictionary" %
                   param)


def substituteParameters(**kwargs):
    '''return a parameter dictionary.

    This method builds a dictionary of parameter values to
    apply for a specific task. The dictionary is built in
    the following order:

    1. take values from the global dictionary (:py:data:`PARAMS`)
    2. substitute values appearing in `kwargs`.
    3. Apply task specific configuration values by looking for the
       presence of ``outfile`` in kwargs.

    The substition of task specific values works by looking for any
    parameter values starting with the value of ``outfile``.  The
    suffix of the parameter value will then be substituted.

    For example::

        PARAMS = {"tophat_threads": 4,
                  "tophat_cutoff": 0.5,
                  "sample1.bam.gz_tophat_threads" : 6}
        outfile = "sample1.bam.gz"
        print substituteParameters(**locals())
        {"tophat_cutoff": 0.5, "tophat_threads": 6}

    Returns
    -------
    params : dict
        Dictionary with parameter values.

    '''

    # build parameter dictionary
    # note the order of addition to make sure that kwargs takes precedence
    local_params = dict(list(PARAMS.items()) + list(kwargs.items()))

    if "outfile" in local_params:
        # replace specific parameters with task (outfile) specific parameters
        outfile = local_params["outfile"]
        for k in list(local_params.keys()):
            if k.startswith(outfile):
                p = k[len(outfile) + 1:]
                if p not in local_params:
                    raise KeyError(
                        "task specific parameter '%s' "
                        "does not exist for '%s' " % (p, k))
                E.debug("substituting task specific parameter "
                        "for %s: %s = %s" %
                        (outfile, p, local_params[k]))
                local_params[p] = local_params[k]

    return local_params


def asList(value):
    '''return a value as a list.

    If the value is a string and contains a ``,``, the string will
    be split at ``,``.

    Returns
    -------
    list

    '''
    if type(value) == str:
        try:
            values = [x.strip() for x in value.strip().split(",")]
        except AttributeError:
            values = [value.strip()]
        return [x for x in values if x != ""]
    elif type(value) in (list, tuple):
        return value
    else:
        return [value]


def isTrue(param, **kwargs):
    '''return True if param has a True value.

    A parameter is False if it is:

    * not set
    * 0
    * the empty string
    * false or False

    Otherwise the value is True.

    Arguments
    ---------
    param : string
        Parameter to be tested
    kwargs : dict
        Dictionary of local configuration values. These will be passed
        to :func:`substituteParameters` before evaluating `param`

    Returns
    -------
    bool

    '''
    if kwargs:
        p = substituteParameters(**kwargs)
    else:
        p = PARAMS
    value = p.get(param, 0)
    return value not in (0, '', 'false', 'False')


def checkParameter(param):
    """check if parameter ``key`` is set"""
    if param not in PARAMS:
        raise ValueError("need `%s` to be set" % param)


def getParams():
    """return handle to global parameter dictionary"""
    return PARAMS

def getCondaEnvironment(name):
    """prepare prefix to call programs in a different conda environment"""
#    return "PATH=" + str(CGATPIPELINES_CONDA_DIR) + "/envs/" + str(name) + "/bin;" + \
#            "CONDA_PREFIX=" + str(CGATPIPELINES_CONDA_DIR) + "/envs/" + str(name)
    return "source " + str(CGATPIPELINES_CONDA_DIR) + "/bin/activate " + str(name)
