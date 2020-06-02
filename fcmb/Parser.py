import numpy as np
from configparser import ConfigParser
import itertools
import pandas as pd
import pickle
import sys
from tqdm import tqdm


class Parser(ConfigParser):
    """parser class.

    Manipulation of configuration parameters. This method allows to read a
    configuration file or to set parameters for a Constrained Causally
    Conected Network (C3Net) model.
    """

    def __init__(self, argv=None, *args, **kwargs):
        """Initialize a parser.

        Parameters
        ----------
            None
        Returns
        -------
            None
        Raises
        ------
            None
        """
        super().__init__()
        self.message = None
        self.check_file(argv)
        self.read_config_file()

        self.load_filenames()
        self.load_config(*args, **kwargs)
        self.check_settings()

    def check_file(self, sys_args=""):
        """Parse paramenters for the simulation from a .ini file.

        Parameters
        ----------
            filename (str): the file name of the map to be read

        Raises
        ------
            None

        Returns
        -------
            None
        """
        from os.path import isfile

        mess = ("Configuration file expected:"
                "\n\t filename or CLI input"
                "\n\t example:  python run_correlation.py"
                "\n\t ../set/experiment.ini"
                "\n\t Using default configuration file")
        if isinstance(sys_args, str):
            if isfile(sys_args):
                msg = f"Loading configuration parameters from {sys_args}"
                self.message = msg
                filename = sys_args
            else:
                self.message = "Input argument is not a valid file\
                                Using default configuration file instead"
                filename = '../set/experiment.ini'

        elif isinstance(sys_args, list):

            if len(sys_args) == 2:
                filename = sys_args[1]

                if isfile(filename):
                    msg = f"Loading configuration parameters from {filename}"
                    self.message = msg
                else:
                    self.message = mess
                    filename = '../set/experiment.ini'
            else:
                self.message = mess
                filename = '../set/experiment.ini'

        else:
            self.message = mess
            filename = '../set/experiment.ini'

        self.filename = filename

    def read_config_file(self):
        """Parse paramenters for the simulation from a .ini file.

        Parameters
        ----------
            None

        Raises
        ------
            None

        Returns
        -------
            None
        """
        self.read(self.filename)

    def load_filenames(self):
        """Make filenames based on info in config file.

        Parameters
        ----------
            None

        Raises
        ------
            None

        Returns
        -------
            list of filenames
        """
        from collections import namedtuple

        # Experiment settings
        exp_id = self['experiment']['exp_id']
        dir_plots = self['output']['dir_plots']
        pars_root = self['output']['pars_root']
        progress_root = self['output']['progress_root']
        dir_output = self['output']['dir_output']
        plot_fname = self['output']['plot_fname']
        plot_ftype = self['output']['plot_ftype']

        fname = dir_plots + plot_fname + '_' + exp_id + plot_ftype

        names = 'exp_id \
                 dir_plots \
                 dir_output \
                 pars_root \
                 progress_root \
                 plot_fname \
                 plot_ftype \
                 fname'

        parset = namedtuple('pars', names)

        res = parset(exp_id,
                     dir_plots,
                     dir_output,
                     pars_root,
                     progress_root,
                     plot_fname,
                     plot_ftype,
                     fname)

        self.filenames = res

    def load_config(self, keys=None, values=None, nran=None,
                    *args, **kwargs):
        """Load parameters from config file.

        Parameters
        ----------
            None

        Raises
        ------
            None

        Returns
        -------
            list of parameters as a named tuple
        """
        if isinstance(keys, list):
            # override configuration file with arguments
            if len(keys) != len(values):
                print('Error overriding parameters (using file values)')
            else:
                for k, v in zip(keys, values):
                    for sec in self.sections():
                        has = self.has_option(sec, k)
                        if has:
                            self[sec][k] = v

        choice = self['UX']['verbose']
        if choice.lower() in 'yesitrue':
            verbose = True
        elif choice.lower() in 'nofalse':
            verbose = False
        else:
            print('warning in .ini file: UX: verbose')
            verbose = False

        if verbose:
            print('loading parameters...')
        from collections import namedtuple

        ghz_inner = float(self['simu']['ghz_inner'])
        ghz_outer = float(self['simu']['ghz_outer'])

        t_max = float(self['simu']['t_max'])

        tau_a_min = float(self['simu']['tau_a_min'])
        tau_a_max = float(self['simu']['tau_a_max'])
        tau_a_nbins = int(self['simu']['tau_a_nbins'])

        tau_s_min = float(self['simu']['tau_s_min'])
        tau_s_max = float(self['simu']['tau_s_max'])
        tau_s_nbins = int(self['simu']['tau_s_nbins'])

        d_max_min = float(self['simu']['d_max_min'])
        d_max_max = float(self['simu']['d_max_max'])
        d_max_nbins = int(self['simu']['d_max_nbins'])

        if nran is None:
            nran = int(self['simu']['nran'])

        choices = self['simu']['run_parallel']
        if choices.lower() in 'yesitrue':
            run_parallel = True
        elif choices.lower() in 'nofalse':
            run_parallel = False
        else:
            run_parallel = False

        # Experiment settings
        exp_id = self['experiment']['exp_id']
        njobs = int(self['simu']['njobs'])
        dir_plots = self['output']['dir_plots']
        dir_output = self['output']['dir_output']
        pars_root = self['output']['pars_root']
        plot_fname = self['output']['plot_fname']
        plot_ftype = self['output']['plot_ftype']
        fname = dir_plots + plot_fname + '_' + exp_id + plot_ftype

        choice = self['UX']['show_progress']
        if choice.lower() in 'yesitrue':
            showp = True
        elif choice.lower() in 'nofalse':
            showp = False
        else:
            print('warning in .ini file: UX: show_progress')
            showp = False

        string_overwrite = self['output']['clobber']
        if string_overwrite.lower() in 'yesitrue':
            overwrite = True
        elif string_overwrite.lower() in 'nofalse':
            overwrite = False
        else:
            print('warning in .ini file: output: clobber')
            overwrite = False

        names = ['ghz_inner',
                 'ghz_outer',
                 't_max',
                 'tau_a_min',
                 'tau_a_max',
                 'tau_a_nbins',
                 'tau_s_min',
                 'tau_s_max',
                 'tau_s_nbins',
                 'd_max_min',
                 'd_max_max',
                 'd_max_nbins',
                 'nran',
                 'run_parallel',
                 'njobs',
                 'exp_id',
                 'dir_plots',
                 'dir_output',
                 'pars_root',
                 'plot_fname',
                 'plot_ftype',
                 'fname',
                 'showp',
                 'overwrite',
                 'verbose']
        names = ' '.join(names)

        parset = namedtuple('pars', names)

        res = parset(ghz_inner,
                     ghz_outer,
                     t_max,
                     tau_a_min,
                     tau_a_max,
                     tau_a_nbins,
                     tau_s_min,
                     tau_s_max,
                     tau_s_nbins,
                     d_max_min,
                     d_max_max,
                     d_max_nbins,
                     nran,
                     run_parallel,
                     njobs,
                     exp_id,
                     dir_plots,
                     dir_output,
                     pars_root,
                     plot_fname,
                     plot_ftype,
                     fname,
                     showp,
                     overwrite,
                     verbose)

        self.p = res

    def check_settings(self):
        """Check if parameters make sense.

        Parameters
        ----------
            None

        Raises
        ------
            None

        Returns
        -------
            Exception if settings have inconsistencies.
        """
        from os import path, makedirs

        if self.p.verbose:
            print(self.message)
            print('Checking settings...')

        # output directory
        if not path.isdir(self.p.dir_output):
            print(f"Directory {self.p.dir_output} does not exist")

            try:
                makedirs(self.p.dir_output)
                if self.p.verbose:
                    print("Directory ", self.p.dir_output,  " Created ")
            except FileExistsError:
                # directory already exists
                pass

        # experiment directory
        ID_dir = self.p.dir_output + self.p.exp_id
        if not path.isdir(ID_dir):
            print(f"Directory {ID_dir} does not exist")

            try:
                makedirs(ID_dir)
                if self.p.verbose:
                    print("Directory ", ID_dir,  " Created ")
            except FileExistsError:
                # directory already exists
                pass

        # plots directory
        if not path.isdir(self.p.dir_plots):
            print(f"Directory {self.p.dir_plots} does not exist")

            try:
                makedirs(self.p.dir_plots)
                if self.p.verbose:
                    print("Directory ", self.p.dir_plots,  " Created ")
            except FileExistsError:
                # directory already exists
                pass

