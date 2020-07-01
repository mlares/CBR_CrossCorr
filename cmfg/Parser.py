import numpy as np
from configparser import ConfigParser
import itertools
import pandas as pd
import pickle
import sys
from tqdm import tqdm
from astropy import units as u

DEFAULT_INI = '../set/set_experiment.ini'
 
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def choice_yn(string, default_choice=None):
    if string.lower() in 'yesitrue':
        choice = True
    elif string.lower() in 'nofalse':
        choice = False
    else:
        if isinstance(default_choice, bool):
            choice = default_choice
        else:
            raise AttributeError('dfasf')
    return choice


# Idea: change named tuple by class:
# class parameters:
# initialize with default parameters.


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
                f"\n\t {DEFAULT_INI}"
                "\n\t Using default configuration file")
        if isinstance(sys_args, str):
            if isfile(sys_args):
                msg = f"Loading configuration parameters from {sys_args}"
                self.message = msg
                filename = sys_args
            else:
                self.message = "Input argument is not a valid file\
                                Using default configuration file instead"
                filename = DEFAULT_INI

        elif isinstance(sys_args, list):

            if len(sys_args) == 2:
                filename = sys_args[1]

                if isfile(filename):
                    msg = f"Loading configuration parameters from {filename}"
                    self.message = msg
                else:
                    self.message = mess
                    filename = DEFAULT_INI
            else:
                self.message = mess
                filename = DEFAULT_INI

        else:
            self.message = mess
            filename = DEFAULT_INI

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

        experiment_id = self['experiment']['experiment_ID']
        datadir_cmb = self['cmb']['datadir_cmb']
        filedata_cmb_mapa = self['cmb']['filedata_cmb_mapa']
        filedata_cmb_mask = self['cmb']['filedata_cmb_mask']
        datadir_glx = self['glx']['datadir_glx']
        filedata_glx = self['glx']['filedata_glx']
        pickle_name_root = self['out']['pickle_name_root']
        pickle_name_ext = self['out']['pickle_name_ext']
        dir_plots = self['out']['dir_plots']
        plot_name_root = self['out']['plot_name_root']
        plot_format = self['out']['plot_format']
        plot_ftype = self['out']['plot_ftype']
        plot_fname = plot_name_root + experiment_id + '.' + plot_format


        fname = dir_plots + plot_fname + '_' + experiment_id + plot_ftype

        names = 'experiment_id \
                 datadir_cmb \
                 filedata_cmb_mapa \
                 filedata_cmb_mask \
                 datadir_glx \
                 filedata_glx \
                 pickle_name_root \
                 pickle_name_ext \
                 dir_plots \
                 plot_name_root \
                 plot_ftype \
                 plot_format'

        parset = namedtuple('pars', names)

        res = parset(experiment_id,
                     datadir_cmb,
                     filedata_cmb_mapa,
                     filedata_cmb_mask,
                     datadir_glx,
                     filedata_glx,
                     pickle_name_root,
                     pickle_name_ext,
                     dir_plots,
                     plot_name_root,
                     plot_ftype,
                     plot_format)

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
        from collections import namedtuple

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

        experiment_id = self['experiment']['experiment_id']

        n_jobs = int(self['run']['n_jobs'])
        
        adaptative_resolution = self['run']['adaptative_resolution']
        
        dir_output = self['out']['dir_output']
        dir_plots = self['out']['dir_plots']

        max_centers = self['glx']['max_centers']
        try:
            int(max_centers)
        except Exception:
            max_centers = None
        else:
            max_centers = int(max_centers)

        choice = self['glx']['control_sample']
        control_sample = choice_yn(choice, default_choice=False)
                               
        choice = self['glx']['control_ranmap']
        control_ranmap = choice_yn(choice, default_choice=False)

        choice = self['glx']['control_angles']
        control_angles = choice_yn(choice, default_choice=False)

        norm_to = False
        r_units_str = self['run']['r_units'].lower()
        if r_units_str == 'arcmin':
            r_units = u.arcmin
        elif r_units_str == 'arcsec':
            r_units = u.arcsec
        elif r_units_str == 'parsec':
            r_units = u.parsec
        elif r_units_str == 'kpc':
            r_units = u.kpc
        elif r_units_str in ['physical']:
            norm_to = 'PHYSICAL'
            r_units = 1.*u.dimensionless_unscaled
        elif r_units_str in ['angular']:
            norm_to = 'ANGULAR'
            r_units = 1.*u.dimensionless_unscaled 
        elif r_units_str in ['cosine']:
            norm_to = 'COS'
            r_units = 1.*u.dimensionless_unscaled 
        else:
            print('Warning: not recognized radial unit or normalization')
            r_units = 1.

        r_start = float(self['run']['r_start'])
        r_stop = float(self['run']['r_stop'])
        r_start = r_start*r_units
        r_stop = r_stop*r_units
        r_n_bins = int(self['run']['r_n_bins'])

        theta_units_str = self['run']['theta_units']
        if theta_units_str == 'arcmin':
            theta_units = u.arcmin
        elif theta_units_str == 'arcsec':
            theta_units = u.arcsec
        elif theta_units_str == 'rad':
            theta_units = u.rad
        elif theta_units_str == 'deg':
            theta_units = u.deg
        else:
            theta_units = 1.
        theta_n_bins = int(self['run']['theta_n_bins']) 
 
        theta_start = self['run']['theta_start']
        if is_number(theta_start):
            theta_start = float(theta_start)
            num = 1.
        elif 'pi' in theta_start:
            n = theta_start.replace('pi','').replace('*','')
            try:
                float(n)
            except:
                num = 1.
            else:
                num = float(n)
            theta_start = num * np.pi
        else:
            print('Error: number not recognized in theta_start')
            exit()
        theta_start = theta_start*theta_units                                 

        theta_stop = self['run']['theta_stop']
        if is_number(theta_stop):
            theta_stop = float(theta_stop)
            num = 1.
        elif 'pi' in theta_stop:
            n = theta_stop.replace('pi','').replace('*','')
            try:
                float(n)
            except:
                num = 1.
            else:
                num = float(n)
            theta_stop = num * np.pi
        else:
            print('Error: number not recognized in theta_stop')
            exit()
        theta_stop = theta_stop*theta_units


        choice = self['run']['disk_align']
        disk_align = choice_yn(choice, default_choice=False)
 
        choice = self['run']['adaptative_resolution']
        adaptative_resolution = choice_yn(choice, default_choice=False)

        choice = self['run']['run_parallel']
        run_parallel = choice_yn(choice, default_choice=False)

        choice = self['UX']['show_progress']
        showp = choice_yn(choice, default_choice=False)

        choice = self['out']['clobber']
        overwrite = choice_yn(choice, default_choice=False)

        choice = self['UX']['verbose']
        verbose = choice_yn(choice, default_choice=True)
        if verbose:
            print('loading parameters...')
 
        galaxy_types = self['run']['galaxy_types']
        galaxy_types = galaxy_types.split(' ')

        redshift_min = float(self['run']['redshift_min'])
        redshift_max = float(self['run']['redshift_max'])
        ellipt_min = float(self['run']['ellipt_min'])
        ellipt_max = float(self['run']['ellipt_max'])


        r_avg_cuts = self['run']['r_avg_cuts'].split(' ')
        if isinstance(r_avg_cuts, str):
            r_avg_cuts = r_n_bins
        else:
            r_avg_cuts = [int(i) for i in r_avg_cuts]
        r_avg_fact = float(self['run']['r_avg_fact'])


        names = ['experiment_id',
                 'n_jobs',
                 'control_sample',
                 'control_ranmap',
                 'control_angles',
                 'r_start',
                 'r_stop',
                 'r_n_bins',
                 'r_units',
                 'r_avg_cuts',
                 'r_avg_fact',
                 'theta_start',
                 'theta_stop',
                 'theta_n_bins',
                 'theta_units',
                 'norm_to',
                 'adaptative_resolution',
                 'disk_align',
                 'galaxy_types',
                 'redshift_min',
                 'redshift_max',
                 'ellipt_min',
                 'ellipt_max',
                 'max_centers',
                 'verbose',
                 'run_parallel',
                 'showp',
                 'overwrite',
                 'dir_output',
                 'dir_plots']

        names = ' '.join(names)

        parset = namedtuple('pars', names)

        res = parset(experiment_id,
                     n_jobs,
                     control_sample,
                     control_ranmap,
                     control_angles,
                     r_start,
                     r_stop,
                     r_n_bins,
                     r_units,
                     r_avg_cuts,
                     r_avg_fact,
                     theta_start,
                     theta_stop,
                     theta_n_bins,
                     theta_units,
                     norm_to,
                     adaptative_resolution,
                     disk_align,
                     galaxy_types,
                     redshift_min,
                     redshift_max,
                     ellipt_min,
                     ellipt_max,
                     max_centers,
                     verbose,
                     run_parallel,
                     showp,
                     overwrite,
                     dir_output,
                     dir_plots)

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

        # Check or create output directory for the current experiment
        if path.isdir(self.p.dir_output):
            dir_exp = (f"{self.p.dir_output}/"
                       f"{self.p.experiment_id}")
            try:
                makedirs(dir_exp)
                if self.p.verbose:
                    print("Directory ", dir_exp,  " Created ")
            except FileExistsError:
                pass
                # directory already exists
        else:
            msg = f"Directory {self.p.dir_output} does not exist!"
            raise NotADirectoryError(msg)
 
        #if path.isdir(self.p.dir_output):
        #    print(f"Directory {self.p.dir_output} does not exist")

        #    try:
        #        makedirs(self.p.dir_output)
        #        if self.p.verbose:
        #            print("Directory ", self.p.dir_output,  " Created ")
        #    except FileExistsError:
        #        # directory already exists
        #        pass


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

        # plots directory for this experiment
        dir_plt = (f"{self.p.dir_plots}/"
                   f"{self.p.experiment_id}")

        if not path.isdir(dir_plt):
            print(f"Directory {dir_plt} does not exist")

            try:
                makedirs(dir_plt)
                if self.p.verbose:
                    print("Directory ", dir_plt,  " Created ")
            except FileExistsError:
                # directory already exists
                pass
 
