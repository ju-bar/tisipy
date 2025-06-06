"""Module providing functions to run Dr. Probe command line programs."""

import os
from datetime import datetime
import subprocess
import numpy as np


def run_job(exename, arguments = None, cwb = None, cwd = None, verbose=False):
    """
    Starts a Dr. Probe command line program. This can be either of the three
    programs CELSLC, MSA or WAVIMG.

    Parameters:
        exename : str
            name of the executable to run (celslc, msa or wavimg)
        arguments : str, default None
            command line argument string
        cwb : str, default None
            executable location
        cwd : str, default None
            execution location
        verbose : boolean, default False
            text output flag

    Returns:
        dict
            number of errors, stdout of the call, error infos
    """
    sexefp = exename
    if cwb is not None:
        sexefp = os.path.join(cwb, exename)
    largs = sexefp
    if arguments is not None:
        largs += (' ' + arguments)
    if verbose:
        print(f"[drprobe_clt] calling {largs}")
        if cwd is not None:
            print(f"   working directory: {cwd}")
    # -----------------------------------------------------
    x = subprocess.run(largs, cwd=cwd, capture_output=True, text=True, check=False)
    # -----------------------------------------------------
    y = x.stdout.replace('\x00','').split('\n')
    l_err = []
    for sout in y:
        if sout.lower().find('error:') >= 0:
            l_err.append(sout)
    num_err = len(l_err)
    return {
        "num_err" : num_err,
        "stdout" : y,
        "errors" : l_err
    }

def check_prm_basics(prm):
    """
    Checks structure of prm for basic required components.
    """
    _ret = {'num_err' : 0, 'stdout' : [], 'errors' : []}
    # instrument section
    if 'instrument' not in prm:
        _ret['num_err'] += 1
        _ret['errors'].append("Instrument settings missing in prm.")
    else:
        prm_instr = prm['instrument']
        if 'beam' not in prm_instr:
            _ret['num_err'] += 1
            _ret['errors'].append("Instrument/beam settings missing in prm.")
        if 'optics' not in prm_instr:
            _ret['num_err'] += 1
            _ret['errors'].append("Instrument/optics settings missing in prm.")
        if 'detectors' not in prm_instr:
            _ret['num_err'] += 1
            _ret['errors'].append("Instrument/detectors settings missing in prm.")
    # sample section
    if 'sample' not in prm:
        _ret['num_err'] += 1
        _ret['errors'].append("Sample settings missing in prm.")
    else:
        prm_smp = prm['sample']
        if 'structure' not in prm_smp:
            _ret['num_err'] += 1
            _ret['errors'].append("Sample/structure settings missing in prm.")
        if 'thickness' not in prm_smp:
            _ret['num_err'] += 1
            _ret['errors'].append("Sample/thickness settings missing in prm.")
    # calculation section
    if 'calculation' not in prm:
        _ret['num_err'] += 1
        _ret['errors'].append("Calculation settings missing in prm.")
    else:
        prm_calc = prm['calculation']
        if 'program' not in prm_calc:
            _ret['num_err'] += 1
            _ret['errors'].append("Calculation/program settings missing in prm.")
        else:
            prm_bin = prm_calc['program']
            if prm_bin['name'] != "drprobe_clt":
                _ret['num_err'] += 1
                _ret['errors'].append("Calculation/program settings are not indicating a call to drprobe_clt.")
    # output section
    if 'output' not in prm:
        _ret['num_err'] += 1
        _ret['errors'].append("Output settings missing in prm.")
    else:
        prm_out = prm['output']
        if 'directory' not in prm_out:
            _ret['num_err'] += 1
            _ret['errors'].append("Missing output/directory in the settings.")

    return _ret


def get_dirnam(prm):
    """
    Returns the output directory and file name prefix
    """
    prm_out = prm['output']
    sdir = prm_out.get('directory', "")
    sout = prm_out.get('prefix', "out")
    return sdir, sout

def addlog(exename, logfile=None, infos=None):
    """
    Writes information to a log file.
    """
    if logfile is None:
        return
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    with open(logfile, 'a', encoding="utf-8") as flog:
        flog.write("\n")
        flog.write(f"{current_time} [tisipy.drprobe_clt] {exename}\n")
        if infos is not None:
            if len(infos) > 0:
                for s in infos:
                    flog.write(f"    {s}\n")
        flog.flush()
        flog.close()
    return


def read_cel(filename):
    """
    Reads structure data from a cel input file.
    """
    info = { 'file' : filename }
    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_lines = len(lines)
    if num_lines < 3:
        info['head'] = "Error reading cel file: unexpected content."
    else:
        info['head'] = lines[0]
        lsa = lines[1].strip().split()
        info['box'] = {
            'a' : float(lsa[1]),
            'b' : float(lsa[2]),
            'c' : float(lsa[3]),
            'alpha' : float(lsa[4]),
            'beta' : float(lsa[5]),
            'gamma' : float(lsa[6])
        }
        info['atoms'] = []
        for line in lines[2:]:
            if "*" in line.strip():
                break
            lsa = line.strip().split()
            info['atoms'].append({
                'symbol' : lsa[0],
                'pos' : [float(lsa[1]), float(lsa[2]), float(lsa[3])],
                'occupancy' : float(lsa[4]),
                'biso' : float(lsa[5]),
                'aniso' : [float(lsa[6]), float(lsa[7]), float(lsa[8])]
            })
    return info


def celslc(prm, verbose=False, logfile=None):
    """
    Execution of a call to the program CELSCL.
    If given, the program is executed in the output
    directory prm['output']['directory'].

    Parameters:
        prm : dict
            Simulation parameters
        verbose : boolean, default False
            Text output flag
        logfile : str, default None
            Log file name
    
    Returns:
        dict
            number of errors, stdout of the call, error infos
    """
    # --------------------
    # init
    _ret = check_prm_basics(prm)
    if _ret['num_err'] > 0:
        _ret['num_err'] += 1
        _ret['errors'].append("celslc initialization failed.")
    prm_calc = prm['calculation']
    prm_bin = prm_calc['program']
    swd, sout = get_dirnam(prm)
    # --------------------
    # generate the call
    _call = "celslc"
    if 'binary_path' in prm_bin:
        if prm_bin['binary_path'] is not None:
            _call = os.path.join(prm_bin['binary_path'], "celslc")
    # --------------------
    # generate the arguments
    _args = ""
    # structure model input
    celfile = prm['sample']['structure']['file']
    celfmt = 'cel' # assume cel format by default
    if celfile.endswith('.cif'):
        celfmt = "cif"
    _args = f"-{celfmt} {celfile}"
    # --------------------
    # slice file output
    _args += f" -slc {sout}"
    # --------------------
    # grid sampling
    _args += f" -nx {prm_calc['sampling']['x']}" + \
             f" -ny {prm_calc['sampling']['y']}" + \
             f" -nz {prm_calc['sampling']['z']}"
    # --------------------
    # beam energy
    _args += f" -ht {prm['instrument']['beam']['energy']}"
    # --------------------
    # phonon handling
    prm_pho = prm_calc['phonon']
    if "QEP" == prm_pho['model'] or "FRFPMS" == prm_pho['model']: # QEP or FRFPMS model
        _args += f" -fl -nv {prm_pho['configurations']}"
        if prm_calc['phonon']['absorption']:
            _args +=  " -abs"
        if "FRFPMS" == prm_pho['model']: # add Debye-Waller factor for FRFPMS model
            _args +=  " -dwf"
    if "absorptive" == prm_pho['model']: # absorptive model
        _args +=  " -dwf -abs"
    if "factor" == prm_pho['model']: # phenomenological model (TEM)
        _args += f" -dwf -abf {prm_pho['fraction']}"
    # --------------------
    # Run celslc
    addlog("celslc", logfile=logfile, infos=['arguments: ' + _args])
    _ret = run_job(_call, _args, cwd=swd, verbose=verbose)
    # --------------------
    # Handle output
    if verbose:
        print(_ret['stdout'])
    if (_ret["num_err"] > 0):
        for serr in _ret["errors"]:
            print("Errors:")
            print(serr)
    addlog("celslc", logfile=logfile, infos=_ret['stdout'])
    return _ret


def get_det_str(prm_det, version=2016021801):
    """
    Returns a STEM detector parameter string for the detector parameter file.
    """
    sprm = ""
    if version==2016021801:
        sprm += f"{prm_det['inner']}, {prm_det['outer']}"
        if 'sector' in prm_det:
            sprm += f", {prm_det['sector']['start']}, {prm_det['sector']['end']}"
        else:
            sprm +=  ", 0.0, 0.0"
        if 'center' in prm_det:
            sprm += f", {prm_det['center']['x']}, {prm_det['center']['y']}"
        else:
            sprm +=  ", 0.0, 0.0"
        sprm += f", '{prm_det['name']}'"
        if 'sensitivity' in prm_det:
            sprm += f", {prm_det['sensitivity']}"
        else:
            sprm +=  ", ''"
    return sprm


def get_periodic_readout_slices(num_slc, stepread):
    """
    Returns a list with the indices of readout slices for a periodic
    structure to identify output thicknesses.
    """
    if stepread == 0:
        return [num_slc]
    l = [0]
    for islc in range(0, num_slc):
        if 0 == (islc+1) % stepread:
            l.append(islc+1)
    if num_slc not in l:
        l.append(num_slc)
    return l


def write_det_prm(filename, prm, verbose=False):
    """
    Writes a detector parameter file for MSA.
    """
    _version = 2016021801
    l_stemdet = []
    if 'instrument' in prm:
        if 'detectors' in prm['instrument']:
            if 'diffraction' in prm['instrument']['detectors']:
                if 'STEM' in prm['instrument']['detectors']['diffraction']:
                    l_stemdet = prm['instrument']['detectors']['diffraction']['STEM']
    num_det = len(l_stemdet)
    with open(filename, 'w') as f: # open / overwrite the detector parameter file
        f.write( "'[Detector Parameters]'  ! parameter block name\n")
        f.write( "2016021801               ! detector parameter file version number\n")
        f.write(f"{num_det}                         ! number of detector definitions in this file\n")
        for det in l_stemdet:
            f.write(f"{get_det_str(det, version=_version)}\n")
        f.flush()
    if verbose:
        stmp = f"{num_det} detector parameter sets written to file [{filename}]."
        print(f"    {stmp}")


def write_msa_prm(filename, prm, verbose=False):
    """
    Writes an MSA parameter file.
    """
    _ret = check_prm_basics(prm)
    if _ret['num_err'] > 0:
        _ret['num_err'] += 1
        _ret['errors'].append("MSA parameter initialization failed.")
        return _ret
    sdir, sout = get_dirnam(prm)
    sdetfile = f"{sout}_detectors.prm"
    if len(sdir) > 0:
        sdetfile = os.path.join(sdir, f"{sout}_detectors.prm")
    write_det_prm(sdetfile, prm, verbose=verbose) # write detector parameter file
    use_fs = 0; fs = 1.0
    use_ss = 0; ss = 0.01
    if 'focus-spread' in prm['instrument']['optics']:
        use_fs = 1
        fs = prm['instrument']['optics']['focus-spread']
    if 'source-size' in prm['instrument']['optics']:
        use_ss = 1
        ss = prm['instrument']['optics']['source-size']
    l_ab = []
    is_stem = (prm['instrument']['optics']['mode'] == "STEM")
    prm_scan = {'size' : {'width' : 0.3905, 'height' : 0.3905},
                'offset' : {'x' : 0.0, 'y' : 0.0}, 'orientation' : 0.0,
                'sampling' : {'x' : 40, 'y' : 40}}
    if is_stem: # get STEM probe aberrations
        l_ab = prm['instrument']['optics'].get('aberrations', [])
        if 'scan' not in prm['instrument']['optics']:
            _ret['num_err'] += 1
            _ret['errors'].append("Missing instrument/optics/scan parameters for STEM simulation.")
            return _ret
        prm_scan = prm['instrument']['optics']['scan']
    tilt = [0.0, 0.0]
    if 'tilt' in prm['sample']: # get sample tilt in mrad
        tilt[0] = prm['sample']['tilt']['x']
        tilt[1] = prm['sample']['tilt']['y']
    # thickness setup (default, take one supercell)
    nzstep = 0
    nzsc = prm['calculation']['sampling']['z'] # number of slices per input cell
    nstack = nzsc
    lstack = np.linspace(0, nzsc-1, nzsc) # stack slices of one input cell
    l_readout = [nstack-1] # only readout final slice
    sc = read_cel(os.path.join(sdir,prm['sample']['structure']['file']))
    if 'box' not in sc:
        _ret['num_err'] += 1
        _ret['errors'].append("Missing sample/structure parameters, error on file reading .")
        return _ret
    sc_z = sc['box']['c']
    zmax = sc_z
    zstep = 0.0
    if 'thickness' in prm['sample']:
        zmax = prm['sample']['thickness'].get('maximum', sc_z)
        zstep = prm['sample']['thickness'].get('step', 0.0)
    nstack = int(np.round(zmax / sc_z * nzsc)) # max. number of slices
    nzstep = int(np.round(zstep / sc_z * nzsc)) # readout step in slices
    lstack = np.linspace(0, nstack-1, nstack, dtype=int) % nzsc
    l_readout = get_periodic_readout_slices(nstack, nzstep)
    prm['output']['readout-slices'] = l_readout
    nv = 1
    npass = 1
    prm_calc = prm['calculation']
    if 'phonon' in prm_calc:
        if "QEP" == prm_calc['phonon']['model'] or "FRFPMS" == prm_calc['phonon']['model']:
            nv = prm_calc['phonon'].get('configurations', 1)
            npass = prm_calc['phonon'].get('passes', 1)
    #
    with open(filename, 'w') as f: # write new MSA parameter file
        f.write( "'[Microscope Parameters]'  ! parameter block name\n")
        f.write(f"{prm['instrument']['beam']['convergence']:.4f}   ! beam convergence semi angle [mrad]\n")
        f.write( "0.0    ! inner detector angle [mrad], replaced by detector parameters\n")
        f.write( "30.0   ! outer detector angle [mrad], replaced by detector parameters\n")
        f.write(f"1, '{sdetfile}'    ! detector parameters (replacing previous two lines)\n")
        f.write(f"{prm['instrument']['beam']['energy']:.4f}    ! beam energy [kev]\n")
        f.write(f"{ss:.4f}    ! source size (HWHM) [nm]\n")
        f.write(f"{fs:.4f}    ! focus spread (1/e half width) [nm]\n")
        f.write( "2.0   ! focus spread kernel relative width\n")
        f.write( "7     ! focus spread kernel number of samples\n")
        f.write(f"{len(l_ab):d}     ! number of aberration definitions for STEM\n")
        for ab in l_ab:
            f.write(f"{ab[0]:d}, {ab[1]:.5E}, {ab[2]:.5E}\n")
        f.write( "'[Multislice Parameters]'  ! parameter block name\n")
        f.write(f"{tilt[0]*0.05729578:.6f}   ! sample tilt x [deg]\n")
        f.write(f"{tilt[1]*0.05729578:.6f}   ! sample tilt y [deg]\n")
        f.write(f"{prm_scan['offset']['x']:.6f}   ! scan offset x [nm] (STEM only)\n")
        f.write(f"{prm_scan['offset']['y']:.6f}   ! scan offset y [nm] (STEM only)\n")
        f.write(f"{prm_scan['size']['width']:.6f}   ! scan size x [nm] (STEM only)\n")
        f.write(f"{prm_scan['size']['height']:.6f}   ! scan size y [nm] (STEM only)\n")
        f.write(f"{prm_scan['orientation']:.6f}   ! scan orientation [deg] (STEM only)\n")
        f.write(f"{int(prm_scan['sampling']['x']):d}   ! scan sampling horizontal (STEM only)\n")
        f.write(f"{int(prm_scan['sampling']['y']):d}   ! scan sampling vertical (STEM only)\n")
        f.write(f"{use_fs:d}   ! flag to apply focus spread convolution (STEM only)\n")
        f.write(f"{use_ss:d}   ! flag to apply source size convolution (STEM only, with specific run)\n")
        f.write( "1   ! supercell repeat x (do not use)\n")
        f.write( "1   ! supercell repeat y (do not use)\n")
        f.write( "1   ! supercell repeat z (do not use)\n")
        f.write(f"{sout}   ! slice file name\n")
        f.write(f"{prm['calculation']['sampling']['z']:d}   ! number of slice file\n")
        f.write(f"{nv:d}   ! number of variants per slice\n")
        f.write(f"{npass:d}   ! number of passes per scan point\n")
        f.write(f"{nzstep:d}   ! readout thickness steps (slices)\n")
        f.write(f"{nstack:d}   ! number of object slices to max. thickness\n")
        for islc in lstack:
            f.write(f"{islc:d}\n")
    if verbose:
        print(f"MSA parameter file: {filename}")
    return _ret


def msa(prm, verbose=False, logfile=None):
    """
    Execution of a call to the program MSA.
    If given, the program is executed in the output
    directory prm['output']['directory'].

    Parameters:
        prm : dict
            Simulation parameters
        verbose : boolean, default False
            Text output flag
        logfile : str, default None
            Log file name
    
    Returns:
        dict
            number of errors, stdout of the call, error infos
    """
    # --------------------
    # init
    _ret = check_prm_basics(prm)
    if _ret['num_err'] > 0:
        _ret['num_err'] += 1
        _ret['errors'].append("msa initialization failed.")
    prm_calc = prm['calculation']
    prm_bin = prm_calc['program']
    swd, sout = get_dirnam(prm)
    # --------------------
    # write parameter files
    sprmfile = os.path.join(swd, sout + '_msa.prm')
    __ret = write_msa_prm(sprmfile, prm, verbose=verbose)
    if __ret['num_err'] > 0:
        _ret['num_err'] += __ret['num_err']
        for serr in __ret['errors']:
            _ret['errors'].append(serr)
        return _ret
    is_stem = (prm['instrument']['optics']['mode'] == "STEM")
    # --------------------
    # generate the call
    _call = "msa"
    if 'binary_path' in prm_bin:
        if prm_bin['binary_path'] is not None:
            _call = os.path.join(prm_bin['binary_path'], "msa")
    # --------------------
    # generate the arguments
    _args = ""
    # parameter file
    _args = f"-prm {sprmfile}"
    # output file
    _args += f" -out {sout}.bin"
    # 3d output data
    _args +=  " /3dout"
    # ctem
    if not is_stem:
        _args +=  " /ctem"
    # verbose
    if verbose:
        _args +=  " /verbose"
    # --------------------
    
    # --------------------
    # Run msa
    addlog("msa", logfile=logfile, infos=['arguments: ' + _args])
    __ret = run_job(_call, _args, cwd=swd, verbose=verbose)
    if __ret['num_err'] > 0:
        _ret['num_err'] += __ret['num_err']
        for serr in __ret['errors']:
            _ret['errors'].append(serr)
        return _ret
    for s in __ret['stdout']:
        _ret['stdout'].append(s)
    # --------------------
    # Handle output
    if verbose:
        print(_ret['stdout'])
    if (_ret["num_err"] > 0):
        for serr in _ret["errors"]:
            print("Errors:")
            print(serr)
    addlog("msa", logfile=logfile, infos=_ret['stdout'])
    return _ret


def run_simulation(prm, verbose=False, logfile=None, do_celslc=True, do_msa=True):
    """
    Runs a simulation defined by parameters prm using the
    Dr. Probe command line tools.
    (currently STEM only)
    """
    _ret = {'num_err' : 0, 'errors' : [], 'stdout' : []}
    if do_celslc:
        __ret = celslc(prm, verbose=verbose, logfile=logfile)
        _ret['num_err'] += __ret['num_err']
        for s in __ret['errors']:
            _ret['errors'].append(s)
        for s in __ret['stdout']:
            _ret['stdout'].append(s)
    if do_msa:
        __ret = msa(prm, verbose=verbose, logfile=logfile)
        _ret['num_err'] += __ret['num_err']
        for s in __ret['errors']:
            _ret['errors'].append(s)
        for s in __ret['stdout']:
            _ret['stdout'].append(s)
    return _ret


def load_scan_image(prm, detname):
    sdir, sout = get_dirnam(prm)
    prm_scan = prm['instrument']['optics']['scan']
    nx = prm_scan['sampling']['x']
    ny = prm_scan['sampling']['y']
    nd = len(prm['output']['readout-slices'])
    ndim = [nd, ny, nx]
    sfile = os.path.join(sdir, sout + f"_{detname}.bin")
    return np.fromfile(sfile, dtype=np.float32).reshape(ndim)