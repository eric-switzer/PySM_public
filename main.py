"""Commandline utility for python sky model
"""
import os
import argparse
import ConfigParser
import pysm_synchrotron, pysm_thermaldust, pysm_cmb, pysm_spinningdust, pysm_noise, pysm_freefree
from pysm import Output, config2list, file_path, write_output_single
import healpy as hp
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Code to simulate galactic foregrounds.')
    parser.add_argument('config_file', help='Main configuration file.')

    # Get the output directory and save the configuration file.
    Config = ConfigParser.ConfigParser()
    Config.read(parser.parse_args().config_file)
    out = Output(Config._sections['GlobalParameters'])

    if not os.path.exists(out.output_dir):
        os.makedirs(out.output_dir)

    cfg_fname = out.output_dir + out.output_prefix + 'main_config.ini'
    with open(cfg_fname, 'w') as configfile:
        Config.write(configfile)

    if out.debug is True:
        print '----------------------------------------------------- \n'
        print ''.join("%s: %s \n" % item for item in vars(out).items())
        print '-----------------------------------------------------'

    model = np.zeros([3, len(out.output_frequency), hp.nside2npix(out.nside)])
    print '----------------------------------------------------- \n'
    if 'synchrotron' in out.components:
        model += pysm_synchrotron.main(parser.parse_args().config_file)

    if 'thermaldust' in out.components:
        model += pysm_thermaldust.main(parser.parse_args().config_file)

    if 'spinningdust' in out.components:
        model += pysm_spinningdust.main(parser.parse_args().config_file)

    if 'freefree' in out.components:
        model += pysm_freefree.main(parser.parse_args().config_file)

    if 'cmb' in out.components:
        model += pysm_cmb.main(parser.parse_args().config_file)

    if out.smoothing:
        print 'Smoothing output maps.'
        print '----------------------------------------------------- \n'
        for i in xrange(len(out.output_frequency)):
            fwhm = (np.pi / 180. / 60.) * out.fwhm[i]

            model[:, i, :] = hp.smoothing(model[:, i, :],
                                          fwhm=fwhm, verbose=False)

    if out.instrument_noise is True:
        model += pysm_noise.instrument_noise(parser.parse_args().config_file)

    if out.instrument_noise:
        out.components.append('noise')

    model = np.swapaxes(model, 0, 1)
    for i in xrange(len(out.output_frequency)):
        write_output_single(model[i, ...], out, Config, i)

    print '-----------------------------------------------------\n'
    print 'PySM completed successfully. \n'
    print '-----------------------------------------------------'


if __name__ == '__main__':
    main()
