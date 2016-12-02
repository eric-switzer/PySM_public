import ConfigParser
import healpy as hp
import numpy as np
import pysm as sm

def main(fname_config):
    config = ConfigParser.ConfigParser()
    config.read(fname_config)
    out = sm.output(config._sections['GlobalParameters'])

    config.read('./ConfigFiles/' + config.get('SpinningDust', 'model') + '_config.ini')
    spdust_general = sm.Component(config._sections['General'], out.nside)
    spdust1 = sm.Component(config._sections['SpinningDust1'], out.nside)
    spdust2 = sm.Component(config._sections['SpinningDust2'], out.nside)

    print('Computing spinning dust map.')
    print '----------------------------------------------------- \n'

    if out.debug is True:
        print ''.join("%s: %s \n" % item for item in vars(spdust1).items())
        print ''.join("%s: %s \n" % item for item in vars(spdust2).items())
        print '----------------------------------------------------- \n'

    cfg_fname = out.output_dir + out.output_prefix + 'spdust_config.ini'
    with open(cfg_fname, 'w') as configfile:
        config.write(configfile)

    # Compute a map of the polarisation angle
    # from the commander dust map polariationn angle.
    pol_angle = np.arctan2(spdust_general.thermaldust_polu,
                           spdust_general.thermaldust_polq)

    # Units to do the scaling in MJysr
    # then bring the result back to the output units.
    conv1 = sm.convert_units(spdust1.template_units,
                             ['u', 'K_RJ'], spdust1.freq_ref)

    conv2 = sm.convert_units(spdust2.template_units,
                             ['u', 'K_RJ'], spdust2.freq_ref)

    conv_end = sm.convert_units(['u', 'K_RJ'], out.output_units,
                                out.output_frequency)

    unit_conversion1 = conv1*conv_end.reshape((len(out.output_frequency), 1))
    unit_conversion2 = conv2*conv_end.reshape((len(out.output_frequency), 1))

    scaled_map_spdust = sm.scale_freqs(spdust1, out, pol=False)*spdust1.em_template*unit_conversion1 + \
                        sm.scale_freqs(spdust2, out, pol=False)*spdust2.em_template*unit_conversion2

    scaled_map_spdust_pol = scaled_map_spdust[np.newaxis, ...]*np.asarray([np.cos(pol_angle), np.sin(pol_angle)])[:, np.newaxis, :]*spdust_general.pol_frac

    if out.debug == True:
        for i in range(0, len(out.output_frequency)):
            hp.write_map(out.output_dir + 'spdust_%d.fits'%(out.output_frequency[i]), scaled_map_spdust[i], coord='G', column_units=out.output_units)

    return np.concatenate([scaled_map_spdust[np.newaxis, ...], scaled_map_spdust_pol])
