import numpy as np
import healpy as hp
import ConfigParser
import pysm as sm

def scale_dust_pop(pop, out, config):
    dust = sm.Component(config._sections[pop], out.nside)
    print('Computing dust maps.')
    print '----------------------------------------------------- \n'
    if out.debug is True:
        print ''.join("%s: %s \n" % item   for item in vars(dust).items())
        print '----------------------------------------------------- \n'

    conv_i = sm.convert_units(dust.template_units, ['u', 'K_RJ'], dust.freq_ref)

    conv_pol = sm.convert_units(dust.template_units, ['u', 'K_RJ'],
                             dust.pol_freq_ref)

    conv2 = sm.convert_units(['u', 'K_RJ'], out.output_units,
                          out.output_frequency)

    unit_conversion_I = conv_i * conv2.reshape((len(out.output_frequency), 1))
    unit_conversion_pol = conv_pol * conv2.reshape((len(out.output_frequency), 1))

    scaled_map_dust = sm.scale_freqs(dust, out, pol=False) * dust.em_template * \
                      unit_conversion_I

    scaled_map_dust_pol = sm.scale_freqs(dust, out, pol=True)[np.newaxis, ...] * \
                          np.array([dust.polq_em_template, dust.polu_em_template])[:, np.newaxis, :]*unit_conversion_pol

    if out.debug is True:
        dus = np.concatenate([scaled_map_dust[np.newaxis, ...],
                              scaled_map_dust_pol])

        for i in range(0, len(out.output_frequency)):
            fname = out.output_dir + out.output_prefix
            fname += 'dust_%d' % (out.output_frequency[i])
            fname += '_' + str(out.nside) + '.fits'

            hp.write_map(fname, dus[:, i, :], coord='G',
                         column_units=out.output_units)

    return np.concatenate([scaled_map_dust[np.newaxis, ...],
                           scaled_map_dust_pol])


def main(fname_config):
    """Read configuration into classes"""
    config = ConfigParser.ConfigParser()
    config_model = ConfigParser.ConfigParser()

    config.read(fname_config)
    out = sm.Output(config._sections['GlobalParameters'])

    config_model.read('./ConfigFiles/' + config.get('ThermalDust', 'model') + '_config.ini')
    pops = config_model.sections()

    cfg_fname = out.output_dir + out.output_prefix + 'thermaldust_config.ini'
    with open(cfg_fname, 'w') as configfile:
        config_model.write(configfile)

    dust_out = 0.

    for p in pops:
        dust_out += scale_dust_pop(p, out, config_model)

    return dust_out
