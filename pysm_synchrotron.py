import numpy as np
import healpy as hp
import pysm as sm
import ConfigParser

def main(fname_config):
    config = ConfigParser.ConfigParser()
    config.read(fname_config)
    out = sm.Output(config._sections['GlobalParameters'])

    cfg_fname = './ConfigFiles/'
    cfg_fname += config.get('Synchrotron', 'model')+'_config.ini'
    config.read(cfg_fname)
    synch = sm.Component(config._sections['Synchrotron'], out.nside)

    cfg_fname = out.output_dir+out.output_prefix+'synchrotron_config.ini'

    with open(cfg_fname, 'w') as configfile:
        config.write(configfile)

    print('Computing synchrotron maps.')
    print '----------------------------------------------------- \n'
    if out.debug is True:
        print ''.join("%s: %s \n" % item for item in vars(synch).items())
        print '----------------------------------------------------- \n'

    # The unit conversion takes care of the scaling being done in uK_RJ.
    # After scaling we convert to whatever the output units are.
    conv_I = sm.convert_units(synch.template_units,
                              ['u', 'K_RJ'], synch.freq_ref)

    conv_pol = sm.convert_units(synch.template_units,
                                ['u', 'K_RJ'], synch.pol_freq_ref)

    conv2 = sm.convert_units(['u', 'K_RJ'], out.output_units,
                             out.output_frequency)

    unit_conversion_i = conv_I * conv2.reshape((len(out.output_frequency), 1))
    unit_conversion_pol = conv_pol * conv2.reshape((len(out.output_frequency), 1))

    # Do the scaling.
    scaled_map_synch = sm.scale_freqs(synch, out, pol=False) * \
                       synch.em_template * unit_conversion_i

    scaled_map_synch_pol = sm.scale_freqs(synch, out, pol=True)[np.newaxis, ...] * \
                           np.array([synch.polq_em_template, synch.polu_em_template])[:, np.newaxis, :] * unit_conversion_pol

    # This section forces P/I<0.75.
    # Same PSM 1.7.8 psm_synchrotron.pro.
    P = np.sqrt(scaled_map_synch_pol[0, :, :] ** 2 + scaled_map_synch_pol[1, :, :] ** 2) / scaled_map_synch
    F = 0.75 * np.tanh(P / 0.75) / P
    scaled_map_synch_pol[0, :, :] = F * scaled_map_synch_pol[0, :, :]
    scaled_map_synch_pol[1, :, :] = F * scaled_map_synch_pol[1, :, :]

    if out.debug is True:
        syn = np.concatenate([scaled_map_synch[np.newaxis, ...],
                              scaled_map_synch_pol])

        for i in range(0, len(out.output_frequency)):
            out_fname = out.output_dir + out.output_prefix
            out_fname += 'synch_%d' % (out.output_frequency[i])
            out_fname += '_'+str(out.nside)+'.fits'

            hp.write_map(out_fname, syn[:, i, :], coord='G',
                         column_units=out.output_units)

    return np.concatenate([scaled_map_synch[np.newaxis, ...],
                           scaled_map_synch_pol])
