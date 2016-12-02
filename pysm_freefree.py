import ConfigParser
import healpy as hp
import numpy as np
import pysm as sm

def main(fname_config):
        config = ConfigParser.ConfigParser()
        config.read(fname_config)
        out = sm.output(config._sections['GlobalParameters'])

        config.read('./ConfigFiles/' + config.get('FreeFree', 'model') + '_config.ini')
        freefree = sm.Component(config._sections['FreeFree'], out.nside)
        with open(out.output_dir + out.output_prefix + 'freefree_config.ini', 'w') as configfile:
            config.write(configfile)

        print('Computing free-free maps.')
        print '----------------------------------------------------- \n'
        if out.debug is True:
                print ''.join("%s: %s \n" % item   for item in vars(freefree).items())
                print '----------------------------------------------------- \n'

        conv_I = sm.convert_units(freefree.template_units, out.output_units, out.output_frequency)

        scaled_map_ff = sm.scale_freqs(freefree, out)*conv_I[..., np.newaxis]*freefree.em_template
        scaled_map_ff_pol = np.zeros((2, np.asarray(out.output_frequency).size, hp.nside2npix(out.nside)))

        if out.debug is True:
            ff = np.concatenate([scaled_map_ff[np.newaxis, ...], scaled_map_ff_pol])
            for i in range(0, len(out.output_frequency)):
                hp.write_map(out.output_dir + out.output_prefix + 'ff_%d'%(out.output_frequency[i]) + '_' + str(out.nside) + '.fits', ff[:, i, :], coord='G', column_units=out.output_units)

        return np.concatenate([scaled_map_ff[np.newaxis, ...], scaled_map_ff_pol])
