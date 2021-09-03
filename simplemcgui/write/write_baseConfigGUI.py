def write_baseConfigGUI(chainsdir, model, analyzer, datasets, general_settings, analyzer_settings):
    """
    Writing the obtained values to the baseConfigGUI.ini

    Args:
        chainsdir [str]: The path of the output directory
        model [str]: Type of the cosmological model
        analyzer [str]: Type of the analyzer/sampler
        datasets [str]: Dataset/Datasets chosen by the user
        general_settings [list]: General settings obtained from the user
        analyzer_settings [list]: Analyzer settings obtained from the user
    """
    with open('baseConfigGUI.ini', 'w') as config_file:
        config_file.write('[custom]\n')
        config_file.write('chainsdir = {}\n'.format(chainsdir))
        config_file.write('model = {}\n'.format(model))
        if general_settings[0] == False:
            config_file.write('prefact = phy\n')
        else:
            config_file.write('prefact = pre\n')
        config_file.write('varys8 = {}\n'.format(general_settings[1]))
        config_file.write('datasets = {}\n'.format(datasets))
        config_file.write('analyzername = {}\n'.format(analyzer))
        config_file.write('addDerived = {}\n'.format(general_settings[3]))
        config_file.write('mcevidence = {}\n'.format(analyzer_settings[0]))
        config_file.write('mcevidence_k = {}\n'.format(analyzer_settings[1]))
        config_file.write('overwrite = {}\n'.format(general_settings[2]))
        config_file.write('getdist = {}\n'.format(analyzer_settings[2]))
        config_file.write('corner = {}\n'.format(analyzer_settings[3]))
        config_file.write('simpleplot = {}\n'.format(analyzer_settings[4]))
        config_file.write('showfig = {}\n'.format(analyzer_settings[5]))
        if analyzer == 'mcmc':
            config_file.write('[mcmc]\n')
            config_file.write('nsamp = {}\n'.format(analyzer_settings[6]))
            config_file.write('skip = {}\n'.format(analyzer_settings[7]))
            config_file.write('temp = {}\n'.format(analyzer_settings[8]))
            config_file.write('GRstop = {}\n'.format(analyzer_settings[9]))
            config_file.write('checkGR = {}\n'.format(analyzer_settings[10]))
            config_file.write('chainno = {}\n'.format(analyzer_settings[11]))
        elif analyzer == 'emcee':
            config_file.write('[emcee]\n')
            config_file.write('walkers = {}\n'.format(analyzer_settings[6]))
            config_file.write('nsamp = {}\n'.format(analyzer_settings[7]))
            config_file.write('burnin = {}\n'.format(analyzer_settings[8]))
            config_file.write('nproc = {}\n'.format(analyzer_settings[9]))
    config_file.close()