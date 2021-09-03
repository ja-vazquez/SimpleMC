import PySimpleGUI as sg


def analyzer_layout(analyzer):
    """
    Creating a GUI layout for a given analyzer

    Args:
        analyzer [str]: Type of the analyzer, given as a string.
    """
    sampler_general_layout = [
                                [sg.Text('mcevidence:', font=('Verdana', 10),
                                tooltip='True to calculate Bayesian evidence with mcevidence'),
                                sg.Input(default_text='True')],

                                [sg.Text('mcevidence_k:', font=('Verdana', 10),
                                tooltip='mcevidence_k is the k number of nearest neighbours in mcevidence'),
                                sg.Input(default_text='2')],

                                [sg.Text('getdist:', font=('Verdana', 10)),
                                sg.Input(default_text='False')],

                                [sg.Text('corner:', font=('Verdana', 10)),
                                sg.Input(default_text='False')],

                                [sg.Text('simpleplot:', font=('Verdana', 10)),
                                sg.Input(default_text='False')],

                                [sg.Text('showfig:', font=('Verdana', 10),
                                tooltip='True to display figures; we recommended false'),
                                sg.Input(default_text='False')],
                            ]
    if analyzer in ['mcmc', 'emcee', 'nested']:
        if analyzer == 'mcmc':
            analyzer_layout = [
                                *sampler_general_layout,

                                [sg.Text('nsamp:', font=('Verdana', 10),
                                tooltip='Number of samples'),
                                sg.Input(default_text='500')],

                                [sg.Text('skip:', font=('Verdana', 10),
                                tooltip='Burn-in'),
                                sg.Input(default_text='0')],

                                [sg.Text('temp:', font=('Verdana', 10),
                                tooltip='Temperature at which to sample'),
                                sg.Input(default_text='2')],

                                [sg.Text('GRstop:', font=('Verdana', 10),
                                tooltip='Gelman-Rubin for convergence'),
                                sg.Input(default_text='0.01')],

                                [sg.Text('checkGR:', font=('Verdana', 10),
                                tooltip='Every number of steps check the GR-criteria'),
                                sg.Input(default_text='50')],

                                [sg.Text('chainno:', font=('Verdana', 10),
                                tooltip='1 if single cpu , otherwise is giving by the nproc-> mpi -np #'),
                                sg.Input(default_text='0')]
                            ]
        elif analyzer == 'emcee':
            analyzer_layout = [
                                *sampler_general_layout,

                                [sg.Text('walkers:', font=('Verdana', 10),
                                tooltip='walkers >= 2*dim'),
                                sg.Input(default_text='10')],

                                [sg.Text('nsamp:', font=('Verdana', 10),
                                tooltip='Number of samples'),
                                sg.Input(default_text='200')],

                                [sg.Text('burnin:', font=('Verdana', 10),
                                tooltip='Burn-in'),
                                sg.Input(default_text='0')],

                                [sg.Text('nproc:', font=('Verdana', 10),
                                tooltip='Number of procesors'),
                                sg.Input(default_text='4')]
                            ]
        return analyzer_layout
