import PySimpleGUI as sg


def sampler_stg_layout(analyzer_name):
    if analyzer_name == 'mcmc':
        layout_sampler_stg = [
                                    [sg.Text('nsamp:', font=('Verdana', 10), tooltip='Number of samples'),
                                    sg.Input(default_text='500')],

                                    [sg.Text('skip:', font=('Verdana', 10), tooltip='Burn-in'),
                                    sg.Input(default_text='0')],

                                    [sg.Text('temp:', font=('Verdana', 10), tooltip='Temperature at which to sample'),
                                    sg.Input(default_text='2')],

                                    [sg.Text('GRstop:', font=('Verdana', 10), tooltip='Gelman-Rubin for convergence'),
                                    sg.Input(default_text='0.01')],

                                    [sg.Text('checkGR:', font=('Verdana', 10), tooltip='Every number of steps check the GR-criteria'),
                                    sg.Input(default_text='50')],

                                    [sg.Text('chainno:', font=('Verdana', 10), tooltip='1 if single cpu , otherwise is giving by the nproc-> mpi -np #'),
                                    sg.Input(default_text='0')],

                                    [sg.Text('mcevidence:', font=('Verdana', 10), tooltip='True to calculate Bayesian evidence with mcevidence'),
                                    sg.Input(default_text='True')],

                                    [sg.Text('mcevidence_k:', font=('Verdana', 10), tooltip='mcevidence_k is the k number of nearest neighbours in mcevidence'),
                                    sg.Input(default_text='2')],

                                    [sg.Frame(layout=[
                                        [sg.Text('getdist:', font=('Verdana', 10)),
                                        sg.Input(default_text='False')],

                                        [sg.Text('corner:', font=('Verdana', 10)),
                                        sg.Input(default_text='False')],

                                        [sg.Text('simpleplot:', font=('Verdana', 10)),
                                        sg.Input(default_text='False')],

                                        [sg.Text('showfig:', font=('Verdana', 10), tooltip='True to display figures; we recommended false'),
                                        sg.Input(default_text='False')]
                                    ], title='Plot Setting', font=('Georgia', 14), tooltip='Options to triangle plots for MCMC, .png files will be saved in chainsdir')]
                                ]
    return layout_sampler_stg
