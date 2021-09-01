#----- IMPORTING MODULES -----#


from pathlib import Path
from subprocess import check_output

import PySimpleGUI as sg

from simplemcgui.secondary_page import *

# The theme of the GUI
sg.ChangeLookAndFeel('SandyBeach')

#----- HELPFUL FUNCTIONS -----#


def find_desktopPATH():
    """
    Finding the Desktop path, which will be used as a default value while setting the output
    directory

    Returns:
        desktopPATH [str]: THe path to the desktop directory
    """
    try:
        desktopPATH = Path.home() / 'Desktop'  # For Windows OS
    except:
        desktopPATH = check_output(
            ['xdg-user-dir', 'DESKTOP'])  # For Linux OS
    return desktopPATH


def find_datasets(values):
    """
    Finding the datasets for a given values in the main page

    Args:
        values [dict]: A dictionary created by the GUI that contains results of the choose values

    Returns:
        [str]: Datasets that is going to use in the program
    """
    avaliable_dataset = ['HD', 'BBAO', 'GBAO', 'GBAO_no6dF', 'CMASS',
                        'LBAO', 'LaBAO', 'LxBAO', 'MGS', 'Planck', 'WRd',
                        'PlDa', 'PlRdx10', 'WMAP', 'PlRd', 'CMBW', 'SN',
                        'SNx10', 'UnionSN', 'RiessH0', '6dFGS', 'dline'
                        ]
    data_range = list(range(4, 26))
    dataset_dict = dict(zip(data_range, avaliable_dataset))
    datasets = ''
    for i in data_range:
        if values[i] == True:
            datasets += dataset_dict[i] + '+'
    return datasets[:-1]

#----- IMPORTANT PARAMETERS -----#


# Simple Model and Simple Cosmological model has been removed
models = ['Anisotropic', 'Binned', 'CPantheon', 'DGP',
          'Decay', 'Decay01', 'Decay05', 'DecayFrac', 'EarlyDE', 'EarlyDE_rd_DE',
          'Grad_Ok', 'JordiCDM', 'LCDM', 'LCDMmasslessnu', 'NeffLCDM',
          'NumuLCDM', 'Phantom', 'PolyCDM', 'PolyOk', 'PolyOkc', 'PolyOkf',
          'Quintess', 'Quintom', 'Quintom_couple', 'Rotation', 'SlowRDE', 'Spline',
          'StepCDM', 'TLight', 'WeirdCDM', 'fPolyCDM', 'noradLCDM',
          'nuLCDM', 'nuoLCDM', 'nuwCDM', 'oLCDM', 'owCDM', 'owaCDM',
          'wCDM', 'waCDM']


sampler = ['mcmc', 'nested', 'emcee', 'maxlike', 'ga_deap']


menu_def = [
    ['&Help', ['&Info', ['&Models', '&Datasets', '&Samplers']]]
]


def initial_stg_gui():
    """
    The main page of the SimpleMC GUI
    """
    data_col1 = [
                    [sg.Checkbox('HD', font=('Verdana', 10))],
                    [sg.Checkbox('BBAO', default=True, font=('Verdana', 10))],
                    [sg.Checkbox('GBAO', font=('Verdana', 10))],
                    [sg.Checkbox('GBAO_no6dF', font=('Verdana', 10))],
                    [sg.Checkbox('CMASS', font=('Verdana', 10))]
                ]

    data_col2 = [
                    [sg.Checkbox('LBAO', font=('Verdana', 10))],
                    [sg.Checkbox('LaBAO', font=('Verdana', 10))],
                    [sg.Checkbox('LxBAO', font=('Verdana', 10))],
                    [sg.Checkbox('MGS', font=('Verdana', 10))],
                    [sg.Checkbox('Planck', default=True, font=('Verdana', 10))]
                ]

    data_col3 = [
                    [sg.Checkbox('WRd', font=('Verdana', 10))],
                    [sg.Checkbox('PlDa', font=('Verdana', 10))],
                    [sg.Checkbox('PlRdx10', font=('Verdana', 10))],
                    [sg.Checkbox('WMAP', font=('Verdana', 10))],
                    [sg.Checkbox('PlRd', font=('Verdana', 10))]
                ]

    data_col4 = [
                    [sg.Checkbox('CMBW', font=('Verdana', 10))],
                    [sg.Checkbox('SN', font=('Verdana', 10))],
                    [sg.Checkbox('SNx10', font=('Verdana', 10))],
                    [sg.Checkbox('UnionSN', default=True, font=('Verdana', 10))],
                    [sg.Checkbox('RiessH0', font=('Verdana', 10))]
                ]

    data_col5 = [
                    [sg.Checkbox('6dFGS', font=('Verdana', 10))],
                    [sg.Checkbox('dline', font=('Verdana', 10))]
                ]

    layout_initial_stg_page = [
                                [sg.Menu(menu_def)],

                                [sg.Text('SimpleMC', size=(44, 1), justification='center',
                                font=('Georgia', 16), relief=sg.RELIEF_SUNKEN)],

                                [sg.Frame(layout=[
                                    [sg.In(find_desktopPATH(), font=('Verdana', 10)),
                                    sg.FolderBrowse(initial_folder=r'simplemc/chains')]
                                ], title='Output Directory', font=('Georgia', 14), tooltip='Directory for chains/output')],

                                [sg.Frame(layout=[
                                    [sg.InputCombo(models, default_value='LCDM', font=('Verdana', 10))]
                                ], title='Models', font=('Georgia', 14), tooltip='Set Model'),
                                sg.Frame(layout=[
                                    [sg.InputCombo(sampler, default_value='mcmc', font=('Verdana', 10))]
                                ], title='Sampler', font=('Georgia', 14))],

                                [sg.Frame(layout=[
                                    [sg.Column(data_col1), sg.Column(data_col2), sg.Column(data_col3), sg.Column(data_col4), sg.Column(data_col5)]
                                ], title='Datasets', font=('Georgia', 14), tooltip='Set datasets used. Ex: UnionSN+BBAO+Planck')],

                                [sg.Button('Continue'), sg.Exit(button_color='red')]
                            ]
    window_initial_stg_page = sg.Window('SimpleMC', layout_initial_stg_page)
    while True:
        event, values = window_initial_stg_page.read()
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
        # obtaining values from GUI
        chainsdir = values[1]
        model = values[2]
        analyzername = values[3]
        datasets = find_datasets(values)
        # analyzing the events
        if event == 'Models':
            sg.popup('Information about the Cosmological Models can be put here')
        elif event == 'Datasets':
            sg.popup('Information about the Datasets can be put here')
        elif event == 'Samplers':
            sg.popup('Information about the Samplers can be put here')
        elif event == 'Continue':
            advanced_stg_gui(chainsdir, model, analyzername, datasets)


