from pathlib import Path
from subprocess import check_output

import PySimpleGUI as sg
from simplemcgui.pages.secondary_page import *

# Theme of the GUI
sg.ChangeLookAndFeel('SandyBeach')

#----- HELPFUL FUNCTIONS -----#


def find_desktopPATH():
    """
    Finding the desktop path, which will be used as a default value while choosing the output directory.

    Returns:
        desktopPATH [str]: Path of the desktop directory
    """
    try:
        desktopPATH = Path.home() / 'Desktop'  # For Windows OS
    except:
        desktopPATH = check_output(['xdg-user-dir', 'DESKTOP'])  # For Linux OS
    return desktopPATH


def find_datasets(values):
    """
    Exctracting the datasets from the given values, obtained from the GUI.

    Args:
        values [dict]: A dictionary created by the GUI that contains results of the selected parameters.

    Returns:
        [str]: Datasets that is going to be used in the program.
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


analyzers = ['mcmc', 'nested', 'emcee', 'maxlike', 'ga_deap']


menu_def = [
    ['&Help', ['&Info', ['&Models', '&Datasets', '&Samplers']]]
]


def main_page_layout():
    """
    Creating the layout of the first page of the SimpleMC GUI.
    """
    dataset_col1 = [
                    [sg.Checkbox('HD', font=('Verdana', 10))],
                    [sg.Checkbox('BBAO', default=True, font=('Verdana', 10))],
                    [sg.Checkbox('GBAO', font=('Verdana', 10))],
                    [sg.Checkbox('GBAO_no6dF', font=('Verdana', 10))],
                    [sg.Checkbox('CMASS', font=('Verdana', 10))]
                ]
    dataset_col2 = [
                    [sg.Checkbox('LBAO', font=('Verdana', 10))],
                    [sg.Checkbox('LaBAO', font=('Verdana', 10))],
                    [sg.Checkbox('LxBAO', font=('Verdana', 10))],
                    [sg.Checkbox('MGS', font=('Verdana', 10))],
                    [sg.Checkbox('Planck', default=True, font=('Verdana', 10))]
                ]
    dataset_col3 = [
                    [sg.Checkbox('WRd', font=('Verdana', 10))],
                    [sg.Checkbox('PlDa', font=('Verdana', 10))],
                    [sg.Checkbox('PlRdx10', font=('Verdana', 10))],
                    [sg.Checkbox('WMAP', font=('Verdana', 10))],
                    [sg.Checkbox('PlRd', font=('Verdana', 10))]
                ]
    dataset_col4 = [
                    [sg.Checkbox('CMBW', font=('Verdana', 10))],
                    [sg.Checkbox('SN', font=('Verdana', 10))],
                    [sg.Checkbox('SNx10', font=('Verdana', 10))],
                    [sg.Checkbox('UnionSN', default=True, font=('Verdana', 10))],
                    [sg.Checkbox('RiessH0', font=('Verdana', 10))]
                ]
    dataset_col5 = [
                    [sg.Checkbox('6dFGS', font=('Verdana', 10))],
                    [sg.Checkbox('dline', font=('Verdana', 10))]
                ]
    main_page_layout = [
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
                                [sg.InputCombo(analyzers, default_value='mcmc', font=('Verdana', 10))]
                            ], title='Analyzers', font=('Georgia', 14))],

                            [sg.Frame(layout=[
                                [sg.Column(dataset_col1),
                                sg.Column(dataset_col2),
                                sg.Column(dataset_col3),
                                sg.Column(dataset_col4),
                                sg.Column(dataset_col5)]
                            ], title='Datasets', font=('Georgia', 14), tooltip='Set datasets used. Ex: UnionSN+BBAO+Planck')],

                            [sg.Button('Continue'), sg.Exit(button_color='red')]
                            ]
    window_main_page = sg.Window('SimpleMC', main_page_layout)
    while True:
        event, values = window_main_page.read()
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
        # obtaining values from GUI
        chainsdir = values[1]
        model = values[2]
        analyzer = values[3]
        datasets = find_datasets(values)
        # analyzing the events
        if event == 'Models':
            sg.popup('Information about the Cosmological Models can be put here')
        elif event == 'Datasets':
            sg.popup('Information about the Datasets can be put here')
        elif event == 'Analyzers':
            sg.popup('Information about the Analyzers can be put here')
        elif event == 'Continue':
            window_main_page.close()
            secondary_page_layout(chainsdir, model, analyzer, datasets)
