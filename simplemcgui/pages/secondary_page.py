import PySimpleGUI as sg
from simplemc.DriverMC import DriverMC

from simplemcgui.settings.analyzer import *
from simplemcgui.settings.base_parameter import *
from simplemcgui.settings.general import *
from simplemcgui.settings.model_parameters import *
from simplemcgui.write.write_baseConfigGUI import *
from simplemcgui.write.write_paramDefs import *


def secondary_page_layout(chainsdir, model, analyzer, datasets):
    """
    Creating the layout of the second page of the SimpleMC GUI.

    Args:
        chainsdir [str]: The path of the output directory
        model [str]: Type of the cosmological model
        analyzer [str]: Type of the analyzer/sampler
        datasets [str]: Dataset/Datasets chosen by the user
    """
    secondary_page_layout = [
                                [sg.Text('SimpleMC', size=(44, 1), justification='center',
                                font=('Georgia', 16), relief=sg.RELIEF_SUNKEN)],

                                [sg.Frame(layout=general_layout(), title='General Settings', font=('Georgia', 14))],

                                [sg.Frame(layout=[
                                    [sg.TabGroup([
                                        [sg.Tab('Base Parameters', base_parameter_layout(), font=('Georgia', 12)),
                                        sg.Tab('Analyzer', analyzer_layout(analyzer), font=('Georgia', 12)),
                                        sg.Tab('Model Parameters', model_parameter_layout(model), font=('Georgia', 12))
                                        ]])
                                    ]
                                ], title='Advanced Settings', font=('Georgia', 14))],

                                [sg.Button('Run SimpleMC')]
                                ]
    window_secondary_page = sg.Window('SimpleMC', secondary_page_layout)
    while True:
        event, values = window_secondary_page.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == 'Run SimpleMC':
            # reading values
            if model == 'LCDM':
                if analyzer == 'mcmc':
                    general_settings = [values[i] for i in range(4)]
                    base_param_value = [values[i] for i in range(16, 28)]
                    base_param_error = [values[i] for i in range(28, 40)]
                    base_param_bound = [values[i] for i in range(40, 52)]
                    analyzer_settings = [values[i] for i in range(52, 64)]
                if analyzer == 'emcee':
                    general_settings = [values[i] for i in range(4)]
                    base_param_value = [values[i] for i in range(16, 28)]
                    base_param_error = [values[i] for i in range(28, 40)]
                    base_param_bound = [values[i] for i in range(40, 52)]
                    analyzer_settings = [values[i] for i in range(52, 62)]
                # writing to baseconfigGUI file
                write_baseConfigGUI(chainsdir, model, analyzer, datasets, general_settings, analyzer_settings)
                #write_base_paramDefs(base_param_value, base_param_error, base_param_bound)

            elif model == 'Grad_Ok':
                if analyzer == 'mcmc':
                    # reading values
                    general_settings = [values[i] for i in range(4)]
                    base_param_value = [values[i] for i in range(16, 28)]
                    base_param_error = [values[i] for i in range(28, 40)]
                    base_param_bound = [values[i] for i in range(40, 52)]
                    analyzer_settings = [values[i] for i in range(52, 64)]
                    model_param_value = [values[i] for i in range(66, 68)]
                    model_param_error = [values[i] for i in range(68, 70)]
                    model_param_bound = [values[i] for i in range(70, 72)]
                if analyzer == 'emcee':
                    general_settings = [values[i] for i in range(4)]
                    base_param_value = [values[i] for i in range(16, 28)]
                    base_param_error = [values[i] for i in range(28, 40)]
                    base_param_bound = [values[i] for i in range(40, 52)]
                    analyzer_settings = [values[i] for i in range(52, 62)]
                    model_param_value = [values[i] for i in range(64, 66)]
                    model_param_error = [values[i] for i in range(66, 68)]
                    model_param_bound = [values[i] for i in range(68, 70)]
                # writing to baseconfigGUI file
                write_baseConfigGUI(chainsdir, model, analyzer, datasets, general_settings, analyzer_settings)
                # writing to paramDefs.py file
                #write_base_paramDefs(base_param_value, base_param_error, base_param_bound)
                #write_model_paramDefs(model, model_param_value, model_param_error, model_param_bound)
            window_secondary_page.close()
            # running the Simple MC
            analyzer = DriverMC(iniFile='baseConfigGUI.ini')
            analyzer.executer()

