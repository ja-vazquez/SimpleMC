import PySimpleGUI as sg

from simplemcgui.general_settings import *
from simplemcgui.parameter_settings import *
from simplemcgui.sampler_settings import *
from simplemcgui.write_baseConfig import *
from simplemcgui.write_paramDefs import *


def advanced_stg_gui(chainsdir, model, analyzername, datasets):
    if model == 'LCDM':
        layout_advanced_stg_page = [
                                    [sg.Text('SimpleMC', size=(44, 1), justification='center',
                                    font=('Georgia', 16), relief=sg.RELIEF_SUNKEN)],

                                    [sg.Frame(layout=general_stg_layout(), title='General Settings', font=('Georgia', 14))],

                                    [sg.Frame(layout=[
                                        [sg.TabGroup([
                                            [sg.Tab('General Parameters', parameter_stg_layout(), font=('Georgia', 12)),
                                            sg.Tab('Sampler', sampler_stg_layout(analyzername), font=('Georgia', 12)),
                                            ]])
                                        ]
                                    ], title='Advanced Settings', font=('Georgia', 14))],

                                    [sg.Button('Run SimpleMC')]
                                    ]
        window_advanced_stg_page = sg.Window(
            'Advanced Settings', layout_advanced_stg_page)
        while True:
            event, values = window_advanced_stg_page.read()
            if event == sg.WIN_CLOSED:
                break
            elif event == 'Run SimpleMC':
                # reading values
                general_settings = values[0], values[2], values[4], values[5]
                param_value_settings = [values[i] for i in range(18, 30)]
                param_error_settings = [values[i] for i in range(30, 42)]
                param_bound_settings = [values[i] for i in range(42, 54)]
                sampler_settings = [values[i] for i in range(54, 66)]
                # writing values
                write_baseConfigGUI(chainsdir, model, analyzername, datasets, general_settings, sampler_settings)
                # write_paramDefs(param_value_settings, param_error_settings, param_bound_settings)
                sg.popup('baseConfigGUI.ini file has been created')
