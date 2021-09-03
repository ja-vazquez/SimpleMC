import PySimpleGUI as sg


def general_layout():
    """
    Creating a GUI layout for the general settings.
    """
    col1 = [
                [sg.Checkbox('Enable Prefactor', font=('Verdana', 10),
                tooltip='Prefactor parameter: c/(r_dH_0)')],
                [sg.Checkbox('Vary Sigma 8', font=('Verdana', 10))]
            ]
    col2 = [
                [sg.Checkbox('Overwrite', default=True, font=('Verdana', 10),
                tooltip='overwrite = True -> overwrite output files with the same name\noverwrite = False -> if the outputname already exist. It sends an error and ends the simplemc execution')],
                [sg.Checkbox('Add derived parameters', font=('Verdana', 10),
                tooltip='i.e. Omega_Lambda, H0, Age of the Universe')]
            ]
    general_layout = [
                        [sg.Column(col1, justification='c', pad=((0,60),(0,0))),
                        sg.Column(col2, justification='c', pad=((60,0),(0,0)))]
                    ]
    return general_layout
