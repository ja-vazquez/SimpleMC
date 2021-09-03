import PySimpleGUI as sg


# Base Parameters
def model_parameter_layout(model):
    """
    Creating a GUI layout for the model parameters.
    """
    if model == 'LCDM':
        return [[sg.Text('No model specific parameters. Please adjust the Base Parameters')]]
    elif model == 'Grad_Ok':
        param_img = [
                        [sg.Text('Parameter', font=('Georgia', 11))],
                        [sg.Image(r'simplemcgui\parameter_images\ggama_par.png', pad=((0, 0), (1, 8.5)))],
                        [sg.Image(r'simplemcgui\parameter_images\glambda_par.png', pad=((0, 0), (1, 8.5)))],
                    ]
        param_values = [
                            [sg.Text('Value', font=('Georgia', 11))],
                            [sg.Input('-1.0', size=(15, 1), font=('Verdana', 10))],
                            [sg.Input('0', size=(15, 1), font=('Verdana', 10))],
                        ]
        param_err = [
                        [sg.Text('Error', font=('Georgia', 11))],
                        [sg.Input('0.1', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.2', size=(15, 1), font=('Verdana', 10))],
                    ]
        param_bounds = [
                            [sg.Text('Bounds', font=('Georgia', 11))],
                            [sg.Input('(-1.3, -0.7)', size=(15, 1), font=('Verdana', 10))],
                            [sg.Input('(-10, 0)', size=(15, 1), font=('Verdana', 10))],
                        ]
        model_parameter_layout = [sg.Column(param_img, element_justification='c', vertical_alignment='c'),
                                sg.Column(param_values, element_justification='c', vertical_alignment='c'),
                                sg.Column(param_err, element_justification='c', vertical_alignment='c'),
                                sg.Column(param_bounds, element_justification='c', vertical_alignment='c')],
        return model_parameter_layout
