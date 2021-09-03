import PySimpleGUI as sg


def base_parameter_layout():
    """
    Creating a GUI layout for the base parameters.
    """
    param_img = [
                    [sg.Text('Parameter', font=('Georgia', 11))],
                    [sg.Image(r'simplemcgui\parameter_images\Om_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\Obh2_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\h_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\mnu_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\Nnu_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\Ok_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\w_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\wa_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\wb_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\wc_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\s8_par.png', pad=((0, 0), (1, 8.5)))],
                    [sg.Image(r'simplemcgui\parameter_images\Pr_par.png', pad=((0, 0), (1, 8.5)))],
                ]
    param_values = [
                        [sg.Text('Value', font=('Georgia', 11))],
                        [sg.Input('0.3038', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.02234', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.6821', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.06', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('3.046', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.0', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('-1', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.0', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.7', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.7', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('0.8', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('28.6', size=(15, 1), font=('Verdana', 10))],
                    ]
    param_err = [
                    [sg.Text('Error', font=('Georgia', 11))],
                    [sg.Input('0.05', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.001', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.05', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.1', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.5', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.01', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.1', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.1', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.2', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.2', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('0.01', size=(15, 1), font=('Verdana', 10))],
                    [sg.Input('4', size=(15, 1), font=('Verdana', 10))],
                ]
    param_bounds = [
                        [sg.Text('Bounds', font=('Georgia', 11))],
                        [sg.Input('(0.1, 0.5)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(0.02, 0.025)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(0.4, 0.9)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(0, 1.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(3.046, 5.046)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(-0.02, 0.02)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(-2.0, 0.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(-2.0, 2.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(-2.0, 3.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(-3.0, 5.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(0.5, 1.0)', size=(15, 1), font=('Verdana', 10))],
                        [sg.Input('(5, 70)', size=(15, 1), font=('Verdana', 10))]
                    ]
    base_parameter_layout = [sg.Column(param_img, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_values, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_err, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_bounds, element_justification='c', vertical_alignment='c')],
    return base_parameter_layout
