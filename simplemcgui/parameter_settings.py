import PySimpleGUI as sg


# General Parameters
def parameter_stg_layout():
    param_img = [
                    [sg.Text('Parameter', font=('Georgia', 11))],
                    [sg.Image(r'simplemcgui\images\Om_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\Obh2_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\h_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\mnu_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\Nnu_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\Ok_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\w_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\wa_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\wb_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\wc_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\s8_par.png', pad=((0, 0), (2, 8)))],
                    [sg.Image(r'simplemcgui\images\Pr_par.png', pad=((0, 0), (2, 8)))],
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

    layout_parameter_stg = [sg.Column(param_img, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_values, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_err, element_justification='c', vertical_alignment='c'),
                            sg.Column(param_bounds, element_justification='c', vertical_alignment='c')],
    return layout_parameter_stg
