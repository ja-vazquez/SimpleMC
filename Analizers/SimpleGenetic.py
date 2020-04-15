from datetime import datetime
from Individual import Individual
from Population import Population

import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import random
import seaborn as sb
import time
import warnings

class SimpleGenetic():
    def __init__(self, funcion_objetivo, n_variables, limites, n_individuos = 50,\
                verbose = False, optimizacion = "maximizar",\
                n_generaciones = 250, metodo_seleccion = "tournament", elitismo = 0.01,\
                prob_mut = 0.1, distribucion = "uniforme", media_distribucion = 1,\
                sd_distribucion = 1, min_distribucion = -1, max_distribucion = 1,\
                parada_temprana = False, rondas_parada = 5, tolerancia_parada = 0.1):

        self.funcion_objetivo = funcion_objetivo
        #Estos limites son una lista donde cada entrada son los limites de un param
        limites = limites

        self.limites_inf = []
        self.limites_sup = []

        for limite in limites:
            self.limites_inf.append(limite[0])
            self.limites_sup.append(limite[1])

        self.n_individuos = n_individuos
        self.n_variables = n_variables
        
        self.distribucion = distribucion
        self.elitismo = elitismo
        self.max_distribucion   = max_distribucion
        self.media_distribucion = media_distribucion
        self.metodo_seleccion = metodo_seleccion
        self.min_distribucion = min_distribucion
        self.n_generaciones = n_generaciones
        self.optimizacion = optimizacion
        self.parada_temprana = parada_temprana
        self.prob_mut = prob_mut
        self.rondas_parada = rondas_parada
        self.sd_distribucion = sd_distribucion
        self.tolerancia_parada  = tolerancia_parada
        self.verbose = verbose

        self.optimizar()

    def optimizar(self):
        poblacion = Population(n_individuos = self.n_individuos,\
                n_variables = self.n_variables,\
                limites_inf = self.limites_inf,\
                limites_sup = self.limites_sup,\
                verbose = self.verbose)
    
        o = poblacion.optimizar(funcion_objetivo = self.funcion_objetivo,\
                            optimizacion = self.optimizacion,\
                            n_generaciones = self.n_generaciones,\
                            metodo_seleccion = self.metodo_seleccion,\
                            elitismo = self.elitismo,\
                            prob_mut = self.prob_mut,\
                            distribucion = self.distribucion,\
                            media_distribucion = self.media_distribucion,\
                            sd_distribucion = self.sd_distribucion ,\
                            min_distribucion = self.min_distribucion,\
                            max_distribucion = self.max_distribucion,\
                            parada_temprana = self.parada_temprana,\
                            rondas_parada = self.rondas_parada,\
                            verbose = self.verbose)
        return o