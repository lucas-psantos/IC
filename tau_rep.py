# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
from glob import glob
import pyhdust as phd
import shutil

 
def calcula_Tau_rep(local, iobs, par_modelos, nmodelos, interpol):

    '''
    Cria uma tabela com os valores de Tau_rep na pasta do projeto
    local -- local onde os modelos estão presentes
    iobs -- qual observador será analisado
    '''

    import pyhdust as phd
    from scipy.integrate import simps
    
    tab = []

    for i in range(len(par_modelos)//nmodelos):
        # Primeiro, preciso que os valores sejam strings
        MOD = i + interpol
        
        if (MOD<10):
            MOD = '0' + str(MOD)
        else:
            MOD = str(MOD)

        amin = str(par_modelos[i*nmodelos,2])
        amax = str(par_modelos[i*nmodelos,3])
        q = str(par_modelos[i*nmodelos,4])

        for j in range(nmodelos):   
            tauv = str(round(par_modelos[i*nmodelos+j,0].astype(float),2))
            rint = str(round(par_modelos[i*nmodelos+j,1].astype(float),2))
            #print(local + f"mod{MOD}/mod{MOD}_Tauv{tauv}_Rint{rint}_amin{amin}_amax{amax}_q{q}*.dust")
            arq = local + f"fullsed/fullsed_mod{MOD}/mod{MOD}_Tauv{tauv}_Rint{rint}_amin{amin}_amax{amax}_q{q}.sed2")) # pega o arquivo sed2
            model = phd.readfullsed2(arq)  # le o arquivo
            tab.append(model)
            print('arquivo lido!')

    tab = np.array(tab)
    
    lbd = tab[:, 0, :, 2]
    flux_all = tab[:, iobs, :, 3]
    #flux_sct = tab[:, iobs, :, 4]
    flux_emit = tab[:, iobs, :, 5]
    #flux_trans = tab[:, iobs, :, 6]
   
   
    #integra cada flux e calcula o tau_rep
    Flux_emit = simps(flux_emit,lbd)
    Flux_all = simps(flux_all,lbd)
    Tau = -np.log(1. - Flux_emit/Flux_all)
                        
    return np.around(Tau ,2) # arredonda para 2 casas decimais

    
def calcula_Tint(local, par_modelos, nmodelos, interpol):
    '''
    Cria uma tabela com os valores de Tint.
    local -- local onde estão presentes os modelos
    '''
        
    ####------------------------------------------------------------####
    ####                    Roda os todos os modelos                ####
    ####------------------------------------------------------------####
    
    Tint = []

    for i in range(len(par_modelos)//nmodelos):
        # Primeiro, preciso que os valores sejam strings
        MOD = i + interpol
        
        if (MOD<10):
            MOD = '0' + str(MOD)
        else:
            MOD = str(MOD)

        amin = str(par_modelos[i*nmodelos,2])
        amax = str(par_modelos[i*nmodelos,3])
        q = str(par_modelos[i*nmodelos,4])

        for j in range(nmodelos):   
            tauv = str(round(par_modelos[i*nmodelos+j,0].astype(float),2))
            rint = str(round(par_modelos[i*nmodelos+j,1].astype(float),2))
            #print(local + f"mod{MOD}/mod{MOD}_Tauv{tauv}_Rint{rint}_amin{amin}_amax{amax}_q{q}*.dust")
            arq = sorted(glob(local + f"mod{MOD}/mod{MOD}_Tauv{tauv}_Rint{rint}_amin{amin}_amax{amax}_q{q}*.dust"))[-1] # pega o último dust da pasta mod
            arq_dust = phd.readdust(arq)
            t_dust = arq_dust[10]  # temperatura da poeira
            aux = t_dust[0, 0, :, :, 0]
            t_int = np.max(aux[0])  #Tint
            Tint.append(t_int)

    return np.array(Tint).astype(int)
    



