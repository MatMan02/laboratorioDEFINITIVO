import re
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import os

#ciao sono mattia, ciao Dr Mattia
#Uagliù, siam proprio degli scunnizzi uaaa
nomeFile = "C240920_0o5OD_1OD_orangeF_16K.txt"

# array con [ nome campione, nome filtro colorato (ignorare) ]
nomeCampioneArr = re.findall("([^_]{7})(?:_)", nomeFile)

# array con [ (lente assorbimento, lente assorbimento decimi), (lente emissione, lente emissione decimi), (temperatura, decimi) ]
match = re.findall("(?:_)(\d+)(?:o)?(\d+)?", nomeFile)
assorbimento = float(match[0][0] + "." + match[0][1])
emissione = float(match[1][0] + "." + match[1][1])
temperatura = float(match[2][0])

w, I = np.loadtxt(nomeFile, delimiter=';', skiprows=7, usecols=(0, 1), unpack=True)

plt.plot(w,I,color='blue', label='signal')
plt.title(nomeCampioneArr[0] + ". Assorbimento: " + str(assorbimento) + " OD; emissione: " + str(emissione) + " OD")
plt.xlabel('wavelenght [nm]')
plt.ylabel('intensity [count]')
plt.legend()
plt.show()

wMin = input("lambda minima: ")
wMax = input("lambda massima: ")

def findIndex(valore, array):
    n = 0
    for el in array:
        if int(el) == int(valore):
            return n
        else:
            n += 1

nMin = findIndex(wMin, w)
nMax = findIndex(wMax, w)

w=w[nMin:nMax]
I=I[nMin:nMax]
#dw=np.abs(w[0]-w[1])/(np.sqrt(12))
dI = np.sqrt(I)

h=6.626e-34
c=3e+8
e=1.6e-19
E=(h*c)/(w*10**(-9)*e)

#od in emissione
I=I*10**(-emissione)
dI=dI*10**(-emissione)

plt.plot(E, I, color='blue', label='signal')
#plt.title(nomeCampioneArr[0] + ". picco da: " + str(E[nMin]) + " e a " + str(E[nMax]) + " nm")
plt.xlabel('energy [eV]')
plt.ylabel('intensity [count]')
plt.legend()
plt.show()

#chi quadro
def chisq(obs, exp, error):
    return np.sum((obs - exp) ** 2 / (error ** 2))


# FUNZIONE DI FIT DOPPIA GAUSSIANA
def gauss(x, A, sigma, x0, B, sigma1, x1):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + B * np.exp(-(x - x1)**2 / (2 * sigma1**2))

(args, var) = curve_fit(gauss, E, I, p0 = [max(I), 0.0005, min(E) + (max(E) + min(E))/4, max(I), 0.0005, max(E) - (max(E) + min(E))/4], sigma=dI)
A = args[0]
sigma = args[1]
Egap = args[2]
B = args[3]
sigma1 = args[4]
Egap_qds = args[5]
dA = np.sqrt(var[0,0])
dsigma = np.sqrt(var[1,1])
dEgap = np.sqrt(var[2,2])
dB = np.sqrt(var[3,3])
dsigma1 = np.sqrt(var[4,4])
dEgap1 = np.sqrt(var[5,5])
FWHM = sigma1*2.355
dFWHM = dsigma1*2.355

plt.errorbar(E, I, marker='', linestyle='-', linewidth=2, color='blue', label = 'signal')
t = np.linspace(np.min(E), np.max(E), 100)
plt.plot(t, gauss(t, *args), linestyle="--", color="red", label="fit gaussiano")
plt.xlabel('energy [eV]')
plt.ylabel('intensity [count]')
plt.legend()
plt.show()

print('Egap = ', Egap, 'Egap Qds = ', Egap_qds)

#salvare nel file esterno
'''energy=round(Egap_qds, 1)
domanda=input('Stai variando le densità ottiche (0) o la temperatura (1)? ')
'''
'''if(int(domanda)==0):
    linea_iniziale='OD \t Egap(eV) \t dEgap \t FWHM \t dFWHM \t Imax \t dImax \n'
    nome=  'dati_' + nomeCampioneArr[0] + '_' + str(temperatura) + "K_picco=" + str(energy) +'.txt'
    # Controlla se il file esiste
    if not os.path.exists(nome):
        # Se non esiste, crealo e aggiungi la linea iniziale
        with open(nome, "w") as output_file:
            output_file.write(linea_iniziale)
    with open (nome, 'a') as output_file :
        output_file.write ( str(assorbimento) + '\t'  + str(Egap) + '\t' + str(dEgap) + '\t' + str(FWHM) + '\t' + str(dFWHM) + '\t' + str(Imax) + '\t' +  str(dImax) + '\n')
else:
    linea_iniziale='T(K) \t Egap(eV) \t dEgap \t FWHM \t dFWHM \t Imax \t dImax \n'
    nome=  'dati_' + nomeCampioneArr[0] + '_' + str(temperatura) + "K_picco=" + str(energy) +'.txt'
    # Controlla se il file esiste
    if not os.path.exists(nome):
        # Se non esiste, crealo e aggiungi la linea iniziale
        with open(nome, "w") as output_file:
            output_file.write(linea_iniziale)
    nome=  'dati_' + nomeCampioneArr[0] + '_' + str(assorbimento) + "K_picco=" + str(energy) +'.txt'
    with open (nome, 'a') as output_file :
        output_file.write ( str(temperatura) + '\t'  + str(Egap) + '\t' + str(dEgap) + '\t' + str(FWHM) + '\t' + str(dFWHM) + '\t' + str(Imax) + '\t' +  str(dImax) + '\n')'''