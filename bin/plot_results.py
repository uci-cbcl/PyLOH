#!/usr/bin/env python

import sys
import scipy
from matplotlib import pyplot

def main():
    infile = open(sys.argv[1])
    
    xTHetA = []
    yTHetA = []
    xPurBayes = []
    yPurBayes = []
    xPyLOH = []
    yPyLOH = []
    
    for line in infile:
        tag, est, true = line.strip('\n').split('\t')
        est = float(est)
        true = float(true)
        
        if tag == 'THetA':
            xTHetA.append(est)
            yTHetA.append(true)
        elif tag == 'PurBayes':
            xPurBayes.append(est)
            yPurBayes.append(true)
        elif tag == 'PyLOH':
            xPyLOH.append(est)
            yPyLOH.append(true)
            
    params = {'legend.fontsize': 6, 'legend.linewidth': 4}
    pyplot.rcParams.update(params)
    pyplot.rcParams['font.size'] = 6
    pyplot.figure(figsize=(5, 5), dpi = 150)
    
    pyplot.plot([0, 1], [0, 1], '--', color = 'black', linewidth=0.5, label = 'Ideal')
    pyplot.plot(xTHetA, yTHetA, 'r^', alpha = 0.6, markersize = 4, label = 'THetA')
    pyplot.plot(xPurBayes, yPurBayes, 'gv', alpha = 0.6, markersize = 4, label = 'PurBayes')
    pyplot.plot(xPyLOH, yPyLOH, 'bo', alpha = 0.6, markersize = 4, label = 'PyLOH')
    pyplot.xlim(0, 1)
    pyplot.ylim(0, 1)
    pyplot.xticks(scipy.linspace(0, 1, 11), ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    pyplot.yticks(scipy.linspace(0, 1, 11), ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    pyplot.xlabel('Estimated tumor cellular frequency')
    pyplot.ylabel('Reported tumor cellular frequency by ABSOLUTE/Ground truth tumor cellular frequency')
    pyplot.grid(True)
    pyplot.legend(loc=4,prop={'size':6})
    pyplot.show()
    
    
    
    infile.close()
    
if __name__ == '__main__':
    main()