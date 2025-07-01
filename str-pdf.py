import numpy as np
#import matplotlib.pyplot as plt
import sys
import os
import subprocess

def sub(exp,calc):
    chi2=0.0
    delta=0.0
    for i in range (calc.shape[0]):
        chi2=chi2+(exp[i,1]-calc[i,1])**2
    return chi2

def main():
    lv=25.24
    cut=630
    g=np.zeros((0,2))
    temp=np.zeros((1,2))
    with open ('exp-pdf', 'rt') as f1:
        for line in f1:
            wordList = []
            for word in line.split():
                wordList.append(word)
            temp[0,0]=float(wordList[0])
            temp[0,1]=float(wordList[1])    
            g=np.append(g,temp,0)
    l=0
    one_pdf=np.zeros((0,2))
    temp=np.zeros((1,2))
    os.chdir("../stored")
    tt_pdf=np.zeros((cut,2))
    if len(os.listdir()) != 0:
        for i in os.listdir():
            l=l+1
            one_pdf=np.zeros((0,2))
            with open (i, 'rt') as f5:
                j=0
                for line in f5:
                    j=j+1
                    wordList = []
                    for word in line.split():
                        wordList.append(word)
                    temp[0,0]=float(wordList[0])
                    temp[0,1]=float(wordList[1])
                    one_pdf=np.append(one_pdf,temp,0)
            tt_pdf=tt_pdf+one_pdf            
    num=l
    os.chdir("../pdf")
    prev_str=open('previous_pdf','w')
    for i in range(cut):
        prev_str.write(str(tt_pdf[i,0])+'  '+str(tt_pdf[i,1])+'  '+'\n')

    if num>0:
        pdf_init_error=1000.0
        step=0.02
        factor_temp=1.0
        scale_factor=0.0
        for i in range(100):
            factor_temp=scale_factor+step*i
            pdf_error=sub(g,factor_temp*tt_pdf/num)
            if pdf_error < pdf_init_error:
                pdf_init_error=pdf_error
                scale_factor=factor_temp
    else:
        scale_factor=1.0
    scale_factor=1.0
    temp=np.zeros((1,1))
    order=np.zeros((0,1))
    with open('order','rt') as f0:
        for line in f0:
            wordList = []
            for word in line.split():
                wordList.append(word)
            temp[0,0]=float(wordList[0])
            order=np.append(order,temp,0)
    f0.close()

    print('going on!')
    np.random.shuffle(order)
    o3=open('order','w')
    for i in range(order.shape[0]):
        o3.write(str(int(order[i,0]))+'\n')
    o4=open('factor', 'w')
    o4.write(str(scale_factor))
main()
