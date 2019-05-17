# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:06:20 2019

@author: Wenqing
"""
import os
import random
import numpy as np
#from Bio.SeqIO import FastaIO


class BenchMark:
    def __init__(self, ICPC, ML, SL, SC, count):
        self.icpc = ICPC #infor per column
        self.ml = ML #motif length
        self.sl = SL #seq length
        self.sc = SC #seq count
        self.count = count
        self.path = "data set"+ str(count)
        
        self.seqs = [] #random sequence
        self.pwm = np.zeros((4,self.ml))
        self.sites = [] #site start from 0
        self.aims = []
        
    
    def getSample(self):
        """generate samples"""
        #create directory
        try:
            os.mkdir(self.path)
        except:
            pass
        
        #write ml in to motiflength.txt
        with open(self.path+'/'+"motiflength.txt", "w") as file:
            file.write(str(self.ml)+'\n')
            
        #generate random sequence
        for i in range(self.sc):
            seq = ''
            for j in range(self.sl):
                seq = seq + random.choice(['A','C','G','T'])
            self.seqs.append(seq)
            
        #generate random motif
        #random.sample(range(8),4)
        if self.icpc ==2:
            #unique possible
            for i in range(self.ml):
                cell = random.choice(range(4))
                self.pwm[cell][i] = 1
        elif self.icpc == 1.5:
            #two possible
            if self.ml%2 ==0:
                #half fix
                fix_col = random.sample(range(self.ml),int(self.ml/2))
                for i in range(self.ml):
                    if i in fix_col:
                        cell = random.choice(range(4))
                        self.pwm[cell][i] = 1
                    else: #random col
                        cells = random.sample(range(4),k=2)
                        self.pwm[cells[0],i] = 0.5
                        self.pwm[cells[1],i] = 0.5
            else: # ignore here
                fix_col = random.sample(range(self.ml),(self.ml-1)/2)
                ran_col = list(set(range(self.ml))-set(fix_col))
                col_15 = random.sample(ran_col,1)[0]
                for i in range(self.ml):
                    if i in fix_col:
                        cell = random.choice(range(4))
                        self.pwm[cell][i] = 1
                    elif i==col_15:
                        
                        cells = random.sample(range(4),k=2)
                        self.pwm[cells[0],i] = 0.110027
                        self.pwm[cells[1],i] = 1-0.110027
                        
                    else: #random col
                        cells = random.sample(range(4),k=2)
                        self.pwm[cells[0],i] = 0.5
                        self.pwm[cells[1],i] = 0.5
                        
        else: #self.icpc==1
            #two possible equal
            for i in range(self.ml):
                temp_ls = list(range(4))
                cell1 = random.choice(temp_ls)
                temp_ls.remove(cell1)
                cell2 = random.choice(temp_ls)
                self.pwm[cell1, i] = 0.5
                self.pwm[cell2, i] = 0.5
        #generate aim string
        for i in range(self.sc):
            aim = ''
            for j in range(self.ml):
                #for each column
                aim = aim + random.choices(['A','C','G','T'], list(self.pwm[:,j]))[0]
            self.aims.append(aim)
        
        #generate sites
        rang = list(range(self.sl-self.ml+1))
        for i in range(self.sc):
            self.sites.append(random.choice(rang)) 
        
        #plant site
        for i in range(self.sc):
            site = self.sites[i]
            self.seqs[i] = self.seqs[i][:site] + self.aims[i] + self.seqs[i][site+self.ml:]
            
        
        #export sequences
        '''
        handle = open(self.path+'/'+"sequences.fa", "w")

        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        for record in self.seqs:
            fasta_out.write_record(record)
        '''
        with open(self.path+'/'+"sequences.txt", "w") as file: 
            for i in range(len(self.seqs)):
                file.write(">Sequence "+str(i)+"\n"+self.seqs[i]+'\n')
                
        #export sites
        with open(self.path+'/'+"sites.txt", "w") as file: 
            for site in self.sites:
                file.write(str(site)+'\n')
                
        #export motif
        with open(self.path+'/'+"motif.txt", "w") as file: 
            file.write(">MOTIF "+ str(self.ml)+'\n') 
            file.flush()
            for i in range(self.ml): #print by column
                temp_ls = list(self.pwm[:,i])
                temp_str = ''
                for ele in temp_ls:
                    temp_str = temp_str+str(ele)+' '
                temp_str = temp_str[:-1]
                file.write(temp_str+'\n')
            file.write(">\n") 
                
                

        
    
if __name__ == '__main__':
    
    icpc_ls = [1, 1.5]
    ml_ls = [6,7]
    sl = 500
    sc_ls = [5, 20]
    count = 1
    loop_c = 0
    
    #default sample
    while loop_c<10:
    #while loop_c<1:
        bm_obj = BenchMark(2, 8, sl, 10, count)
        bm_obj.getSample()
        count +=1
        loop_c +=1
   
    loop_c = 0
    for icpc in icpc_ls:
        while loop_c<10:
            bm_obj = BenchMark(icpc, 8, sl, 10, count)
            bm_obj.getSample()
            count +=1
            loop_c +=1
        loop_c = 0
          
    loop_c = 0    
    
    for ml in ml_ls:
        while loop_c<10:
            bm_obj = BenchMark(2, ml, sl, 10, count)
            bm_obj.getSample()  
            count +=1
            loop_c +=1
        loop_c = 0
    
    loop_c = 0
    for sc in sc_ls:
        while loop_c<10:
            bm_obj = BenchMark(2, 8, sl, sc, count)
            bm_obj.getSample() 
            count +=1
            loop_c +=1
        loop_c = 0
    
        
