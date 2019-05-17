# -*- coding: utf-8 -*-
"""
Created on Wed May  1 22:38:55 2019

@author: Wenqing
"""

import numpy as np
import random
import math
import time


class Motif:
    def __init__(self, ml, seqs, path):
        
        self.ml =ml
        self.seqs = seqs.copy()
        #self.pwm = np.zeros((4, ml))
        self.cand_dic = {}
        self.chars = ['A','C','G','T']
        self.path = path
        self.sc = len(seqs)
        self.sl = 500
        self.max_icpc = 2
        self.loop_num = 5
        self.start_pair = 100
        
        
    def getPWM(self, loc_ls):
        """given sites,
        calculate pwm"""
        pwm = np.zeros((4, ml))
        for ind in range(self.ml):
            col = list(self.seqs[i][loc_ls[i]+ind] for i in range(self.sc))
            for char_id in range(4):
                pwm[char_id, ind] = sum(list(1 for ele in col if ele==self.chars[char_id]))/self.sc
        return(pwm)
    
    def getPWM2(self, loc_ls_seqs, loc_ls_sites):
        """Given index of seqs and sites
        calculate pwm and content info"""
        pwm = np.zeros((4, ml))
        cnt_info = 0
        for ind in range(self.ml):
            col = list(self.seqs[loc_ls_seqs[i]][loc_ls_sites[i]+ind] for i in range(len(loc_ls_seqs)))
            for char_id in range(4):
                pwm[char_id, ind] = sum(list(1 for ele in col if ele==self.chars[char_id]))/len(loc_ls_seqs)
            
            cnt_info +=  sum(list(p*math.log(p/0.25,2) for p in pwm[:,ind] if p!=0))

        return(pwm, cnt_info)
    
    def output(self,pwm, loc_ls):
        #export motif
        with open(self.path+'/'+"predictedmotif.txt", "w") as file: 
            file.write(">MOTIF "+ str(self.ml)+'\n') 
            for i in range(self.ml): #print by column
                temp_ls = list(pwm[:,i])
                temp_str = ''
                for ele in temp_ls:
                    temp_str = temp_str+str(ele)+' '
                temp_str = temp_str[:-1]
                file.write(temp_str+'\n')
            file.write(">\n") 
            
        #export sites
        with open(self.path+'/'+"predictedsites.txt", "w") as file: 
            for site in loc_ls:
                file.write(str(site)+'\n')        

    def findFollowMtf(self,pair, seq_ord):
        """for each tied pair, find following sites
        return pwm, sites and ci"""
        sites = pair.copy()
        for i in range(2,self.sc):
            #get sub seq with min difference
            temp_info_flw = 0
            temp_site = -1
            for loc in range(self.sl-self.ml+1):
                #get info
                pwm, cnt_info = self.getPWM2(seq_ord[:i+1],sites+[loc])
                #print("debug: cnt_info is",cnt_info )
                if cnt_info > temp_info_flw:
                    temp_info_flw = cnt_info
                    temp_site = loc
                    if i == self.sc-1: #last sequence
                        temp_pwm = pwm
            sites.append(temp_site)
        return([temp_info_flw, sites, temp_pwm])
        
        
    def findMtf2(self):
        '''greedy algorithm
        find motif seqs and pwm
        statistics for each column''' 
        cur_info = 0
        cur_sites = []
        cur_pwm = {}
        cur_order = []
        count = 0
        while count< self.loop_num:
            #disordder seqs
            seq_ord = random.sample(range(self.sc), k=self.sc)
            print("order :", seq_ord)
            pairs = []
            pairs_ci = []
            #find first pair: min difference/ max info
            for i in range(self.sl-self.ml+1):
                for j in range(self.sl-self.ml+1):
                    pwm, cnt_info = self.getPWM2(seq_ord[:2],[i,j])
                    if len(pairs)<self.start_pair or (len(pairs)>=self.start_pair and cnt_info==self.max_icpc*self.ml and min(pairs_ci)==self.max_icpc*self.ml):
                        pairs.append([i,j])
                        pairs_ci.append(cnt_info)
                    else:
                        if cnt_info > min(pairs_ci): #need to replace
                            idx = pairs_ci.index(min(pairs_ci))
                            pairs[idx] = [i,j]
                            pairs_ci[idx] = cnt_info
            split_res = []
            for pair in pairs:
                #for each pair, find following sites
                split_res.append(self.findFollowMtf(pair, seq_ord))
                if split_res[-1][0] == self.max_icpc*self.ml:
                    break
            #find split with highest ci
            chosen_split = max(split_res, key = lambda splt: splt[0])
            #check if larger ci found
            if chosen_split[0] > cur_info:
                cur_info = chosen_split[0]
                cur_sites = chosen_split[1].copy()
                cur_pwm = chosen_split[2]
                cur_order = seq_ord
            print("bebug:",cur_info)
            if cur_info == self.max_icpc*self.ml:
                break
            count += 1
            
        #export sites and pwm 
        #reorder sites
        order_site = []
        for i in range(len(self.seqs)):
            site_loc = cur_sites[cur_order.index(i)]
            order_site.append(site_loc)
        #order_site = sorted(cur_sites, key=lambda site: cur_order[cur_sites.index(site)])
        print("final result:", cur_info, '\n', cur_pwm, '\n', order_site)
        
        self.output(cur_pwm, order_site)
        
    
                
        
    def findMtf(self):
        '''find motif seqs and pwm'''
        for ind in range(len(self.seqs[0])-self.ml+1):
            #check all candidate
            cand = self.seqs[0][ind:ind+self.ml]
            loc_ls = []
            for seq in self.seqs:
                loc = seq.find(cand)
                if loc != -1:
                    loc_ls.append(loc)
            if len(loc_ls) == self.sc:
                #if cand appear in all seqs
                self.cand_dic[cand] = loc_ls
                
        if len(self.cand_dic) != 0:
            print("attemp1")
            cand = list(self.cand_dic.keys())[0]
            pwm = self.getPWM(self.cand_dic[cand])           
            #export motif
            self.output(pwm, loc_ls)
           
        else:
            print("attemp2")
            self.findMtf2()
            
                


if __name__ == '__main__':
    
    
    for count in range(11,21):
        print("working on",count)
        t1 = time.time()
        #path 
        path = "data set"+str(count)
        #sequences
        seqs = []
        #read in length
        with open(path+'/'+"motiflength.txt", "r") as file: 
            ml = int(file.readline())
            
        with open(path+'/'+"sequences.txt", "r") as file: 
            line = file.readline()
            while line:
                if line[0]=='>': #name line
                    line = file.readline()
                    continue
                seqs.append(line[:-1]) #remove \n
                line = file.readline()
        motif_obj = Motif(ml, seqs, path)
        motif_obj.findMtf2()
        t2 = time.time()
        with open(path + '/' + "evalruntime.txt", "w") as file:
            file.write(str(t2-t1) + '\n')
        print(count,"time used:", t2-t1)
