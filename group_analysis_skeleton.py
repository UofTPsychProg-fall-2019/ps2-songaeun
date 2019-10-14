#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scene-cat problem set for PSY 1210 - Fall 2018

@author: Michael Mack
"""

#%% import block 
import numpy as np
import scipy as sp
import scipy.stats
import os
import shutil


#%%
# copy files from testing room folders to raw data, rename files to include
# testing room letter in the filename
#
baseDir = os.getcwd() # If your cwd is not 'ps2-songaeun', please change the path.
testingrooms = ['A','B','C']
for room in testingrooms:
    oPath = os.path.join(baseDir, 'testingroom'+room, 'experiment_data.csv') # get original data path
    nPath = os.path.join(baseDir, 'rawdata', 'experiment_data_'+room+'.csv') # get new data path
    shutil.copyfile(oPath, nPath) # gather data from three room in a new data folder


#%%
# read in all the data files in rawdata directory using a for loop
# columns: subject, stimulus, pairing, accuracy, median RT
#
data = np.empty((0,5))
for room in testingrooms:
    nPath  = os.path.join(baseDir, 'rawdata', 'experiment_data_'+room+'.csv') # get new data path to load
    tmp = sp.loadtxt(nPath,delimiter=',') # load data
    data = np.vstack([data,tmp]) # stack data in a column


#%%
# calculate overall average accuracy and average median RT
#
acc_avg = np.mean(data[:,3])   # 91.48%
mrt_avg = np.mean(data[:,4])   # 477.3ms
print(f'{np.round(100*acc_avg,2)}%, {np.round(mrt_avg,1)}ms')

#%%
# calculate averages (accuracy & RT) split by stimulus using a for loop and an 
# if statement. (i.e., loop through the data to make a sum for each condition, 
# then divide by the number of data points going into the sum)
#

# Gaeun: I splitted data into sjbdata (word, face) using for loop and if statement,
#        and directly calculated averages.
subdata_stim = [0,0]
subdata_stim[0] = np.empty((0,5))
subdata_stim[1] = np.empty((0,5))
for d in np.arange(0, np.size(data,axis=0), 1):
    if data[d,1]==1:
        tmp = data[d,:]
        subdata_stim[0] = np.vstack([subdata_stim[0],tmp])
    else:
        tmp = data[d,:]
        subdata_stim[1] = np.vstack([subdata_stim[1],tmp])

acc_word = np.mean(subdata_stim[0][:,3])
mrt_word = np.mean(subdata_stim[0][:,4])
acc_face = np.mean(subdata_stim[1][:,3])
mrt_face = np.mean(subdata_stim[1][:,4])

print(f'word: {np.round(100*acc_word,1)}% {np.round(mrt_word,1)}ms, ace: {np.round(100*acc_face,1)}% {np.round(mrt_face,1)}ms')
# words: 88.6%, 489.4ms   faces: 94.4%, 465.3ms


#%%
# calculate averages (accuracy & RT) split by congruency using indexing, 
# slicing, and numpy's mean function 
# wp - white/pleasant, bp - black/pleasant
# (hint: only one line of code is needed per average)
#

pair1 = data[data[:,2]==1]
pair2 = data[data[:,2]==2]

acc_wp = np.mean(pair1[:,3])  # 94.0%
acc_bp = np.mean(pair2[:,3])  # 88.9%
mrt_wp = np.mean(pair1[:,4])  # 469.6ms
mrt_bp = np.mean(pair2[:,4])  # 485.1ms


#%% 
# calculate average median RT for each of the four conditions
# use for loops, indexing/slicing, or both!
# (hint: might be easier to slice data into separate words and faces datasets)
#
words = data[data[:,1]==1]
faces = data[data[:,1]==2]
mrt_words_wp = np.mean(words[words[:,2]==1,4])
mrt_words_bp = np.mean(words[words[:,2]==2,4])
mrt_faces_wp = np.mean(faces[faces[:,2]==1,4])
mrt_faces_bp = np.mean(faces[faces[:,2]==2,4])

# words - white/pleasant: 478.4ms
# words - black/pleasant: 500.3ms
# faces - white/pleasant: 460.8ms
# faces - black/pleasant: 469.9ms


#%%        
# compare pairing conditions' effect on RT within stimulus using scipy's 
# paired-sample t-test: scipy.stats.ttest_rel()
#
import scipy.stats
words_tTest = scipy.stats.ttest_rel(words[words[:,2]==1,4], words[words[:,2]==2,4])
faces_tTest = scipy.stats.ttest_rel(faces[faces[:,2]==1,4], faces[faces[:,2]==2,4])
# words: t=-5.36, p=2.19e-5
# faces: t=-2.84, p=0.0096


#%%
# print out averages and t-test results
# (hint: use the ''.format() method to create formatted strings)
#
acc_words_wp = np.mean(words[words[:,2]==1,3])
acc_words_bp = np.mean(words[words[:,2]==2,3])
acc_faces_wp = np.mean(faces[faces[:,2]==1,3])
acc_faces_bp = np.mean(faces[faces[:,2]==2,3])

print('Words >> wp: {:.2f}%, {:.1f}ms, bp: {:.2f}%, {:.1f}ms, t={:.2f} p={:.4f}'\
      '\nFaces >> wp: {:.2f}%, {:.1f}ms, bp: {:.2f}%, {:.1f}ms, t={:.2f} p={:.4f}'\
      .format(100*acc_words_wp, mrt_words_wp, 100*acc_words_bp, mrt_words_bp, words_tTest[0], words_tTest[1],\
              100*acc_faces_wp, mrt_faces_wp, 100*acc_faces_bp, mrt_faces_bp, faces_tTest[0], faces_tTest[1]))
