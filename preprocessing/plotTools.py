#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:57:59 2019

@author: ch199899
"""

import numpy as np
import os

import matplotlib
matplotlib.use("Agg")
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.backend_bases
import matplotlib.animation as manimation



#%% plot single slice over time
#
#
def plotImagesWithMaskVideo( data, times, slices, masks=None, frameRateX=4, titles=None, maxVal=None, vidName=None ):
    
    if not vidName is None:
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='', artist='')
        writer = FFMpegWriter(fps=15, metadata=metadata)
        
        
    # determine max frame rate
    difTime = np.Inf
    maxTime = 0
    for rm in range(len(data)):
        tmpDif = np.min(np.abs(np.diff(times[rm])))
        if tmpDif < difTime:
            difTime = tmpDif
            
        if maxTime < np.max(times[rm]):
            maxTime = np.max(times[rm])
        
    timeVec = np.arange(0,maxTime,difTime)
    
    fig, ax = plt.subplots(1,len(data))
    
    if not ax is list:
        ax = [ax]
    
    with writer.saving(fig, vidName, 125):
       
        # set each frame
        for sl in slices[rm]:
             ## for every time instance
             for tt in range(len(timeVec)):
             
                 for rm in range(len(data)):
                    
                    # determine the appropriate time
                    timeIX = np.argmin( np.abs(timeVec[tt]  - times[rm]) )
                    
                    
                    ax[rm].cla()
                    image = data[rm][:,:,sl,timeIX]
                    
                    if maxVal is None:
                        maxValscale = np.max(image.ravel())
                    else:
                        maxValscale = maxVal[rm]
                    
                    # plot image
                    ax[rm].imshow( image.T  \
                          ,clim=(-1,maxValscale), cmap='gray',alpha=1)
                    
                    if not titles is None:
                        ax[rm].set_title( titles[rm] )
                    
                    # plot mask contours
                    if not masks is None:
                        for km in masks.keys():
                            ax[rm].contour(masks[km][:,:,sl].T, levels=[1],colors='r')
                    
                    # set axis
                    ax[rm].set_xticks( [] )
                    ax[rm].set_xticklabels( [] )
                    ax[rm].set_yticklabels( [] )
                    ax[rm].set_yticks( [] )
                    
                 if not vidName is None:
                    for scale in range(frameRateX): # multiply the number of frames per plot
                        writer.grab_frame()
            
                 plt.savefig('./'+vidName[:-4]+'_sl'+str(sl)+ '_tt'+str(tt)+ '.png')



#%% plot single slice over time
#
#
def plotMIPVideo( data, times=None, frameRateX=4, titles=None, maxVal=None, vidName=None ):
    
    if not vidName is None:
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='', artist='')
        writer = FFMpegWriter(fps=15, metadata=metadata)
        
    if times is None:
        times = []
        for rm in range(len(data)):
            times.append( range(data[rm].shape[-1]) )
        
    # determine max frame rate
    difTime = np.Inf
    maxTime = 0
    for rm in range(len(data)):
        tmpDif = np.min(np.abs(np.diff(times[rm])))
        if tmpDif < difTime:
            difTime = tmpDif
            
        if maxTime < np.max(times[rm]):
            maxTime = np.max(times[rm])
        
    timeVec = np.arange(0,maxTime,difTime)
    
    fig, ax = plt.subplots(1,len(data))
    
    if not ax is list:
        ax = [ax]
    
    figureNames = []
    with writer.saving(fig, vidName, 125):
       
        ## for every time instance
        for tt in range(len(timeVec)):
             
            for rm in range(len(data)):
                    
                # determine the appropriate time
                timeIX = np.argmin( np.abs(timeVec[tt]  - times[rm]) )
                
                
                ax[rm].cla()
                image = np.max( data[rm][:,:,:,timeIX], axis=2)
                
                if maxVal is None:
                    maxValscale = np.max(image.ravel())
                else:
                    maxValscale = maxVal[rm]
                
                # plot image
                ax[rm].imshow( image.T  \
                        ,clim=(-1,maxValscale), cmap='gray',alpha=1)
                
                if not titles is None:
                    ax[rm].set_title( titles[rm] )
                else:
                    ax[rm].set_title( 'Frame: ' + str(tt) )
                
                # set axis
                ax[rm].set_xticks( [] )
                ax[rm].set_xticklabels( [] )
                ax[rm].set_yticklabels( [] )
                ax[rm].set_yticks( [] )
                
                if not vidName is None:
                    for scale in range(frameRateX): # multiply the number of frames per plot
                        writer.grab_frame()
                
                figureNames.append( vidName[:-4] + '_MIP_tt'+str(tt)+ '.png' )
                plt.savefig( figureNames[-1])

    # remove png files
    for ff in figureNames:
        os.remove(ff)

#%%
#
#   This class is designed to plot a desired image and wait for the user to manually select a number of points on it.
#   Once it is done, it returns the temporal data of the pixels (or averages around it) selected.
#
class selectDataOnImg:

    ## TODO: create a list of rectangles and get rid of so much redundancy in the code
    rect = []
    press = None
    cidpress = []
    cidrelease = []
    cidmotion = []
    numPts = 0
    currPts = 0
    tt = 0
    
    results = []
    data = None
    data2 = None
    timer = 300

    fig = None
    ax = None
    width = 1

    def __init__( self, data, data2=[], numPts=1, width=(1,1,1), pos=None, title='' ):
        
        self.press = None
        self.numPts = numPts

        self.data = data
        self.data2 = data2
        volShape = self.data.shape

        self.results = list(range(self.numPts))

        if isinstance(width,(tuple,list)):
            self.width = width
        else:
            self.width = [ width for dd in range(3)]

        if pos is None:
            pos = [0,0,0,0]
            pos[3] = int(np.round( volShape[3] / 2))
            pos[2] = int(np.round( volShape[2] / 2))
            pos[1] = int(np.round( volShape[1] / 2))
            pos[0] = int(np.round( volShape[0] / 2))
        self.tt = pos[3]

        # create figure
        self.fig, ax= plt.subplots(2,2)
        self.fig.subplots_adjust( left=0.1, bottom=0.1,  wspace=0.05, hspace=0.1)
        self.fig.suptitle(title)
        self.ax = [ ax[0,0], ax[0,1], ax[1,0], ax[1,1] ]
        

        self.ax[0].imshow( self.data[:,:,pos[2],pos[3]].T, cmap='gray' )
        self.ax[1].imshow( self.data[pos[0],:,:,pos[3]], cmap='gray' )
        self.ax[2].imshow( self.data[:,pos[1],:,pos[3]].T, cmap='gray' )

        self.ax[3].cla()
        for rr in range(3):
            self.ax[rr].set_xticks([])
            self.ax[rr].set_yticks([])
            self.ax[rr].set_aspect('auto')
        

        # initialize rectangles
        self.rect = []
        for rr in range(self.numPts):

            tmpRec = [  matplotlib.patches.Rectangle( xy=( 0,0 ), width=self.width[0], height=self.width[1], alpha=0.25, color='r' ), \
                        matplotlib.patches.Rectangle( xy=( 0,0 ), width=self.width[2], height=self.width[1], alpha=0.25, color='r' ), \
                        matplotlib.patches.Rectangle( xy=( 0,0 ), width=self.width[0], height=self.width[2], alpha=0.25, color='r' )]

            self.rect.append( tmpRec )

            self.ax[0].add_patch(self.rect[rr][0])
            self.ax[1].add_patch(self.rect[rr][1])
            self.ax[2].add_patch(self.rect[rr][2])

        self.values = [ np.array((0,)) for rr in range(numPts) ]

        self.connect()

        print('Operation paused to select points. If no point is selected, will exit in %0.1fmin.' %(self.timer/60))
        plt.show()
        plt.pause(self.timer)

    def updatePlots(self):

        x0, y0 = self.rect[self.currPts][0].get_xy()
        z0 = self.rect[self.currPts][1].get_xy()[0]

        self.ax[0].imshow( self.data[:,:,int(z0),self.tt].T, cmap='gray' )
        self.ax[1].imshow( self.data[int(x0),:,:,self.tt], cmap='gray' )
        self.ax[2].imshow( self.data[:,int(y0),:,self.tt].T, cmap='gray' )
        
        self.ax[3].cla()
        for rr in range(3):
            self.ax[rr].set_xticks([])
            self.ax[rr].set_yticks([])
            self.ax[rr].set_aspect('auto')

        self.results[self.currPts] = {'position': [x0,y0,z0], 'value':self.data[int(x0),int(y0),int(z0),:], \
                                        'value2':[ self.data2[rr][int(x0),int(y0),int(z0),:] for rr in range(len(self.data2))]} 

        self.ax[3].cla()
        for pt in range(self.currPts+1):
            self.ax[3].plot( self.results[pt]['value'] )
            for rr in range(len(self.data2)):
                self.ax[3].plot( self.results[pt]['value2'][rr] )
            

    def connect(self):
        'connect to all the events we need'
        self.cidpress = [list(range(3)) for rr in range(self.numPts) ]
        self.cidrelease = [list(range(3)) for rr in range(self.numPts) ]
        self.cidmotion = [list(range(3)) for rr in range(self.numPts) ]
        for pp in range(self.numPts):
            for rr in range(3):
                self.cidpress[pp][rr] = self.rect[pp][rr].figure.canvas.mpl_connect(
                    'button_press_event', self.on_press)
                self.cidrelease[pp][rr] = self.rect[pp][rr].figure.canvas.mpl_connect(
                    'button_release_event', self.on_release)
                self.cidmotion[pp][rr] = self.rect[pp][rr].figure.canvas.mpl_connect(
                    'motion_notify_event', self.on_motion)


    def on_press(self, event):

        if event.button == 3:
            self.nextPoint()

        if event.button == 1:
            for rr in range(3):

                if event.inaxes == self.ax[rr]:
                    
                    if rr == 0:
                        width = [np.round(self.width[0]/2.), np.round(self.width[1]/2.)]
                    if rr == 1:
                        width = [np.round(self.width[2]/2.), np.round(self.width[1]/2.)]
                    if rr == 2:
                        width = [np.round(self.width[0]/2.), np.round(self.width[2]/2.)]

                    print(width)
                    x0 = event.xdata - width[0]
                    y0 = event.ydata - width[1]
                    
                    self.rect[self.currPts][rr].set_x(x0)
                    self.rect[self.currPts][rr].set_y(y0)

                    if rr == 0:
                        self.rect[self.currPts][2].set_x(x0)
                        self.rect[self.currPts][1].set_y(y0)
                    elif rr == 1:
                        self.rect[self.currPts][0].set_y(y0)
                        self.rect[self.currPts][2].set_x(x0)
                    elif rr == 2:
                        self.rect[self.currPts][0].set_x(x0)
                        self.rect[self.currPts][1].set_x(y0)
                    

                    x0, y0 = self.rect[self.currPts][0].get_xy()
                    z0 = self.rect[self.currPts][1].get_xy()[0]
                    print('Current point: %d. Position %d,%d,%d' %(self.currPts,x0,y0,z0) )

                    x0, y0 = self.rect[self.currPts][rr].xy
                    self.press = x0, y0, event.xdata, event.ydata

            self.updatePlots()

    def on_motion(self, event):

        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return

        x0, y0, xpress, ypress = self.press

        if event.button == 1:
            for rr in range(3):
                if event.inaxes == self.ax[rr]:

                    if rr == 0:
                        width = [np.round(self.width[0]/2.), np.round(self.width[1]/2.)]
                    if rr == 1:
                        width = [np.round(self.width[2]/2.), np.round(self.width[1]/2.)]
                    if rr == 2:
                        width = [np.round(self.width[0]/2.), np.round(self.width[2]/2.)]

                    dx = event.xdata - xpress - width[0]
                    dy = event.ydata - ypress - width[1]
                    
                    self.rect[self.currPts][rr].set_x(x0+dx)
                    self.rect[self.currPts][rr].set_y(y0+dy)

                    if rr == 0:
                        self.rect[self.currPts][2].set_x(x0+dx)
                        self.rect[self.currPts][1].set_y(y0+dy)
                    elif rr == 1:
                        self.rect[self.currPts][0].set_y(y0+dy)
                        self.rect[self.currPts][2].set_x(x0+dx)
                    elif rr == 2:
                        self.rect[self.currPts][0].set_x(x0+dx)
                        self.rect[self.currPts][1].set_x(y0+dy)
                    
        else:
            return

        for pp in range(self.numPts):
            for rr in range(3):
                self.rect[pp][rr].figure.canvas.draw()

        self.updatePlots()
        

    def on_release(self, event):

        'on release we reset the press data'
        self.press = None
        for pp in range(self.numPts):
            for rr in range(3):
                self.rect[pp][rr].figure.canvas.draw()
        
        self.updatePlots()

    def disconnect(self):
        'disconnect all the stored connection ids'
        for pp in range(self.numPts):
            for rr in range(3):
                self.rect[pp][rr].figure.canvas.mpl_disconnect(self.cidpress[pp][rr])
                self.rect[pp][rr].figure.canvas.mpl_disconnect(self.cidrelease[pp][rr])
                self.rect[pp][rr].figure.canvas.mpl_disconnect(self.cidmotion[pp][rr])


    def nextPoint(self):

        x0, y0 = self.rect[self.currPts][0].get_xy()
        z0 = self.rect[self.currPts][1].get_xy()[0]
        print('Button %d selected at position (%d,%d,%d). Value: %0.3f' %(self.currPts, x0, y0,z0, np.mean(self.values[-1])))

        self.currPts += 1
        if self.currPts < self.numPts:
            print('Next point')
        else:
            self.disconnect()
            self.fig.canvas.stop_event_loop()