import ROOT
import cv2 
import numpy as np
import os
import math
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def filtercircle(inputpic,tag):
    img = inputpic.copy()
    nrows, ncols, rgb = img.shape
    cx     = ncols/2
    cy     = nrows/2
    r      = nrows/2
    for i in range(0, nrows):         
      for j in range(0, ncols):
           distance = math.sqrt(  (j-cx)*(j-cx) + (i-cy)*(i-cy))
           if distance < 0.75*r and distance>0.12*r: continue
           for k in range(0, rgb): img[i,j,k] = 0
    cv2.imwrite("plots_multishot/test_circle_%s.png"%tag,img)
    return img

def convertcameracenter(coordinates,ncols,nrows):
    k = 0.04589621286 #mm per pixel
    print "[INFO] Converting corner information"
    print "   - Corner w.r.t photo origin: (%.3f,%.3f)"%(coordinates[0],coordinates[1])
    xp = coordinates[0]-(ncols/2)
    yp = coordinates[1]-(nrows/2)
    print "   - Corner w.r.t photo center: (%.3f,%.3f)"%(xp,yp)
    cx = round(xp*k,3)
    cy = round(yp*-k,3)
    print "   - Corner [mm]: (%.3f,%.3f)"%(cx,cy)
    return cx,cy

def preselect_averagey(pixels_xy):
    pixels =[]
    x=-1
    for k in pixels_xy:
        temp_x=k[0]
        temp_y=k[1]
        if temp_x==x: continue
        sumy=0
        ny  =0
        avg =0
        done=False
        for j in pixels_xy:
           if j[0] == temp_x:
              sumy+=j[1]
              ny +=1
              done=True  
        if done==True: 
            avg = float(sumy)/float(ny)
            x = temp_x
            pixels.append([x,avg])
    return pixels

def makegraph(inputpic):
    img = inputpic.copy()
    nrows, ncols, rgb = img.shape
    pixels_xy = []
    #Find all white pixel in x-axis
    for i in range(0, nrows):         
      for j in range(0, ncols):         
           if img[i,j,0] != 0: 
              pixels_xy.append([j,i])
    #order them from low x to high x          
    pixels_xy = sorted(pixels_xy, key=lambda x: x[0], reverse=False)
    #Preselect pixels by computing average in y
    pixels_xy = preselect_averagey(pixels_xy)
    x, y = array( 'f' ), array( 'f' )
    for i in range(0, len(pixels_xy)):
       x.append( pixels_xy[i][0] )
       y.append( pixels_xy[i][1] )
    gr = ROOT.TGraph( len(pixels_xy), x, y )
    return gr

def makepoint(a,b):
    x, y = array( 'f' ), array( 'f' )
    pixels_xy = []
    x.append(a)
    y.append(b)
    gr = ROOT.TGraph(1,x,y)
    return gr

def findintersection(fit1,fit2,low,high):
    max_dk = high
    x = 0.0
    y = 0.0
    steps = (high - low)*1000 
    points = [ 0 + k*0.001 for k in range(0,steps+1)]
    for k in points:
         k1 = fit1.Eval(k)
         k2 = fit2.Eval(k)
         if abs(k2-k1)<max_dk:
            max_dk = abs(k2-k1)
            x = k
            y = fit2.Eval(k)
    return float(x),float(y)

def makefits(inputh,inputv,tag):

    #Get points for fits
    nrows, ncols, rgb = inputh.shape
    gr1 = makegraph(inputh)
    gr2 = makegraph(inputv)

    #Get min/max for fit
    n_1    = gr1.GetN()      
    maxx_1 = gr1.GetX()[n_1-1]
    minx_1 = gr1.GetX()[0]
    n_2    = gr2.GetN()      
    maxx_2 = gr2.GetX()[n_2-1]
    minx_2 = gr2.GetX()[0]

    #Fit to graphs with linear function 
    print "[INFO] Starting fits in corner %s"%tag
    gr1.Fit("pol1","RQ","",minx_1-1,maxx_1+1);
    gr2.Fit("pol1","RQ","",minx_2-0.1,maxx_2+0.1);
    gr1.GetFunction("pol1").SetLineColor(ROOT.kRed)
    gr2.GetFunction("pol1").SetLineColor(ROOT.kRed)

    #Find the intersection
    fit1 = gr1.GetFunction("pol1")
    fit2 = gr2.GetFunction("pol1")
    x,y  = findintersection(fit1,fit2,0,ncols)
    print "   - The fitted corner positions is (x,y)=(%.3f,%.3f)"%(x,y)
    gr0  = makepoint(x,y)

    #Create Canvas
    c1 = ROOT.TCanvas("c1", "c1", 1600, 1000)
    c1.SetFrameLineWidth(3)
    c1.SetBottomMargin (0.15)
    c1.SetRightMargin (0.05)
    c1.SetLeftMargin (0.10)

    #First just a frame to give the format to the plot
    histoframe = ROOT.TH2F("","",ncols*10,0,ncols,nrows*10,0,nrows)
    histoframe.GetYaxis().SetTitleSize(0.050)
    histoframe.GetXaxis().SetTitleSize(0.055)
    histoframe.GetYaxis().SetLabelSize(0.05)
    histoframe.GetXaxis().SetLabelSize(0.05)
    histoframe.GetXaxis().SetLabelOffset(0.010)
    histoframe.GetYaxis().SetTitleOffset(1.0)
    histoframe.GetXaxis().SetTitleOffset(1.1)
    histoframe.GetXaxis().SetTitle("X pixel coordinate")
    histoframe.GetYaxis().SetTitle("Y pixel coordinate")
    histoframe.Draw()

    #format to tgraphs
    gr1.SetLineColor( 2 )
    gr1.SetLineWidth( 1 )
    gr1.SetMarkerColor( ROOT.kBlack )
    gr1.SetMarkerSize( 1.5 )
    gr1.SetMarkerStyle( 8 )
    gr1.Draw( 'P same' )
    gr2.SetLineColor( 2 )
    gr2.SetLineWidth( 1 )
    gr2.SetMarkerColor(ROOT.kBlack )
    gr2.SetMarkerSize( 1.5 )
    gr2.SetMarkerStyle( 8 )
    gr2.Draw( 'P same' )
    gr0.SetMarkerColor(ROOT.kBlue)
    gr0.SetMarkerSize( 2)
    gr0.SetMarkerStyle(8)
    gr0.Draw( 'P same' )

    #Draw the intersection point
    pt1 = ROOT.TPaveText(0.1863218,0.886316,0.3045977,0.978947,"brNDC")
    pt1.SetBorderSize(0)
    pt1.SetTextAlign(12)
    pt1.SetTextFont(62)
    pt1.SetTextSize(0.05)
    pt1.SetFillColor(0)
    pt1.SetFillStyle(0)
    pt1.AddText("CMS HGCAL SiPM-on-tile")
    pt1.Draw("SAME")
    
    #Draw legends
    leg_1 = ROOT.TLegend(0.18,0.75,0.45,0.90)
    leg_1.SetNColumns(1)
    leg_1.SetBorderSize(0)
    leg_1.SetTextSize(0.030)
    leg_1.SetTextFont(42)
    leg_1.SetLineColor(1)
    leg_1.SetLineWidth(10)
    leg_1.SetFillColor(0)
    leg_1.SetFillStyle(0)
    leg_1.Draw()
    leg_1.AddEntry(fit1,"Fits", "l")
    leg_1.AddEntry(gr1,"Points", "p")
    leg_1.AddEntry(gr0,"Corner %s: (%.1f,%.1f)"%(tag,x,y), "p")

    c1.Update()
    c1.SaveAs("results_multishot/fit_%s.pdf"%tag)

    return [float(x),float(y)], [fit1.GetParameter(0),fit1.GetParameter(1)], [fit2.GetParameter(0),fit2.GetParameter(1)]


def filterhorizontalandvertical(inputpic,res,tag):
    print "[INFO] Pixel analysis on corner %s"%tag
    imgh = inputpic.copy()
    imgv = inputpic.copy()
    nrows, ncols, rgb = imgh.shape
    w_h_xmin  = 0
    w_h_xmax  = 0
    w_h_ymax  = 0
    w_h_ymin  = 0
    w_v_xmin  = 0
    w_v_xmax  = 0
    w_v_ymax  = 0
    w_v_ymin  = 0
    pixels_xy = []
    #Find all white pixel in x-axis
    for i in range(0, nrows):         
      for j in range(0, ncols):         
           if imgh[i,j,0] != 0: 
              pixels_xy.append([j,i])
    pixels_xy_byx= sorted(pixels_xy, key=lambda x: x[0], reverse=True)
    pixels_xy_byy= sorted(pixels_xy, key=lambda x: x[1], reverse=True)
    #Find key points for edge recognition
    minx = pixels_xy_byx[len(pixels_xy_byx)-1]
    maxx = pixels_xy_byx[0]
    miny = pixels_xy_byy[len(pixels_xy_byy)-1]
    maxy = pixels_xy_byy[0]
    print "   - Coordinates of min x pixel: ",minx
    print "   - Coordinates of max x pixel: ",maxx
    print "   - Coordinates of min y pixel: ",miny
    print "   - Coordinates of max y pixel: ",maxy

    if minx[0]<ncols*0.4 and maxy[1]>nrows*0.6: 
        #Windows for horizontal
        w_h_xmin = minx[0]
        w_h_xmax = maxx[0]
        w_h_ymax = minx[1]+res
        w_h_ymin = minx[1]-res
        #Windows for vertical
        w_v_xmin = maxx[0]-res
        w_v_xmax = maxx[0]+res
        w_v_ymax = maxy[1]
        w_v_ymin = miny[1]
        print "   - It is top right corner (a.k.a. 1)!"
    elif maxx[0]>ncols*0.6 and maxy[1]>nrows*0.6:
        #Windows for horizontal
        w_h_xmin = minx[0]
        w_h_xmax = maxx[0]
        w_h_ymax = maxx[1]+res
        w_h_ymin = maxx[1]-res
        #Windows for vertical
        w_v_xmin = minx[0]-res
        w_v_xmax = minx[0]+res
        w_v_ymax = maxy[1]
        w_v_ymin = miny[1]        
        print "   - It is top left corner (a.k.a. 2)!"
    elif minx[0]<ncols*0.4 and miny[1]<nrows*0.4: 
        #Windows for horizontal
        w_h_xmin = minx[0]
        w_h_xmax = maxx[0]
        w_h_ymax = minx[1]+res
        w_h_ymin = minx[1]-res
        #Windows for vertical
        w_v_xmin = maxx[0]-res
        w_v_xmax = maxx[0]+res
        w_v_ymax = maxy[1]
        w_v_ymin = miny[1]
        print "   - It is bottom right corner (a.k.a. 4)!"
    elif maxx[0]>ncols*0.6 and miny[1]<nrows*0.4:
        #Windows for horizontal
        w_h_xmin = minx[0]
        w_h_xmax = maxx[0]
        w_h_ymax = maxx[1]+res
        w_h_ymin = maxx[1]-res
        #Windows for vertical
        w_v_xmin = minx[0]-res
        w_v_xmax = minx[0]+res
        w_v_ymax = maxy[1]
        w_v_ymin = miny[1]        
        print "   - It is bottom left corner (a.k.a. 3)!"
    else:
        print "[WARNING] Exceptional case: Not covered by algo, so filter will not make sense"

    #Filter only horizontal
    for i in range(0, nrows):         
      for j in range(0, ncols):
               if i>=w_h_ymin and i<=w_h_ymax and j>=w_h_xmin and j<=w_h_xmax: continue
               for k in range(0, rgb): imgh[i,j,k] = 0
    #Filter only vertical
    for i in range(0, nrows):         
      for j in range(0, ncols):      
               if i>=w_v_ymin and i<=w_v_ymax and j>=w_v_xmin and j<=w_v_xmax: continue
               for k in range(0, rgb): imgv[i,j,k] = 0
    cv2.imwrite("plots_multishot/horizontal_%s.png"%tag,imgh)
    cv2.imwrite("plots_multishot/vertical_%s.png"%tag,imgv)

    return imgh,imgv

def drawresults(inputpic,p,fit1,fit2,tag):
    img = inputpic.copy()
    nrows, ncols, rgb = inputpic.shape
    #Draw all fitted corners
    for i in range(0, nrows):         
      for j in range(0, ncols):         
           if math.sqrt(  (j-p[0])*(j-p[0]) + (i-p[1])*(i-p[1]) ) < 3:
               img[i,j,0]=0
               img[i,j,1]=0               
               img[i,j,2]=255            
    cv2.imwrite("results_multishot/results_%s.png"%tag,img)

#make directories
os.system("mkdir plots_multishot")
os.system("mkdir results_multishot")

#original shots and contours (in RGB)
img_original_1   = cv2.imread('images/thickness_1/extrashots/contour_1.png')
img_original_2   = cv2.imread('images/thickness_1/extrashots/contour_2.png')
img_original_3   = cv2.imread('images/thickness_1/extrashots/contour_3.png')
img_original_4   = cv2.imread('images/thickness_1/extrashots/contour_4.png')
img_shot_1       = cv2.imread('images/thickness_1/extrashots/shot_1.png')
img_shot_2       = cv2.imread('images/thickness_1/extrashots/shot_2.png')
img_shot_3       = cv2.imread('images/thickness_1/extrashots/shot_3.png')
img_shot_4       = cv2.imread('images/thickness_1/extrashots/shot_4.png')

#filter circle from contours (in RGB)
img_nocircle_1   = filtercircle(img_original_1,"1")
img_nocircle_2   = filtercircle(img_original_2,"2")
img_nocircle_3   = filtercircle(img_original_3,"3")
img_nocircle_4   = filtercircle(img_original_4,"4")

#filter horizontal and vertical lines (in RGB)
img_horizontal_1,img_vertical_1 = filterhorizontalandvertical(img_nocircle_1,10,"1")
img_horizontal_2,img_vertical_2 = filterhorizontalandvertical(img_nocircle_2,10,"2")
img_horizontal_3,img_vertical_3 = filterhorizontalandvertical(img_nocircle_3,10,"3")
img_horizontal_4,img_vertical_4 = filterhorizontalandvertical(img_nocircle_4,10,"4")

#fit points
xy_1,fit1_1,fit2_1 = makefits(img_horizontal_1,img_vertical_1,"1")
xy_2,fit1_2,fit2_2 = makefits(img_horizontal_2,img_vertical_2,"2")
xy_3,fit1_3,fit2_3 = makefits(img_horizontal_3,img_vertical_3,"3")
xy_4,fit1_4,fit2_4 = makefits(img_horizontal_4,img_vertical_4,"4")

#Draw fitted corners in shots and contours
drawresults(img_original_1,xy_1,fit1_1,fit2_1,"contour_1")
drawresults(img_original_2,xy_2,fit1_2,fit2_2,"contour_2")
drawresults(img_original_3,xy_3,fit1_3,fit2_3,"contour_3")
drawresults(img_original_4,xy_4,fit1_4,fit2_4,"contour_4")
drawresults(img_shot_1    ,xy_1,fit1_1,fit2_1,"shot_1")
drawresults(img_shot_2    ,xy_2,fit1_2,fit2_2,"shot_2")
drawresults(img_shot_3    ,xy_3,fit1_3,fit2_3,"shot_3")
drawresults(img_shot_4    ,xy_4,fit1_4,fit2_4,"shot_4")

#Convert fitted corners in mm (w.r.t to camera center)
c_1 = convertcameracenter(xy_1,img_original_1.shape[1],img_original_1.shape[0])
c_2 = convertcameracenter(xy_2,img_original_2.shape[1],img_original_2.shape[0])
c_3 = convertcameracenter(xy_3,img_original_3.shape[1],img_original_3.shape[0])
c_4 = convertcameracenter(xy_4,img_original_4.shape[1],img_original_4.shape[0])

#Convert to global coordinates
#Noozle position in extra shots (for thickness 1 folder)
n_1 = [-5.837610, 95.795117]
n_2 = [27.162390, 95.795117] 
n_3 = [27.162390, 129.795117]
n_4 = [-5.837610, 129.795117]

##Noozle position in extra shots (for thickness 2 folder)
#n_1 = [-5.793762, 95.692207]
#n_2 = [27.206238, 95.692207] 
#n_3 = [27.206238, 129.692207]
#n_4 = [-5.793762, 129.692207]

p_1 = [  n_3[0] + c_1[0],  n_3[1] + c_1[1] ]
p_2 = [  n_4[0] + c_2[0],  n_4[1] + c_2[1] ]
p_3 = [  n_1[0] + c_3[0],  n_1[1] + c_3[1] ]
p_4 = [  n_2[0] + c_4[0],  n_2[1] + c_4[1] ]

#Print out the results
print "[RESULTS] Corner points in global coordinate system (mm):"
print "    - Top right    (x,y): [%.3f,%.3f]"%(p_1[0],p_1[1])
print "    - Top left     (x,y): [%.3f,%.3f]"%(p_2[0],p_2[1])
print "    - Bottom left  (x,y): [%.3f,%.3f]"%(p_3[0],p_3[1])
print "    - Bottom right (x,y): [%.3f,%.3f]"%(p_4[0],p_4[1])
top = math.sqrt(  (p_1[0]-p_2[0])*(p_1[0]-p_2[0]) + (p_1[1]-p_2[1])*(p_1[1]-p_2[1]) )
lef = math.sqrt(  (p_2[0]-p_3[0])*(p_2[0]-p_3[0]) + (p_2[1]-p_3[1])*(p_2[1]-p_3[1]) )
bot = math.sqrt(  (p_3[0]-p_4[0])*(p_3[0]-p_4[0]) + (p_3[1]-p_4[1])*(p_3[1]-p_4[1]) )
rig = math.sqrt(  (p_4[0]-p_1[0])*(p_4[0]-p_1[0]) + (p_4[1]-p_1[1])*(p_4[1]-p_1[1]) )
print "[RESULTS] Trapezoid sizes (mm):"
print "    - Top      (33.83): %.3f"%top
print "    - Bottom   (33.05): %.3f"%bot
print "    - Left     (33.81): %.3f"%lef
print "    - Right    (33.81): %.3f"%rig

