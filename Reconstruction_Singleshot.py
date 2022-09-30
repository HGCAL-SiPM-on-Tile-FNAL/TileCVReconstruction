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
           if distance < 0.85*r: continue
           for k in range(0, rgb): img[i,j,k] = 0
    cv2.imwrite("plots_singleshot/test_circle_%s.png"%tag,img)
    return img

def converttomm(coordinates,ncols,nrows):
    k = 0.04589621286 #mm per pixel
    print "[INFO] Converting corner information to mm"
    cx = []
    cy = []
    for i in range(0,len(coordinates)):
            cx.append( round(coordinates[i][0]*k,3)) 
            cy.append( round(coordinates[i][1]*k,3)) 
            print "   - Corner %i [mm]: (%.3f,%.3f)"%(i+1,cx[i],cy[i])
    return [cx[0],cy[0]],[cx[1],cy[1]],[cx[2],cy[2]],[cx[3],cy[3]],

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
    #Preselect pixels
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

def makefits(inputt,inputb,inputl,inputr,tag):
    #Get points for fits
    nrows, ncols, rgb = inputt.shape
    grt = makegraph(inputt)
    grb = makegraph(inputb)
    grl = makegraph(inputl)
    grr = makegraph(inputr)
    #Get min/max for fit
    n_t    = grt.GetN()      
    maxx_t = grt.GetX()[n_t-1]
    minx_t = grt.GetX()[0]
    n_b    = grb.GetN()      
    maxx_b = grb.GetX()[n_b-1]
    minx_b = grb.GetX()[0]
    n_l    = grl.GetN()      
    maxx_l = grl.GetX()[n_l-1]
    minx_l = grl.GetX()[0]
    n_r    = grr.GetN()      
    maxx_r = grr.GetX()[n_r-1]
    minx_r = grr.GetX()[0]

    #Fit to graphs with linear function 
    print "[INFO] Starting fits in all corners"
    grt.Fit("pol1","RQ","",minx_t-1,maxx_t+1);
    grb.Fit("pol1","RQ","",minx_b-1,maxx_b+1);
    grl.Fit("pol1","RQ","",minx_l-0.5,maxx_l+0.1);
    grr.Fit("pol1","RQ","",minx_r-0.1,maxx_r+0.5);
    grt.GetFunction("pol1").SetLineColor(ROOT.kRed)
    grb.GetFunction("pol1").SetLineColor(ROOT.kRed)
    grl.GetFunction("pol1").SetLineColor(ROOT.kRed)
    grr.GetFunction("pol1").SetLineColor(ROOT.kRed)

    #Find the intersection
    fitt = grt.GetFunction("pol1");
    fitb = grb.GetFunction("pol1");
    fitl = grl.GetFunction("pol1");
    fitr = grr.GetFunction("pol1");
    x_1,y_1  = findintersection(fitt,fitr,0,ncols)
    x_2,y_2  = findintersection(fitt,fitl,0,ncols)
    x_3,y_3  = findintersection(fitb,fitl,0,ncols)
    x_4,y_4  = findintersection(fitb,fitr,0,ncols)
    print "   - The fitted corner positions:"
    print "           -  Corner 1 (x,y): [%.3f,%.3f]"%(x_1,y_1)
    print "           -  Corner 2 (x,y): [%.3f,%.3f]"%(x_2,y_2)
    print "           -  Corner 3 (x,y): [%.3f,%.3f]"%(x_3,y_3)
    print "           -  Corner 4 (x,y): [%.3f,%.3f]"%(x_4,y_4)
    gr1  = makepoint(x_1,y_1)
    gr2  = makepoint(x_2,y_2)
    gr3  = makepoint(x_3,y_3)
    gr4  = makepoint(x_4,y_4)
 
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
    grt.SetLineColor( 2 )
    grt.SetLineWidth( 1 )
    grt.SetMarkerColor( ROOT.kBlack )
    grt.SetMarkerSize( 1.5 )
    grt.SetMarkerStyle( 8 )
    grt.Draw( 'P same' )
    grb.SetLineColor( 2 )
    grb.SetLineWidth( 1 )
    grb.SetMarkerColor( ROOT.kBlack )
    grb.SetMarkerSize( 1.5 )
    grb.SetMarkerStyle( 8 )
    grb.Draw( 'P same' )
    grl.SetLineColor( 2 )
    grl.SetLineWidth( 1 )
    grl.SetMarkerColor( ROOT.kBlack )
    grl.SetMarkerSize( 1.5 )
    grl.SetMarkerStyle( 8 )
    grl.Draw( 'P same' )
    grr.SetLineColor( 2 )
    grr.SetLineWidth( 1 )
    grr.SetMarkerColor( ROOT.kBlack )
    grr.SetMarkerSize( 1.5 )
    grr.SetMarkerStyle( 8 )
    grr.Draw( 'P same' )

    gr1.SetMarkerColor(ROOT.kBlue)
    gr1.SetMarkerSize( 2)
    gr1.SetMarkerStyle(8)
    gr1.Draw( 'P same' )
    gr2.SetMarkerColor(ROOT.kBlue)
    gr2.SetMarkerSize( 2)
    gr2.SetMarkerStyle(8)
    gr2.Draw( 'P same' )
    gr3.SetMarkerColor(ROOT.kBlue)
    gr3.SetMarkerSize( 2)
    gr3.SetMarkerStyle(8)
    gr3.Draw( 'P same' )
    gr4.SetMarkerColor(ROOT.kBlue)
    gr4.SetMarkerSize( 2)
    gr4.SetMarkerStyle(8)
    gr4.Draw( 'P same' )

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
    leg_1 = ROOT.TLegend(0.12,0.3,0.25,0.70)
    leg_1.SetNColumns(1)
    leg_1.SetBorderSize(0)
    leg_1.SetTextSize(0.030)
    leg_1.SetTextFont(42)
    leg_1.SetLineColor(1)
    leg_1.SetLineWidth(10)
    leg_1.SetFillColor(0)
    leg_1.SetFillStyle(0)
    leg_1.Draw()
    leg_1.AddEntry(fitt,"Fits", "l")
    leg_1.AddEntry(grt,"Points", "p")
    leg_1.AddEntry(gr1,"Corner 1: (%.1f,%.1f)"%(x_1,y_1), "p")
    leg_1.AddEntry(gr2,"Corner 2: (%.1f,%.1f)"%(x_2,y_2), "p")
    leg_1.AddEntry(gr3,"Corner 3: (%.1f,%.1f)"%(x_3,y_3), "p")
    leg_1.AddEntry(gr4,"Corner 4: (%.1f,%.1f)"%(x_4,y_4), "p")

    c1.Update()
    c1.SaveAs("results_singleshot/fit_%s.pdf"%tag)

    return [[x_1,y_1],[x_2,y_2],[x_3,y_3],[x_4,y_4]], [fitt,fitb,fitr,fitl]

def filtersides(inputpic,res,tag):
    print "[INFO] Pixel analysis all corners"
    imgt = inputpic.copy()
    imgb = inputpic.copy()
    imgl = inputpic.copy()
    imgr = inputpic.copy()

    nrows, ncols, rgb = inputpic.shape
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
           if inputpic[i,j,0] != 0: 
              pixels_xy.append([j,i])
    pixels_xy_byx= sorted(pixels_xy, key=lambda x: x[0], reverse=True)
    pixels_xy_byy= sorted(pixels_xy, key=lambda x: x[1], reverse=True)
    #Find key points
    minx = pixels_xy_byx[len(pixels_xy_byx)-1]
    maxx = pixels_xy_byx[0]
    miny = pixels_xy_byy[len(pixels_xy_byy)-1]
    maxy = pixels_xy_byy[0]
    print "   - Coordinates of min x pixel: ",minx
    print "   - Coordinates of max x pixel: ",maxx
    print "   - Coordinates of min y pixel: ",miny
    print "   - Coordinates of max y pixel: ",maxy

    #Windows for top
    w_t_xmin = minx[0]
    w_t_xmax = maxx[0]
    w_t_ymax = miny[1]+res
    w_t_ymin = miny[1]-res
    #Windows for bot
    w_b_xmin = minx[0]
    w_b_xmax = maxx[0]
    w_b_ymax = maxy[1]+res
    w_b_ymin = maxy[1]-res
    #Windows for lef
    w_l_xmin = minx[0]-res
    w_l_xmax = minx[0]+res
    w_l_ymax = maxy[1]
    w_l_ymin = miny[1] 
    #Windows for lef
    w_r_xmin = maxx[0]-res
    w_r_xmax = maxx[0]+res
    w_r_ymax = maxy[1]
    w_r_ymin = miny[1] 
    #Filter only top
    for i in range(0, nrows):         
      for j in range(0, ncols):
               if i>=w_t_ymin and i<=w_t_ymax and j>=w_t_xmin and j<=w_t_xmax: continue
               for k in range(0, rgb): imgt[i,j,k] = 0
    #Filter only bot
    for i in range(0, nrows):         
      for j in range(0, ncols):
               if i>=w_b_ymin and i<=w_b_ymax and j>=w_b_xmin and j<=w_b_xmax: continue
               for k in range(0, rgb): imgb[i,j,k] = 0
    #Filter lef top
    for i in range(0, nrows):         
      for j in range(0, ncols):
               if i>=w_l_ymin and i<=w_l_ymax and j>=w_l_xmin and j<=w_l_xmax: continue
               for k in range(0, rgb): imgl[i,j,k] = 0
    #Filter rig bot
    for i in range(0, nrows):         
      for j in range(0, ncols):
               if i>=w_r_ymin and i<=w_r_ymax and j>=w_r_xmin and j<=w_r_xmax: continue
               for k in range(0, rgb): imgr[i,j,k] = 0
    cv2.imwrite("plots_singleshot/top_%s.png"%tag,imgt)
    cv2.imwrite("plots_singleshot/bot_%s.png"%tag,imgb)
    cv2.imwrite("plots_singleshot/lef_%s.png"%tag,imgl)
    cv2.imwrite("plots_singleshot/rig_%s.png"%tag,imgr)
    return imgt,imgb,imgl,imgr

def drawresults(inputpic,xy_1234,fit_1234,tag):
    img = inputpic.copy()
    nrows, ncols, rgb = inputpic.shape
    p1 = xy_1234[0]
    p2 = xy_1234[1]
    p3 = xy_1234[2]
    p4 = xy_1234[3]
    #Draw all fitted corners
    for i in range(0, nrows):         
      for j in range(0, ncols):         
           if math.sqrt(  (j-p1[0])*(j-p1[0]) + (i-p1[1])*(i-p1[1]) ) < 3: 
               img[i,j,0]=0
               img[i,j,1]=0               
               img[i,j,2]=255 
           if math.sqrt(  (j-p2[0])*(j-p2[0]) + (i-p2[1])*(i-p2[1]) ) < 3:
               img[i,j,0]=0
               img[i,j,1]=0               
               img[i,j,2]=255
           if math.sqrt(  (j-p3[0])*(j-p3[0]) + (i-p3[1])*(i-p3[1]) ) < 3:
               img[i,j,0]=0
               img[i,j,1]=0               
               img[i,j,2]=255
           if math.sqrt(  (j-p4[0])*(j-p4[0]) + (i-p4[1])*(i-p4[1]) ) < 3:
               img[i,j,0]=0
               img[i,j,1]=0               
               img[i,j,2]=255             
    cv2.imwrite("results_singleshot/results_%s.png"%tag,img)

#make directories
os.system("mkdir plots_singleshot")
os.system("mkdir results_singleshot")

#original shot and contour (in RGB)
img_contour  = cv2.imread('images/thickness_1/singleshot/contour.png')
img_shot     = cv2.imread('images/thickness_1/singleshot/shot.png')

#filter circle from contours (in RGB)
img_nocircle  = filtercircle(img_contour,"singleshot")

#filter horizontal and vertical lines (in RGB)
img_top,img_bottom,img_left,img_right = filtersides(img_nocircle,10,"singleshot")

#fit points
xy_1234,fit_1234 = makefits(img_top,img_bottom,img_left,img_right,"singleshot")

#Draw fitted corners in shots and contours
drawresults(img_contour,xy_1234,fit_1234,"contour")
drawresults(img_shot   ,xy_1234,fit_1234,"shot")

#Convert corners to mm 
p_1,p_2,p_3,p_4 = converttomm(xy_1234,img_contour.shape[1],img_contour.shape[0])
top = math.sqrt(  (p_1[0]-p_2[0])*(p_1[0]-p_2[0]) + (p_1[1]-p_2[1])*(p_1[1]-p_2[1]) )
lef = math.sqrt(  (p_2[0]-p_3[0])*(p_2[0]-p_3[0]) + (p_2[1]-p_3[1])*(p_2[1]-p_3[1]) )
bot = math.sqrt(  (p_3[0]-p_4[0])*(p_3[0]-p_4[0]) + (p_3[1]-p_4[1])*(p_3[1]-p_4[1]) )
rig = math.sqrt(  (p_4[0]-p_1[0])*(p_4[0]-p_1[0]) + (p_4[1]-p_1[1])*(p_4[1]-p_1[1]) )

#Print out the results
print "[RESULTS] Trapezoid sizes (mm):"
print "    - Top      (33.83): %.3f"%top
print "    - Bottom   (33.05): %.3f"%bot
print "    - Left     (33.81): %.3f"%lef
print "    - Right    (33.81): %.3f"%rig
