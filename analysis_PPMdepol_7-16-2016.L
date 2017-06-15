/* Analysis for PPMdepol2009 Runs */
/* A. T. Holley */

/* Created: July 16, 2016
   Modified: August 3, 2016           */

/*
Changelog:
[DAD]: Added plots in Runsummary() to show all field values. also added in semi working version of holdcount() to list field values and counts of hold interval.
[DAD]: Fixed holdcount(), works properly

[DAD]: Added holdfitamp() to plot hold amplitude versus run number, color codes for field value and lines seperate guide types.

[DAD]: added unldfitamp() to plot unload amplitude versus run number, color codes for field value and lines seperate guide types

[DAD]: added holdfittau() to plot hold tau versus run number, color codes for field value and lines seperate guide types

[DAD]: added Holdfitbg() and unldfitbg() to plot the backgrounds versus run number, same color coding and guide lines

[DAD]: split runsummary() into seperate functions and will use run summary to select which plots to view along with the load counts.

[DAD]: added in more instances of the useg() function in drawallusegs() in order to draw multiple sets of guide types on one plot

[DAD]:added in guidesep() and sepselect() which added a function which filters by guide type and field and plots norm unld vs. hold interval, and makes a menu for using it
*/

#include <stdint.h>
#include <math.h>
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "PLOT_PPMdepol_7-16-2016.h"
#include "getobs_PPMdepol2.L"
#include "getobs.h"

#define ODIR "/home/damien/Desktop/PPMdepol"
#define RUMAX 100

int runs[];

//Plots
void runsummary(); //Plots of various quantities and ratios over the entire set of runs.
void getloadcount(); //Plots of the load counts versus run number
void getholdcount(); //Plots of the normalized hold counts versus run number
void getunldcount(); //plots of normalized unload counts.
void getholdfitamp(); //Plots of the normalized hold amp counts versus run number
void getunldfitamp(); //Plots of the normalized unload amp counts versus run number
void getholdfittau(); //Plots of normalized hold tau counts versus run number
void getunldfittau(); //Plots of normalized unld tau counts versus run number
void getholdfitbg(); //plots of normalized hold background counts versus run number
void getunldfitbg(); //plots of normalized unld background counts versus run number
void holdint(); //plots hold interval versus normalized unload
void guidesep(); //plots norm unload counts versus hold time seperating by guide and field
double dpb(int); //calcultes depolarization per bounce
void plotObs(); //draws a plot of different observables v. run number
void plotdpb(); // calculates dpb for all runs and plots it by guide type and field strength
void vlinedraw(double, double, int); 
double getTauhat(int);// gets tauhat for each guide type


//Retrieving extracted observables lemmings
TTree* xdload(const char*,const char*);

//General lemmings
void guidelabel(TCanvas*, const char*,double, double);//makes the labels for the guide types in the plots
void edivide(double,double,double,double,double*,double*);
void vlinedraw(double,double,int);
void wsum(int,double*,double*,double*,double*);
int RunUnits(const char*,int[]); //Get runs in a run unit for a particular run type.
void useg(const char*,int); //Set off a particular run unit.
void drawallusegs();
void sepselect();


int _xdtinit = 0; //Global varible to keep from re-initializing the data tree.

  /*-------------------------------Plot List--------------------------------*/
  /*plt1 = load counts
   *plt2 = unload counts
   *plt3 = hold counts
   *plt4 = hold amplitude
   *plt5 = unload amplitude
   *plt6 = hold tau
   *plt7 = unload tau
   *plt8 = hold background
   *plt9 = unload background
  
  /* ------------------------- Data Analysis Plots ----------------------- */
void runsummary() {
  int pnum;
  
  cout << "Which Observable Plot would you like to view?\n\n1. Normalized Hold\n2. Normalized Unload\n3. Normalized Hold Amp\n4. Normalized Unload Amp\n5. Normalized Hold Tau\n6. Normalized Unload Tau\n7. Normalized Hold Background\n8. Normalized Unload Background\n9. Hold Interval v. Normalized Unload\n10. Normalized unload v. Hold interval(filtered)" << endl;
  
  cin >> pnum;
  
  if(pnum == 1) {
     getholdcount();
     getloadcount();
    return;
  }
  
  if(pnum == 2) {
     getunldcount();
     getloadcount();
    return;
  }
  
  if(pnum == 3) {
     getholdfitamp();
     getloadcount();
    return;
  }
  
  if(pnum == 4) {
     getunldfitamp();
     getloadcount();
    return;
  }
  
 if(pnum == 5) {
     getholdfittau();
     getloadcount();
    return;
  }
  
  if(pnum == 6) {
     getunldfittau();
     getloadcount();
    return;
  }
  
  if(pnum == 7) {
     getholdfitbg();
     getloadcount();
    return;
  }
  
  if(pnum == 8) {
     getunldfitbg();
     getloadcount();
    return;
  }
  
  if(pnum == 9) {
     holdint(); 
     getloadcount();
    return;
  }
  
   if(pnum == 10) {
     sepselect();
    return;
  }
  
  return;
}

void getloadcount() {
    int i;
    
   char name[500];
   double normload,enormload;
   
   
  //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //Initialize plots
  PLOT plt1a;  plt1a.init();
  PLOT plt1b;  plt1b.init();
  PLOT plt1c;  plt1c.init();
  PLOT plt1d;  plt1d.init();
  PLOT plt1e;  plt1e.init();
  PLOT plt1f;  plt1f.init();
  PLOT plt1g;  plt1g.init();
  PLOT plt1h;  plt1h.init();
  
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    normload = xd.loadcounts;
    if(xd.Tfield == 250) plt1a.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 100) plt1b.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 260) plt1c.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 175) plt1d.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 50) plt1e.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 30) plt1f.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 25) plt1g.loadpnt(xd.rnum,0,normload,enormload);
    if(xd.Tfield == 10) plt1h.loadpnt(xd.rnum,0,normload,enormload);
  }  
  
  //Draw load counts plot.
  TCanvas *can1 = new TCanvas("can1","Load Counts",735,5,700,500);
  can1->SetFrameLineColor(0);
  can1->cd();
  plt1a.make("plt1a");
  plt1a.haxis->SetTitle("Load Counts");
  plt1a.haxis->GetXaxis()->SetTitle("run#");
  plt1a.haxis->GetYaxis()->SetTitle("UCN Counts");
  plt1a.haxis->GetXaxis()->SetRangeUser(0,1000);
  plt1a.haxis->GetYaxis()->SetRangeUser(0,40000);
  plt1a.haxis->Draw();
  plt1a.plot->SetMarkerColor(4);
  plt1a.plot->SetLineColor(4);
  plt1a.plot->SetMarkerStyle(20);
  plt1a.plot->Draw("PSame");
  
  plt1b.make("plt1b");
  plt1b.plot->SetMarkerColor(2);
  plt1b.plot->SetLineColor(2);
  plt1b.plot->SetMarkerStyle(20);
  plt1b.plot->SetMarkerSize(1.3);
  plt1b.plot->Draw("PSame");
  
  plt1c.make("plt1c");
  plt1c.plot->SetMarkerColor(1);
  plt1c.plot->SetLineColor(1);
  plt1c.plot->SetMarkerStyle(20);
  plt1c.plot->SetMarkerSize(1.3);
  plt1c.plot->Draw("PSame");
  
  plt1d.make("plt1d");
  plt1d.plot->SetMarkerColor(3);
  plt1d.plot->SetLineColor(3);
  plt1d.plot->SetMarkerStyle(20);
  plt1d.plot->SetMarkerSize(1.3);
  plt1d.plot->Draw("PSame");
  
  plt1e.make("plt1e");
  plt1e.plot->SetMarkerColor(5);
  plt1e.plot->SetLineColor(5);
  plt1e.plot->SetMarkerStyle(20);
  plt1e.plot->SetMarkerSize(1.3);
  plt1e.plot->Draw("PSame");
  
  plt1f.make("plt1f");
  plt1f.plot->SetMarkerColor(6);
  plt1f.plot->SetLineColor(6);
  plt1f.plot->SetMarkerStyle(20);
  plt1f.plot->SetMarkerSize(1.3);
  plt1f.plot->Draw("PSame");
  
  plt1g.make("plt1g");
  plt1g.plot->SetMarkerColor(7);
  plt1g.plot->SetLineColor(7);
  plt1g.plot->SetMarkerStyle(20);
  plt1g.plot->SetMarkerSize(1.3);
  plt1g.plot->Draw("PSame");
  
  plt1h.make("plt1h");
  plt1h.plot->SetMarkerColor(14);
  plt1h.plot->SetLineColor(14);
  plt1h.plot->SetMarkerStyle(20);
  plt1h.plot->SetMarkerSize(1.3);
  plt1h.plot->Draw("PSame");
  
  TLegend *leg1 = new TLegend(0.125,0.635021,0.317529,0.92616);
  leg1->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg1->AddEntry(plt1a.plot,"B = 250 G","LP");
  leg1->AddEntry(plt1b.plot,"B = 100 G","LP");
  leg1->AddEntry(plt1c.plot,"B = 260 G", "LP");
  leg1->AddEntry(plt1d.plot,"B = 175 G", "LP");
  leg1->AddEntry(plt1e.plot,"B = 50 G", "LP");
  leg1->AddEntry(plt1f.plot,"B = 30 G", "LP");
  leg1->AddEntry(plt1g.plot,"B = 25 G", "LP");
  leg1->AddEntry(plt1h.plot,"B = 10 G", "LP");
  leg1->Draw("Same");

  drawallusegs();
  
 return; 
}

void getholdcount() {
    int i;
   char name[500];
   double normhold,enormhold;

  //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt3a;  plt3a.init(); //Normalized hold counts as a function of run number, B=250G.
  PLOT plt3b;  plt3b.init(); //Normalized hold counts as a function of run number, B=100G.
  PLOT plt3c;  plt3c.init(); //Normalized hold counts as a function of run number, B=260G.
  PLOT plt3d;  plt3d.init(); //Normalized hold counts as a function of run number, B=175G.
  PLOT plt3e;  plt3e.init(); //Normalized hold counts as a function of run number, B=50G.
  PLOT plt3f;  plt3f.init(); //Normalized hold counts as a function of run number, B=30G.
  PLOT plt3g;  plt3g.init(); //Normalized hold counts as a function of run number, B=25G.
  PLOT plt3h;  plt3h.init(); //Normalized hold counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.holdcounts,xd.eholdcounts,xd.loadcounts,xd.eloadcounts,&normhold,&enormhold);
    if(xd.Tfield == 250) plt3a.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 100) plt3b.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 260) plt3c.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 175) plt3d.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 50) plt3e.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 30) plt3f.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 25) plt3g.loadpnt(xd.rnum,0,normhold,enormhold);
    if(xd.Tfield == 10) plt3h.loadpnt(xd.rnum,0,normhold,enormhold);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can3 = new TCanvas("can3","Normalized Hold",735,5,700,500);
  can3->SetFrameLineColor(0);
  can3->cd();
  
  plt3a.make("plt3a");
  plt3a.haxis->SetTitle("Total Normalized Hold Counts");
  plt3a.haxis->GetXaxis()->SetTitle("run#");
  plt3a.haxis->GetYaxis()->SetTitle("hold counts/load counts");
  plt3a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt3a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt3a.haxis->Draw();
  plt3a.plot->SetMarkerColor(4);
  plt3a.plot->SetLineColor(4);
  plt3a.plot->SetMarkerStyle(20);
  plt3a.plot->Draw("PSame");
  
  plt3b.make("plt3b");
  plt3b.plot->SetMarkerColor(2);
  plt3b.plot->SetLineColor(2);
  plt3b.plot->SetMarkerStyle(20);
  plt3b.plot->SetMarkerSize(1.3);
  plt3b.plot->Draw("PSame");
  
  plt3c.make("plt3c");
  plt3c.plot->SetMarkerColor(1);
  plt3c.plot->SetLineColor(1);
  plt3c.plot->SetMarkerStyle(20);
  plt3c.plot->SetMarkerSize(1.3);
  plt3c.plot->Draw("PSame");
  
  plt3d.make("plt3d");
  plt3d.plot->SetMarkerColor(3);
  plt3d.plot->SetLineColor(3);
  plt3d.plot->SetMarkerStyle(20);
  plt3d.plot->SetMarkerSize(1.3);
  plt3d.plot->Draw("PSame");
  
  plt3e.make("plt3e");
  plt3e.plot->SetMarkerColor(5);
  plt3e.plot->SetLineColor(5);
  plt3e.plot->SetMarkerStyle(20);
  plt3e.plot->SetMarkerSize(1.3);
  plt3e.plot->Draw("PSame");
  
  plt3f.make("plt3f");
  plt3f.plot->SetMarkerColor(6);
  plt3f.plot->SetLineColor(6);
  plt3f.plot->SetMarkerStyle(20);
  plt3f.plot->SetMarkerSize(1.3);
  plt3f.plot->Draw("PSame");
  
  plt3g.make("plt3g");
  plt3g.plot->SetMarkerColor(7);
  plt3g.plot->SetLineColor(7);
  plt3g.plot->SetMarkerStyle(20);
  plt3g.plot->SetMarkerSize(1.3);
  plt3g.plot->Draw("PSame");
  
  plt3h.make("plt3h");
  plt3h.plot->SetMarkerColor(14);
  plt3h.plot->SetLineColor(14);
  plt3h.plot->SetMarkerStyle(20);
  plt3h.plot->SetMarkerSize(1.3);
  plt3h.plot->Draw("PSame");
  
  TLegend *leg3 = new TLegend(0.899,0.364,0.997,0.917);
  leg3->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg3->AddEntry(plt3a.plot,"B = 250 G","LP");
  leg3->AddEntry(plt3b.plot,"B = 100 G","LP");
  leg3->AddEntry(plt3c.plot,"B = 260 G", "LP");
  leg3->AddEntry(plt3d.plot,"B = 175 G", "LP");
  leg3->AddEntry(plt3e.plot,"B = 50 G", "LP");
  leg3->AddEntry(plt3f.plot,"B = 30 G", "LP");
  leg3->AddEntry(plt3g.plot,"B = 25 G", "LP");
  leg3->AddEntry(plt3h.plot,"B = 10 G", "LP");
  leg3->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can3, "Non Coated Copper",.093, .004);
  guidelabel(can3,"DLC1",.254,.004);
  guidelabel(can3,"EP",.484,.004);
  guidelabel(can3,"DLC2",.613,.004);
  guidelabel(can3,"SS",.699,.004);
  guidelabel(can3,"MP",.788,.004);
  guidelabel(can3,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can3, "DLC PSN 1a",.198, .850);
  guidelabel(can3,"DLC PSN 1b",.308,.850);
  guidelabel(can3,"EP PSN 1",.497,.850);
  guidelabel(can3,"DLC PSN 2",.591,.850);
  guidelabel(can3,"MP PSN",.777,.850);
  guidelabel(can3,"EP PSN 2",.836,.907);
  
  
  
 return; 
}

void getunldcount() {
  int i;
   char name[500];
   double normunld,enormunld;
   
   //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;

  
  PLOT plt2a;  plt2a.init(); //Normalized unload counts as a function of run number, B=250G.
  PLOT plt2b;  plt2b.init(); //Normalized unload counts as a function of run number, B=100G.
  PLOT plt2c;  plt2c.init(); //Normalized unload counts as a function of run number, B=260G.
  PLOT plt2d;  plt2d.init(); //Normalized unload counts as a function of run number, B=175G.
  PLOT plt2e;  plt2e.init(); //Normalized unload counts as a function of run number, B=50G.
  PLOT plt2f;  plt2f.init(); //Normalized unload counts as a function of run number, B=30G.
  PLOT plt2g;  plt2g.init(); //Normalized unload counts as a function of run number, B=25G.
  PLOT plt2h;  plt2h.init(); //Normalized unload counts as a function of run number, B=10G.
  
  
  //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.unldcounts,xd.eunldcounts,xd.loadcounts,xd.eloadcounts,&normunld,&enormunld);
    if(xd.Tfield == 250) plt2a.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 100) plt2b.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 260) plt2c.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 175) plt2d.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 50) plt2e.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 30) plt2f.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 25) plt2g.loadpnt(xd.rnum,0,normunld,enormunld);
    if(xd.Tfield == 10) plt2h.loadpnt(xd.rnum,0,normunld,enormunld);
    
  }

  

  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized Unload",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt2a.make("plt2a");
  plt2a.haxis->SetTitle("Total Normalized Unload Counts");
  plt2a.haxis->GetXaxis()->SetTitle("run#");
  plt2a.haxis->GetYaxis()->SetTitle("unload counts/load counts");
  plt2a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt2a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt2a.haxis->Draw();
  plt2a.plot->SetMarkerColor(4);
  plt2a.plot->SetLineColor(4);
  plt2a.plot->SetMarkerStyle(20);
  plt2a.plot->Draw("PSame");
  
  plt2b.make("plt2b");
  plt2b.plot->SetMarkerColor(2);
  plt2b.plot->SetLineColor(2);
  plt2b.plot->SetMarkerStyle(20);
  plt2b.plot->SetMarkerSize(1.3);
  plt2b.plot->Draw("PSame");
  
  plt2c.make("plt2c");
  plt2c.plot->SetMarkerColor(1);
  plt2c.plot->SetLineColor(1);
  plt2c.plot->SetMarkerStyle(20);
  plt2c.plot->SetMarkerSize(1.3);
  plt2c.plot->Draw("PSame");
  
  plt2d.make("plt2d");
  plt2d.plot->SetMarkerColor(3);
  plt2d.plot->SetLineColor(3);
  plt2d.plot->SetMarkerStyle(20);
  plt2d.plot->SetMarkerSize(1.3);
  plt2d.plot->Draw("PSame");
  
  plt2e.make("plt2e");
  plt2e.plot->SetMarkerColor(5);
  plt2e.plot->SetLineColor(5);
  plt2e.plot->SetMarkerStyle(20);
  plt2e.plot->SetMarkerSize(1.3);
  plt2e.plot->Draw("PSame");
  
  plt2f.make("plt2f");
  plt2f.plot->SetMarkerColor(6);
  plt2f.plot->SetLineColor(6);
  plt2f.plot->SetMarkerStyle(20);
  plt2f.plot->SetMarkerSize(1.3);
  plt2f.plot->Draw("PSame");
  
  plt2g.make("plt2g");
  plt2g.plot->SetMarkerColor(7);
  plt2g.plot->SetLineColor(7);
  plt2g.plot->SetMarkerStyle(20);
  plt2g.plot->SetMarkerSize(1.3);
  plt2g.plot->Draw("PSame");
  
  plt2h.make("plt2h");
  plt2h.plot->SetMarkerColor(14);
  plt2h.plot->SetLineColor(14);
  plt2h.plot->SetMarkerStyle(20);
  plt2h.plot->SetMarkerSize(1.3);
  plt2h.plot->Draw("PSame");
  
  TLegend *leg2 = new TLegend(0.899,0.364,0.997,0.917);
  leg2->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg2->AddEntry(plt2a.plot,"B = 250 G","LP");
  leg2->AddEntry(plt2b.plot,"B = 100 G","LP");
  leg2->AddEntry(plt2c.plot,"B = 260 G", "LP");
  leg2->AddEntry(plt2d.plot,"B = 175 G", "LP");
  leg2->AddEntry(plt2e.plot,"B = 50 G", "LP");
  leg2->AddEntry(plt2f.plot,"B = 30 G", "LP");
  leg2->AddEntry(plt2g.plot,"B = 25 G", "LP");
  leg2->AddEntry(plt2h.plot,"B = 10 G", "LP");
  leg2->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  return;
}

void getholdfitamp() {
   int i;
   char name[500];
   double normholdamp,enormholdamp;

  //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt4a;  plt4a.init(); //Normalized hold amp counts as a function of run number, B=250G.
  PLOT plt4b;  plt4b.init(); //Normalized hold amp counts as a function of run number, B=100G.
  PLOT plt4c;  plt4c.init(); //Normalized hold amp counts as a function of run number, B=260G.
  PLOT plt4d;  plt4d.init(); //Normalized hold amp counts as a function of run number, B=175G.
  PLOT plt4e;  plt4e.init(); //Normalized hold amp counts as a function of run number, B=50G.
  PLOT plt4f;  plt4f.init(); //Normalized hold amp counts as a function of run number, B=30G.
  PLOT plt4g;  plt4g.init(); //Normalized hold amp counts as a function of run number, B=25G.
  PLOT plt4h;  plt4h.init(); //Normalized hold amp counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.holdfitamp,xd.eholdfitamp,xd.loadcounts,xd.eloadcounts,&normholdamp,&enormholdamp);
    if(xd.Tfield == 250) plt4a.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 100) plt4b.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 260) plt4c.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 175) plt4d.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 50) plt4e.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 30) plt4f.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 25) plt4g.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    if(xd.Tfield == 10) plt4h.loadpnt(xd.rnum,0,normholdamp,enormholdamp);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized Hold Amp",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt4a.make("plt4a");
  plt4a.haxis->SetTitle("Total Normalized Hold Amp");
  plt4a.haxis->GetXaxis()->SetTitle("run#");
  plt4a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt4a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt4a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt4a.haxis->Draw();
  plt4a.plot->SetMarkerColor(4);
  plt4a.plot->SetLineColor(4);
  plt4a.plot->SetMarkerStyle(20);
  plt4a.plot->Draw("PSame");
  
  plt4b.make("plt4b");
  plt4b.plot->SetMarkerColor(2);
  plt4b.plot->SetLineColor(2);
  plt4b.plot->SetMarkerStyle(20);
  plt4b.plot->SetMarkerSize(1.3);
  plt4b.plot->Draw("PSame");
  
  plt4c.make("plt4c");
  plt4c.plot->SetMarkerColor(1);
  plt4c.plot->SetLineColor(1);
  plt4c.plot->SetMarkerStyle(20);
  plt4c.plot->SetMarkerSize(1.3);
  plt4c.plot->Draw("PSame");
  
  plt4d.make("plt4d");
  plt4d.plot->SetMarkerColor(3);
  plt4d.plot->SetLineColor(3);
  plt4d.plot->SetMarkerStyle(20);
  plt4d.plot->SetMarkerSize(1.3);
  plt4d.plot->Draw("PSame");
  
  plt4e.make("plt4e");
  plt4e.plot->SetMarkerColor(5);
  plt4e.plot->SetLineColor(5);
  plt4e.plot->SetMarkerStyle(20);
  plt4e.plot->SetMarkerSize(1.3);
  plt4e.plot->Draw("PSame");
  
  plt4f.make("plt4f");
  plt4f.plot->SetMarkerColor(6);
  plt4f.plot->SetLineColor(6);
  plt4f.plot->SetMarkerStyle(20);
  plt4f.plot->SetMarkerSize(1.3);
  plt4f.plot->Draw("PSame");
  
  plt4g.make("plt4g");
  plt4g.plot->SetMarkerColor(7);
  plt4g.plot->SetLineColor(7);
  plt4g.plot->SetMarkerStyle(20);
  plt4g.plot->SetMarkerSize(1.3);
  plt4g.plot->Draw("PSame");
  
  plt4h.make("plt4h");
  plt4h.plot->SetMarkerColor(14);
  plt4h.plot->SetLineColor(14);
  plt4h.plot->SetMarkerStyle(20);
  plt4h.plot->SetMarkerSize(1.3);
  plt4h.plot->Draw("PSame");
  
  TLegend *leg4 = new TLegend(0.899,0.364,0.997,0.917);
  leg4->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg4->AddEntry(plt4a.plot,"B = 250 G","LP");
  leg4->AddEntry(plt4b.plot,"B = 100 G","LP");
  leg4->AddEntry(plt4c.plot,"B = 260 G", "LP");
  leg4->AddEntry(plt4d.plot,"B = 175 G", "LP");
  leg4->AddEntry(plt4e.plot,"B = 50 G", "LP");
  leg4->AddEntry(plt4f.plot,"B = 30 G", "LP");
  leg4->AddEntry(plt4g.plot,"B = 25 G", "LP");
  leg4->AddEntry(plt4h.plot,"B = 10 G", "LP");
  leg4->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return; 
}

void getunldfitamp() {
  int i;
  char name[500];
  double normunldamp,enormunldamp;
  
  //Load the extracted data tree
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt5a;  plt5a.init(); //Normalized unload amp counts as a function of run number, B=250G.
  PLOT plt5b;  plt5b.init(); //Normalized unload amp counts as a function of run number, B=100G.
  PLOT plt5c;  plt5c.init(); //Normalized unload amp counts as a function of run number, B=260G.
  PLOT plt5d;  plt5d.init(); //Normalized unload amp counts as a function of run number, B=175G.
  PLOT plt5e;  plt5e.init(); //Normalized unload amp counts as a function of run number, B=50G.
  PLOT plt5f;  plt5f.init(); //Normalized unload amp counts as a function of run number, B=30G.
  PLOT plt5g;  plt5g.init(); //Normalized unload amp counts as a function of run number, B=25G.
  PLOT plt5h;  plt5h.init(); //Normalized unload amp counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.unldfitamp,xd.eunldfitamp,xd.loadcounts,xd.eloadcounts,&normunldamp,&enormunldamp);
    if(xd.Tfield == 250) plt5a.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 100) plt5b.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 260) plt5c.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 175) plt5d.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 50) plt5e.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 30) plt5f.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 25) plt5g.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    if(xd.Tfield == 10) plt5h.loadpnt(xd.rnum,0,normunldamp,enormunldamp);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized unload Amp",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt5a.make("plt5a");
  plt5a.haxis->SetTitle("Total Normalized unload Amp");
  plt5a.haxis->GetXaxis()->SetTitle("run#");
  plt5a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt5a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt5a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt5a.haxis->Draw();
  plt5a.plot->SetMarkerColor(4);
  plt5a.plot->SetLineColor(4);
  plt5a.plot->SetMarkerStyle(20);
  plt5a.plot->Draw("PSame");
  
  plt5b.make("plt5b");
  plt5b.plot->SetMarkerColor(2);
  plt5b.plot->SetLineColor(2);
  plt5b.plot->SetMarkerStyle(20);
  plt5b.plot->SetMarkerSize(1.3);
  plt5b.plot->Draw("PSame");
  
  plt5c.make("plt5c");
  plt5c.plot->SetMarkerColor(1);
  plt5c.plot->SetLineColor(1);
  plt5c.plot->SetMarkerStyle(20);
  plt5c.plot->SetMarkerSize(1.3);
  plt5c.plot->Draw("PSame");
  
  plt5d.make("plt5d");
  plt5d.plot->SetMarkerColor(3);
  plt5d.plot->SetLineColor(3);
  plt5d.plot->SetMarkerStyle(20);
  plt5d.plot->SetMarkerSize(1.3);
  plt5d.plot->Draw("PSame");
  
  plt5e.make("plt5e");
  plt5e.plot->SetMarkerColor(5);
  plt5e.plot->SetLineColor(5);
  plt5e.plot->SetMarkerStyle(20);
  plt5e.plot->SetMarkerSize(1.3);
  plt5e.plot->Draw("PSame");
  
  plt5f.make("plt5f");
  plt5f.plot->SetMarkerColor(6);
  plt5f.plot->SetLineColor(6);
  plt5f.plot->SetMarkerStyle(20);
  plt5f.plot->SetMarkerSize(1.3);
  plt5f.plot->Draw("PSame");
  
  plt5g.make("plt5g");
  plt5g.plot->SetMarkerColor(7);
  plt5g.plot->SetLineColor(7);
  plt5g.plot->SetMarkerStyle(20);
  plt5g.plot->SetMarkerSize(1.3);
  plt5g.plot->Draw("PSame");
  
  plt5h.make("plt5h");
  plt5h.plot->SetMarkerColor(14);
  plt5h.plot->SetLineColor(14);
  plt5h.plot->SetMarkerStyle(20);
  plt5h.plot->SetMarkerSize(1.3);
  plt5h.plot->Draw("PSame");
  
  TLegend *leg5 = new TLegend(0.899,0.364,0.997,0.917);
  leg5->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg5->AddEntry(plt5a.plot,"B = 250 G","LP");
  leg5->AddEntry(plt5b.plot,"B = 100 G","LP");
  leg5->AddEntry(plt5c.plot,"B = 260 G", "LP");
  leg5->AddEntry(plt5d.plot,"B = 175 G", "LP");
  leg5->AddEntry(plt5e.plot,"B = 50 G", "LP");
  leg5->AddEntry(plt5f.plot,"B = 30 G", "LP");
  leg5->AddEntry(plt5g.plot,"B = 25 G", "LP");
  leg5->AddEntry(plt5h.plot,"B = 10 G", "LP");
  leg5->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return;  
}

void getholdfittau() {
  int i;
  char name[500];
  double normholdtau,enormholdtau;
  
  //Load the extracted data tree
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt6a;  plt6a.init(); //Normalized hold tau counts as a function of run number, B=250G.
  PLOT plt6b;  plt6b.init(); //Normalized hold tau counts as a function of run number, B=100G.
  PLOT plt6c;  plt6c.init(); //Normalized hold tau counts as a function of run number, B=260G.
  PLOT plt6d;  plt6d.init(); //Normalized hold tau counts as a function of run number, B=175G.
  PLOT plt6e;  plt6e.init(); //Normalized hold tau counts as a function of run number, B=50G.
  PLOT plt6f;  plt6f.init(); //Normalized hold tau counts as a function of run number, B=30G.
  PLOT plt6g;  plt6g.init(); //Normalized hold tau counts as a function of run number, B=25G.
  PLOT plt6h;  plt6h.init(); //Normalized hold tau counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.holdfittau,xd.eholdfittau,xd.loadcounts,xd.eloadcounts,&normholdtau,&enormholdtau);
    if(xd.Tfield == 250) plt6a.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 100) plt6b.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 260) plt6c.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 175) plt6d.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 50) plt6e.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 30) plt6f.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 25) plt6g.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    if(xd.Tfield == 10) plt6h.loadpnt(xd.rnum,0,normholdtau,enormholdtau);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized Hold tau",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt6a.make("plt6a");
  plt6a.haxis->SetTitle("Total Normalized Hold tau");
  plt6a.haxis->GetXaxis()->SetTitle("run#");
  plt6a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt6a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt6a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt6a.haxis->Draw();
  plt6a.plot->SetMarkerColor(4);
  plt6a.plot->SetLineColor(4);
  plt6a.plot->SetMarkerStyle(20);
  plt6a.plot->Draw("PSame");
  
  plt6b.make("plt6b");
  plt6b.plot->SetMarkerColor(2);
  plt6b.plot->SetLineColor(2);
  plt6b.plot->SetMarkerStyle(20);
  plt6b.plot->SetMarkerSize(1.3);
  plt6b.plot->Draw("PSame");
  
  plt6c.make("plt6c");
  plt6c.plot->SetMarkerColor(1);
  plt6c.plot->SetLineColor(1);
  plt6c.plot->SetMarkerStyle(20);
  plt6c.plot->SetMarkerSize(1.3);
  plt6c.plot->Draw("PSame");
  
  plt6d.make("plt6d");
  plt6d.plot->SetMarkerColor(3);
  plt6d.plot->SetLineColor(3);
  plt6d.plot->SetMarkerStyle(20);
  plt6d.plot->SetMarkerSize(1.3);
  plt6d.plot->Draw("PSame");
  
  plt6e.make("plt6e");
  plt6e.plot->SetMarkerColor(5);
  plt6e.plot->SetLineColor(5);
  plt6e.plot->SetMarkerStyle(20);
  plt6e.plot->SetMarkerSize(1.3);
  plt6e.plot->Draw("PSame");
  
  plt6f.make("plt6f");
  plt6f.plot->SetMarkerColor(6);
  plt6f.plot->SetLineColor(6);
  plt6f.plot->SetMarkerStyle(20);
  plt6f.plot->SetMarkerSize(1.3);
  plt6f.plot->Draw("PSame");
  
  plt6g.make("plt6g");
  plt6g.plot->SetMarkerColor(7);
  plt6g.plot->SetLineColor(7);
  plt6g.plot->SetMarkerStyle(20);
  plt6g.plot->SetMarkerSize(1.3);
  plt6g.plot->Draw("PSame");
  
  plt6h.make("plt6h");
  plt6h.plot->SetMarkerColor(14);
  plt6h.plot->SetLineColor(14);
  plt6h.plot->SetMarkerStyle(20);
  plt6h.plot->SetMarkerSize(1.3);
  plt6h.plot->Draw("PSame");
  
  TLegend *leg6 = new TLegend(0.899,0.364,0.997,0.917);
  leg6->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg6->AddEntry(plt6a.plot,"B = 250 G","LP");
  leg6->AddEntry(plt6b.plot,"B = 100 G","LP");
  leg6->AddEntry(plt6c.plot,"B = 260 G", "LP");
  leg6->AddEntry(plt6d.plot,"B = 175 G", "LP");
  leg6->AddEntry(plt6e.plot,"B = 50 G", "LP");
  leg6->AddEntry(plt6f.plot,"B = 30 G", "LP");
  leg6->AddEntry(plt6g.plot,"B = 25 G", "LP");
  leg6->AddEntry(plt6h.plot,"B = 10 G", "LP");
  leg6->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return;  
}

void getunldfittau() {
  int i;
  char name[500];
  double normunldtau,enormunldtau;
  
  //Load the extracted data tree
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt7a;  plt7a.init(); //Normalized unload tau counts as a function of run number, B=250G.
  PLOT plt7b;  plt7b.init(); //Normalized unload tau counts as a function of run number, B=100G.
  PLOT plt7c;  plt7c.init(); //Normalized unload tau counts as a function of run number, B=260G.
  PLOT plt7d;  plt7d.init(); //Normalized unload tau counts as a function of run number, B=175G.
  PLOT plt7e;  plt7e.init(); //Normalized unload tau counts as a function of run number, B=50G.
  PLOT plt7f;  plt7f.init(); //Normalized unload tau counts as a function of run number, B=30G.
  PLOT plt7g;  plt7g.init(); //Normalized unload tau counts as a function of run number, B=25G.
  PLOT plt7h;  plt7h.init(); //Normalized unload tau counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.unldfittau,xd.eunldfittau,xd.loadcounts,xd.eloadcounts,&normunldtau,&enormunldtau);
    if(xd.Tfield == 250) plt7a.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 100) plt7b.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 260) plt7c.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 175) plt7d.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 50) plt7e.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 30) plt7f.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 25) plt7g.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    if(xd.Tfield == 10) plt7h.loadpnt(xd.rnum,0,normunldtau,enormunldtau);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized unload tau",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt7a.make("plt7a");
  plt7a.haxis->SetTitle("Total Normalized unload tau");
  plt7a.haxis->GetXaxis()->SetTitle("run#");
  plt7a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt7a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt7a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt7a.haxis->Draw();
  plt7a.plot->SetMarkerColor(4);
  plt7a.plot->SetLineColor(4);
  plt7a.plot->SetMarkerStyle(20);
  plt7a.plot->Draw("PSame");
  
  plt7b.make("plt7b");
  plt7b.plot->SetMarkerColor(2);
  plt7b.plot->SetLineColor(2);
  plt7b.plot->SetMarkerStyle(20);
  plt7b.plot->SetMarkerSize(1.3);
  plt7b.plot->Draw("PSame");
  
  plt7c.make("plt7c");
  plt7c.plot->SetMarkerColor(1);
  plt7c.plot->SetLineColor(1);
  plt7c.plot->SetMarkerStyle(20);
  plt7c.plot->SetMarkerSize(1.3);
  plt7c.plot->Draw("PSame");
  
  plt7d.make("plt7d");
  plt7d.plot->SetMarkerColor(3);
  plt7d.plot->SetLineColor(3);
  plt7d.plot->SetMarkerStyle(20);
  plt7d.plot->SetMarkerSize(1.3);
  plt7d.plot->Draw("PSame");
  
  plt7e.make("plt7e");
  plt7e.plot->SetMarkerColor(5);
  plt7e.plot->SetLineColor(5);
  plt7e.plot->SetMarkerStyle(20);
  plt7e.plot->SetMarkerSize(1.3);
  plt7e.plot->Draw("PSame");
  
  plt7f.make("plt7f");
  plt7f.plot->SetMarkerColor(6);
  plt7f.plot->SetLineColor(6);
  plt7f.plot->SetMarkerStyle(20);
  plt7f.plot->SetMarkerSize(1.3);
  plt7f.plot->Draw("PSame");
  
  plt7g.make("plt7g");
  plt7g.plot->SetMarkerColor(7);
  plt7g.plot->SetLineColor(7);
  plt7g.plot->SetMarkerStyle(20);
  plt7g.plot->SetMarkerSize(1.3);
  plt7g.plot->Draw("PSame");
  
  plt7h.make("plt7h");
  plt7h.plot->SetMarkerColor(14);
  plt7h.plot->SetLineColor(14);
  plt7h.plot->SetMarkerStyle(20);
  plt7h.plot->SetMarkerSize(1.3);
  plt7h.plot->Draw("PSame");
  
  TLegend *leg7 = new TLegend(0.899,0.364,0.997,0.917);
  leg7->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg7->AddEntry(plt7a.plot,"B = 250 G","LP");
  leg7->AddEntry(plt7b.plot,"B = 100 G","LP");
  leg7->AddEntry(plt7c.plot,"B = 260 G", "LP");
  leg7->AddEntry(plt7d.plot,"B = 175 G", "LP");
  leg7->AddEntry(plt7e.plot,"B = 50 G", "LP");
  leg7->AddEntry(plt7f.plot,"B = 30 G", "LP");
  leg7->AddEntry(plt7g.plot,"B = 25 G", "LP");
  leg7->AddEntry(plt7h.plot,"B = 10 G", "LP");
  leg7->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return;     
}

void getholdfitbg() {
  int i;
   char name[500];
   double normholdbg,enormholdbg;

  //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt8a;  plt8a.init(); //Normalized hold amp counts as a function of run number, B=250G.
  PLOT plt8b;  plt8b.init(); //Normalized hold amp counts as a function of run number, B=100G.
  PLOT plt8c;  plt8c.init(); //Normalized hold amp counts as a function of run number, B=260G.
  PLOT plt8d;  plt8d.init(); //Normalized hold amp counts as a function of run number, B=175G.
  PLOT plt8e;  plt8e.init(); //Normalized hold amp counts as a function of run number, B=50G.
  PLOT plt8f;  plt8f.init(); //Normalized hold amp counts as a function of run number, B=30G.
  PLOT plt8g;  plt8g.init(); //Normalized hold amp counts as a function of run number, B=25G.
  PLOT plt8h;  plt8h.init(); //Normalized hold amp counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.holdfitbg,xd.eholdfitbg,xd.loadcounts,xd.eloadcounts,&normholdbg,&enormholdbg);
    if(xd.Tfield == 250) plt8a.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 100) plt8b.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 260) plt8c.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 175) plt8d.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 50) plt8e.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 30) plt8f.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 25) plt8g.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    if(xd.Tfield == 10) plt8h.loadpnt(xd.rnum,0,normholdbg,enormholdbg);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized Hold background",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt8a.make("plt8a");
  plt8a.haxis->SetTitle("Total Normalized Hold background");
  plt8a.haxis->GetXaxis()->SetTitle("run#");
  plt8a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt8a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt8a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt8a.haxis->Draw();
  plt8a.plot->SetMarkerColor(4);
  plt8a.plot->SetLineColor(4);
  plt8a.plot->SetMarkerStyle(20);
  plt8a.plot->Draw("PSame");
  
  plt8b.make("plt8b");
  plt8b.plot->SetMarkerColor(2);
  plt8b.plot->SetLineColor(2);
  plt8b.plot->SetMarkerStyle(20);
  plt8b.plot->SetMarkerSize(1.3);
  plt8b.plot->Draw("PSame");
  
  plt8c.make("plt8c");
  plt8c.plot->SetMarkerColor(1);
  plt8c.plot->SetLineColor(1);
  plt8c.plot->SetMarkerStyle(20);
  plt8c.plot->SetMarkerSize(1.3);
  plt8c.plot->Draw("PSame");
  
  plt8d.make("plt8d");
  plt8d.plot->SetMarkerColor(3);
  plt8d.plot->SetLineColor(3);
  plt8d.plot->SetMarkerStyle(20);
  plt8d.plot->SetMarkerSize(1.3);
  plt8d.plot->Draw("PSame");
  
  plt8e.make("plt8e");
  plt8e.plot->SetMarkerColor(5);
  plt8e.plot->SetLineColor(5);
  plt8e.plot->SetMarkerStyle(20);
  plt8e.plot->SetMarkerSize(1.3);
  plt8e.plot->Draw("PSame");
  
  plt8f.make("plt8f");
  plt8f.plot->SetMarkerColor(6);
  plt8f.plot->SetLineColor(6);
  plt8f.plot->SetMarkerStyle(20);
  plt8f.plot->SetMarkerSize(1.3);
  plt8f.plot->Draw("PSame");
  
  plt8g.make("plt8g");
  plt8g.plot->SetMarkerColor(7);
  plt8g.plot->SetLineColor(7);
  plt8g.plot->SetMarkerStyle(20);
  plt8g.plot->SetMarkerSize(1.3);
  plt8g.plot->Draw("PSame");
  
  plt8h.make("plt8h");
  plt8h.plot->SetMarkerColor(14);
  plt8h.plot->SetLineColor(14);
  plt8h.plot->SetMarkerStyle(20);
  plt8h.plot->SetMarkerSize(1.3);
  plt8h.plot->Draw("PSame");
  
  TLegend *leg8 = new TLegend(0.899,0.364,0.997,0.917);
  leg8->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg8->AddEntry(plt8a.plot,"B = 250 G","LP");
  leg8->AddEntry(plt8b.plot,"B = 100 G","LP");
  leg8->AddEntry(plt8c.plot,"B = 260 G", "LP");
  leg8->AddEntry(plt8d.plot,"B = 175 G", "LP");
  leg8->AddEntry(plt8e.plot,"B = 50 G", "LP");
  leg8->AddEntry(plt8f.plot,"B = 30 G", "LP");
  leg8->AddEntry(plt8g.plot,"B = 25 G", "LP");
  leg8->AddEntry(plt8h.plot,"B = 10 G", "LP");
  leg8->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return;   
}

void getunldfitbg() {
 int i;
  char name[500];
  double normunldbg,enormunldbg;
  
  //Load the extracted data tree
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //Create plot containers.
  PLOT xd;
  
  //initialize plots
  PLOT plt9a;  plt9a.init(); //Normalized unload tau counts as a function of run number, B=250G.
  PLOT plt9b;  plt9b.init(); //Normalized unload tau counts as a function of run number, B=100G.
  PLOT plt9c;  plt9c.init(); //Normalized unload tau counts as a function of run number, B=260G.
  PLOT plt9d;  plt9d.init(); //Normalized unload tau counts as a function of run number, B=175G.
  PLOT plt9e;  plt9e.init(); //Normalized unload tau counts as a function of run number, B=50G.
  PLOT plt9f;  plt9f.init(); //Normalized unload tau counts as a function of run number, B=30G.
  PLOT plt9g;  plt9g.init(); //Normalized unload tau counts as a function of run number, B=25G.
  PLOT plt9h;  plt9h.init(); //Normalized unload tau counts as a function of run number, B=10G.
  
   //Fill plot containers...
  int entries = tree->GetEntries();
  for(i=0 ; i < entries ; i++) {
    xd.get_entry(tree,i);
    edivide(xd.unldfitbg,xd.eunldfitbg,xd.loadcounts,xd.eloadcounts,&normunldbg,&enormunldbg);
    if(xd.Tfield == 250) plt9a.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 100) plt9b.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 260) plt9c.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 175) plt9d.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 50) plt9e.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 30) plt9f.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 25) plt9g.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    if(xd.Tfield == 10) plt9h.loadpnt(xd.rnum,0,normunldbg,enormunldbg);
    
  }
  
  
  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized unload background",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt9a.make("plt9a");
  plt9a.haxis->SetTitle("Total Normalized unload Background");
  plt9a.haxis->GetXaxis()->SetTitle("run#");
  plt9a.haxis->GetYaxis()->SetTitle("hold amp/load counts");
  plt9a.haxis->GetXaxis()->SetRangeUser(65,1000);
  plt9a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt9a.haxis->Draw();
  plt9a.plot->SetMarkerColor(4);
  plt9a.plot->SetLineColor(4);
  plt9a.plot->SetMarkerStyle(20);
  plt9a.plot->Draw("PSame");
  
  plt9b.make("plt9b");
  plt9b.plot->SetMarkerColor(2);
  plt9b.plot->SetLineColor(2);
  plt9b.plot->SetMarkerStyle(20);
  plt9b.plot->SetMarkerSize(1.3);
  plt9b.plot->Draw("PSame");
  
  plt9c.make("plt9c");
  plt9c.plot->SetMarkerColor(1);
  plt9c.plot->SetLineColor(1);
  plt9c.plot->SetMarkerStyle(20);
  plt9c.plot->SetMarkerSize(1.3);
  plt9c.plot->Draw("PSame");
  
  plt9d.make("plt9d");
  plt9d.plot->SetMarkerColor(3);
  plt9d.plot->SetLineColor(3);
  plt9d.plot->SetMarkerStyle(20);
  plt9d.plot->SetMarkerSize(1.3);
  plt9d.plot->Draw("PSame");
  
  plt9e.make("plt9e");
  plt9e.plot->SetMarkerColor(5);
  plt9e.plot->SetLineColor(5);
  plt9e.plot->SetMarkerStyle(20);
  plt9e.plot->SetMarkerSize(1.3);
  plt9e.plot->Draw("PSame");
  
  plt9f.make("plt9f");
  plt9f.plot->SetMarkerColor(6);
  plt9f.plot->SetLineColor(6);
  plt9f.plot->SetMarkerStyle(20);
  plt9f.plot->SetMarkerSize(1.3);
  plt9f.plot->Draw("PSame");
  
  plt9g.make("plt9g");
  plt9g.plot->SetMarkerColor(7);
  plt9g.plot->SetLineColor(7);
  plt9g.plot->SetMarkerStyle(20);
  plt9g.plot->SetMarkerSize(1.3);
  plt9g.plot->Draw("PSame");
  
  plt9h.make("plt9h");
  plt9h.plot->SetMarkerColor(14);
  plt9h.plot->SetLineColor(14);
  plt9h.plot->SetMarkerStyle(20);
  plt9h.plot->SetMarkerSize(1.3);
  plt9h.plot->Draw("PSame");
  
  TLegend *leg9 = new TLegend(0.899,0.364,0.997,0.917);
  leg9->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg9->AddEntry(plt9a.plot,"B = 250 G","LP");
  leg9->AddEntry(plt9b.plot,"B = 100 G","LP");
  leg9->AddEntry(plt9c.plot,"B = 260 G", "LP");
  leg9->AddEntry(plt9d.plot,"B = 175 G", "LP");
  leg9->AddEntry(plt9e.plot,"B = 50 G", "LP");
  leg9->AddEntry(plt9f.plot,"B = 30 G", "LP");
  leg9->AddEntry(plt9g.plot,"B = 25 G", "LP");
  leg9->AddEntry(plt9h.plot,"B = 10 G", "LP");
  leg9->Draw("Same");

  drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  
 return;   
  
  
  
}

void holdint() {
  int i;
  char name[500];
  double normunld,enormunld;
 
  //Load the extracted data tree.
 sprintf(name,"%s/obsvbls.root",ODIR);
 TTree *tree = xdload(name,"tree1");

 //Create plot containers.
 PLOT xd;

 PLOT plt10a;  plt10a.init(); //Normalized unload counts as a function of run number, B=250G.
 PLOT plt10b;  plt10b.init(); //Normalized unload counts as a function of run number, B=100G.
 PLOT plt10c;  plt10c.init(); //Normalized unload counts as a function of run number, B=260G.
 PLOT plt10d;  plt10d.init(); //Normalized unload counts as a function of run number, B=175G.
 PLOT plt10e;  plt10e.init(); //Normalized unload counts as a function of run number, B=50G.
 PLOT plt10f;  plt10f.init(); //Normalized unload counts as a function of run number, B=30G.
 PLOT plt10g;  plt10g.init(); //Normalized unload counts as a function of run number, B=25G.
 PLOT plt10h;  plt10h.init(); //Normalized unload counts as a function of run number, B=10G.
  
  
  //Fill plot containers...
 int entries = tree->GetEntries();
 for(i=0 ; i < entries ; i++) {
   xd.get_entry(tree,i);
    edivide(xd.unldcounts,xd.eunldcounts,xd.loadcounts,xd.eloadcounts,&normunld,&enormunld);
    if(xd.Tfield == 250) plt10a.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 100) plt10b.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 260) plt10c.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 175) plt10d.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 50) plt10e.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 30) plt10f.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 25) plt10g.loadpnt(xd.Thold,0,normunld,enormunld);
    if(xd.Tfield == 10) plt10h.loadpnt(xd.Thold,0,normunld,enormunld);
    
  }

  

  //Draw normalized unload counts plot.
  TCanvas *can2 = new TCanvas("can2","Normalized Unload",735,5,700,500);
  can2->SetFrameLineColor(0);
  can2->cd();
  
  plt10a.make("plt10a");
  plt10a.haxis->SetTitle("Total Normalized Unload Counts");
  plt10a.haxis->GetXaxis()->SetTitle("Hold interval duation (s)");
  plt10a.haxis->GetYaxis()->SetTitle("unload counts/load counts");
  plt10a.haxis->GetXaxis()->SetRangeUser(0,100);
  plt10a.haxis->GetYaxis()->SetRangeUser(0,1);
  plt10a.haxis->Draw();
  plt10a.plot->SetMarkerColor(4);
  plt10a.plot->SetLineColor(4);
  plt10a.plot->SetMarkerStyle(20);
  plt10a.plot->Draw("PSame");
  
  plt10b.make("plt10b");
  plt10b.plot->SetMarkerColor(2);
  plt10b.plot->SetLineColor(2);
  plt10b.plot->SetMarkerStyle(20);
  plt10b.plot->SetMarkerSize(1.3);
  plt10b.plot->Draw("PSame");
  
  plt10c.make("plt10c");
  plt10c.plot->SetMarkerColor(1);
  plt10c.plot->SetLineColor(1);
  plt10c.plot->SetMarkerStyle(20);
  plt10c.plot->SetMarkerSize(1.3);
  plt10c.plot->Draw("PSame");
  
  plt10d.make("plt10d");
  plt10d.plot->SetMarkerColor(3);
  plt10d.plot->SetLineColor(3);
  plt10d.plot->SetMarkerStyle(20);
  plt10d.plot->SetMarkerSize(1.3);
  plt10d.plot->Draw("PSame");
  
  plt10e.make("plt10e");
  plt10e.plot->SetMarkerColor(5);
  plt10e.plot->SetLineColor(5);
  plt10e.plot->SetMarkerStyle(20);
  plt10e.plot->SetMarkerSize(1.3);
  plt10e.plot->Draw("PSame");
  
  plt10f.make("plt10f");
  plt10f.plot->SetMarkerColor(6);
  plt10f.plot->SetLineColor(6);
  plt10f.plot->SetMarkerStyle(20);
  plt10f.plot->SetMarkerSize(1.3);
  plt10f.plot->Draw("PSame");
  
  plt10g.make("plt10g");
  plt10g.plot->SetMarkerColor(7);
  plt10g.plot->SetLineColor(7);
  plt10g.plot->SetMarkerStyle(20);
  plt10g.plot->SetMarkerSize(1.3);
  plt10g.plot->Draw("PSame");
  
  plt10h.make("plt10h");
  plt10h.plot->SetMarkerColor(14);
  plt10h.plot->SetLineColor(14);
  plt10h.plot->SetMarkerStyle(20);
  plt10h.plot->SetMarkerSize(1.3);
  plt10h.plot->Draw("PSame");
  
  TLegend *leg10 = new TLegend(0.899,0.364,0.997,0.917);
  leg10->SetFillColor(0);
  //leg2->SetLineColor(0);
  //leg2->SetTextSize(0.04);
  leg10->AddEntry(plt10a.plot,"B = 250 G","LP");
  leg10->AddEntry(plt10b.plot,"B = 100 G","LP");
  leg10->AddEntry(plt10c.plot,"B = 260 G", "LP");
  leg10->AddEntry(plt10d.plot,"B = 175 G", "LP");
  leg10->AddEntry(plt10e.plot,"B = 50 G", "LP");
  leg10->AddEntry(plt10f.plot,"B = 30 G", "LP");
  leg10->AddEntry(plt10g.plot,"B = 25 G", "LP");
  leg10->AddEntry(plt10h.plot,"B = 10 G", "LP");
  leg10->Draw("Same");

  /*drawallusegs();
  //standard guides
  guidelabel(can2, "Non Coated Copper",.093, .004);
  guidelabel(can2,"DLC1",.254,.004);
  guidelabel(can2,"EP",.484,.004);
  guidelabel(can2,"DLC2",.613,.004);
  guidelabel(can2,"SS",.699,.004);
  guidelabel(can2,"MP",.788,.004);
  guidelabel(can2,"EP2",.850,.004);
  
  //poison guides
  guidelabel(can2, "DLC PSN 1a",.198, .850);
  guidelabel(can2,"DLC PSN 1b",.308,.850);
  guidelabel(can2,"EP PSN 1",.497,.850);
  guidelabel(can2,"DLC PSN 2",.591,.850);
  guidelabel(can2,"MP PSN",.777,.850);
  guidelabel(can2,"EP PSN 2",.836,.907);
  
  
  char name[100];
  int fitmax;
 
  sprintf(name,"[0]*exp(-(x-T)/[1])+[2]",Tloadend);
  //[0] is amplitude
  //[1] is Tau
  //[2] is background
  TF1 *fit = new TF1("fit",name,7,41);
 
  fit->SetParameter(0,.834);
  fit->SetParameter(1,3);
  fit->SetParameter(2,1);
 
  can2->Fit(fit,"LMN","",Tloadend,Tunldstart);
  fit->Draw("same");
  
  holdfitamp = fit->GetParameter(0);
  eholdfitamp = fit->GetParError(0);
 
  holdfittau = fit->GetParameter(1);
  eholdfittau = fit->GetParError(1);
 
  holdfitbg = fit->GetParameter(2);
  eholdfitbg = fit->GetParError(2);
 
 
 //reporting
  printf("Hold interval Fit amplitude value for Run %d: %lf +/- %lf\n\n",RNUM,holdfitamp,eholdfitamp);
  printf("Hold interval Fit Tau value for Run %d: %lf +/- %lf\n\n",RNUM,holdfittau,eholdfittau);
  printf("Hold interval Fit Background value for Run %d: %lf +/- %lf\n\n",RNUM,holdfitbg,eholdfitbg);
  eholdfitbg*/
  return;  
}

  /* ------------------------- Lemmings ----------------------- */
  
void guidelabel(TCanvas *can, const char* text, double horz, double vert) {
  
  
  TPaveText *noncoat = new TPaveText(horz, vert,horz+.07, vert+.06 ,"NDC");
  noncoat->AddText(text);
  noncoat->SetTextSize(.025);
  noncoat->Draw("Same");
  
  
 return; 
}

TTree* xdload(const char *xdfname, const char *treename) {  //Load an extracted data tree.
  TDirectory *work = gDirectory;
  TFile *tempdir = new TFile(xdfname,"READ");
  tempdir->cd();
  TTree *tree = (TTree*)gDirectory->Get(treename);
  if(tree == NULL) printf("\n  Requested TTree %s not found in %s!\n",treename,xdfname);
  tree->GetEntry(0);  /*  For some reason a .root file produced by the new version of root
                                           compiled on the Mac won't transfer directories correctly unless it is accessed first. However, .root files produced by an older version of Root don't seem to have this problem!?  */
  tree->SetDirectory(work);
  work->cd();
  
  return tree;
}


void edivide(double numr, double enumr, double denom, double edenom, double *result, double *eresult) {

  double penumr = 0;
  if(numr != 0) penumr = enumr/numr;
  double pedenom = 0;
  if(denom != 0) pedenom = edenom/denom;
  double div = numr/denom;
  double pediv = sqrt(pow(penumr,2) + pow(pedenom,2));
  double ediv = pediv * div;
  *result = div;
  *eresult = fabs(ediv);
  
  return;
}


void vlinedraw(double low, double high,int color) {
  double cutlinelx[2], cutlinehx[2];
  double cutlinely[2], cutlinehy[2];
  cutlinely[0]=0;
  cutlinehy[0]=0;
  cutlinely[1]=1e10;
  cutlinehy[1]=1e10;
  cutlinelx[0] = low;
  cutlinelx[1] = low;
  cutlinehx[0] = high;
  cutlinehx[1] = high;
  TGraph *glow = new TGraph(2,cutlinelx,cutlinely);
  TGraph *ghigh = new TGraph(2,cutlinehx,cutlinehy);
  glow->SetLineColor(color);
  ghigh->SetLineColor(color);
  glow->Draw("same");
  ghigh->Draw("same");
  
  return;
}


void wsum(int num, double *in, double *ein, double *out, double *eout) {
  int i,j;
  double sumw,sumxw;
  
  int cnt = 0;
  
  for(i=0 ; i < 5000 ; i++) {
    out[i] = -1;
    eout[i] = -1;
  }
  
  for(i=0 ; i < 5000 ; i++) if(in[i] == -1) break;
  int incnt = i;
  
  for(i=0 ; i < incnt ; i = i + num) {
    sumw = 0;
    sumxw = 0;
    for(j=0 ; j < num ; j++) {
      if(ein[i+j] < 0) return;
      sumw = sumw + 1/pow(ein[i+j],2);
    }
    for(j=0 ; j < num ; j++) {
      if(in[i+j] < 0) return;
      sumxw = sumxw + in[i+j]/pow(ein[i+j],2);
    }
    out[cnt]  = sumxw/sumw;
    eout[cnt] = sqrt(1/sumw);
    cnt++;
  }
}


int RunUnits(const char *type, int runs[RUMAX]) {
  int i,j,index;
  int fflag = 0;
  char pars[500];
  
  for(i=0; i < RUMAX ; i++) runs[i] = -1;
  
  FILE *rufp = fopen("RunUnits.txt","r");
  if(rufp == NULL) {
    printf("\n  'RunUnits.txt' not found!\n\n");
    return -1;
  }
  
  for(i=0 ; i < 2 ; i++) fgets(pars,500,rufp);
  index = 0;
  while(feof(rufp) == 0 && fflag == 0) {
    fscanf(rufp,"%s",pars);
    if(strcmp(pars,type) == 0) {
      fflag = 1;
      for(j=0 ; j < RUMAX ; j++) {
        fscanf(rufp,"%s",pars);
        if((int)atof(pars) != 0) {
          runs[index] = (int)atof(pars);
          index++;
        }
        if(strcmp(pars,"----") == 0) break;
      }
    }
  }
    
    if(fflag == 0 || runs[0] == -1) {
      printf("\n  Unknown run type in 'RunUnits'!\n\n");
      return -1;
    }
    else {
      //printf("%d runs from run unit: %s\n",index,type);
      for(i=0 ; i < RUMAX ; i++) {
        if(runs[i] == -1) break;
        //printf("    %d\n",runs[i]);
      }
    }

  fclose(rufp);
  return index;
}


void useg(const char *type, int color) {
  int i,start,stop;
  
  int runs[RUMAX];
  RunUnits(type,runs);
  
  start = runs[0];
  
  for(i=0 ; i < RUMAX ; i++) {
    if(runs[i] == -1) break;
    stop = runs[i];
  }
  
  vlinedraw(start,stop,color);
  

  return;
}


void drawallusegs() {
  
  useg("NON_COAT_CU",1); //black
  useg("DLC_1",2);//red
  useg("DLC_PSN_1a",3);//green
  useg("DLC_PSN_1b",3);
  useg("EP_CU_1",4);//blue
  useg("EP_CU_PSN_1",5);//yellow
  useg("DLC_2",2);//red
  useg("DLC_PSN_2",3);//green
  useg("SS",8);//dark green
  useg("MP_CU",7);//light blue
  useg("MP_CU_PSN",6);//pink
  useg("EP_CU_2",4);//blue
  useg("EP_CU_PSN_2",5);//yellow
 
  return;
}

void guidesep(const char* type, double field) {
  int i,j;
  char name[500];
  double normunld,enormunld;
  double guide;
  
  int runs[RUMAX];
  guide = RunUnits(type,runs);
  
  //Load the extracted data tree.
  sprintf(name,"%s/obsvbls.root",ODIR);
  TTree *tree = xdload(name,"tree1");

  //create plot containers
  PLOT xd;
  
  //initialize plots
  PLOT plt11;  plt11.init();
  
  int entries = tree->GetEntries(); 
  int flag = 0;
 for(i=0 ; i < entries ; i++) { 
   flag = 0;
   xd.get_entry(tree,i);
   for(j=0; j < RUMAX; j++) {
     if(runs[j] == xd.rnum) {
      flag = 1;
      break;
     } 
  }
  if(flag == 0) continue;
  
    edivide(xd.unldcounts,xd.eunldcounts,xd.loadcounts,xd.eloadcounts,&normunld,&enormunld);
    if(xd.Tfield == field) plt11.loadpnt(xd.Thold,0,normunld,enormunld);
  }
    
    //Draw normalized unload counts plot.
  TCanvas *can4 = new TCanvas("can4","Normalized Unload",735,5,700,500);
  can4->SetFrameLineColor(0);
  can4->cd();
  
  plt11.make("plt11");
  plt11.haxis->SetTitle("Total Normalized Unload Counts");
  plt11.haxis->GetXaxis()->SetTitle("Hold");
  plt11.haxis->GetYaxis()->SetTitle("unload counts/load counts");
  plt11.haxis->GetXaxis()->SetRangeUser(0,100);
  plt11.haxis->GetYaxis()->SetRangeUser(0.0001,1);
  plt11.haxis->Draw();
  can4->SetLogy();
  plt11.plot->SetMarkerColor(4);
  plt11.plot->SetLineColor(4);
  plt11.plot->SetMarkerStyle(20);
  plt11.plot->Draw("PSame");
  
  double valx, valy, tau;
  tau = 13;
  
  plt11.plot->GetPoint(0,valx,valy);
  sprintf(name,"[0]*exp(-(x-%f)/[1])",valy);
  
 
  TF1 *fit = new TF1("fit",name,7,80);
 
  fit->SetParameter(0,valy);
  fit->SetParameter(1,tau);
  
  plt11.plot->Fit(fit,"MNQ","",valx,80);
  fit->Draw("Same");
  
  char report[100];
  sprintf(report,"Amp for this plot is: %lf\n\nTau for the plot is: %lf\n\n",valy,tau);
  TPaveText *vals = new TPaveText(.5,.5,.5,.5,"NDC");
  vals->AddText(report);
  vals->Draw("Same");
  
   printf("Unload Count Fit amplitude value: %lf \n\n",valy);
   printf("Unload Count Fit Tau value: %lf \n\n",tau);
  
  return; 
}

void sepselect() {
  char type[100];
  double field;
  //int numSelect;
  
  
  printf("\n Type the guide you wish to see:\n\n 1. NON_COAT_CU\n 2. DLC_1\n 3. DLC_PSN_1a\n 4. DLC_PSN_1b\n 5. EP_CU_1\n 6. EP_CU_PSN_1\n 7. DLC_2\n 8. DLC_PSN_2\n 9. SS\n 10. MP_CU\n 11. MP_CU_PSN\n 12. EP_CU_2\n 13. EP_CU_PSN_2\n");
  scanf("%s",type);
  
  /*switch(numSelect){
      
      
  }*/
  
  printf("Type the field value you wish to see:\n\n260\n250\n175\n100\n50\n30\n25\n10\n\n");
  scanf("%lf",&field);
  
  //sprintf(" %s %lf ", &type, &field);
  
  guidesep(type,field);
  
  return;
}

double dpb(int rnum) {
    
   int i; 
   char name[500];  
   
   double sumclean, sumunload;
   double a = 3.14*.03515625;
   double A = 3.14*pow((2.75/2), 2);
   double beta = 45;
   double E = exp(1);
   double delta = 20;
   double Delta;
   
    
   double Thold,Tunload,tautwidle, etautwidle, tauhat, etauhat, tau1, etau1, tau2, etau2;
   
   //Load the extracted data tree.
   sprintf(name,"%s/obsvbls.root",ODIR);
   TTree *tree = xdload(name,"tree1");

   PLOT xd;
   int entries = tree->GetEntries();
   
   for(i=0 ; i < entries ; i++) {
     xd.get_entry(tree,i); 
  
     if(rnum == xd.rnum) { 
       Thold = xd.Thold;
       Tunload = xd.Tunload;
       sumclean = xd.bsnholdcounts;
       sumunload = xd.bsnunldcounts;
       tautwidle = xd.tautwidle;
       etautwidle = xd.etautwidle;
       tauhat = xd.tauhat;
       etauhat = xd.etauhat;
       tau1 = xd.tau1;
       etau1 = xd.etau1;
       tau2 = xd.tau2;
       etau2 = xd.etau2;
     }
   }
    /*
    printf("a: %f\n",a);
    printf("A: %f\n",A);
    printf("Thold: %f\n",Thold);
    printf("Tunload: %f\n",Tunload);
    printf("sumclean: %f\n",sumclean);
    printf("sumunload: %f\n",sumunload);
    printf("tautwidle: %lf\n",tautwidle);
    printf("etautwidle: %lf\n",etautwidle);
    printf("tauhat: %lf\n",tauhat);
    printf("etauhat: %lf\n",etauhat);
    printf("tau1: %lf\n",tau1);
    printf("etau1: %lf\n",etau1);
    printf("tau2: %lf\n",tau2);
    printf("etau2: %lf\n",etau2);
    */
    double ntest1 = exp(Thold/tauhat - Thold/tautwidle);
    double ntest2 = -1 + pow(E,Tunload/tau1);
    double ntest3 = sumclean*tau1 - a*pow(E,Tunload/tau1);
    double ntest4 = -1 + pow(E,Thold/tautwidle);
    double ntest5 = sumunload*tautwidle;
    
    double dtest1 = (-1 + pow(E,Tunload/tau2));
    double dtest2 = sumclean*tau2;
    double dtest3 = a*(pow(E,Tunload/tau2));
    double dtest4 = (-1 + pow(E,Thold/tauhat));
    double dtest5 = sumunload*tauhat;
    
    //printf("\n\nntest1: %lf\nntest2: %lf\nntest3: %lf\nntest4: %lf\nntest5: %lf\n",ntest1,ntest2,ntest3,ntest4,ntest5);
    //printf("dtest1: %lf\ndtest2: %lf\ndtest3: %lf\ndtest4: %lf\ndtest5: %lf\n",dtest1,dtest2,dtest3,dtest4,dtest5);
    
    double dT = ntest1*(A*(ntest2)*ntest3*(ntest4)*ntest5);
    double dB = delta*beta*(-(A*(dtest1)*dtest2)+dtest3*(dtest4)*dtest5);
    
    //printf("\n\ndT: ntest1*(A(ntest2)*ntest3*(ntest4)*ntest5) = %lf", dT);
    //printf("\ndB: delta*beta*(-(A(dtest1)*dtest2)+dtest3(dtest4)dtest5 = %lf)\n\n", dB);
    
    Delta = (dT/dB);

  /*
   Delta =  (pow(E,Thold/tauhat - Thold/tautwidle - Tunload/tau1 + Tunload/tau2)*(A*(-1 + pow(E,Tunload/tau1))*sumclean*tau1 - a*pow(E,Tunload/tau1)*(-1 + pow(E,Thold/tautwidle))*sumunload*tautwidle)) 
   /
   (beta*delta*(-(A*(-1 + pow(E,Tunload/tau2))*sumclean*tau2) + a*pow(E,Tunload/tau2)*(-1 + pow(E,Thold/tauhat))*sumunload*tauhat));
   */
  
   if(Delta!=Delta || Delta > 40000 || Delta < -40000){
       Delta = 0;
       return Delta;   
   }
   
   //printf("\n\nDelta %d: %f\n\n", rnum,Delta);
    
   
  
   
 return Delta; 
}

void plotObs(){
   int i; 
   char name[500];
   int choice;
   int rnum;
   double sumclean, sumunload, Thold, Tunload, tautwidle, etautwidle, tauhat, etauhat, tau1, etau1, tau2, etau2;
   
    //Load the extracted data tree.
   sprintf(name,"%s/obsvbls.root",ODIR);
   TTree *tree = xdload(name,"tree1");

   PLOT xd;
   
   PLOT pltSumclean; pltSumclean.init();
   PLOT pltSumunload; pltSumunload.init();
   PLOT pltTautwidle; pltTautwidle.init();
   PLOT pltTauhat; pltTauhat.init();
   PLOT pltTau1; pltTau1.init();
   PLOT pltTau2; pltTau2.init();
   
   
   int entries = tree->GetEntries();
   
   for(i=0 ; i < entries ; i++) {
     xd.get_entry(tree,i); 
  
     if(rnum == xd.rnum) { 
       Thold = xd.Thold;
       Tunload = xd.Tunload;
       sumclean = xd.holdcounts;    //
       sumunload = xd.unldcounts;   //
       tautwidle = xd.tautwidle;    //
       etautwidle = xd.etautwidle;
       tauhat = xd.tauhat;          //
       etauhat = xd.etauhat;
       tau1 = xd.tau1;              //
       etau1 = xd.etau1;
       tau2 = xd.tau2;              //
       etau2 = xd.etau2;
     }
     pltSumclean.loadpnt(xd.rnum,0,xd.holdcounts,0);
     pltSumunload.loadpnt(xd.rnum,0,xd.unldcounts,0);
     pltTautwidle.loadpnt(xd.rnum,0,xd.tautwidle,0);
     pltTauhat.loadpnt(xd.rnum,0,xd.tauhat,0);
     pltTau1.loadpnt(xd.rnum,0,xd.tau1,0);
     pltTau2.loadpnt(xd.rnum,0,xd.tau2,0);
   }
   
   cout << "\n1.sumclean\n2.sumunload\n3.tautwidle\n4.tauhat\n5.tau1/tau2\nSelect a Observable plot: ";
   cin >> choice;
   
   switch(choice){
       case 1: //Draw sumclean plot.
       {
           TCanvas *canClean = new TCanvas("canClean","Sumclean v. Rnum",735,5,700,500);
           canClean->SetFrameLineColor(0);
           canClean->cd();
           pltSumclean.make("pltSumclean");
           pltSumclean.haxis->SetTitle("Sumclean v. Run");
           pltSumclean.haxis->GetXaxis()->SetTitle("run#");
           pltSumclean.haxis->GetYaxis()->SetTitle("UCN Counts");
           pltSumclean.haxis->GetXaxis()->SetRangeUser(0,400);
           pltSumclean.haxis->GetYaxis()->SetRangeUser(0,6000);
           pltSumclean.haxis->Draw();
           pltSumclean.plot->SetMarkerColor(4);
           pltSumclean.plot->SetLineColor(4);
           pltSumclean.plot->SetMarkerStyle(20);
           pltSumclean.plot->Draw("PSame");
           break;
       }
       case 2:
       {
           TCanvas *canUnload = new TCanvas("canUnload","Sumunload v. Rnum",735,5,700,500);
           canUnload->SetFrameLineColor(0);
           canUnload->cd();
           pltSumunload.make("pltSumunload");
           pltSumunload.haxis->SetTitle("Sumunld v. Run");
           pltSumunload.haxis->GetXaxis()->SetTitle("run#");
           pltSumunload.haxis->GetYaxis()->SetTitle("UCN Counts");
           pltSumunload.haxis->GetXaxis()->SetRangeUser(0,400);
           pltSumunload.haxis->GetYaxis()->SetRangeUser(0,6000);
           pltSumunload.haxis->Draw();
           pltSumunload.plot->SetMarkerColor(4);
           pltSumunload.plot->SetLineColor(4);
           pltSumunload.plot->SetMarkerStyle(20);
           pltSumunload.plot->Draw("PSame");
           break;
       }
       case 3:
       {
           TCanvas *canTautwidle = new TCanvas("canTautwidle","Tautwidle v. Rnum",735,5,700,500);
           canTautwidle->SetFrameLineColor(0);
           canTautwidle->cd();
           pltTautwidle.make("pltTautwidle");
           pltTautwidle.haxis->SetTitle("Tautwidle v. Run");
           pltTautwidle.haxis->GetXaxis()->SetTitle("run#");
           pltTautwidle.haxis->GetYaxis()->SetTitle("Tautwidle value");
           pltTautwidle.haxis->GetXaxis()->SetRangeUser(0,400);
           pltTautwidle.haxis->GetYaxis()->SetRangeUser(0,15);
           pltTautwidle.haxis->Draw();
           pltTautwidle.plot->SetMarkerColor(4);
           pltTautwidle.plot->SetLineColor(4);
           pltTautwidle.plot->SetMarkerStyle(20);
           pltTautwidle.plot->Draw("PSame");
           break;
       }
       case 4: 
       {
           TCanvas *canTauhat = new TCanvas("canTauhat","Tauhat v. Rnum",735,5,700,500);
           canTauhat->SetFrameLineColor(0);
           canTauhat->cd();
           pltTauhat.make("pltTauhat");
           pltTauhat.haxis->SetTitle("tauhat v. Run");
           pltTauhat.haxis->GetXaxis()->SetTitle("run#");
           pltTauhat.haxis->GetYaxis()->SetTitle("Tauhat value");
           pltTauhat.haxis->GetXaxis()->SetRangeUser(0,400);
           pltTauhat.haxis->GetYaxis()->SetRangeUser(0,4000);
           pltTauhat.haxis->Draw();
           pltTauhat.plot->SetMarkerColor(4);
           pltTauhat.plot->SetLineColor(4);
           pltTauhat.plot->SetMarkerStyle(20);
           pltTauhat.plot->Draw("PSame");
           break;
       }
       case 5:
       {
           TCanvas *canTau1 = new TCanvas("canTau1","Tau1/tau2 v. Rnum",735,5,700,500);
           canTau1->SetFrameLineColor(0);
           canTau1->cd();
           pltTau1.make("pltTau1");
           pltTau1.haxis->SetTitle("tau1/tau2 v. Run");
           pltTau1.haxis->GetXaxis()->SetTitle("run#");
           pltTau1.haxis->GetYaxis()->SetTitle("Tau1/tau2 value");
           pltTau1.haxis->GetXaxis()->SetRangeUser(0,400);
           pltTau1.haxis->GetYaxis()->SetRangeUser(0,20);
           pltTau1.haxis->Draw();
           pltTau1.plot->SetMarkerColor(4);
           pltTau1.plot->SetLineColor(4);
           pltTau1.plot->SetMarkerStyle(20);
           pltTau1.plot->Draw("PSame");
           break;
       }
   }
    
    return;
}

void plotdpb(){
    int i; 
    char name[500];
    double depol;
    
    //Load the extracted data tree.
    sprintf(name,"%s/obsvbls.root",ODIR);
    TTree *tree = xdload(name,"tree1");

    PLOT xd;
    
    int entries = tree->GetEntries();
    
    
    
    PLOT pltDelta10;  pltDelta10.init();
    PLOT pltDelta25;  pltDelta25.init();
    PLOT pltDelta30;  pltDelta30.init();
    PLOT pltDelta50;  pltDelta50.init();
    PLOT pltDelta100;  pltDelta100.init();
    PLOT pltDelta175;  pltDelta175.init();
    PLOT pltDelta250;  pltDelta250.init();
    PLOT pltDelta260;  pltDelta260.init();
    
    for(i=0 ; i < entries; i++) {
        xd.get_entry(tree,i);
        
        depol = dpb(xd.rnum); 
        printf("Depol %d: %f\n\n",i,depol); 
        //pltDelta.loadpnt(i,0,depol,0);
        
        if(xd.Tfield == 260) pltDelta260.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 250) pltDelta250.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 175) pltDelta175.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 100) pltDelta100.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 50) pltDelta50.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 30) pltDelta30.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 25) pltDelta25.loadpnt(xd.rnum,0,depol,0);
        if(xd.Tfield == 10) pltDelta10.loadpnt(xd.rnum,0,depol,0);
    }
   
   TCanvas *canDPB = new TCanvas("canDPB","Depol v. Rnum",735,5,900,700);
   canDPB->SetFrameLineColor(0);
   canDPB->cd();
           
           pltDelta250.make("pltDelta250");
           pltDelta250.haxis->SetTitle("Depol v. Run");
           pltDelta250.haxis->GetXaxis()->SetTitle("run#");
           pltDelta250.haxis->GetYaxis()->SetTitle("DPB");
           pltDelta250.haxis->GetXaxis()->SetRangeUser(70,400);
           pltDelta250.haxis->GetYaxis()->SetRangeUser(-40000,40000);
           pltDelta250.haxis->Draw();
           pltDelta250.plot->SetMarkerColor(12);
           pltDelta250.plot->SetLineColor(4);
           pltDelta250.plot->SetMarkerStyle(20);
           pltDelta250.plot->Draw("PSame");
           
           pltDelta100.make("pltDelta100");
           pltDelta100.plot->SetMarkerColor(5);
           pltDelta100.plot->SetLineColor(4);
           pltDelta100.plot->SetMarkerStyle(20);
           pltDelta100.plot->Draw("PSame");
           
           pltDelta260.make("pltDelta260");
           pltDelta260.plot->SetMarkerColor(8);
           pltDelta260.plot->SetLineColor(4);
           pltDelta260.plot->SetMarkerStyle(20);
           pltDelta260.plot->Draw("PSame");
           
           pltDelta175.make("pltDelta175");
           pltDelta175.plot->SetMarkerColor(7);
           pltDelta175.plot->SetLineColor(4);
           pltDelta175.plot->SetMarkerStyle(20);
           pltDelta175.plot->Draw("PSame");
           
           pltDelta50.make("pltDelta50");
           pltDelta50.plot->SetMarkerColor(6);
           pltDelta50.plot->SetLineColor(4);
           pltDelta50.plot->SetMarkerStyle(20);
           pltDelta50.plot->Draw("PSame");
           
           pltDelta30.make("pltDelta30");
           pltDelta30.plot->SetMarkerColor(2);
           pltDelta30.plot->SetLineColor(4);
           pltDelta30.plot->SetMarkerStyle(20);
           pltDelta30.plot->Draw("PSame");
           
           pltDelta25.make("pltDelta25");
           pltDelta25.plot->SetMarkerColor(3);
           pltDelta25.plot->SetLineColor(4);
           pltDelta25.plot->SetMarkerStyle(20);
           pltDelta25.plot->Draw("PSame");
           
           pltDelta10.make("pltDelta10");
           pltDelta10.plot->SetMarkerColor(4);
           pltDelta10.plot->SetLineColor(4);
           pltDelta10.plot->SetMarkerStyle(20);
           pltDelta10.plot->Draw("PSame");
           
           TLegend *legdpb = new TLegend(0.769,0.669,0.997,0.996);
           legdpb->SetFillColor(0);
           legdpb->AddEntry( pltDelta260.plot,"B = 250 G","LP");
           legdpb->AddEntry( pltDelta250.plot,"B = 100 G","LP");
           legdpb->AddEntry( pltDelta175.plot,"B = 260 G", "LP");
           legdpb->AddEntry( pltDelta100.plot,"B = 175 G", "LP");
           legdpb->AddEntry( pltDelta50.plot,"B = 50 G", "LP");
           legdpb->AddEntry( pltDelta30.plot,"B = 30 G", "LP");
           legdpb->AddEntry( pltDelta25.plot,"B = 25 G", "LP");
           legdpb->AddEntry( pltDelta10.plot,"B = 10 G", "LP");
           legdpb->Draw("Same");
           
           vlinedraw(70,78,2); //CU
           vlinedraw(79,103,3);//Non coat Cu
           vlinedraw(104,198,4);//DLC
           vlinedraw(199,273,5);//EP Cu
           vlinedraw(304,338,6);//SS
           vlinedraw(339,373,7);//MP Cu
           vlinedraw(374,386,5);//EP Cu
    return;
}
void vlinedraw(double low, double high, int color) {
  double cutlinelx[2], cutlinehx[2];
  double cutlinely[2], cutlinehy[2];
  cutlinely[0]=0;
  cutlinehy[0]=0;
  cutlinely[1]=1e10;
  cutlinehy[1]=1e10;
  cutlinelx[0] = low;
  cutlinelx[1] = low;
  cutlinehx[0] = high;
  cutlinehx[1] = high;
  TGraph *glow = new TGraph(2,cutlinelx,cutlinely);
  TGraph *ghigh = new TGraph(2,cutlinehx,cutlinehy);
  glow->SetLineColor(color);
  ghigh->SetLineColor(color);
  glow->Draw("same");
  ghigh->Draw("same");
  
  return;
}

double getTauhat(int guide){
   char name[500];
   double tauhat = 0;
   double temp[200];
   int counter = 0;
   int i,j,k;
   double test;
   
   for(k = 0;k < 200; k++) temp[i] = 0;
   
   
   //Load the extracted data tree.
   sprintf(name,"%s/obsvbls.root",ODIR);
   TTree *tree = xdload(name,"tree1");

   PLOT xd;
   int entries = tree->GetEntries();
  
   for(i = 0; i < entries; i++){
       xd.get_entry(tree,i);
            
       //if((xd.Tguide == guide) && (xd.rtyp == 1)) {
           test = gettauhat(xd.rnum);
           
           printf("tau %f",test);
       //}
       
       
   }
    
   /* for(j = 0; j < 200; j++){
        if(temp[i] != 0){
            tauhat = tauhat + temp[i];
            counter++;
        }
            
    } 
    
    tauhat = tauhat/counter;
    
    printf("Average Tauhat: %f\n\n", tauhat); */
    
   return tauhat;
}

