//Online Analysis for PPM Depol Runs [A.T. Holley]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Group Numbers:
    0 -> high field, SS guide
    1 -> low field, SS guide
    2 -> high field, EP Cu 
    3 -> low field, EP Cu  
    4 -> high field, DLC Cu
    5 -> low field, DLC Cu 
    6 -> field zero, DLC Cu  */

TCanvas *allcan;
void ppmdepol();
/* extract -> (run#, Tclean, group#, global Tstart, global Tend) */
void extract(int,double,int,double,double);
void plot1(TGraphErrors*,int);
void plot2(TGraphErrors*,int);
void plot3(TgraphErrors*,int);
void plot4(TgarphErrors*,int);
double calc_dpb(double,double,double,double,double,double,double,double,double,double,double);
double calc_edpb(double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double);
int runcnt[10];
double Thld[100][10];
double eThld[100][10];
double nunld[100][10];
double enunld[100][10];
double taudrain[100][10];
double etaudrain[100][10];
double tau1[100][10];
double etau1[100][10];
double nr[100][10];
double enr[100][10];
double dpb[100][10];
double edpb[100][10];
double x[100], y[100], ex[100], ey[100];
double wavg_tau1 = 0;
double w_tau1 = 0;
double wavg_bgrate = 0;
double w_bgrate = 0;

void ppmdepol() {

  int i,j;
  double w;
  for(i=0 ; i < 10 ; i++) runcnt[i]=0;
  allcan = new TCanvas();
  TH1F *hisum = new TH1F("hisum","hisum",110,0,110);

//SS bottle guide, high field [group 0]
  //extract(11510,40,0,0,5000);
  //extract(11512,40,0,0,5000);
  //extract(11513,30,0,0,5000);
  //extract(11515,20,0,0,5000);
  //extract(11516,0.1,0,0,5000);
  //extract(11517,5,0,0,5000);
  //extract(11518,10,0,0,5000);
  
//SS bottle guide, low field [group 1]
  //extract(11519,10,1,0,5000);
  //extract(11520,5,1,0,5000);
  //extract(11521,20,1,0,5000);
  //extract(11522,30,1,0,5000);
  //extract(11525,40,1,0,5000);
  //extract(11526,0.1,1,0,5000);
  
//EP Cu bottle guide, low field [group 3]
  extract(11531,20,3,18.75,3763);
  hisum->Add(swt11531);
  extract(11532,40,3,14.5,4258);
  hisum->Add(swt11532);
  extract(11542,60,3,11.4,3444.5);
  hisum->Add(swt11542);
  //extract(11543,60,3,17,1390.5);
  //hisum->Add(swt11543);
  extract(11549,30,3,92.3,4253);
  hisum->Add(swt11549);
  extract(11548,0.1,3,21.1638,541.018);

//EP Cu bottle guide, high field [group 2]
  extract(11535,20,2,446.703,3505.17);
  hisum->Add(swt11535);
  extract(11537,40,2,15,3946);
  hisum->Add(swt11537);
  extract(11538,60,2,18,3451);
  hisum->Add(swt11538);
  //extract(11540,60,2,10.8,3902);
  //hisum->Add(swt11540);
  //extract(11541,40,2,32,2279);
  //hisum->Add(swt11541);

//DLC Cu bottle guide, high field [group 4]
  extract(11558,20,4,15.3,3658);
  hisum->Add(swt11558);
  extract(11559,40,4,181.563,3830.50);
  hisum->Add(swt11559);
  extract(11560,25,4,42.2773,3786);
  hisum->Add(swt11560);
  extract(11565,0.1,4,47.8,567.877);

//DLC Cu bottle guide, low field [group 5]
  extract(11561,25,5,24.6767,3691.38);
  hisum->Add(swt11561);
  extract(11562,40,5,30.8,3869);
  hisum->Add(swt11562);
  extract(11563,20,5,14.26,3728);
  hisum->Add(swt11563);
  
//DLC Cu bottle guide, field zero [group 6]
  extract(11570,25,6,14.2,794);
  hisum->Add(swt11570);
  extract(11571,40,6,8.37,944.5);
  hisum->Add(swt11571);
  extract(11572,20,6,17.0,599.4);
  hisum->Add(swt11572);
  
//Fit the summed clean/hold spectrum
  new TCanvas(); hisum->SetLineColor(4); hisum->Draw();
  TF1 *fitsum = new TF1("fitsum","exp([0]-x/[1])+exp([2]-x/[3])+[4]",20,40);
  fitsum->SetLineColor(2);  fitsum->SetLineWidth(2);
  fitsum->SetParameter(0,26.5);
  fitsum->SetParameter(1,1.15);
  fitsum->SetParameter(2,8.2);
  fitsum->SetParameter(3,10);
  fitsum->FixParameter(4,0);
  hisum->Fit("fitsum","QLM+","",21,40);
  double cltc1 = fitsum->GetParameter(1);
  double ecltc1 = fitsum->GetParError(1);
  double cltc2 = fitsum->GetParameter(3);
  double ecltc2 = fitsum->GetParError(3);
  double fit_const = fitsum->GetParameter(4);
  double efit_const = fitsum->GetParError(4);
  printf("\nFrom summed clean/hold spectrum...\n");
  printf("  Clean tau1: %f +/- %f\n",cltc1,ecltc1);
  printf("  Clean tau2: %f +/- %f\n",cltc2,ecltc2);
  printf("  const: %f +/- %f\n",fit_const,efit_const);
  printf("  Chi2/dof = %f\n\n",fitsum->GetChisquare()/fitsum->GetNDF());

//Legend for the following plots
  TLegend *leg = new TLegend(0.5,0.67,0.88,0.88);

//Plot group 2
  for(i=0 ; i < runcnt[2] ; i++) {
    x[i] = Thld[i][2];
    ex[i] = eThld[i][2];
    y[i] = nunld[i][2];
    ey[i] = enunld[i][2];
  }
  //----------------------------------//
  /* w = 1/pow(ey[2],2) + 1/pow(ey[3],2);
  y[3] = y[2]/pow(ey[2],2) + y[3]/pow(ey[3],2);
  y[3] = y[3]/w;
  ey[3] = sqrt(1/w);
  y[2] = 0;  ey[2] = 0; */
  //----------------------------------//
  /* w = 1/pow(ey[1],2) + 1/pow(ey[3],2);
  y[3] = y[1]/pow(ey[1],2) + y[3]/pow(ey[3],2);
  y[3] = y[3]/w;
  ey[3] = sqrt(1/w);
  y[1] = 0;  ey[1] = 0; */
  //----------------------------------//
  TGraphErrors *graph1 = new TGraphErrors(runcnt[2],x,y,ex,ey);
  TH1F *haxis1 = new TH1F("haxis1","",150,0,80);
  haxis1->SetTitle("Fits for #hat{#tau}");
  new TCanvas();
  haxis1->GetXaxis()->SetTitle("Holding Time [s]");
  haxis1->GetYaxis()->SetTitle("BS Unload Cnts/Load Cnts");
  haxis1->GetYaxis()->SetTitleOffset(1.2);
  haxis1->GetYaxis()->SetRangeUser(1e-5,1.0);
  haxis1->SetStats(kFALSE);
  haxis1->Draw();
  graph1->SetMarkerColor(4);
  graph1->SetMarkerStyle(21);
  graph1->SetMarkerSize(2);
  graph1->Draw("P");
  leg->AddEntry(graph1,"EP Cu, B0 high","LP");
  
//Plot group 3
  for(i=0 ; i < runcnt[3] ; i++) {
    x[i] = Thld[i][3];
    ex[i] = eThld[i][3];
    y[i] = nunld[i][3];
    ey[i] = enunld[i][3];
  }
  //------------------------------------//
  /* w = 1/pow(ey[2],2) + 1/pow(ey[3],2);
  y[3] = y[2]/pow(ey[2],2) + y[3]/pow(ey[3],2);
  y[3] = y[3]/w;
  ey[3] = sqrt(1/w);
  y[2] = 0;  ey[2] = 0; */
  //------------------------------------//
  TGraphErrors *graph2 = new TGraphErrors(runcnt[3],x,y,ex,ey);
  graph2->SetMarkerColor(42);
  graph2->SetMarkerStyle(22);
  graph2->SetMarkerSize(2);
  graph2->Draw("PSame");
  leg->AddEntry(graph2,"EP Cu, B0 low","LP");
  
//Plot group 4
  for(i=0 ; i < runcnt[4] ; i++) {
    x[i] = Thld[i][4];
    ex[i] = eThld[i][4];
    y[i] = nunld[i][4];
    ey[i] = enunld[i][4];
  }
  TGraphErrors *graph5 = new TGraphErrors(runcnt[4],x,y,ex,ey);
  graph5->SetMarkerColor(2);
  graph5->SetMarkerStyle(20);
  graph5->SetMarkerSize(2);
  graph5->Draw("PSame");
  leg->AddEntry(graph5,"DLC Cu, B0 high","LP");
  
  //Plot group 5
  for(i=0 ; i < runcnt[5] ; i++) {
    x[i] = Thld[i][5];
    ex[i] = eThld[i][5];
    y[i] = nunld[i][5];
    ey[i] = enunld[i][5];
  }
  TGraphErrors *graph6 = new TGraphErrors(runcnt[5],x,y,ex,ey);
  graph6->SetMarkerColor(3);
  graph6->SetMarkerStyle(23);
  graph6->SetMarkerSize(2);
  graph6->Draw("PSame");
  leg->AddEntry(graph6,"DLC Cu, B0 low","LP");
  
  //Plot group 6
  for(i=0 ; i < runcnt[6] ; i++) {
    x[i] = Thld[i][6];
    ex[i] = eThld[i][6];
    y[i] = nunld[i][6];
    ey[i] = enunld[i][6];
  }
  TGraphErrors *graph7 = new TGraphErrors(runcnt[6],x,y,ex,ey);
  graph7->SetMarkerColor(kCyan);
  graph7->SetMarkerStyle(29);
  graph7->SetMarkerSize(2);
  graph7->Draw("PSame");
  leg->AddEntry(graph7,"DLC Cu, field zero","LP");

  leg->Draw();
  
//Fit group 4
  TF1 *fit1 = new TF1("fit1","exp([0]+[1]*x)+[2]",20,40);
  fit1->SetParameter(0,1);
  fit1->SetParameter(1,-1/30.);
  fit1->SetParameter(2,0);
  fit1->FixParameter(2,0);
  fit1->SetLineColor(2);  fit1->SetLineWidth(2);
  graph5->Fit("fit1","QLM","",20,40);
  double trapT1 = -1/fit1->GetParameter(1);
  double etrapT1 = trapT1**2 * fit1->GetParError(1);
  printf("\nTrapped Lifetime (group 4): %f +/- %f\n",trapT1,etrapT1);
  printf("Chi2/dof = %f\n",fit1->GetChisquare()/fit1->GetNDF());
  
  //Fit group 2
  TF1 *fit2 = new TF1("fit2","exp([0]+[1]*x)+[2]",20,60);
  fit2->SetParameter(0,1);
  fit2->SetParameter(1,-1/30.);
  fit2->SetParameter(2,0);
  fit2->FixParameter(2,0);
  fit2->SetLineColor(4);  fit2->SetLineWidth(2);
  graph1->Fit("fit2","QLM","",20,60);
  //graph1->Fit("fit2","QLM","",20,40);
  double trapT2 = -1/fit2->GetParameter(1);
  double etrapT2 = trapT2**2 * fit2->GetParError(1);
  printf("\nTrapped Lifetime (group 2): %f +/- %f\n",trapT2,etrapT2);
  printf("Chi2/dof = %f\n",fit2->GetChisquare()/fit2->GetNDF());
  
  //Fit group 3
  TF1 *fit3 = new TF1("fit3","exp([0]+[1]*x)+[2]",20,60);
  fit3->SetParameter(0,1);
  fit3->SetParameter(1,-1/30.);
  fit3->SetParameter(2,0);
  fit3->FixParameter(2,0);
  fit3->SetLineColor(42);  fit3->SetLineWidth(2);
  graph2->Fit("fit3","QLM","",20,60);
  //graph2->Fit("fit3","QLM","",20,40);
  double trapT3 = -1/fit3->GetParameter(1);
  double etrapT3 = trapT3**2 * fit3->GetParError(1);
  printf("\nTrapped Lifetime (group 3): %f +/- %f\n",trapT3,etrapT3);
  printf("Chi2/dof = %f\n",fit3->GetChisquare()/fit3->GetNDF());
  
  //Fit group 5
  TF1 *fit4 = new TF1("fit4","exp([0]+[1]*x)+[2]",20,40);
  fit4->SetParameter(0,1);
  fit4->SetParameter(1,-1/30.);
  fit4->SetParameter(2,0);
  fit4->FixParameter(2,0);
  fit4->SetLineColor(3);  fit4->SetLineWidth(2);
  graph6->Fit("fit4","QLM","",20,40);
  double trapT4 = -1/fit4->GetParameter(1);
  double etrapT4 = trapT4**2 * fit4->GetParError(1);
  printf("\nTrapped Lifetime (group 5): %f +/- %f\n",trapT4,etrapT4);
  printf("Chi2/dof = %f\n",fit4->GetChisquare()/fit4->GetNDF());
  
  //Fit group 6
  TF1 *fit5 = new TF1("fit5","exp([0]+[1]*x)+[2]",20,40);
  fit5->SetParameter(0,1);
  fit5->SetParameter(1,-1/30.);
  fit5->SetParameter(2,0);
  fit5->FixParameter(2,0);
  fit5->SetLineColor(kCyan);  fit5->SetLineWidth(2);
  graph7->Fit("fit5","QLM","",20,40);
  double trapT5 = -1/fit5->GetParameter(1);
  double etrapT5 = trapT5**2 * fit5->GetParError(1);
  printf("\nTrapped Lifetime (group 6): %f +/- %f\n",trapT5,etrapT5);
  printf("Chi2/dof = %f\n",fit5->GetChisquare()/fit5->GetNDF());

  //Weighted average of selected fit results from above...
   w = 1/pow(etrapT4,2) + 1/pow(etrapT1,2) + 1/pow(etrapT3,2);
   w = w + 1/pow(etrapT2,2) + 1/pow(etrapT5,2);
  double trapT = trapT4/pow(etrapT4,2) + trapT1/pow(etrapT1,2);
  trapT = trapT + trapT2/pow(etrapT2,2) + trapT5/pow(etrapT5,2);
  trapT = trapT + trapT3/pow(etrapT3,2);
  trapT = trapT/w;
  double etrapT = sqrt(1/w);
  printf("\nWeighted average of g2, g3, g4, g5, and g6: ");
  printf("%f +/- %f\n",trapT,etrapT);
  
  //Report weighted average of tau1...
  printf("\nWeighted average of selected runs' tau1: ");
  printf("%f +/- %f\n",wavg_tau1/w_tau1,sqrt(1/w_tau1));
  
  //Report weighted average of bgrates...
  printf("\nWeighted average of selected runs' bgrate: ");
  printf("%f +/- %f\n",wavg_bgrate/w_bgrate,sqrt(1/w_bgrate));
  printf("  **No cycle normalization!\n\n");

//Plot clean/hold fit results
  new TCanvas("cdrt");
  TH1F *haxis2 = new TH1F("haxis2","haxis2",150,0,80);
  haxis2->GetYaxis()->SetTitle("Clean/Hold Drain Time [s]");
  haxis2->GetYaxis()->SetRangeUser(0,5);
  haxis2->GetYaxis()->SetTitleOffset(1.2);
  haxis2->GetXaxis()->SetTitle("Holding Time [s]");
  haxis2->SetStats(kFALSE);
  haxis2->Draw();
  TGraphErrors *plot0;
  TGraphErrors *plot1;
  TGraphErrors *plot2;
  TGraphErrors *plot3;
  TGraphErrors *plot4;
  TGraphErrors *plot5;
  TGraphErrors *plot6;
  TLegend *leg2 = new TLegend(0.5,0.67,0.88,0.88);
  //plot1(plot0,0);  leg2->AddEntry(plot0,"group 0","LP");
  //plot1(plot1,1);  leg2->AddEntry(plot1,"group 1","LP");
  plot1(plot2,2);  leg2->AddEntry("plot2","EP Cu, High Field (g2)","LP");
  plot1(plot3,3);  leg2->AddEntry("plot3","EP Cu, Low Field (g3)","LP");
  plot1(plot4,4);  leg2->AddEntry("plot4","DLC Cu, High Field (g4)","LP");
  plot1(plot5,5);  leg2->AddEntry("plot5","DLC Cu, Low Field (g5)","LP");
  plot1(plot6,6);  leg2->AddEntry("plot6","DLC Cu, Field Zero (g6)","LP");
  leg2->Draw();
  
//Plot tau1 (unload) fit results)
  new TCanvas("ctau1");
  TH1F *haxis3 = new TH1F("haxis3","haxis3",150,0,80);
  haxis3->GetYaxis()->SetTitle("Unload Drain Time [s]");
  haxis3->GetYaxis()->SetRangeUser(0,5);
  haxis3->GetYaxis()->SetTitleOffset(1.2);
  haxis3->GetXaxis()->SetTitle("Holding Time [s]");
  haxis3->SetStats(kFALSE);
  haxis3->Draw();
  TLegend *leg3 = new TLegend(0.5,0.67,0.88,0.88);
  //plot2(plot0,0);  leg3->AddEntry(plot0,"group 0","LP");
  //plot2(plot1,1);  leg3->AddEntry(plot1,"group 1","LP");
  plot2(plot2,2);  leg3->AddEntry("plot2","EP Cu, High Field (g2)","LP");
  plot2(plot3,3);  leg3->AddEntry("plot3","EP Cu, Low Field (g3)","LP");
  plot2(plot4,4);  leg3->AddEntry("plot4","DLC Cu, High Field (g4)","LP");
  plot2(plot5,5);  leg3->AddEntry("plot5","DLC Cu, Low Field (g5)","LP");
  plot2(plot6,6);  leg3->AddEntry("plot6","DLC Cu, Field Zero (g6)","LP");
  leg3->Draw();
  
  //Plot nr (clean/hold) fit results
  new TCanvas("cnr");
  TH1F *haxis4 = new TH1F("haxis4","haxis4",150,0,80);
  haxis4->GetYaxis()->SetTitle("nr");
  haxis4->GetYaxis()->SetRangeUser(0,1000);
  haxis4->GetYaxis()->SetTitleOffset(1.3);
  haxis4->GetXaxis()->SetTitle("Holding Time [s]");
  haxis4->SetStats(kFALSE);
  haxis4->Draw();
  TLegend *leg4 = new TLegend(0.5,0.67,0.88,0.88);
  //plot3(plot0,0);  leg4->AddEntry(plot0,"group 0","LP");
  //plot3(plot1,1);  leg4->AddEntry(plot1,"group 1","LP");
  plot3(plot2,2);  leg4->AddEntry("plot2","EP Cu, High Field (g2)","LP");
  plot3(plot3,3);  leg4->AddEntry("plot3","EP Cu, Low Field (g3)","LP");
  plot3(plot4,4);  leg4->AddEntry("plot4","DLC Cu, High Field (g4)","LP");
  plot3(plot5,5);  leg4->AddEntry("plot5","DLC Cu, Low Field (g5)","LP");
  plot3(plot6,6);  leg4->AddEntry("plot6","DLC Cu, Field Zero (g6)","LP");
  leg4->Draw();
  
  //Plot extracted dpb
  new TCanvas("cdpb");
  TH1F *haxis5 = new TH1F("haxis5","",150,0,80);
  haxis5->SetTitle("Estimated Depolarization per Bounce");
  haxis5->SetLineColor(0);
  haxis5->GetYaxis()->SetTitle("dpb");
  haxis5->GetYaxis()->SetRangeUser(-0.4e-6,0.03e-3);
  haxis5->GetYaxis()->SetTitleOffset(1.2);
  haxis5->GetXaxis()->SetTitle("Holding Time [s]");
  haxis5->SetStats(kFALSE);
  haxis5->Draw();
  TLegend *leg5 = new TLegend(0.5,0.67,0.88,0.88);
  //plot4(plot0,0);  leg5->AddEntry(plot0,"group 0","LP");
  //plot4(plot1,1);  leg5->AddEntry(plot1,"group 1","LP");
  plot4(plot2,2);  leg5->AddEntry("plot2","EP Cu, High Field (g2)","LP");
  plot4(plot3,3);  leg5->AddEntry("plot3","EP Cu, Low Field (g3)","LP");
  plot4(plot4,4);  leg5->AddEntry("plot4","DLC Cu, High Field (g4)","LP");
  plot4(plot5,5);  leg5->AddEntry("plot5","DLC Cu, Low Field (g5)","LP");
  //plot4(plot6,6);  leg5->AddEntry("plot6","DLC Cu, Field Zero (g6)","LP");
  leg5->Draw();
  //cdpb->SetLogy();
  
  
  return;
  
}


void plot1(TGraphErrors* plot, int group) {
  int i;
  char string[50];
  
  for(i=0 ; i < runcnt[group] ; i++) {
    x[i] = Thld[i][group];
    ex[i] = eThld[i][group];
    y[i] = taudrain[i][group];
    ey[i] = etaudrain[i][group];
  }
  
  plot = new TGraphErrors(runcnt[group],x,y,ex,ey);
  sprintf(string,"plot%d",group);
  plot->SetName(string);
  if(group == 2) {plot->SetMarkerColor(2); plot->SetMarkerStyle(20);}
  if(group == 3) {plot->SetMarkerColor(3); plot->SetMarkerStyle(21);}
  if(group == 4) {plot->SetMarkerColor(4); plot->SetMarkerStyle(22);}
  if(group == 5) {plot->SetMarkerColor(6); plot->SetMarkerStyle(23);}
  if(group == 6) {plot->SetMarkerColor(7); plot->SetMarkerStyle(29);}
  plot->SetMarkerSize(1.5);
  plot->Draw("PSame");

  return;
}


void plot2(TGraphErrors* plot, int group) {
  int i;
  char string[50];
  
  for(i=0 ; i < runcnt[group] ; i++) {
    x[i] = Thld[i][group];
    ex[i] = eThld[i][group];
    y[i] = tau1[i][group];
    ey[i] = etau1[i][group];
  }
  
  plot = new TGraphErrors(runcnt[group],x,y,ex,ey);
  sprintf(string,"plot%d",group);
  plot->SetName(string);
  if(group == 2) {plot->SetMarkerColor(2); plot->SetMarkerStyle(20);}
  if(group == 3) {plot->SetMarkerColor(3); plot->SetMarkerStyle(21);}
  if(group == 4) {plot->SetMarkerColor(4); plot->SetMarkerStyle(22);}
  if(group == 5) {plot->SetMarkerColor(6); plot->SetMarkerStyle(23);}
  if(group == 6) {plot->SetMarkerColor(7); plot->SetMarkerStyle(29);}
  plot->SetMarkerSize(1.5);
  plot->Draw("PSame");

  return;
}

void plot3(TGraphErrors* plot, int group) {
  int i;
  char string[50];
  
  for(i=0 ; i < runcnt[group] ; i++) {
    x[i] = Thld[i][group];
    ex[i] = eThld[i][group];
    y[i] = nr[i][group];
    ey[i] = enr[i][group];
  }
  
  plot = new TGraphErrors(runcnt[group],x,y,ex,ey);
  sprintf(string,"plot%d",group);
  plot->SetName(string);
  if(group == 2) {plot->SetMarkerColor(2); plot->SetMarkerStyle(20);}
  if(group == 3) {plot->SetMarkerColor(3); plot->SetMarkerStyle(21);}
  if(group == 4) {plot->SetMarkerColor(4); plot->SetMarkerStyle(22);}
  if(group == 5) {plot->SetMarkerColor(6); plot->SetMarkerStyle(23);}
  if(group == 6) {plot->SetMarkerColor(7); plot->SetMarkerStyle(29);}
  plot->SetMarkerSize(1.5);
  plot->Draw("PSame");

  return;
}


void plot4(TGraphErrors* plot, int group) {
  int i;
  char string[50];
  
  for(i=0 ; i < runcnt[group] ; i++) {
    x[i] = Thld[i][group];
    ex[i] = eThld[i][group];
    y[i] = dpb[i][group];
    ey[i] = edpb[i][group];
  }

  /* if(group == 3) {
    printf("\nAveraging dpb group# %d, clean/hold time ",group);
    printf("%fs and %fs...\n",x[2],x[3]);
    double w = 1/pow(ey[2],2) + 1/pow(ey[3],2);
    double avg = y[2]/pow(ey[2],2) + y[3]/pow(ey[3],2);
    avg = avg/w;
    double eavg = sqrt(1/w);
    printf("  Weighted average: %e +/- %e\n\n",avg,eavg);
    y[2] = avg; ey[2] = eavg; y[3] = avg; ey[3] = eavg;
  } */
  
  /* if(group == 2) {
    printf("\nAveraging dpb group# %d, clean/hold time ",group);
    printf("%fs and %fs...\n",x[1],x[3]);
    double w = 1/pow(ey[1],2) + 1/pow(ey[3],2);
    double avg = y[1]/pow(ey[1],2) + y[3]/pow(ey[3],2);
    avg = avg/w;
    double eavg = sqrt(1/w);
    printf("  Weighted average: %e +/- %e\n\n",avg,eavg);
    y[1] = avg; ey[1] = eavg; y[3] = avg; ey[3] = eavg;
    
    printf("\nAveraging dpb group# %d, clean/hold time ",group);
    printf("%fs and %fs...\n",x[2],x[3]);
    double w = 1/pow(ey[2],2) + 1/pow(ey[3],2);
    double avg = y[2]/pow(ey[2],2) + y[3]/pow(ey[3],2);
    avg = avg/w;
    double eavg = sqrt(1/w);
    printf("  Weighted average: %e +/- %e\n\n",avg,eavg);
    y[2] = avg; ey[2] = eavg; y[3] = avg; ey[3] = eavg;
  } */

  plot = new TGraphErrors(runcnt[group],x,y,ex,ey);
  sprintf(string,"plot%d",group);
  plot->SetName(string);
  plot->SetMarkerColor(1+group);
  if(group == 2) {plot->SetMarkerColor(2); plot->SetMarkerStyle(20);}
  if(group == 3) {plot->SetMarkerColor(3); plot->SetMarkerStyle(21);}
  if(group == 4) {plot->SetMarkerColor(4); plot->SetMarkerStyle(22);}
  if(group == 5) {plot->SetMarkerColor(6); plot->SetMarkerStyle(23);}
  if(group == 6) {plot->SetMarkerColor(7); plot->SetMarkerStyle(29);}
  plot->SetMarkerSize(1.5);
  plot->Draw("PSame");

  return;
}


int extract(int rnum, double hold, int grpnum, double Tstart, double Tend) {

  char string[100];
  char stringe[200];
  double par1,epar1;
  int index, flag;
  
  double Tclnhldi = 20.0;
  double Tclnhldf = 20.0 + hold - 1.0;
  double Tunldi = 20.0 + hold;
  double Tunldf = 20.0 + hold + 10.0; //End of unload for getting counts.
  double Tunldf_fit = 20.0 + hold + 15.0; //End of unload for fitting.
  double spbin = 1.0;
  
  runcnt[grpnum]++;
  index = runcnt[grpnum] - 1;
  Thld[index][grpnum] = hold;
  eThld[index][grpnum] = 0.0;

  printf("\n--------------------------------------\n");
  printf("Run#: %d | Group#: %d | Clean/Hold: %f\n",rnum,grpnum,hold);
  
  sprintf(string,"cswt%d",rnum);
  TCanvas *canswt = new TCanvas(string,string);

  sprintf(string,"full%d.root",rnum);
  TFile *rfile = new TFile(string);
  canswt->cd();
  
  sprintf(string,"swt%d",rnum);
  TH1F *swt = new TH1F(string,string,110,0,110);
  swt->SetLineColor(index+3);
  swt->SetTitle("ppm_depol");
  swt->GetXaxis()->SetTitle("Time [s]");
  swt->GetYaxis()->CenterTitle();
  swt->GetYaxis()->SetTitle("Main Counts");
  swt->GetYaxis()->SetTitleOffset(1.2);
  //swt->SetStats(kFALSE);
  sprintf(string,"S8200/1e6>>swt%d",rnum);
  sprintf(stringe,"Sis00==2052 && Pdc311 > 400 && Pdc311 < 1000 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  if(rnum == 11541) sprintf(stringe,"Sis00==2052 && Pdc311 > 200 && Pdc311 < 600 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  if(rnum == 11532) sprintf(stringe,"Sis00==2052 && Pdc311 > 400 && Pdc311 < 1000 && S83028/1e6 > %f && S83029 < 207",Tstart);
  if(rnum == 11542 || rnum == 11543 ) sprintf(stringe,"Sis00==2052 && Pdc311 > 100 && Pdc311 < 330 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  if(rnum == 11549) sprintf(stringe,"Sis00==2052 && Pdc311 > 500 && Pdc311 < 1200 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  if(rnum >= 11558) sprintf(stringe,"Sis00==2052 && Pdc311 > 155 && Pdc311 < 350 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  if(rnum >= 11570) sprintf(stringe,"Sis00==2052 && Pdc311 > 450 && Pdc311 < 1200 && S83028/1e6 > %f && S83028/1e6 < %f",Tstart,Tend);
  canswt->cd();
  h1->Draw(string,stringe);

//Fitting clean/hold for taudrain (tau_twidle)
  sprintf(string,"fit%d",rnum);
  TF1 *fit = new TF1(string,"exp([0]+[1]*x)+[2]",Tclnhldi,Tclnhldf);
  fit->SetParameter(0,10);
  fit->SetParameter(1,-1/2.0);
  fit->SetParameter(2,0.1);
  fit->SetLineColor(3);  fit->SetLineWidth(2);
  swt->Fit(string,"LMQ+","",Tclnhldi,Tclnhldf);
  par1 = fit->GetParameter(1);
  if(par1 !=0) double drainT = -1/par1;
  else double drainT = 0;
  double edrainT = drainT**2 * fit->GetParError(1);
  printf("tau_twidle: %f +/- %f\n",drainT,edrainT);
  taudrain[index][grpnum] = drainT;
  etaudrain[index][grpnum] = edrainT;

//Fitting unload curve for tau1
  sprintf(string,"fit_tau1%d",rnum);
  fit = new TF1(string,"exp([0]-x/[1])+[2]",Tunldi,Tunldf_fit);
  double parest = swt->Integral(Tunldi/spbin - 1.0,Tunldi/spbin + 1.0)/2;
  parest = log(parest) + Tunldi;
  fit->SetParameter(0,parest);
  fit->SetParameter(1,1.0);
  fit->SetParameter(2,0);
  fit->SetLineColor(3);  fit->SetLineWidth(2);
  swt->Fit(string,"LMQ+","",Tunldi+1,Tunldf_fit);
  par1 = fit->GetParameter(1);
  if(par1 != 0) drainT = par1;
  else drainT = 0;
  edrainT = fit->GetParError(1);
  printf("tau1: ");
  printf("%f +/- %f [const: %f]\n",drainT,edrainT,fit->GetParameter(2));
  tau1[index][grpnum] = drainT;
  etau1[index][grpnum] = edrainT;
  flag = 1;
  if(rnum == 11562 || rnum == 11559 || rnum == 11541) flag = 0;
  if(rnum == 11540 || rnum == 11538 || rnum == 11535) flag = 0;
  if(rnum == 11549) flag = 0;
  if(flag == 1) {
    w_tau1 = w_tau1 + 1/pow(edrainT,2);
    wavg_tau1 = wavg_tau1 + drainT/pow(edrainT,2);
  }
    
/*------------------------*/
  double tau_twidle = taudrain[index][grpnum];  tau_twidle = 1.113051;
  double etau_twidle = etaudrain[index][grpnum];  etau_twidle = 0.041307;
  double tau_hat = 19.190029;  double etau_hat = 1.286095;
  double tau_1 = tau1[index][grpnum];  tau_1 = 1.052982;
  double etau_1 = etau1[index][grpnum];  etau_1 = 0.011784;
  double tau_2 = tau1[index][grpnum];  tau_2 = 1.052982;
  double etau_2 = etau1[index][grpnum];  etau_2 = 0.011784;
  double beta = 45;
  double delta = 20;
/*------------------------*/

//Fitting clean/hold for n_r
  sprintf(string,"fit_nr%d",rnum);
  sprintf(stringe,"[0]*exp([1]*(x-%f)) + [2]*exp([3]*(x-%f))+[4]",Tclnhldi,Tclnhldi);
  fit = new TF1(string,stringe,Tclnhldi,Tclnhldf);
  fit->SetParameter(0,1000);
  fit->FixParameter(1,-1/tau_twidle);
  fit->SetParameter(2,0);
  fit->FixParameter(3,-1/tau_hat);
  fit->SetParameter(4,0.1);
  fit->SetLineColor(50);  fit->SetLineWidth(2);
  swt->Fit(string,"LMQ+","",Tclnhldi,Tclnhldf);
  par1 = fit->GetParameter(0);
  epar1 = fit->GetParError(0);
  printf("nr: %f +/- %f\n",par1,epar1);
  nr[index][grpnum] = par1;
  enr[index][grpnum] = epar1;

/*------------------------*/
  double n_r = nr[index][grpnum];
/*------------------------*/

//Fitting for n_w, extraction of dpb
  sprintf(string,"fit_nw%d",rnum);
  sprintf(stringe,"[0]*exp(-%f/%f)*exp([1]*(x-%f)) + [2]*exp(-%f/%f)*exp([3]*(x-%f)) + [4]",hold,tau_twidle,Tunldi+1,hold,tau_hat,Tunldi+1);
  fit = new TF1(string,stringe,Tunldi,Tunldf_fit);
  fit->SetParameter(0,n_r);
  fit->FixParameter(1,-1/tau_1);
  fit->SetParameter(2,0.002*n_r);
  fit->FixParameter(3,-1/tau_2);
  fit->FixParameter(4,0.1);  fit->SetParLimits(4,0,5);
  fit->SetLineColor(50);  fit->SetLineWidth(2);
  swt->Fit(string,"LMQ+","",Tunldi+1,Tunldf_fit);
  par1 = fit->GetParameter(2);
  epar1 = fit->GetParError(2);
  printf("nw: %f +/- %f\n",par1,epar1);
  double nw = par1/exp(-hold/tau_hat);
  double enw = 0;
  double a = pow(3.0/16.0,2)*3.1415927;
  double A = pow(2.75/2.0,2)*3.1415927;
  double temp = (a/A) * (nw/nr[index][grpnum]) * 1/(45*20);
  printf("dpb from fit values ~ %e\n",temp);

//Extraction of load, clean/hold, and unload counts
  double Tbg = 40.0 + hold - Tunldf;
  int bgbini = (int)(40.0 + hold - Tbg)/spbin + 1;
  int bgbinf = (int)(40.0 + hold)/spbin;
  int clnhldbini = (int)Tclnhldi/spbin + 1;
  int clnhldbinf = (int)Tclnhldf/spbin;
  double clnhldDelta = Tclnhldf - Tclnhldi;
  int unldbini = (int)Tunldi/spbin + 1;
  int unldbinf = (int)Tunldf/spbin;
  double unldDelta = Tunldf - Tunldi;
  double bgcnts = swt->Integral(bgbini,bgbinf);
  double ebgcnts = sqrt(bgcnts);
  if(bgcnts != 0) double pebgcnts = ebgcnts/bgcnts;
  else double pebgcnts = 0;
  double bgrate = bgcnts/(double)Tbg;
  double ebgrate = ebgcnts/Tbg;
  /* Calculating the weighted average of the backgrounds */
  flag = 1;
  if(rnum == 11565 || rnum == 11548 || rnum == 11540) flag = 0;
  if(rnum == 11541) flag = 0;
  if(flag == 1) {
    w_bgrate = w_bgrate + 1/pow(ebgrate,2);
    wavg_bgrate = wavg_bgrate + bgrate/pow(ebgrate,2);
  }
  /* --------------------------------------------------- */
  double BGrate = bgrate;    //The actual background rate used below for
  double eBGrate = ebgrate; //background subtraction.
  //BGrate = 0.260870;  eBGrate = 0.039173;
  //if(rnum == 11541) {BGrate = 0; eBGrate = 0;}
  //if(rnum == 11541) printf("\n  *** BG rate manually set to zero ***  \n");
  double clnhldcnts = swt->Integral(clnhldbini,clnhldbinf);
  double eclnhldcnts = sqrt(clnhldcnts);
  double bsclnhldcnts = clnhldcnts - BGrate*clnhldDelta;
  double ebsclnhldcnts = sqrt(pow(eclnhldcnts,2) + pow(eBGrate*clnhldDelta,2));
  double unldcnts = swt->Integral(unldbini,unldbinf);
  double eunldcnts = sqrt(unldcnts);
  if(unldcnts != 0) double peunldcnts = eunldcnts/unldcnts;
  else double peunldcnts = 0;
  double bsunldcnts = unldcnts - BGrate*unldDelta;
  if(bsunldcnts == 0) bsunldcnts = 1e-6;
  double ebsunldcnts = sqrt(pow(eunldcnts,2) + pow(eBGrate*unldDelta,2));
  double ldcnts = swt->Integral(0/spbin+1,20.0/spbin);
  double eldcnts = sqrt(ldcnts);
  double peldcnts = eldcnts/ldcnts;
  nunld[index][grpnum]=bsunldcnts/ldcnts;
  double penunld = sqrt(peunldcnts**2 + peldcnts**2);
  enunld[index][grpnum] = penunld*nunld[index][grpnum];
  
  double var[18];
  var[0]  =  clnhldDelta;
  var[1]  =  unldDelta;
  var[2]  =  tau_twidle;    var[3]  =  etau_twidle;
  var[4]  =  tau_hat;       var[5]  =  etau_hat;
  var[6]  =  tau_1;         var[7]  =  etau_1;
  var[8]  =  tau_2;         var[9]  =  etau_2;
  var[10] =  A;             var[11] =  a;
  var[12] =  bsclnhldcnts;  var[13] =  ebsclnhldcnts;
  var[14] =  bsunldcnts;    var[15] =  ebsunldcnts;
  var[16] =  beta;          var[17] =  delta;
  
  dpb[index][grpnum] = calc_dpb(var[0],var[1],var[2],var[4],var[6],var[8],var[10],var[11],var[12],var[14],var[16],var[17]);
  edpb[index][grpnum] = calc_edpb(var[0],var[1],var[2],var[3],var[4],var[5],var[6],var[7],var[8],var[9],var[10],var[11],var[12],var[13],var[14],var[15],var[16],var[17]);
  
  printf("\nBG Counts over %fs (bins %d -> %d): %f +/- %f",Tbg,bgbini,bgbinf,bgcnts,ebgcnts);
  printf("\nLoad Counts: %f +/- %f",ldcnts,eldcnts);
  printf("\nClean/Hold Counts over %fs (bins %d -> %d): %f +/- %f",clnhldDelta,clnhldbini,clnhldbinf,clnhldcnts,eclnhldcnts);
  printf("\nBS Clean/Hold Counts: %f +/- %f",bsclnhldcnts,ebsclnhldcnts);
  printf("\nUnload Counts over %fs (bins %d -> %d): %f +/- %f",unldDelta,unldbini,unldbinf,unldcnts,eunldcnts);
  printf("\nBS Unload Counts: %f +/- %f",bsunldcnts,ebsunldcnts);
  printf("\nNormed BS Unload Counts: %f +/- %f",nunld[index][grpnum],enunld[index][grpnum]);
  printf("\n\n  --> depol per bounce:");
  printf(" %e +/- %e\n\n",dpb[index][grpnum],edpb[index][grpnum]);
  
  
  allcan->cd();  swt->Draw("Same");
  canswt->cd();
  
  return;
}


double calc_dpb(double Tclean, double Tunload, double tautwidle, double tauhat, double tau1, double tau2, double A, double a, double sumclean, double sumunload, double beta, double delta) {

  double E = exp(1);

  return (pow(E,Tclean/tauhat - Tclean/tautwidle - 
       Tunload/tau1 + Tunload/tau2)*
     (A*(-1 + pow(E,Tunload/tau1))*sumclean*
        tau1 - a*pow(E,Tunload/tau1)*
        (-1 + pow(E,Tclean/tautwidle))*
        sumunload*tautwidle))/
   (beta*delta*(-(A*
          (-1 + pow(E,Tunload/tau2))*
          sumclean*tau2) + 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*
        sumunload*tauhat));
}


double calc_edpb(double Tclean, double Tunload, double tautwid, double etautwid, double tauhat, double etauhat, double tau1, double etau1, double tau2, double etau2, double A, double a, double sumclean, double esumclean, double sumunload, double esumunload, double beta, double delta) {

  double E = exp(1);

  double e1 = (pow(E,(-2*Tunload)/tau1 + (2*Tunload)/tau2 + 
       (2*Tclean)/tauhat - (2*Tclean)/tautwid)*
     pow(etautwid,2)*
     pow(A*(-1 + pow(E,Tunload/tau1))*
        Tclean*sumclean*tau1 + 
       a*pow(E,Tunload/tau1)*sumunload*tautwid*
        (Tclean + tautwid - 
          pow(E,Tclean/tautwid)*tautwid),2))/
   (pow(beta,2)*pow(delta,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,2)*pow(tautwid,4));

  double e2 = (pow(E,(-2*Tunload)/tau1 + (2*Tunload)/tau2 + 
       (2*Tclean)/tauhat - (2*Tclean)/tautwid)*
     pow(etauhat,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        Tclean*sumclean*tau2 + 
       a*pow(E,Tunload/tau2)*sumunload*tauhat*
        (Tclean + tauhat - 
          pow(E,Tclean/tauhat)*tauhat),2)*
     pow(A*(-1 + pow(E,Tunload/tau1))*
        sumclean*tau1 - 
       a*pow(E,Tunload/tau1)*
        (-1 + pow(E,Tclean/tautwid))*sumunload*
        tautwid,2))/
   (pow(beta,2)*pow(delta,2)*pow(tauhat,4)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,4));

  double e3 = (pow(A,2)*pow(E,
      (-2*Tunload)/tau1 + (2*Tunload)/tau2 + 
       (2*Tclean)/tauhat - (2*Tclean)/tautwid)*
     pow(etau1,2)*pow(sumclean,2)*
     pow(Tunload + tau1 - 
       pow(E,Tunload/tau1)*tau1,2))/
   (pow(beta,2)*pow(delta,2)*pow(tau1,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,2));

  double e4 = (pow(A,2)*pow(E,
      (-2*Tunload)/tau1 + (2*Tunload)/tau2 + 
       (2*Tclean)/tauhat - (2*Tclean)/tautwid)*
     pow(etau2,2)*pow(sumclean,2)*
     pow(Tunload + tau2 - 
       pow(E,Tunload/tau2)*tau2,2)*
     pow(A*(-1 + pow(E,Tunload/tau1))*
        sumclean*tau1 - 
       a*pow(E,Tunload/tau1)*
        (-1 + pow(E,Tclean/tautwid))*sumunload*
        tautwid,2))/
   (pow(beta,2)*pow(delta,2)*pow(tau2,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,4));

  double e5 = (pow(a,2)*pow(A,2)*
     pow(E,(-2*Tunload)/tau1 + 
       (2*Tunload)/tau2 + (2*Tclean)/tauhat - 
       (2*Tclean)/tautwid)*pow(esumclean,2)*
     pow(sumunload,2)*
     pow(pow(E,Tunload*(1/tau1 + 1/tau2))*tau1*
        tauhat - pow(E,Tunload/tau2)*tau1*tauhat - 
       pow(E,Tunload*(1/tau1 + 1/tau2) + 
          Tclean/tauhat)*tau1*tauhat + 
       pow(E,Tunload/tau2 + Tclean/tauhat)*tau1*
        tauhat + pow(E,Tunload/tau1)*tau2*tautwid - 
       pow(E,Tunload*(1/tau1 + 1/tau2))*tau2*
        tautwid + pow(E,
         Tunload/tau1 + Tclean/tautwid)*
        (-1 + pow(E,Tunload/tau2))*tau2*tautwid,2)
     )/(pow(beta,2)*pow(delta,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,4));

  double e6 = (pow(a,2)*pow(A,2)*
     pow(E,(-2*Tunload)/tau1 + 
       (2*Tunload)/tau2 + (2*Tclean)/tauhat - 
       (2*Tclean)/tautwid)*pow(esumunload,2)*
     pow(sumclean,2)*
     pow(pow(E,Tunload*(1/tau1 + 1/tau2))*tau1*
        tauhat - pow(E,Tunload/tau2)*tau1*tauhat - 
       pow(E,Tunload*(1/tau1 + 1/tau2) + 
          Tclean/tauhat)*tau1*tauhat + 
       pow(E,Tunload/tau2 + Tclean/tauhat)*tau1*
        tauhat + pow(E,Tunload/tau1)*tau2*tautwid - 
       pow(E,Tunload*(1/tau1 + 1/tau2))*tau2*
        tautwid + pow(E,
         Tunload/tau1 + Tclean/tautwid)*
        (-1 + pow(E,Tunload/tau2))*tau2*tautwid,2)
     )/(pow(beta,2)*pow(delta,2)*
     pow(A*(-1 + pow(E,Tunload/tau2))*
        sumclean*tau2 - 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*sumunload*
        tauhat,4));

  return sqrt(e1 + e2 + e3 + e4 + e5 + e6);
  
}
