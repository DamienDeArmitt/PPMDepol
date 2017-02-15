/* Code to extract data from PPM depol runs from LANL 2009 */
/* A. T. Holley and D. DeArmitt */
/* Created: June 1, 2016 */
/* Modified: June 15, 2016 */
/* Changelog:
June 4, 2016 [ath] : The ACLiC compiler, it turns out, didn't like the 'gROOT' because I had hadn't included 'TROOT.h', which is now included. 'name' character array was not used in 'gethists', so I removed it. In 'data2hists', we weren't using the sumok to check that the sums had been correctly completed, so I added a check for that and non-empty sum histograms before moving them to the directory for saving. Extraneous variables 'rfname', 'hfname', and 'dfp' removed from 'data2hists'. Extraneous variable 'name' removed from 'sumcycles'. Now the code compiles via ACLiC without any warnings. However, it fails to run in interpreted mode because in 'findcyclestart' when I fixed the indexing bug I neglected to ensure that the new index didn't get set to -1 for the first cycle. I also noticed that 'TSTARTCYCLE' wasn't being initialized, so I fixed that too. I also found that some drawing code had somehow made its way into the parts of 'data2hists' where the histograms are saved. I removed it.

June 4, 2016 [ath] : Added a feature where the non-cut extracted cycle start times are written to a text file called rls<run number>.txt in a 'cycles' directory. If the code finds that such a file exists, it uses those values rather than running the finding algorithm again. You can just delete the files if you want the finding algorithm to run. You can also manually remove specific start time lines if they are detected incorrectly, which should allow for better run reconstruction.

June 15, 2016 [ath] : Added a function to allow runs to be examined without exiting from root and cycle start finding parameters to be changed on the fly. Also added a plot of the raw histogram with the cut-out data actually cut out.

June 22, 2016 [dad] : Added the path for the cycles directory since i moved files around.
*/

#include <stdint.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TROOT.h>

#define DDIR "/home/damien/PPMdepol/DATA" //Path to the data files.
#define HISDIR "./" //Path to which to write the histograms file.
#define HISFN "hists_PPMdepol.root" //File that contains the extracted histograms and cycle sum histograms for each run.
#define MAXCYCLES 1000 //Maximum number of cycles expected.
#define MAXCUTS 20 //Maximum number if cut values.
#define CYPATH ".."

//Function declarations.
void gethists(const char *); /* Enter: "####,####" for a range of runs, or "file" to extract runs from the 'PPMdepolruns.txt' file, or "force" to  */
int data2hists(int, const char *); //Read text files for each detector in a run and convert them to histograms, then save the histograms in a file.
int sumcycles(int); //Sum the raw data histograms into a single cycle for each detector and save the histograms into the histograms file.
int getdata(int,const char*,TH1F*); //Fill a histogram from the specificed text file.
void cleanup(); //Clean up histograms.
int FindCycleStart(int rnum); //Find the start times of all the cycles.
void vlinedraw(double,double,int);
void runcheck(); //A qol function for checking runs.

//Define global variables
int TSTARTBINS[MAXCYCLES];
TH1F *gvmon,*main,*a2,*cnmon; //Raw histograms.
TH1F *gvmonwc,*mainwc,*a2wc,*cnmonwc; //Raw histograms with cuts removed.
TH1F *sgvmon,*smain,*sa2,*scnmon; //Summed histograms.
TH1F *usgvmon,*usmain,*usa2,*uscnmon; //Uncut summed histograms.
int win = 10; //Number of bins to integrate over for deciding if a cycle has begun.
int thresh = 5; //Threshold for total counts in integration window for deciding if a cycle has begun.
int bindelay = 550; //Number of bins that must pass before a new cycle start will be sought.
double tres = 0.1; //Desired histogram bin resolution.

void runcheck() {
  int rnum;
  char cmd[50];
  char response[100];

	cleanup();

  printf("\n  Enter run number: ");
  scanf("%d",&rnum);
 
  printf("\n  Current bin resolution: %f s per bin\n",tres);
  printf("  Current cycle detection parameters:\n");
  printf("         win: %d\n",win);
  printf("      thresh: %d\n",thresh);
  printf("    bindelay: %d\n\n",bindelay);
  printf("Change rl parameters and delete the rls file for run# %d (y/n)? ",rnum);
    scanf("%s",response);
    if(strcmp(response,"y")==0) {
      printf("Enter a new set of parmeters: tres,win,thresh,bindelay: ");
      scanf("%lf,%d,%d,%d",&tres,&win,&thresh,&bindelay);
      printf("\n  New bin resolution: %f s per bin\n",tres);
      printf("  New cycle detection parameters:\n");
      printf("         win: %d\n",win);
      printf("      thresh: %d\n",thresh);
      printf("    bindelay: %d\n\n",bindelay);
      sprintf(response,".!rm %s/cycles/rls%d.txt",CYPATH,rnum);
      gROOT->ProcessLine(response);
  }
  
  sprintf(cmd,"%d,%d",rnum,rnum);
  gethists(cmd);
  

  return;
}

void gethists(const char *command) {
  int i=0;
  int runlow,runhigh,prun;
  int doflag = 0;
  FILE *fp;
  
  //Set range of run numbers to search
  int scanrangelow = 1;
  int scanrangehigh = 397;
  
  if(atof(command) != 0) {
    sscanf(command,"%d,%d",&runlow,&runhigh);
    doflag = 1;
  }
  if(strcmp(command,"force") == 0) doflag = 2;
  if(strcmp(command,"") == 0) doflag = 3;
  if(strcmp(command,"file") == 0) doflag = 4;
  
  if(doflag == 0) {
    printf("\n\n  Conversion command not recognized.\n\n");
    return;
  }
  
  if(doflag == 4) { //We want to read runs from the PPMdepolruns.txt file.
      fp = fopen("./PPMdepolruns.txt","r");
      if(fp == NULL) {
	printf("\n  'PPMdepolruns.txt' file not found!\n\n");
	return;
      }
      while(feof(fp) == 0) {
        fscanf(fp,"%d",&prun);
        printf("\n ------ Getting histograms from run# %d ------\n",prun);
        data2hists(prun,"F"); //This function should get histograms for a particular run.
      }
      fclose(fp);
      return;
    }
  
  for(i=scanrangelow ; i <= scanrangehigh ; i++) {
    if(doflag == 2) {
      printf("\n ------ Getting histograms from run# %d ------\n",i);
      data2hists(i,"F");
    }
    if(doflag == 3) {
      printf("\n ------ Getting histograms from from run# %d ------\n",i);
      data2hists(i,"");
    }
    if(doflag == 1 && i >= runlow && i <= runhigh) {
      printf("\n ------ Getting histograms from run# %d ------\n",i);
      if(runlow == runhigh) {
	data2hists(i,"D");
      }
      else {
	data2hists(i,"");
      }
    }
  }

  return;
}


int data2hists(int rnum, const char *drawop) {  //Function to get histograms from a run and save them into a root file.
  int i;
  double Tcyclestart;
  char hisname[100];
  char name[100];

  //Don't draw canvases etc. if we aren't drawing plots.
  if(strcmp(drawop,"D") != 0 && strcmp(drawop,"D+") != 0) gROOT->SetBatch(kTRUE);

           /* Setup histogram containers */

  //Binning
  double tlow = 0.0;
  double thigh = 6000.0;
  int bins = (int)ceil((thigh-tlow)/tres);
  double tlowsums = 0.0;
  double thighsums = 600.0;
  int sumbins = (int)ceil((thighsums-tlowsums)/tres);
  if((tlow-thigh)/bins != (tlowsums-thighsums)/sumbins) {
    printf("  !!!! Raw <--> Summed Binning ERROR !!!!\n");
    return -1;
  }
 
  
  //Setup a2 histograms.
  sprintf(hisname,"a%d",rnum);
  a2 = new TH1F(hisname,hisname,bins,tlow,thigh);
  a2->GetXaxis()->SetTitle("Time [s]");
  a2->GetYaxis()->SetTitle("Rate [s^{-1}]");
  a2->GetYaxis()->CenterTitle();
  a2->SetStats(kFALSE);
  a2->SetLineColor(15);
  
  sprintf(hisname,"awc%d",rnum);
  a2wc = new TH1F(hisname,hisname,bins,tlow,thigh);
  a2wc->GetXaxis()->SetTitle("Time [s]");
  a2wc->GetYaxis()->SetTitle("Rate [s^{-1}]");
  a2wc->GetYaxis()->CenterTitle();
  a2wc->SetStats(kFALSE);
  a2wc->SetLineColor(6);
  
  sprintf(hisname,"sa%d",rnum);
  sa2 = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  sa2->GetXaxis()->SetTitle("Time [s]");
  sa2->GetYaxis()->SetTitle("Rate [s^{-1}]");
  sa2->GetYaxis()->CenterTitle();
  sa2->SetStats(kFALSE);
  sa2->SetLineColor(6);
  
  sprintf(hisname,"usa%d",rnum);
  usa2 = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  usa2->GetXaxis()->SetTitle("Time [s]");
  usa2->GetYaxis()->SetTitle("Rate [s^{-1}]");
  usa2->GetYaxis()->CenterTitle();
  usa2->SetStats(kFALSE);
  usa2->SetLineColor(15);
  
  //Setup cnmon histograms.
  sprintf(hisname,"cnmon%d",rnum);
  cnmon = new TH1F(hisname,hisname,bins,tlow,thigh);
  cnmon->GetXaxis()->SetTitle("Time [s]");
  cnmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  cnmon->GetYaxis()->CenterTitle();
  cnmon->SetStats(kFALSE);
  cnmon->SetLineColor(15);
  
  sprintf(hisname,"cnmonwc%d",rnum);
  cnmonwc = new TH1F(hisname,hisname,bins,tlow,thigh);
  cnmonwc->GetXaxis()->SetTitle("Time [s]");
  cnmonwc->GetYaxis()->SetTitle("Rate [s^{-1}]");
  cnmonwc->GetYaxis()->CenterTitle();
  cnmonwc->SetStats(kFALSE);
  cnmonwc->SetLineColor(2);
  
  sprintf(hisname,"scnmon%d",rnum);
  scnmon = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  scnmon->GetXaxis()->SetTitle("Time [s]");
  scnmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  scnmon->GetYaxis()->CenterTitle();
  scnmon->SetStats(kFALSE);
  scnmon->SetLineColor(2);
  
  sprintf(hisname,"uscnmon%d",rnum);
  uscnmon = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  uscnmon->GetXaxis()->SetTitle("Time [s]");
  uscnmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  uscnmon->GetYaxis()->CenterTitle();
  uscnmon->SetStats(kFALSE);
  uscnmon->SetLineColor(15);
  
  //Setup gvmon histograms.
  sprintf(hisname,"gvmon%d",rnum);
  gvmon = new TH1F(hisname,hisname,bins,tlow,thigh);
  gvmon->GetXaxis()->SetTitle("Time [s]");
  gvmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  gvmon->GetYaxis()->CenterTitle();
  gvmon->SetStats(kFALSE);
  gvmon->SetLineColor(15);
  
  sprintf(hisname,"gvmonwc%d",rnum);
  gvmonwc = new TH1F(hisname,hisname,bins,tlow,thigh);
  gvmonwc->GetXaxis()->SetTitle("Time [s]");
  gvmonwc->GetYaxis()->SetTitle("Rate [s^{-1}]");
  gvmonwc->GetYaxis()->CenterTitle();
  gvmonwc->SetStats(kFALSE);
  gvmonwc->SetLineColor(3);
  
  sprintf(hisname,"sgvmon%d",rnum);
  sgvmon = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  sgvmon->GetXaxis()->SetTitle("Time [s]");
  sgvmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  sgvmon->GetYaxis()->CenterTitle();
  sgvmon->SetStats(kFALSE);
  sgvmon->SetLineColor(3);
  
  sprintf(hisname,"usgvmon%d",rnum);
  usgvmon = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  usgvmon->GetXaxis()->SetTitle("Time [s]");
  usgvmon->GetYaxis()->SetTitle("Rate [s^{-1}]");
  usgvmon->GetYaxis()->CenterTitle();
  usgvmon->SetStats(kFALSE);
  usgvmon->SetLineColor(15);
  
  //Setup main histograms.
  sprintf(hisname,"main%d",rnum);
  main = new TH1F(hisname,hisname,bins,tlow,thigh);
  main->GetXaxis()->SetTitle("Time [s]");
  main->GetYaxis()->SetTitle("Rate [s^{-1}]");
  main->GetYaxis()->CenterTitle();
  main->SetStats(kFALSE);
  main->SetLineColor(15);
  
  sprintf(hisname,"mainwc%d",rnum);
  mainwc = new TH1F(hisname,hisname,bins,tlow,thigh);
  mainwc->GetXaxis()->SetTitle("Time [s]");
  mainwc->GetYaxis()->SetTitle("Rate [s^{-1}]");
  mainwc->GetYaxis()->CenterTitle();
  mainwc->SetStats(kFALSE);
  mainwc->SetLineColor(4);
  
  sprintf(hisname,"smain%d",rnum);
  smain = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  smain->GetXaxis()->SetTitle("Time [s]");
  smain->GetYaxis()->SetTitle("Rate [s^{-1}]");
  smain->GetYaxis()->CenterTitle();
  smain->SetStats(kFALSE);
  smain->SetLineColor(4);

  sprintf(hisname,"usmain%d",rnum);
  usmain = new TH1F(hisname,hisname,sumbins,tlowsums,thighsums);
  usmain->GetXaxis()->SetTitle("Time [s]");
  usmain->GetYaxis()->SetTitle("Rate [s^{-1}]");
  usmain->GetYaxis()->CenterTitle();
  usmain->SetStats(kFALSE);
  usmain->SetLineColor(15);
  
       /* --------------------------------------------------------------------- */
     
  //Get the histograms.
  int a2ok = getdata(rnum,"A2run",a2);
  int cnok = getdata(rnum,"CNrun",cnmon);
  int gvok = getdata(rnum,"GVmon-run",gvmon);
  int mainok = getdata(rnum,"main-run",main);
  
  int sumok = sumcycles(rnum); //Sum the cycles in the histograms.

  //Save and clean up or draw.
  if(strcmp(drawop,"D") != 0) {  //If drawing is not requested, save the histograms.
    sprintf(name,"%s/%s",HISDIR,HISFN);
    TDirectory *wrkdir = gDirectory;
    TFile *rfp = new TFile(name,"UPDATE");
    
    //Check to see if histograms for the requested run already exist
    int hisexistflag = 1;
    sprintf(name,"main%d",rnum);
    TH1F *hisexistck = (TH1F*)gDirectory->Get(name);
    if(!hisexistck) hisexistflag = 0;
    delete hisexistck;  hisexistck = NULL;
    if(hisexistflag == 1) {
      if(strcmp(drawop,"F") != 0) {
        printf("\n\n Histograms for run# %d already exist and",rnum);
        printf(" will not be overwritten.\n");
        rfp->Close();  delete rfp;  rfp = NULL;
        return 0;
      }
    }

    wrkdir->cd();
    if(a2ok == 0 && a2->GetEntries() > 5) a2->SetDirectory(rfp);
    if(cnok == 0 && cnmon->GetEntries() > 5) cnmon->SetDirectory(rfp);
    if(gvok == 0 && gvmon->GetEntries() > 5) gvmon->SetDirectory(rfp);
    if(mainok == 0 && main->GetEntries() > 5) main->SetDirectory(rfp);
    if(sumok == 0) {
      if(sa2->GetEntries() > 5) sa2->SetDirectory(rfp);
      if(scnmon->GetEntries() > 5) scnmon->SetDirectory(rfp);
      if(sgvmon->GetEntries() > 5) sgvmon->SetDirectory(rfp);
      if(smain->GetEntries() > 5) smain->SetDirectory(rfp);
    }
    rfp->cd();
    printf("\n\n  Done... saving histograms\n\n");
    rfp->Write("",TObject::kOverwrite);
    cleanup();
    rfp->Close();
  }
  else {  //Otherwise draw the histograms
    printf("\n\n  Done... displaying select histograms\n");
    if(gvok == 0) {
      TCanvas *gvcan = new TCanvas("gvcan");
      gvcan->Divide(1,2);
      gvcan->cd(1);
      gvmon->Draw();
      gvmonwc->Draw("Same");
      gvcan->cd(2);
      sgvmon->Draw();
      usgvmon->Draw("Same");
      gvcan->cd(0);
    }
    else printf("\n  GV histogram not plotted-- data file not found!\n");
    if(cnok == 0) {
      TCanvas *cncan = new TCanvas("cncan");
      cncan->Divide(1,2);
      cncan->cd(1);
      cnmon->Draw();
      cnmonwc->Draw("Same");
      cncan->cd(2);
      scnmon->Draw();
      uscnmon->Draw("Same");
      cncan->cd(0);
    }
    else printf("\n  CN histogram not plotted-- data file not found!\n");
    if(a2ok == 0) {
      TCanvas *a2can = new TCanvas("a2can");
      a2can->Divide(1,2);
      a2can->cd(1);
      a2->Draw();
      a2wc->Draw("Same");
      a2can->cd(2);
      sa2->Draw();
      usa2->Draw("Same");
      a2can->cd(0);
     }
     else printf("\n  a2 histogram not plotted-- data file not found!\n");
     if(mainok == 0) {
       TCanvas *maincan = new TCanvas("maincan");
       maincan->Divide(1,2);
       maincan->cd(1);
       main->Draw(); 
       for(i=0 ; i < MAXCYCLES ; i++) {
         if(TSTARTBINS[i] == -1) break;
	     Tcyclestart = main->GetXaxis()->GetBinLowEdge(TSTARTBINS[i]);
	     vlinedraw(Tcyclestart,Tcyclestart,2);
       }
       
       mainwc->Draw("Same");
       maincan->cd(2);
       //smain->Draw();
       usmain->Draw();
       smain->Draw("Same");
       maincan->cd(0);
     }
     else printf("\n  main histogram not plotted-- data file not found!\n"); 
  }
  
  return 0;
}


int getdata(int rnum, const char *runame, TH1F *his) {  //Here we fill the histograms from the correct raw data file.
  int i;
  int val1,val2,val3;
  char rname[100];
  char pars[500];
  FILE *dfp;
  
  if(rnum < 10) sprintf(rname,"%s/%s00%d.csv",DDIR,runame,rnum);
  if(rnum > 10 && rnum < 100) sprintf(rname,"%s/%s0%d.csv",DDIR,runame,rnum);
  if(rnum >= 100) sprintf(rname,"%s/%s%d.csv",DDIR,runame,rnum);
  
  dfp = fopen(rname,"r");
  if(dfp == NULL) {
   printf("\n  Run file '%s' was not found!\n",rname);
   return -1;
  }
  
  for(i=0 ; i < 13 ; i++) fgets(pars,500,dfp); //Skip over the thirteen-line header.
  while(feof(dfp) == 0) {  //Read in the data
    fscanf(dfp,"%d %d %d",&val1,&val2,&val3);
    his->SetBinContent(val1,val3);
  }
  
  fclose(dfp);
  return 0;
}


int sumcycles(int rnum) {
  int i,j,k,rlbin_low,rlbin_high,rawcnts,sumcnts,sumbin;
  int runcheck,val;
  int cuts[MAXCUTS];
  int numofbindeltas = 0;
  double avgbindelta = 0;
  char pars[500];
  char vals[10];
  FILE *cutsfp;
  
  for(i=0 ; i < MAXCYCLES ; i++) TSTARTBINS[i] = -1; //Initialize array.
  
   //Read in cuts for the current run and put them into the 'cuts' array.
  for(i=0 ; i < MAXCUTS ; i++) cuts[i] = -1;
  cutsfp = fopen("RawDataCuts.txt","r");
  if(cutsfp == NULL) {
    printf("  'RawDataCuts' file not found!\n");
    return -1;
  }
  for(i=0 ; i < 2 ; i++) fgets(pars,500,cutsfp);
  while(feof(cutsfp) == 0) {
    fscanf(cutsfp,"%d",&runcheck);
    if(runcheck == rnum) {
      for(j=0 ; j < MAXCUTS ; j++) {
        fscanf(cutsfp,"%s",vals);
        val = (int)atof(vals);
        if(val == 0) break;
        else cuts[j] = val;
      }
    }
    else fgets(pars,500,cutsfp);
  }
  fclose(cutsfp);
  
  printf("  Finding cycle start times...\n");
  FindCycleStart(rnum);
  
  printf("  Summing the cycles...\n");
  int nbins = main->GetXaxis()->GetNbins(); //All the histograms should have the same number of bins.
  for(i=0 ; i < MAXCYCLES ; i++) {
   if(TSTARTBINS[i] == -1) break; 
  }
  int cycles = i;
  
  int cutflag = 0;
  int cutflag4avg = 0;
  for(i = 0 ; i < cycles ; i++) {  //i here indexes the red lines.
    rlbin_low = TSTARTBINS[i];
    rlbin_high = TSTARTBINS[i+1];
    for(k = 0 ; k < MAXCUTS ; k = k + 2) {
	  if(cuts[k] == -1) break;
	  if(rlbin_low >= cuts[k] && rlbin_high <= cuts[k+1]) cutflag4avg = 1;
    }
    if(cutflag4avg == 0) {
      numofbindeltas++;
      avgbindelta = avgbindelta + rlbin_high - rlbin_low;
    }
    else cutflag4avg = 0;
    if(rlbin_high == -1) rlbin_high = nbins;
    for(j = rlbin_low ; j < rlbin_high ; j++) {
      //Don't don any summing if j is between a pair of cut bins.
      cutflag = 0;
      for(k = 0 ; k < MAXCUTS ; k = k + 2) {
	    if(cuts[k] == -1) break;
	    if(j >= cuts[k] && j <= cuts[k+1]) cutflag = 1;
      }
      
      sumbin = j - rlbin_low;
      
      rawcnts = main->GetBinContent(j);
      sumcnts = usmain->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      usmain->SetBinContent(sumbin,sumcnts);
      
      rawcnts = gvmon->GetBinContent(j);
      sumcnts = usgvmon->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      usgvmon->SetBinContent(sumbin,sumcnts);
      
      rawcnts = cnmon->GetBinContent(j);
      sumcnts = uscnmon->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      uscnmon->SetBinContent(sumbin,sumcnts);
      
      rawcnts = a2->GetBinContent(j);
      sumcnts = usa2->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      usa2->SetBinContent(sumbin,sumcnts);
      
      if(cutflag == 1) continue; //If we are in a cut region, don't do any summing.

      rawcnts = main->GetBinContent(j);
      sumcnts = smain->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      smain->SetBinContent(sumbin,sumcnts);
      mainwc->SetBinContent(j,rawcnts);
      
      rawcnts = gvmon->GetBinContent(j);
      sumcnts = sgvmon->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      sgvmon->SetBinContent(sumbin,sumcnts);
      gvmonwc->SetBinContent(j,rawcnts);
      
      rawcnts = cnmon->GetBinContent(j);
      sumcnts = scnmon->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      scnmon->SetBinContent(sumbin,sumcnts);
      cnmonwc->SetBinContent(j,rawcnts);
      
      rawcnts = a2->GetBinContent(j);
      sumcnts = sa2->GetBinContent(sumbin);
      sumcnts = sumcnts + rawcnts;
      sa2->SetBinContent(sumbin,sumcnts);
      a2wc->SetBinContent(j,rawcnts);
    }
  }
  
  printf("      %d different bin deltas, with an average delta of %f\n",numofbindeltas,avgbindelta/numofbindeltas);
  
  return 0;
}


int FindCycleStart(int rnum) {
  int i;
  int skipflag = 0;
  int cycletimecnt = 0;
  char pars[50];
  char name[100];
  FILE *cfp; //Pointer to the cycles file.
  
  for(i=0 ; i < MAXCYCLES ; i++) TSTARTBINS[i] = -1; //Initialize the start bins array.
  
  //Check to see if a cycles file exists.
  sprintf(name,"%s/cycles/rls%d.txt",CYPATH,rnum);
  cfp = fopen(name,"r");
  if(cfp != NULL) {  //The file exists, so read the start times into TSTARTBINS and return.
	i=0;
    while(feof(cfp) == 0) {
      if(i == MAXCYCLES) {
        printf("  'rls%d' has too many cycle times... ignoring any past MAXCYCLES = %d!\n",rnum,MAXCYCLES);
        fclose(cfp);
        return -1;
      }
      fscanf(cfp,"%s",pars);
      if(strcmp(pars,"end") == 0) {
        printf("  %d cycles times read from 'rls%d'\n",i,rnum);
        cycletimecnt = i;
        skipflag = 1;
        fclose(cfp);
        break;
      }
      else {
        TSTARTBINS[i] = (int)atof(pars);
        i++;
      }
    }
    if(strcmp(pars,"end") != 0) {
      printf("  ERROR: No 'end' designation in 'rls%d' file!\n",rnum);
      fclose(cfp);
      return -1;
    }
  }
  //for(i=0 ; i < MAXCYCLES ; i++) printf("%d\n",TSTARTBINS[i]);
  if(skipflag == 0) {  //If there was no red lines file, find the cycle start times and save them to a red lines file.
    TH1F *his = (TH1F*)main->Clone(); //Make a copy of the histogram that we want to use for determining the cycle start times.

    //Scan through the histogram and find cycle start times.
    int bins = his->GetXaxis()->GetNbins();
    int lastcyclebin = 0;
    double cnts = 0;
    cycletimecnt = 0; //Keeps track of the mail slot number in Tstartbin.
    for(i=1 ; i <= bins ; i++) {  //Scanning all bins in the histogram.
      cnts = his->Integral(i,i+win); //Getting number of counts for some period after the current bin.
      if(cycletimecnt == 0) lastcyclebin = 0; //Make sure we don't try to access the '-1' mail slot.
      else lastcyclebin = TSTARTBINS[cycletimecnt-1]; //Get bin number of previous cycle start time.
      if(cnts > thresh && i > lastcyclebin + bindelay) {
        TSTARTBINS[cycletimecnt] = i;
        cycletimecnt++;
	    if(cycletimecnt == MAXCYCLES) {
	      printf("\n  Max cycles reached!\n");
	      return -1;
	    }
      }
    }
    //Save the cycle start times to a red line file.
   cfp = fopen(name,"w");
   for(i=0 ; i < MAXCYCLES ; i++) {
     if(TSTARTBINS[i] == -1) break;
     fprintf(cfp,"%d\n",TSTARTBINS[i]);
   }
   fprintf(cfp,"end");
   fclose(cfp);
  }
    
  printf("  Detected cycles: %d\n",cycletimecnt);
  printf("  Bin Delta(s):\n");
  int tempbindelta;
  int maxbindelta = 0;
  int bindeltacnt = 1;
  int bindelta = TSTARTBINS[1]-TSTARTBINS[0];
  double avgbindelta = 0;
  //printf("    #0: %d\n",bindelta);
  for(i=1 ; i < MAXCYCLES ; i++) {
    if(TSTARTBINS[i] == -1) break;
    tempbindelta = TSTARTBINS[i]-TSTARTBINS[i-1];
    if(tempbindelta > maxbindelta) maxbindelta = tempbindelta;
    if(tempbindelta != bindelta) {
      //printf("    #%d: %d\n",bindeltacnt,tempbindelta);
      avgbindelta = avgbindelta + tempbindelta;
	  bindeltacnt++;
    }
  }
  printf("      %d different bin deltas, with average %f and max delta: %d\n",bindeltacnt,(double)avgbindelta/(double)bindeltacnt,maxbindelta);

  return 0;
}


void cleanup() {
 
  delete a2;       a2 = NULL;
  delete gvmon;    gvmon = NULL;
  delete cnmon;    cnmon = NULL;
  delete main;     main = NULL;
  
  delete a2wc;       a2wc = NULL;
  delete gvmonwc;    gvmonwc = NULL;
  delete cnmonwc;    cnmonwc = NULL;
  delete mainwc;     mainwc = NULL;
  
  TCanvas *maincan = (TCanvas*)gROOT->FindObject("maincan");
  if(maincan != NULL) {delete maincan; maincan = NULL;}
  TCanvas *gvcan = (TCanvas*)gROOT->FindObject("gvcan");
  if(gvcan != NULL) {delete gvcan; gvcan = NULL;}
  TCanvas *a2can = (TCanvas*)gROOT->FindObject("a2can");
  if(a2can != NULL) {delete a2can; a2can = NULL;}
  TCanvas *cncan = (TCanvas*)gROOT->FindObject("cncan");
  if(cncan != NULL) {delete cncan; cncan = NULL;}
  
  
  delete sa2;     sa2 = NULL;
  delete sgvmon;  sgvmon = NULL;
  delete scnmon;  scnmon = NULL;
  delete smain;   smain = NULL;
  
  delete usa2;     usa2 = NULL;
  delete usgvmon;  usgvmon = NULL;
  delete uscnmon;  uscnmon = NULL;
  delete usmain;   usmain = NULL;
  
  gDirectory->GetList()->Delete(); //This gets rid of any named histograms still hanging around because ROOT keeps track of names.
  
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