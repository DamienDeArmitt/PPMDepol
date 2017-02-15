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
#include <TROOT.h>

#define HISDIR "./Sum" //Hist. file dir.
#define HISFN "hists_PPMdepol.root" //Histograms file name.
#define EXOFN "obsvbls.root" //Extracted observables file name.


int dpb(){
  double Tclean, Tunload, tautwidle, tauhat, tau1, tau2, sumclean, sumunload;
  double a = pow(3.14*(3/24), 2);
  double A = pow(3.14*(2.75/2), 2);
  double beta = 45;
  double E = exp(1);
  double delta, delta2;
  
  printf("What is the loading phase duration?");
  cin >> delta2;
  
  printf("What is your THold value?");
  cin >> Tclean;
  
  printf("\n\nWhat is your Tunload value?");
  cin >> Tunload;
  
  printf("\n\nWhat is the Tau for the normal holding interval?");
  cin >> tautwidle;
  
  printf("\n\nWhat is the Tau for the poison holding interval?");
  cin >> tauhat;
  
  printf("\n\nWhat is the Tau for the normal unloading interval?");
  cin >> tau1;
  
  printf("\n\nWhat is the Tau for poison unloading interval?");
  cin >> tau2;
  
  printf("\n\nWhat is the total number of hold counts minus background?");
  cin >> sumclean;
  
  printf("\n\nWhat is the total number of unload counts minus background?");
  cin >> sumunload;
  
  delta = 
  (pow(E,Tclean/tauhat - Tclean/tautwidle - 
       Tunload/tau1 + Tunload/tau2)*
     (A*(-1 + pow(E,Tunload/tau1))*sumclean*
        tau1 - a*pow(E,Tunload/tau1)*
        (-1 + pow(E,Tclean/tautwidle))*
        sumunload*tautwidle))/
   (beta*delta2*(-(A*
          (-1 + pow(E,Tunload/tau2))*
          sumclean*tau2) + 
       a*pow(E,Tunload/tau2)*
        (-1 + pow(E,Tclean/tauhat))*
        sumunload*tauhat));
   
   printf(" The depolariztion per bounce for this set is: %f", &delta);
  
  
  
  
  
  return 0;
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