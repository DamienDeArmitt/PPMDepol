/* Class for making plots for PPMdepol data analysis... */
/* A. T. Holley */

//Created: July 16, 2016
//Modified: --

/*
Changelog:
*/

class PLOT {
public:
  void get_entry(TTree*,int);
  void loadpnt(double,double,double,double);
  void make(const char*);
  void init();

  TGraphErrors *plot;
  TH1F *haxis;
  TF1 *fit;

  //Timing and other run parameters
  int rnum,rtyp;
  double Tload,Thold,Tunload,Tpoison,Tfield;
  int Tloadstartbin,Tloadendbin,Tholdstartbin,Tholdendbin,Tunldstartbin,Tunldendbin;
  double Tloadstart, Tloadend,Tholdstart,Tholdend,Tunldstart,Tunldend;
  double Tbackgroundstart, Tbackgroundend;
  int Tbackgroundstartbin,Tbackgroundendbin;
  
  //Analysis variables
  double loadcounts,eloadcounts;
  double holdcounts,eholdcounts;
  double unldcounts,eunldcounts;
  
  double holdfitamp, eholdfitamp;
  double holdfittau, eholdfittau;
  double unldfitamp,eunldfitamp,unldfitbg,eunldfitbg;
  double unldfittau, eunldfittau;
  
  
  double holdfitbg, eholdfitbg;
  
  double bgrate, ebgrate;
  double bsnholdcounts, ebsnholdcounts, bsnunldcounts, ebsnunldcounts;
  double tauhat, etauhat, tautwidle, etautwidle; //tautwidle = short, tauhat = long
  double tau1, etau1, tau2, etau2;
  double Tguide;

  
  
  //Plotting range variables
  double xmax;
  double xmin;

private:

  int index;
  double x[1000];
  double ex[1000];
  double y[1000];
  double ey[1000];
  double ymax;
  double eymax;
  
};


void PLOT::init() {
  index = 0;
  xmax = 0;
  ymax = 0;
  eymax = 0;
  xmin = 0;
}


void PLOT::get_entry(TTree *tree, int i) {
  //Timing and other run parameters
  tree->SetBranchAddress("rnum",&rnum);
  tree->SetBranchAddress("rtyp",&rtyp);
  
  tree->SetBranchAddress("Tload",&Tload);
  tree->SetBranchAddress("Thold",&Thold);
  tree->SetBranchAddress("Tunload",&Tunload);
  tree->SetBranchAddress("Tpoison",&Tpoison);
  
  tree->SetBranchAddress("Tfield",&Tfield);
  
  tree->SetBranchAddress("Tloadstart",&Tloadstart);
  tree->SetBranchAddress("Tloadend",&Tloadend);
  tree->SetBranchAddress("Tholdstart",&Tholdstart);
  tree->SetBranchAddress("Tholdend",&Tholdend);
  tree->SetBranchAddress("Tunldstart",&Tunldstart);
  tree->SetBranchAddress("Tunldend",&Tunldend);
  tree->SetBranchAddress("Tbackgroundstart",&Tbackgroundstart);
  tree->SetBranchAddress("Tbackgroundend",&Tbackgroundend);
  
  tree->SetBranchAddress("Tloadstartbin",&Tloadstartbin);
  tree->SetBranchAddress("Tloadendbin",&Tloadendbin);                       
  tree->SetBranchAddress("Tholdstartbin",&Tholdstartbin);
  tree->SetBranchAddress("Tholdendbin",&Tholdendbin);
  tree->SetBranchAddress("Tunldstartbin",&Tunldstartbin);
  tree->SetBranchAddress("Tunldendbin",&Tunldendbin);
  tree->SetBranchAddress("Tbackgroundstartbin",&Tbackgroundstartbin);
  tree->SetBranchAddress("Tbackgroundendbin",&Tbackgroundendbin);
  
  //Analysis variables
  tree->SetBranchAddress("loadcounts",&loadcounts);
  tree->SetBranchAddress("eloadcounts",&eloadcounts);
  tree->SetBranchAddress("holdcounts",&holdcounts);
  tree->SetBranchAddress("eholdcounts",&eholdcounts);
  tree->SetBranchAddress("unldcounts",&unldcounts);
  tree->SetBranchAddress("eunldcounts",&eunldcounts);
  
  tree->SetBranchAddress("holdfitamp",&holdfitamp);
  tree->SetBranchAddress("eholdfitamp",&eholdfitamp);
  tree->SetBranchAddress("holdfittau",&holdfitamp);
  tree->SetBranchAddress("eholdfittau",&eholdfitamp);
  tree->SetBranchAddress("holdfitbg",&holdfitbg);
  tree->SetBranchAddress("eholdfitbg",&eholdfitbg);
  
  tree->SetBranchAddress("unldfitamp",&unldfitamp);
  tree->SetBranchAddress("eunldfitamp",&eunldfitamp);
  tree->SetBranchAddress("unldfitbg",&unldfitbg);
  tree->SetBranchAddress("eunldfitbg",&eunldfitbg);
  tree->SetBranchAddress("tau1",&tau1);
  tree->SetBranchAddress("etau1",&etau1);
  tree->SetBranchAddress("tau2",&tau2);
  tree->SetBranchAddress("etau2",&etau2);
  
  tree->SetBranchAddress("bgrate", &bgrate);
  tree->SetBranchAddress("ebgrate", &ebgrate);
  tree->SetBranchAddress("bsnholdcounts", &bsnholdcounts);
  tree->SetBranchAddress("ebsnholdcounts", &ebsnholdcounts);
  tree->SetBranchAddress("bsnunldcounts", &bsnunldcounts);
  tree->SetBranchAddress("ebsnunldcounts", &ebsnunldcounts);
  
  tree->SetBranchAddress("Tguide", &Tguide);
  
  
  tree->GetEntry(i);
  
  return;
}


void PLOT::loadpnt(double lx, double lex, double ly, double ley) {
  x[index] = lx;
  if(x[index] > xmax) xmax = x[index];
  y[index] = ly;
  if(y[index] > ymax) {
    ymax = y[index];
    eymax = ley;
  }
  ex[index] = lex;
  ey[index] = ley;
  index++;
  return;
}


void PLOT::make(const char *plotname) {
  char hisname[100];
  sprintf(hisname,"%s_hisname",plotname);
  if(xmin != 0)
    haxis = new TH1F(hisname,"",(int)ceil(xmax)*100,xmin-0.1*xmin,xmax+0.1*xmax);
  else
    haxis = new TH1F(hisname,"",(int)ceil(xmax)*100,-10,xmax+0.1*xmax);
  plot = new TGraphErrors(index,x,y,ex,ey);
  haxis->GetYaxis()->SetRangeUser(0,ymax + 2*eymax);
  haxis->SetStats(kFALSE);
  haxis->GetYaxis()->CenterTitle();
  haxis->GetYaxis()->SetTitleOffset(1.2);
  haxis->SetLineColor(0);
  plot->SetName(plotname);
  return;
}
