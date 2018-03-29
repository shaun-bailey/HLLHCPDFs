#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "assert.h"
#include <vector>
#include "TLatex.h"
using namespace std;

void write_nplotter_Z_only(double, double);
vector<double> bin_limits(vector<double>,double,double);
void m116_150_13TeV()
{bool GENERATE_MCFM=0;
double p9030_d40x1y6_xval[] = { 1.0, 3.0, 5.0, 7.0, 9.0, 11.5, 14.5, 18.0, 22.5,
  27.5, 33.5, 41.0, 50.0, 60.0, 70.0, 80.0, 95.0, 127.5, 175.0,
  550.0 };
double p9030_d40x1y6_xerrminus[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 2.0, 2.5,
  2.5, 3.5, 4.0, 5.0, 5.0, 5.0, 5.0, 10.0, 22.5, 25.0,
  350.0 };
double p9030_d40x1y6_xerrplus[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 2.0, 2.5,
  2.5, 3.5, 4.0, 5.0, 5.0, 5.0, 5.0, 10.0, 22.5, 25.0,
  350.0 };
double p9030_d40x1y6_yval[] = { 0.1207, 0.2567, 0.2719, 0.2426, 0.2161, 0.1721, 0.1398, 0.1078, 0.08209,
  0.06068, 0.04451, 0.03176, 0.02176, 0.01481, 0.01004, 0.007255, 0.004879, 0.002125, 6.599E-4,
  2.57E-5 };
double p9030_d40x1y6_yerrminus[] = { 0.006721375416534922, 0.005762900707109224, 0.004717288692670822, 0.0042089526915849275, 0.004134236430104112, 0.002985823405695655, 0.0026266035330822202, 0.0019463796751918675, 0.0015335677575183952,
  0.0013155110520250296, 9.649538056300935E-4, 7.414446655010743E-4, 5.971127973842128E-4, 5.002624212350954E-4, 4.3218672446062017E-4, 3.701472275662753E-4, 2.0532418077518294E-4, 8.155682451518082E-5, 3.2543174404473824E-5,
  1.4144342367179892E-6 };
double p9030_d40x1y6_yerrplus[] = { 0.006721375416534922, 0.005762900707109224, 0.004717288692670822, 0.0042089526915849275, 0.004134236430104112, 0.002985823405695655, 0.0026266035330822202, 0.0019463796751918675, 0.0015335677575183952,
  0.0013155110520250296, 9.649538056300935E-4, 7.414446655010743E-4, 5.971127973842128E-4, 5.002624212350954E-4, 4.3218672446062017E-4, 3.701472275662753E-4, 2.0532418077518294E-4, 8.155682451518082E-5, 3.2543174404473824E-5,
  1.4144342367179892E-6 };
double p9030_d40x1y6_ystatminus[] = { 0.002414, 0.0030803999999999996, 0.0032627999999999993, 0.0029112, 0.0030253999999999993, 0.0020652, 0.0019572, 0.0014014, 0.0010671699999999999,
  9.102E-4, 6.676500000000001E-4, 5.0816E-4, 3.6992E-4, 3.2582E-4, 2.8112E-4, 2.53925E-4, 1.31733E-4, 5.3125000000000004E-5, 2.6396000000000003E-5,
  1.2336000000000001E-6 };
double p9030_d40x1y6_ystatplus[] = { 0.002414, 0.0030803999999999996, 0.0032627999999999993, 0.0029112, 0.0030253999999999993, 0.0020652, 0.0019572, 0.0014014, 0.0010671699999999999,
  9.102E-4, 6.676500000000001E-4, 5.0816E-4, 3.6992E-4, 3.2582E-4, 2.8112E-4, 2.53925E-4, 1.31733E-4, 5.3125000000000004E-5, 2.6396000000000003E-5,
  1.2336000000000001E-6 };
int p9030_d40x1y6_numpoints = 20;
TGraphAsymmErrors * p9030_d40x1y6;
p9030_d40x1y6 = new TGraphAsymmErrors(p9030_d40x1y6_numpoints, p9030_d40x1y6_xval, p9030_d40x1y6_yval, p9030_d40x1y6_xerrminus, p9030_d40x1y6_xerrplus, p9030_d40x1y6_yerrminus, p9030_d40x1y6_yerrplus);


p9030_d40x1y6->SetTitle("ATLAS_Z;p^{T}_{ll} [GeV];d\\sigma/dp^{T}_{ll} [pb]");

vector<double> bin_centrals;
int nbins = sizeof(p9030_d40x1y6_xval ) / sizeof(p9030_d40x1y6_xval[0]);
for (int i=0;i<nbins;i++)
{
  bin_centrals.push_back(p9030_d40x1y6_xval[i]);
}
//DEFINING BINNING ===============================
//double bin_limit[] = {10, 13, 16, 20, 25, 30, 37, 45, 55, 65, 75, 85, 105, 150, 200, 900,2000,5000,13000};
vector<double> bin_limit;
 equal log steps
int n=30;
double pt_max=20000.;//GeV
double pt_min=10.;//GeV
double delta_pt = 0.25;
double k=0;
cout<<"[";
while(true)
{
    double xbin;
    xbin=pt_min*pow(10,k);
    k+=delta_pt;
    bin_limit.push_back(xbin);
    if(xbin>=pt_max) break;
    cout<<xbin<<", ";
  }
cout<<"]"<<endl;

/* log steps
int i=11;
while(true)
{
    double xbin;
    if(i%10!=0 || i%10!=9)
    {xbin=(i%10)*pow(10,(i-i%10)/10);i++;}
    if(xbin==0) {continue;i++;}
    if(xbin>=13*pow(10,3)) break;
    bin_limit.push_back(xbin);
    cout<<xbin<<endl;
}
*/
//
if(GENERATE_MCFM)
{
for(int i=0; i<nbins; i++)
{
write_nplotter_Z_only(bin_limit[i],bin_limit[i+1]);

gSystem->Exec("(cd ../ && make -j8)");
gSystem->Exec("./mcfm input.DAT");
gSystem->Exec(Form("mv Z_only_nlo_CT14.NN_80___80___13TeV.C 116_150/%d.C",i));
gSystem->Exec("rm Z_only_nlo_CT14.NN_80___80___13TeV*");
gSystem->Exec("find 116_150/. -type f -name '*.C' -exec sed -i '' s/NaN/0/ {} +");
gSystem->Exec(Form("root -l -q 116_150/%d.C",i));
gSystem->Exec(Form("mv Z_only_nlo_CT14.NN_80___80___13TeV.root 116_150/%d.root",i));

}
}
TCanvas * c = new TCanvas();
p9030_d40x1y6->SetMaximum(10);
p9030_d40x1y6->Draw("AP");
TFile *_file0[nbins];
TH1 * h[nbins];
for(int i=0; i<nbins; i++)
{
_file0[i] = TFile::Open(Form("116_150/%d.root",i));
h[i] = (TH1 *)_file0[i]->Get("id8");
if(h[i]==NULL) continue;
h[i]->Scale(1./1000);
h[i]->SetLineColor(6);
h[i]->Draw("same");

gPad->SetLogx();
gPad->SetLogy();
}
c->Update();
c->Modified();

TLatex * t; t = new TLatex(.2,.26,"116< m_{ll} < 150 GeV");
TLatex * t2; t2 = new TLatex(.2,.2,"p^{T}_{lepton} > 20 GeV");
TLatex * t3; t3 = new TLatex(.2,.14,"|#eta|_{lepton} < 2.4");
t->SetNDC(kTRUE);t2->SetNDC(kTRUE);t3->SetNDC(kTRUE);
t->Draw("SAME");t2->Draw("SAME");t3->Draw("SAME");

c->SaveAs("xsection_13TeV.pdf");
c->SaveAs("xsection_13TeV.root");

//========================== Events plot =====================================
double lumi_HLLHC = 3.*pow(10,18); // 3 ato barn-1
double lumi_ATLAS = 20.3*pow(10,15); // 20.3 femto barn-1

TCanvas * cc = new TCanvas("cc","cc");
p9030_d40x1y6->SetTitle("ATLAS_Z;p^{T}_{ll} [GeV];dN/dp^{T}_{ll} [pb]");
double N_MAX=0;
for(int i=0; i<nbins; i++)
{
  _file0[i] = TFile::Open(Form("116_150/%d.root",i));
  h[i] = (TH1 *)_file0[i]->Get("id8");
  if(h[i]==NULL) continue;
  h[i]->Scale(pow(10,-12)*lumi_HLLHC);
  if(h[i]->GetMaximum()>N_MAX) N_MAX = h[i]->GetMaximum();
}
p9030_d40x1y6->SetMaximum(N_MAX*10);
auto temp1 = p9030_d40x1y6;

for (int i=0;i<temp1->GetN();i++) temp1->GetY()[i] = temp1->GetY()[i]*pow(10,-12)*lumi_ATLAS;

temp1->Draw("AP");
temp1->GetXaxis()->SetLimits(10,10000);
for(int i=0; i<nbins; i++)
{
//  _file0[i] = TFile::Open(Form("116_150/%d.root",i));
//  h[i] = (TH1 *)_file0[i]->Get("id8");
  if(h[i]==NULL) continue;
  //h[i]->Scale(pow(10,-12)*lumi_HLLHC);
  h[i]->SetLineColor(6);
  h[i]->Draw("same");
}
gPad->SetLogx();
gPad->SetLogy();
cc->Update();
cc->Modified();
cc->SaveAs("events_13TeV.pdf");
cc->SaveAs("events_13TeV.root");


}

void write_nplotter_Z_only(double ptmin, double ptmax)
{
gSystem->Exec(Form("sed '140s/ptmin/%lf/' /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only_template.f > /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only_temporary1.f",ptmin));
gSystem->Exec(Form("sed '140s/ptmax/%lf/' /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only_temporary1.f > /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only_temporary2.f",ptmax));
gSystem->Exec(Form("sed '140s/ptbin/%lf/' /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only_temporary2.f > /Users/rabah/Documents/mcfm/src/User/nplotter_Z_only.f",ptmax-ptmin));
}

vector<double> bin_limits(vector<double> bin_centrals,double min, double max)
{
  int nbins = bin_centrals.capacity();
  vector<double> result;
  cout<<"nbins = "<<nbins<<endl;

  for(int i=0; i<nbins;i++)
  {
    if(i==0)
    {result.push_back(min);result.push_back(bin_centrals[i]+(bin_centrals[i+1]-bin_centrals[i])/2);}
    else
    {
      if(i==nbins-1)
      {
        cout<<" i = "<<i<<endl;
        result.push_back(max);
      }
      else
      { cout<<" i = "<<i<<endl;
        result.push_back(bin_centrals[i]+(bin_centrals[i+1]-bin_centrals[i])/2);
      }
    }
  }
return result;

}
