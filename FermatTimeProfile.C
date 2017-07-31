#include "TMath.h"
#include "Riostream.h"
#include "TF1.h"

Double_t index(Double_t *x, Double_t *par)
{
  return (par[0] + par[1] * TMath::Power(TMath::E(),(par[2]*x[0])));
}

Double_t planeDistance(Double_t x1, Double_t x2, Double_t y1, Double_t y2)
{
  return TMath::Sqrt(TMath::Power((x2 - x1),2.0) + TMath::Power((y2 - y1),2.0));
}

Double_t thetaF1(Double_t *x, Double_t *par)
{
  Double_t result = 0.0;
  TF1 *f1 = new TF1("f1", "TMath::Tan(TMath::ASin(TMath::Power(TMath::E(),"
  "-([0]*[1]*TMath::Power(TMath::E(),[2]*x)+"
  "(TMath::Power([1],2)/2)*TMath::Power(TMath::E(),2*[2]*x)))) + [3] + [4])");

  f1->SetParameters(par[0],par[1],par[2],x[0],par[6]);

  result = TMath::Abs(par[5] - (f1->Integral(par[3],par[4])));
  delete f1;
  return result;
}

Double_t path(Double_t *x, Double_t *par)
{
  TF1 *f7 = new TF1("f1", "TMath::Tan(TMath::ASin(TMath::Power(TMath::E(),"
  "-([0]*[1]*TMath::Power(TMath::E(),[2]*x)+"
  "(TMath::Power([1],2)/2)*TMath::Power(TMath::E(),2*[2]*x)))) + [3] + [4])");

  f7->SetParameters(par[0],par[1],par[2],par[4],par[5]);
  Double_t result = f7->Integral(par[3],x[0]);
  delete f7;
  return result;
}

Double_t timeF(Double_t *x, Double_t *par)
{
  return (par[0] + par[1] * TMath::Power(TMath::E(),(par[2]*x[0]))) *
  TMath::Sqrt(TMath::Power(TMath::Tan(TMath::ASin(TMath::Power(TMath::E(),
   -(par[0]*par[1]*TMath::Power(TMath::E(),par[2]*x[0])+
   (TMath::Power(par[1],2)/2)*TMath::Power(TMath::E(),2*par[2]*x[0]))))
   + par[3] + par[4]), 2) + 1);
}

Double_t findTheta(Double_t dx, Double_t z1, Double_t z2, Double_t a, Double_t c, Double_t n)
{
  Double_t realAngle = 0.0;

  TF1 *f3 = new TF1("f3",thetaF1,0,TMath::Pi(),7);
  f3->SetParameters(a,c,n,z1,z2,dx,0);
  f3->SetNpx(10);

  realAngle = f3->GetMinimumX(0,TMath::Pi());

  delete f3;
  return realAngle;
}
/*Double_t integral(Double_t *x, Double_t *par)
{
  TF1 *integral = new TF1("integral","([0] + [1] * TMath::Power(TMath::E(),([2]*x))) *"
  "TMath::Sqrt(TMath::Power(TMath::Tan(TMath::ASin(TMath::Power(TMath::E(),"
   "-([0]*[1]*TMath::Power(TMath::E(),[2]*x)+"
   "(TMath::Power([1],2)/2)*TMath::Power(TMath::E(),2*[2]*x))))"
   "+ [3] + [4]), 2) + 1)");

   integral->SetParameters(par[0],par[1],par[2],par[3],par[4]);

   return integral->Integral(x[0],0);
}*/
Double_t findTime(Double_t angle, Double_t z1, Double_t z2, Double_t a, Double_t c, Double_t n)
  {
    Double_t result = 0.0;
    TF1 *f5 = new TF1("f5", timeF,-1000, 100, 5);
    f5->SetParameters(a,c,n,angle,0);

    result = (1/TMath::C()) * f5->Integral(z1, z2);

    /*TF1 *f6 = new TF1("f6", integral ,-1000, 100, 5);
    f6->SetParameters(a,c,n,angle,0);
    f6->Draw();*/
    delete f5;

    return result;
  }

void FermatTimeProfile(Double_t a, Double_t c, Double_t n, const char *filename, int npoints)
{
  ifstream in(filename);
  ofstream out("Predicted_Time.txt");
  out << "Direct Time:\n";

  if (!in)
  {
    std::cout << "File could not be read\n";
    return;
  }

  Double_t x1, x2, y1, y2, z1, z2;
  Double_t dx = 0.0;
  Double_t theta = 0.0;
  Double_t directTime = 0.0;

  TF1 *indexProfile = new TF1("indexProfile", index , -500, 10,3);
  indexProfile->SetParameters(a,c,n);
  //indexProfile->Draw();


  for(int i = 0; i < npoints; ++i)
  {
    in >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
    z1 = -z1;
    dx = planeDistance(x1, x2, y1, y2);
    theta = findTheta(dx, z1, z2, a, c, n);
    directTime = findTime(theta, z1, z2, a, c, n);
    TF1 *path1 = new TF1("path1",path,-500,100,6);
    path1->SetParameters(a,c,n,z1,theta,0);
    path1->Draw();

    out << directTime << "\n";
  }



}
