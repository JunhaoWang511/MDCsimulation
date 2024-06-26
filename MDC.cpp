
#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>
#include <string.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH2.h>
#include <TH1.h>
#include <TPolyLine.h>
#include <TStyle.h>
#include <TF1.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/Track.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;
using namespace std;
using namespace TMath;

#define cout cout << "[UserInfo]: "

void start_run(string filename = "Gasfile/He_C3H8_60.40_1T.gas", int fileid = 0)
{
  // ------- set plotting option -------

  bool drawOption = false;
  bool saveOption = false;
  if (!drawOption)
    saveOption = false;

  // ------- creat root file -------

  size_t lastslash = filename.find_last_of('/') + 1;
  size_t lastdot = filename.find_last_of('.');
  if (lastslash > 1e10)
    lastslash = 0;
  string gasname = filename.substr(lastslash, lastdot - lastslash);
  TFile *outfile = new TFile(Form("output/MDC_%s_%i.root", gasname.c_str(), fileid), "RECREATE");
  TTree *Data = new TTree("Time", "Time");
  vector<Double_t> reachtime;
  double firstTime;
  double momentum;
  double distance;
  double angle;
  double xstart;
  Data->Branch("reachTime", &reachtime);
  Data->Branch("firstTime", &firstTime);
  Data->Branch("momentum", &momentum);
  Data->Branch("distance", &distance);
  Data->Branch("angle", &angle);
  Data->Branch("xstart", &xstart);

  // ------- set gas transport properties -------

  MediumMagboltz *gas = new MediumMagboltz();
  gas->LoadGasFile(filename);
  double vx, vy, vz;
  gas->ElectronVelocity(2000, 0, 0, 0, 0, 1, vx, vy, vz);
  double Velocity = sqrt(pow(vx, 2) + pow(vy, 2));
  double Ldiffusion, Tdiffusion;
  gas->ElectronDiffusion(760, 0, 0, 0, 0, 1, Ldiffusion, Tdiffusion);
  cout << "drift velocity at 2000 v/cm: " << Velocity * 1000 << " cm/us" << endl;
  cout << "transverse diffusion coefficient at 760 v/cm: " << Tdiffusion * 10000 << " um/sqrt(cm)" << endl;
  // Set the Penning transfer efficiency.
  // constexpr double rPenning = 0.51;
  // constexpr double lambdaPenning = 0.;
  // gas->EnablePenningTransfer(rPenning, lambdaPenning, "He");
  const std::string path = std::getenv("GARFIELD_INSTALL");
  // ions cannot be tracked microscopically in Garfield++, so we set ion transport properties manually.
  gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_He+_He.txt");
  gas->PrintGas();
  ViewMedium mediumView;
  TCanvas *medium_can = nullptr;
  mediumView.SetMedium(gas);
  if (drawOption)
  {
    medium_can = new TCanvas("medium", "", 960, 680);
    mediumView.SetCanvas(medium_can);
    mediumView.PlotElectronVelocity('e');
  }
  if (saveOption)
    medium_can->SaveAs(Form("drift_velocity-%s.png", gasname.c_str()));

  // ------- define geometry and component field -------

  const double rSense = 0.00105; // Wire radius [cm]
  const double lbox = 100.;      // Half-length of the box[cm]
  const double halfcell = 0.6;
  const double hh = halfcell;
  const double width = 2 * halfcell;
  const double ex = 0.;

  // 3:1
  const double p = 0.;
  const double sh = 1.;
  const double di = 1.;
  const double rField = 0.00505;

  GeometrySimple *geo = new GeometrySimple();
  SolidBox *box = new SolidBox(0., 0., 0., width / 2 * di, width / 2, lbox);
  geo->AddSolid(box, gas);

  ComponentAnalyticField *cmp = new ComponentAnalyticField();
  cmp->SetGeometry(geo);
  double bx, by, bz;
  int st;
  cmp->SetMagneticField(0, 0, 1.0);

  cmp->EnablePeriodicityX();
  cmp->EnablePeriodicityY();
  cmp->SetPeriodicityX(2 * halfcell * di);
  cmp->SetPeriodicityY(2 * halfcell);
  string cellType;
  cellType = cmp->GetCellType();
  cout << "cell type: " << cellType << endl;

  const double vSense = 2000;
  const double vField = 0.;

  cmp->AddWire(0., 0., 2 * rSense, vSense, "s");
  cmp->AddWire(0., halfcell, 2 * rField, vField, "t");
  cmp->AddWire(halfcell * di, p * halfcell, 2 * rField, vField, "t");
  cmp->AddWire(0., -halfcell, 2 * rField, vField, "t");
  cmp->AddWire(-halfcell * di, p * halfcell, 2 * rField, vField, "t");
  cmp->AddWire(-halfcell * di - ex, halfcell, 2 * rField, vField, "t");
  cmp->AddWire(-halfcell * di - ex, -halfcell, 2 * rField, vField, "t");
  cmp->AddWire(halfcell * di - ex, halfcell, 2 * rField, vField, "t");
  cmp->AddWire(halfcell * di - ex, -halfcell, 2 * rField, vField, "t");
  cout << "wire number: " << cmp->GetNumberOfWires() << endl;
  cmp->AddReadout("s");
  // get electric and magnetic field strength at specific position
  cmp->MagneticField(0.3, 0.3, 0, bx, by, bz, st);
  // cout << "magnetic field at pos(0.3,0.3,0): " << bx << "  " << by << "  " << bz << endl;
  std::array<double, 3> electric;
  electric = cmp->ElectricField(0.01, 0.01, 0);
  // cout << "electric field at pos(0.01,0.01,0): " << electric[0] << "  " << electric[1] << "  " << electric[2] << endl;

  // ------- electric field visualization -------

  plottingEngine.SetDefaultStyle();
  ViewField *fieldView = new ViewField();
  fieldView->SetComponent(cmp);
  fieldView->SetPlane(0., 0., 1., 0., 0., 0.);
  fieldView->SetArea(-width * 2, -width * 2, width * 2, width * 2);
  fieldView->SetNumberOfSamples2d(1000, 1000);
  fieldView->SetNumberOfContours(20);
  fieldView->SetElectricFieldRange(0, 20000);
  TCanvas *field_can = nullptr;
  if (drawOption)
  {
    field_can = new TCanvas("electric_field", "", 920, 680);
    fieldView->SetCanvas(field_can);
    fieldView->PlotContour("e");
  }
  if (saveOption)
    field_can->SaveAs(Form("field_%s.png", gasname.c_str()));

  // ------- sensor for induced current signal -------

  Sensor *sensor = new Sensor();
  sensor->AddComponent(cmp);
  sensor->AddElectrode(cmp, "s");
  // sensor->ConvoluteSignals();
  // sensor->GetIonSignal();
  // sensor->IntegrateSignal();
  const double tMin = 0.;
  const double tMax = 550.;
  const double tStep = 0.01;
  const int nTimeBins = int((tMax - tMin) / tStep);
  sensor->SetTimeWindow(0., tStep, nTimeBins);

  // ------- primary track definition && drift process simulation -------

  gRandom->SetSeed(time(NULL));
  TrackHeed *track = new TrackHeed();
  track->SetSensor(sensor);
  track->SetParticle("pion");
  track->EnableMagneticField();
  track->EnableDeltaElectronTransport();

  // AvalancheMC* aval = new AvalancheMC();
  // AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  // aval->SetSensor(sensor);
  // aval->EnableSignalCalculation();
  // aval->SetTimeWindow(0.05); //test
  // aval->EnableMagneticField();

  AvalancheMC *drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->EnableSignalCalculation();
  drift->SetCollisionSteps(300);

  // drift->SetDistanceSteps(1.e-3);
  // drift->EnableMagneticField();
  // drift->SetTimeSteps(5);

  // DriftLineRKF *driftRKF = new DriftLineRKF();
  // driftRKF->SetSensor(sensor);
  // driftRKF->EnableSignalCalculation();
  // driftRKF->EnableAvalanche(false);

  // ------- geometry, track and signal visualization -------
  ViewDrift driftView;
  ViewSignal signalView;
  TCanvas *drift_can = nullptr;
  TCanvas *signal_can = nullptr;
  if (drawOption)
  {
    drift_can = new TCanvas("drift", "", 960, 920);
    driftView.SetCanvas(drift_can);
    // drift->EnablePlotting(&driftView);
    drift->EnablePlotting(&driftView);
    track->EnablePlotting(&driftView);
    ViewCell cellView;
    // cellView.EnableWireMarkers(false);
    cellView.SetCanvas(drift_can);
    cellView.SetComponent(cmp);
    cellView.SetArea(-width / 2 * di, -width / 2, width / 2 * di, width / 2);
    cellView.Plot2d();
  }
  if (drawOption)
  {
    signalView.SetSensor(sensor);
    signal_can = new TCanvas("signal", "", 920, 680);
    signalView.SetCanvas(signal_can);
  }
  // track level loop
  int trackNubmer = 1;
  double xe1, ye1, ze1, te1, e1;
  double xe2, ye2, ze2, te2, e2;
  int status;

  double x, y, z, te, e, xe, ye, ze, ax, ay, az;
  double xc = 0., yc = 0., zc = 0., tc = 0.;
  int nc = 0;
  double ec = 0.;
  double d;

  double extra = 0.;
  double min = 1e6;
  for (int trackloop = 0; trackloop < trackNubmer; trackloop++)
  {
    cout << "track id: " << trackloop << endl;

    // d = gRandom->Uniform(-width / 2, width / 2);
    d = (sqrt(2) - 1) / 2 * width;
    xstart = d;
    // momentum = gRandom->Uniform(1.e8, 1.5e8);
    momentum = 1e8;
    // angle = gRandom->Uniform(M_PI / 18., M_PI / 2.);
    angle = Garfield::Pi / 4;
    // distance between track and origin? what's the difference between d > or < 0;
    // if (d < 0)
    // {
    //   distance = (-d - 0.6 / tan(angle)) * sin(angle);
    // }
    // else
    // {
    //   distance = (d + 0.6 / tan(angle)) * sin(angle);
    // }
    distance = fabs(0.6 + d * tan(angle)) * cos(angle);
    track->SetMomentum(momentum);
    track->NewTrack(d, -width / 2, 0, 0, cos(angle), sin(angle), 0);
    // cout << " density    " << track->GetClusterDensity() << endl;
    // cout << " energy loss /cm   " << track->GetStoppingPower() << endl;

    // ------- drift electrons in clusters of primary track -------

    const vector<TrackHeed::Cluster> &clusters = track->GetClusters();
    cout << "number of clusters: " << clusters.size() << endl;
    // cluster level loop
    for (vector<TrackHeed::Cluster>::const_iterator ite_cluster = clusters.begin(); ite_cluster != clusters.end(); ite_cluster++)
    {
      cout << " cluster id: " << ite_cluster - clusters.begin() << endl;
      xc = ite_cluster->x;
      yc = ite_cluster->y;
      zc = ite_cluster->z;
      tc = ite_cluster->t;
      nc = ite_cluster->electrons.size();
      ec = ite_cluster->energy;
      extra = ite_cluster->extra;
      // cout << "number of electrons: " << nc << endl;
      // electron level loop
      for (vector<TrackHeed::Electron>::const_iterator ite_electron = ite_cluster->electrons.begin(); ite_electron != ite_cluster->electrons.end(); ite_electron++)
      {
        // drift->DriftElectron(ite_electron->x, ite_electron->y, ite_electron->z, ite_electron->t);
        // cout << "drift end point number: " << drift->GetElectrons().size() << endl;
        // drift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);
        drift->AvalancheElectron(ite_electron->x, ite_electron->y, ite_electron->z, ite_electron->t);
        drift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);
        if (fabs(xe2) < 0.05 && fabs(ye2) < 0.05)
        {
          reachtime.push_back(te2 - te1);
          if (te2 - te1 < min)
            min = te2 - te1;
        }
      }
    }
    firstTime = min;
    Data->Fill();
    reachtime.clear();
    min = 1e6;

    // ------- plot geometry, track, drift and signal -------

    if (drawOption)
    {
      constexpr bool twod = true;
      constexpr bool drawaxis = false;
      driftView.Plot(twod, drawaxis);
      if (saveOption)
        drift_can->SaveAs(Form("drift_%s.png", gasname.c_str()));
      // sleep(5);
      driftView.Clear();
      signalView.PlotSignal("s");
      if (saveOption)
        signal_can->SaveAs(Form("signal_%s.png", gasname.c_str()));
    }
    sensor->ClearSignal();
  }

  outfile->Write();
  outfile->Close();
}

int main(int argc, char *argv[])
{
  int argnum = argc;
  string filename;
  int fileid;
  if (argnum == 2)
    filename = argv[1];
  if (argnum == 3)
  {
    filename = argv[1];
    fileid = std::stoi(argv[2]);
  }
  TApplication app("app", &argc, argv);
  if (argnum == 1)
  {
    start_run();
  }
  else if (argnum == 2)
  {
    start_run(filename);
  }
  else if (argnum == 3)
  {
    start_run(filename, fileid);
  }
  else
    cout << "input parameter number error!" << endl;
  cout << "program ends with time consumption: " << clock() / CLOCKS_PER_SEC << "s." << endl;
  app.Run();
  return 0;
}
