#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TLorentzVector.h>

#include "T49Run.h"
#include "T49EventRoot.h"
#include "T49EventMCRoot.h"
#include "T49VetoRoot.h"
#include "T49VertexRoot.h"
#include "T49VertexMCRoot.h"
#include "T49ParticleRoot.h"
#include "T49ParticleMCRoot.h"
#include "T49Vertex.h"
#include "T49Ring.h"

#include "DataTreeEvent.h"
#include "productionMap.h"
#include "particleCodes.h"

using namespace std;

bool isSim = false;
T49Run *run;
TTree *tree;
DataTreeEvent *fDTEvent;
TTree *myTree;
T49EventRoot *event;
T49ParticleRoot *particle;
T49RingRoot *ring;
T49VertexRoot *vertex;
T49VertexMCRoot *MCvertex;
T49VetoRoot *veto;
TObjArray *particles;
TLorentzVector TrackPar;
T49EventMCRoot *MCevent;
T49ParticleMCRoot *MCparticle;
TObjArray *MCparticles;

void ReadEvent ();
void ReadMCEvent ();
int T49_to_DT(string inputFileList, string productionTag, string outputPath = ".", int maxFileNumber = 10000);

int main(int argc, char *argv[])
{
  TStopwatch timer;
  timer.Reset();
  timer.Start();

  switch (argc)
  {
    case 3:
    return T49_to_DT(argv[1], argv[2]);
    break;
    case 4:
    return T49_to_DT(argv[1], argv[2], argv[3]);
    break;
    case 5:
    return T49_to_DT(argv[1], argv[2], argv[3], atoi(argv[4]));
    break;
    default:
    cout << "Specify input filelist and production tag!" << endl;
    cout << argv[0] << " inputFileList productionTag [outputPath] [maxFileNumber]" << endl << endl;
    return 0;
  }

  timer.Stop();
  printf("Real time: %f\n",timer.RealTime());
  printf("CPU time: %f\n",timer.CpuTime());
  return 1;
}

int T49_to_DT(string inputFileList, string productionTag, string outputPath, int maxFileNumber)
{
  if (prodMap.count(productionTag) == 0)
    {
    cout << "There is no such production tag as " << productionTag << endl;
//    return 0;
    }

  string runType = prodMap [productionTag];
  outputPath = outputPath + "/" + runType;
  string command = "mkdir -p " + outputPath;
  gSystem -> Exec (command.c_str ());
  cout << "inputFileList: " << inputFileList << endl;
  cout << "outputPath: " << outputPath << endl;

  cout << runType << endl;

  myTree = new TTree("DataTree", "NA49 Pb-Pb DataTree");
  fDTEvent = new DataTreeEvent();
  myTree -> Branch("DTEvent", "DataTree", fDTEvent);

  run = new T49Run();
  run->SetVerbose(0);

  int fileNumber = 0;

  ifstream inputFileListStream (inputFileList);
  if (!inputFileListStream) cout << "list.txt not found!\n";
  string inputFileName;
  while (getline(inputFileListStream, inputFileName) && fileNumber < maxFileNumber)
  {
    string outputFileName = inputFileName.substr(inputFileName.rfind("t49run" + 1));
    outputFileName = outputPath + "/DT" + outputFileName;

    TFile* outputFile = new TFile((Char_t*) outputFileName.c_str(),"recreate");
    tree = myTree->CloneTree();
    fileNumber++;
    cout << "Opening file " << fileNumber << ": " << inputFileName.c_str() << endl;
    if(!run -> Open ((Char_t*) (inputFileName).c_str())) continue;
    if (runType.find ("VENUS") != string::npos)
    {
        isSim = true;
        cout << "isSim = true!!!\n";
    }

    if (runType.find ("40") != string::npos)
    {
        cout << "Veto calibration for 40A GeV\n";
        run->SetEvetoCal (3); // set 1 for other energies?, 0 for no time-dependent calibration
    }
    else if (runType.find ("160") != string::npos)
    {
        cout << "Veto calibration for 160A GeV\n"; // set 1 for other energies?, 0 for no time-dependent calibration
        run->SetEvetoCal (2);
    }
    else
    {
        cout << "Veto calibration for energies other than 40A and 160A GeV\n"; // set 1 for other energies?, 0 for no time-dependent calibration
        run->SetEvetoCal (1);
    }
run->SetEvetoCal (1); // patch
    run -> SetRunType ((Char_t*) (runType.c_str ()));
    event = (T49EventRoot*) run -> GetNextEvent ();
    run -> SetRunNumber (event -> GetNRun ());
    run -> CheckEvetoCal ();
    cout << "Run number: " << run -> GetRunNumber () << endl;
    cout << "Centrality of first event: " << event -> GetCentrality () << endl << endl;

    int nEvents = run -> GetMaxEvents();
    for (int nEvent = 0; nEvent < nEvents; nEvent++)
    {
        cout << "\rEvent " << nEvent;
        fDTEvent -> ClearEvent ();
        ReadEvent ();
        if (isSim) ReadMCEvent ();
        tree -> Fill ();
        run -> GetNextEvent ();
    }
    outputFile -> cd ();
  	tree -> Write ();
 	outputFile -> Close ();
    run -> Close ();
    cout << endl;
  }

  return 1;
}

void ReadEvent ()
{
    event = (T49EventRoot*) run -> GetEvent ();

    fDTEvent -> SetRunId( event->GetNRun() );
    fDTEvent -> SetEventId( event->GetNEvent() );
    fDTEvent -> SetVertexPosition(event->GetVertexX(), event->GetVertexY(), event->GetVertexZ(), EnumVertexType::kReconstructedVertex);

    vertex = ((T49VertexRoot*)event -> GetPrimaryVertex());
		fDTEvent -> SetHasVertex ();
    fDTEvent -> SetVertexQuality ( vertex->GetPchi2(), 0 );
    fDTEvent -> SetVertexQuality ( vertex->GetIflag(), 1 );
    fDTEvent -> AddTrigger();
    if (vertex->GetIflag())
      fDTEvent -> GetTrigger(0)->SetIsFired(1);
    else
      fDTEvent -> GetTrigger(0)->SetIsFired(0);
    fDTEvent -> SetPsdPosition(0., 0., 1800.);
    fDTEvent -> SetPsdEnergy(event -> GetTDCalEveto());
//    fDTEvent -> SetRPAngle ( event -> GetCentralityClass () );
    fDTEvent -> SetImpactParameter ( event -> GetCentrality () );

//      fTriggerMask =(event->GetTriggerMask()); // Trigger Mask
//      fDate=(event->GetDate()); // Event date
//      fTime =(event->GetTime()); // Event time
//      fWfaNbeam=(event->GetWfaNbeam()); // returns the number of the beam particles (including the interacting one) arriving in the ±25µs around the interaction time
//      fWfaNinter=(event->GetWfaNinter()); // returns the number of the interacting particles (including the one that triggered the event) arriving in the ±25µs around the interaction time
//      fWfaBeamTime=(event->GetWfaBeamTime() ); // returns the array containing the times in nanoseconds for all beam particles in the ±25µs around the interaction time. The length of this array is given by the method above
//      fWfaInterTime=(event->GetWfaInterTime()); //returns the array containing the times in nanoseconds for all interacting particles in the ±25µs around the interaction time. The length of this array is given by the method above

    float vetoModuleAngle [4] = {1.25, 1.75, 0.25, 0.75};
    float psdModuleEnergy;
    veto = (T49VetoRoot*)event->GetVeto();
    for (unsigned int iVeto=0; iVeto<4; ++iVeto)
    {
        fDTEvent -> AddPSDModule(0);
        fDTEvent -> GetLastPSDModule() -> SetPosition(25. * cos(vetoModuleAngle [iVeto] * TMath::Pi()),
                                                      25. * sin(vetoModuleAngle [iVeto] * TMath::Pi()),
                                                      2600.);
        psdModuleEnergy = veto->GetADChadron(iVeto)+veto->GetADCphoton(iVeto);
        //psdModuleEnergy = veto->GetADChadron(iVeto);
        if (psdModuleEnergy < 0.) psdModuleEnergy = 0;
        fDTEvent -> GetLastPSDModule() -> SetEnergy(psdModuleEnergy);
	//cout << iVeto << "\t" << psdModuleEnergy << endl;
    }

    ring = (T49RingRoot *)event->GetRing();

    int PSDModuleId;

    float ring_middle_radius[10]={30.5,36,42.5,50.2,59.4,70.3,83,98.1,116,137};

//    for (unsigned int iRing=0; iRing<240; ++iRing)
//    {
//        fDTEvent -> AddPSDModule(0);
////        fDTEvent -> GetLastPSDModule() -> SetPosition( cos( (iRing+0.5)*TMath::Pi()/12 ), sin( (iRing+0.5)*TMath::Pi()/12 ), 1800 );
//        fDTEvent -> GetLastPSDModule() -> SetPosition(-ring_middle_radius[iRing%10]*sin(((iRing/10)+0.5)*TMath::Pi()/12),
//                                                       ring_middle_radius[iRing%10]*cos(((iRing/10)+0.5)*TMath::Pi()/12), 1800);
//        fDTEvent -> GetLastPSDModule() -> SetEnergy(ring->GetADChadron(iRing));
//    }

    for (unsigned int iRing=0; iRing<240; ++iRing)
    {
        fDTEvent -> AddPSDModule(0);
    }

    for (unsigned int iRing=0; iRing<240; ++iRing)
    {
        PSDModuleId = 4 + iRing / 10 + iRing % 10 * 24;
        fDTEvent -> GetPSDModule (PSDModuleId) -> SetPosition (-ring_middle_radius [iRing%10] * sin (((iRing/10)+0.5) * TMath::Pi() / 12.),
                                                                ring_middle_radius [iRing%10] * cos (((iRing/10)+0.5) * TMath::Pi() / 12.),
                                                                1800);                                                              
        //psdModuleEnergy = ring->GetADChadron(iRing);
        psdModuleEnergy = ring->GetADChadron(iRing)+ring->GetADCphoton(iRing);
        if (psdModuleEnergy < 0.) psdModuleEnergy = 0;
        fDTEvent -> GetPSDModule (PSDModuleId) -> SetEnergy(psdModuleEnergy);
	//cout << iRing << "\t" << psdModuleEnergy << endl;
    }

    particles = (TObjArray*) event -> GetPrimaryParticles();
    TIter particles_iter(particles);
    for (int iTrack=0; particle = (T49ParticleRoot*) particles_iter.Next(); iTrack++)
    {
        fDTEvent -> AddVertexTrack();
        fDTEvent -> GetLastVertexTrack() -> SetId(iTrack);
        TrackPar.SetPtEtaPhiM (particle->GetPt(), particle->GetEta(), particle->GetPhi(), 0.14);
        fDTEvent -> GetLastVertexTrack() -> SetMomentum(TrackPar);

        fDTEvent -> GetLastVertexTrack() -> SetDCA(particle->GetBx(), particle->GetBy(), 0);
        fDTEvent -> GetLastVertexTrack() -> SetChi2(particle->GetPchi2());
        fDTEvent -> GetLastVertexTrack() -> SetNDF(1);
        fDTEvent -> GetLastVertexTrack() -> SetCharge(particle->GetCharge());

        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHits(particle->GetNFitPoint(0), EnumTPC::kVTPC1);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHits(particle->GetNFitPoint(1), EnumTPC::kVTPC2);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHits(particle->GetNFitPoint(2), EnumTPC::kMTPC);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHits(particle->GetNFitPoint(), EnumTPC::kTPCAll);

        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHitsPotential(particle->GetNMaxPoint(0), EnumTPC::kVTPC1);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHitsPotential(particle->GetNMaxPoint(1), EnumTPC::kVTPC2);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHitsPotential(particle->GetNMaxPoint(2), EnumTPC::kMTPC);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfHitsPotential(particle->GetNMaxPoint(), EnumTPC::kTPCAll);

        fDTEvent -> GetLastVertexTrack() -> SetdEdx(particle->GetTmeanCharge(0), EnumTPC::kVTPC1);
        fDTEvent -> GetLastVertexTrack() -> SetdEdx(particle->GetTmeanCharge(1), EnumTPC::kVTPC2);
        fDTEvent -> GetLastVertexTrack() -> SetdEdx(particle->GetTmeanCharge(2), EnumTPC::kMTPC);
        fDTEvent -> GetLastVertexTrack() -> SetdEdx(particle->GetTmeanCharge(), EnumTPC::kTPCAll);

//cout << "dEdx = " << particle->GetTmeanCharge(2) << endl; 
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfdEdxClusters(0.001 * particle-> GetNDedxPoint(0), EnumTPC::kVTPC1);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfdEdxClusters(0.001 * particle-> GetNDedxPoint(1), EnumTPC::kVTPC2);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfdEdxClusters(0.001 * particle-> GetNDedxPoint(2), EnumTPC::kMTPC);
        fDTEvent -> GetLastVertexTrack() -> SetNumberOfdEdxClusters(0.001 * particle-> GetNDedxPoint(), EnumTPC::kTPCAll);

        fDTEvent -> AddTOFHit();
        fDTEvent -> GetLastTOFHit() -> SetPathLength( particle->GetTofPathl() );
        fDTEvent -> GetLastTOFHit() -> SetSquaredMass( particle->GetTofMass2() );
        fDTEvent -> GetLastTOFHit() -> SetSquaredMassError( particle->GetTofSigMass2() );
        fDTEvent -> GetLastTOFHit() -> AddRecoTrackId( iTrack );
        fDTEvent -> GetLastTOFHit() -> SetCharge( particle->GetTofCharge() );
        fDTEvent -> GetLastTOFHit() -> SetPosition( particle->GetTofX(), particle->GetTofY(), 500);

        /*
        fPrimaryParticles_fIdDet[k] =(particle->GetIdDet());
              for(unsigned int j=0;j<3;j++)
              {
                fPrimaryParticles_fXFirst[k][j]=(particle->GetXFirst(j));
                fPrimaryParticles_fYFirst[k][j]=(particle->GetYFirst(j));
                fPrimaryParticles_fZFirst[k][j]=(particle->GetZFirst(j));
                fPrimaryParticles_fXLast[k][j] =(particle->GetXLast(j));
                fPrimaryParticles_fYLast[k][j] =(particle->GetYLast(j));
                fPrimaryParticles_fZLast[k][j] =(particle->GetZLast(j));
              }
          fPrimaryParticles_fXFirst[k][3]=(particle->GetXFirst());
          fPrimaryParticles_fYFirst[k][3]=(particle->GetYFirst());
          fPrimaryParticles_fZFirst[k][3]=(particle->GetZFirst());
          fPrimaryParticles_fXLast[k][3] =(particle->GetXLast());
          fPrimaryParticles_fYLast[k][3] =(particle->GetYLast());
          fPrimaryParticles_fZLast[k][3] =(particle->GetZLast());
          fPrimaryParticles_fLabel[k]=(particle->GetLabel());
          fPrimaryParticles_fTofIflag[k]= (particle->GetTofIflag());
          fPrimaryParticles_fTofId[k]=(particle->GetTofId());
          fPrimaryParticles_fIflag[k]=(particle->GetIflag());
          fPrimaryParticles_fSigPx[k]=(particle->GetSigPx());
          fPrimaryParticles_fSigPy[k]=(particle->GetSigPy());
          */
    }
}

void ReadMCEvent ()
{
    MCevent = (T49EventMCRoot*) run -> GetEvent ();
    MCvertex = (T49VertexMCRoot*) (MCevent -> GetVertices () -> At (0));
    fDTEvent -> SetMCVertexPosition (MCvertex->GetX(), MCvertex->GetY(), MCvertex->GetZ());
    fDTEvent -> SetPsdEnergy(event -> GetEveto());
    fDTEvent -> SetRunId( run -> GetRunNumber () );

    MCparticles = (TObjArray*) MCevent -> GetMCParticles();
    TIter MCparticles_iter (MCparticles);
    int priMatchedTrackIndex;
    int nPriMatchedTracks;
    int nPriMatchPoints;
    int nPoints;
    double trackPurity;

    for (int iTrack=0; MCparticle = (T49ParticleMCRoot*) MCparticles_iter.Next(); iTrack++)
    {
        if (MCparticle -> GetStartVertex() != 0) continue;
        fDTEvent -> AddMCTrack();
        TrackPar.SetXYZM (MCparticle -> GetPx(), MCparticle -> GetPy(), MCparticle -> GetPz(), MCparticle -> GetMass());
        fDTEvent -> GetLastMCTrack() -> SetMomentum(TrackPar);
        fDTEvent -> GetLastMCTrack() -> SetCharge(MCparticle->GetCharge());
        fDTEvent -> GetLastMCTrack() -> SetPdgId(GEANT_to_PDG [MCparticle -> GetPid()]);
        fDTEvent -> GetLastMCTrack() -> SetNumberOfHits(MCparticle->GetNPoint(0), EnumTPC::kVTPC1);
        fDTEvent -> GetLastMCTrack() -> SetNumberOfHits(MCparticle->GetNPoint(1), EnumTPC::kVTPC2);
        fDTEvent -> GetLastMCTrack() -> SetNumberOfHits(MCparticle->GetNPoint(2), EnumTPC::kMTPC);
        fDTEvent -> GetLastMCTrack() -> SetNumberOfHits(MCparticle->GetNPoint(), EnumTPC::kTPCAll);

        int nPriPart = particles->GetSize();
        nPriMatchedTracks = MCparticle -> GetNPriMatched ();
        for (int iMatch = 0; iMatch < nPriMatchedTracks; iMatch++)
        {
            priMatchedTrackIndex = MCparticle -> GetPriMatched (iMatch);
            nPriMatchPoints = MCparticle -> GetNPriMatchPoint (iMatch);
            nPoints = MCparticle -> GetNPoint ();
            trackPurity = (double) nPriMatchPoints / nPoints;
            fDTEvent -> AddTrackMatch (fDTEvent -> GetVertexTrack(priMatchedTrackIndex), fDTEvent -> GetLastMCTrack(), nPriMatchPoints, trackPurity);
        }
    }
}

