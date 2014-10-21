// -*- C++ -*-
//
// Package:    H2TestBeamAnalyzer
// Class:      H2TestBeamAnalyzer
// 
/**\class H2TestBeamAnalyzer H2TestBeamAnalyzer.cc UserCode/H2TestBeamAnalyzer/src/H2TestBeamAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Viktor Khristenko,510 1-004,+41227672815,
//         Created:  Tue Sep 16 15:47:09 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBBeamCounters.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBEventPosition.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBParticleId.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTiming.h"

#include "RecoTBCalo/HcalTBObjectUnpacker/interface/HcalTBTriggerDataUnpacker.h"
#include "RecoTBCalo/HcalTBObjectUnpacker/interface/HcalTBSlowDataUnpacker.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define NUMCHS 2000
#define NUMADCS 128
#define NUMTS 50
double adc2fC[NUMADCS]={
	-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5, 10.5,11.5,12.5,
	13.5,15.,17.,19.,21.,23.,25.,27.,29.5,32.5,35.5,38.5,42.,46.,50.,54.5,59.5,
	64.5,59.5,64.5,69.5,74.5,79.5,84.5,89.5,94.5,99.5,104.5,109.5,114.5,119.5,
	124.5,129.5,137.,147.,157.,167.,177.,187.,197.,209.5,224.5,239.5,254.5,272.,
	292.,312.,334.5,359.5,384.5,359.5,384.5,409.5,434.5,459.5,484.5,509.5,534.5,
	559.5,584.5,609.5,634.5,659.5,684.5,709.5,747.,797.,847.,897.,947.,997.,
	1047.,1109.5,1184.5,1259.5,1334.5,1422.,1522.,1622.,1734.5,1859.5,1984.5,
	1859.5,1984.5,2109.5,2234.5,2359.5,2484.5,2609.5,2734.5,2859.5,2984.5,
	3109.5,3234.5,3359.5,3484.5,3609.5,3797.,4047.,4297.,4547.,4797.,5047.,
	5297.,5609.5,5984.5,6359.5,6734.5,7172.,7672.,8172.,8734.5,9359.5,9984.5};

struct TCalibLedInfo
{
	int numChs;
	int iphi[50];
	int ieta[50];
	int cBoxChannel[50];
	vector<string> cBoxString;
	int nTS[50];
	double pulse[50][10];
};

struct THFInfo
{
	int numChs;
	int numTS;
	int iphi[NUMCHS];
	int ieta[NUMCHS];
	int depth[NUMCHS];
	double pulse[NUMCHS][NUMTS];
};

struct H2Triggers
{
	//
	//	Standard Triggers
	//
	int ped;
	int led;
	int laser;
	int beam;
	string str;

	//
	//	Added for completeness
	//
	int fakeTrg;
	int inSpillTrg;
};

struct H2BeamCounters
{
	double cer1adc;
	double cer2adc;
	double cer3adc;
	double s1adc;
	double s2adc;
	double s3adc;
	double s4adc;
};

struct H2Timing
{
	int s1Count;
	int s2Count;
	int s3Count;
	int s4Count;

	double triggerTime;
	double ttcL1Atime;
};

//
// class declaration
//

class H2TestBeamAnalyzer : public edm::EDAnalyzer {
   public:
      explicit H2TestBeamAnalyzer(const edm::ParameterSet&);
      ~H2TestBeamAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
	  void getData(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

		TFile *_file;
		TTree *_tree;
		TTree *_treeHF;
		TTree *_treeTriggers;
		TTree *_treeWC;
		TTree *_treeBC;
		TTree *_treeTiming;


		string _outFileName;
		int _verbosity;
		TCalibLedInfo _calibInfo;
		THFInfo _hfInfo;

		H2Triggers _triggers;
		H2BeamCounters _BCData;
		H2Timing _timing;

		vector<double> wcX[5];
		vector<double> wcY[5];

		TH1D *x[5];
		TH1D *y[5];
		TH1D *s1, *s2, *s3, *s4;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
H2TestBeamAnalyzer::H2TestBeamAnalyzer(const edm::ParameterSet& iConfig) :
	_outFileName(iConfig.getUntrackedParameter<string>("OutFileName")),
	_verbosity(iConfig.getUntrackedParameter<int>("Verbosity"))
{
   //now do what ever initialization is needed

	_file = new TFile(_outFileName.c_str(), "recreate");
	_file->mkdir("HFData");
	_file->mkdir("Triggers");
	_file->mkdir("WCData");
	_file->mkdir("Timing");
	_file->mkdir("BeamCounters");
	_file->mkdir("Histos");

/*
	_file->cd("DiodeInfo");
	_tree = new TTree("Events", "Events");
	_tree->Branch("numChs", &_calibInfo.numChs, "numChs/I");
	_tree->Branch("iphi", _calibInfo.iphi, "iphi[numChs]/I");
	_tree->Branch("ieta", _calibInfo.ieta, "ieta[numChs]/I");
	_tree->Branch("cBoxChannel", _calibInfo.cBoxChannel, "cBoxChannel[numChs]/I");
	_tree->Branch("cBoxString", &_calibInfo.cBoxString);
	_tree->Branch("nTS", _calibInfo.nTS, "nTS[numChs]/I");
	_tree->Branch("pulse", _calibInfo.pulse, "pulse[numChs][10]/D");
*/

	_file->cd("HFData");
	_treeHF = new TTree("Events", "Events");
	_treeHF->SetAutoSave(10000000);
	_treeHF->Branch("numChs", &_hfInfo.numChs, "numChs/I");
	_treeHF->Branch("numTS", &_hfInfo.numTS, "numTS/I");
	_treeHF->Branch("iphi", _hfInfo.iphi, "iphi[numChs]/I");
	_treeHF->Branch("ieta", _hfInfo.ieta, "ieta[numChs]/I");
	_treeHF->Branch("depth", _hfInfo.depth, "depth[numChs]/I");
	_treeHF->Branch("pulse", _hfInfo.pulse, "pulse[numChs][50]/D");

	_file->cd("Triggers");
	_treeTriggers = new TTree("Events", "Events");
	_treeTriggers->Branch("ped", &_triggers.ped, "ped/I");
	_treeTriggers->Branch("led", &_triggers.led, "led/I");
	_treeTriggers->Branch("laser", &_triggers.laser, "laser/I");
	_treeTriggers->Branch("beam", &_triggers.beam, "beam/I");
//	_treeTriggers->Branch("str", &_triggers.str, "char[200]")
	_treeTriggers->Branch("fakeTrg", &_triggers.fakeTrg, "fakeTrg/I");
	_treeTriggers->Branch("inSpillTrg", &_triggers.inSpillTrg, "inSpillTrg/I");

	_file->cd("WCData");
	_treeWC = new TTree("Events", "Events");
	_treeWC->Branch("xA", &(wcX[0]));
	_treeWC->Branch("yA", &(wcY[0]));
	_treeWC->Branch("xB", &(wcX[1]));
	_treeWC->Branch("yB", &(wcY[1]));
	_treeWC->Branch("xC", &(wcX[2]));
	_treeWC->Branch("yC", &(wcY[2]));
	_treeWC->Branch("xD", &(wcX[3]));
	_treeWC->Branch("yD", &(wcY[3]));
	_treeWC->Branch("xE", &(wcX[4]));
	_treeWC->Branch("yE", &(wcY[4]));

	_file->cd("Timing");
	_treeTiming = new TTree("Events", "Events");
	_treeTiming->Branch("s1Count", &_timing.s1Count, "s1Count/I");
	_treeTiming->Branch("s2Count", &_timing.s2Count, "s2Count/I");
	_treeTiming->Branch("s3Count", &_timing.s3Count, "s3Count/I");
	_treeTiming->Branch("s4Count", &_timing.s4Count, "s4Count/I");
	_treeTiming->Branch("triggerTime", &_timing.triggerTime, "triggerTime/D");
	_treeTiming->Branch("ttcL1Atime", &_timing.ttcL1Atime, "ttcL1Atime/D");


	_file->cd("Histos");
	s1 = new TH1D("s1", "s1", 10000, 0, 1000);
	s2 = new TH1D("s2", "s2", 10000, 0, 1000);
	s3 = new TH1D("s3", "s3", 10000, 0, 1000);
	s4 = new TH1D("s4", "s4", 10000, 0, 1000);
	x[0] = new TH1D("xA", "xA", 10000, -100, 100);
	x[1] = new TH1D("xB", "xB", 10000, -100, 100);
	x[2] = new TH1D("xC", "xC", 10000, -100, 100);
	x[3] = new TH1D("xD", "xD", 10000, -100, 100);
	x[4] = new TH1D("xE", "xE", 10000, -100, 100);
	y[0] = new TH1D("yA", "yA", 10000, -100, 100);
	y[1] = new TH1D("yB", "yB", 10000, -100, 100);
	y[2] = new TH1D("yC", "yC", 10000, -100, 100);
	y[3] = new TH1D("yD", "yD", 10000, -100, 100);
	y[4] = new TH1D("yE", "yE", 10000, -100, 100);

	_file->cd("BeamCounters");
	_treeBC = new TTree("Events", "Events");
	_treeBC->Branch("s1adc", &_BCData.s1adc, "s1adc/D");
	_treeBC->Branch("s2adc", &_BCData.s2adc, "s2adc/D");
	_treeBC->Branch("s3adc", &_BCData.s3adc, "s3adc/D");
	_treeBC->Branch("s4adc", &_BCData.s4adc, "s4adc/D");
}


H2TestBeamAnalyzer::~H2TestBeamAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

	_file->cd("Histos");
	s1->Write();
	s2->Write();
	s3->Write();
	s4->Write();
	x[0]->Write();
	x[1]->Write();
	x[2]->Write();
	x[3]->Write();
	x[4]->Write();
	y[0]->Write();
	y[1]->Write();
	y[2]->Write();
	y[3]->Write();
	y[4]->Write();

	_file->Write();
	_file->Close();
}

void H2TestBeamAnalyzer::getData(const edm::Event &iEvent, 
		const edm::EventSetup &iSetup)
{
	using namespace edm;

	//
	//	Extracting All the Collections containing useful Info
	//
	edm::Handle<HFDigiCollection> hfDigiCollection;
	edm::Handle<HBHEDigiCollection> hbheDigiCollection;
	edm::Handle<HODigiCollection> hoDigiCollection;
	edm::Handle<HcalTBTriggerData> trigData;
	edm::Handle<HcalTBEventPosition> eventPos;
	edm::Handle<HcalTBBeamCounters> beamCounters;
	edm::Handle<HcalTBParticleId> pId;
	edm::Handle<HcalTBTiming> timing;

	iEvent.getByType(hfDigiCollection);	
	iEvent.getByType(hbheDigiCollection);
	iEvent.getByType(hoDigiCollection);
	iEvent.getByType(trigData);
	iEvent.getByLabel("tbunpack", eventPos);
	iEvent.getByType(beamCounters);
	iEvent.getByType(pId);
	iEvent.getByType(timing);

	//
	//	Extract Trigger Info
	//
	_triggers.ped = 0;
	_triggers.led = 0;
	_triggers.laser = 0;
	_triggers.beam = 0;
	_triggers.fakeTrg = 0;
	_triggers.inSpillTrg = 0;
	_triggers.str = trigData->runNumberSequenceId();
	if (trigData->wasInSpillPedestalTrigger() || 
			trigData->wasOutSpillPedestalTrigger() ||
			trigData->wasSpillIgnorantPedestalTrigger())
		_triggers.ped = 1;
	if (trigData->wasLEDTrigger())
		_triggers.led = 1;
	if (trigData->wasLaserTrigger())
		_triggers.laser = 1;
	if (trigData->wasBeamTrigger())
		_triggers.beam = 1;
	if (trigData->wasFakeTrigger())
		_triggers.fakeTrg = 1;
	if (trigData->wasInSpill())
		_triggers.inSpillTrg = 1;

	//
	//	Extract Event Position
	//
//	double tableX = eventPos->hfTableX();
//	double tableY = eventPos->hfTableY();
//	double tableV = eventPos->hfTableV();
	vector<double> xxx, yyy;
	eventPos->getChamberHits('A', xxx, yyy);
	wcX[0] = xxx;
	wcY[0] = yyy;
	xxx.clear();
	yyy.clear();

	eventPos->getChamberHits('B', xxx, yyy);
	wcX[1] = xxx;
	wcY[1] = yyy;
	xxx.clear(); yyy.clear();

	eventPos->getChamberHits('C', xxx, yyy);
	wcX[2] = xxx;
	wcY[2] = yyy;
	xxx.clear();
	yyy.clear();

	eventPos->getChamberHits('D', xxx, yyy);
	wcX[3] = xxx;
	wcY[3] = yyy;
	xxx.clear();
	yyy.clear();

	eventPos->getChamberHits('E', xxx, yyy);
	wcX[4] = xxx;
	wcY[4] = yyy;
	xxx.clear();
	yyy.clear();

	for (int i=0; i<5; i++)
	{
		for (vector<double>::iterator it=wcX[i].begin(); 
				it!=wcX[i].end(); ++it)
			x[i]->Fill(*it);
		for (vector<double>::iterator it=wcY[i].begin(); 
				it!=wcY[i].end(); ++it)
			y[i]->Fill(*it);
	}


	//
	//	Extract Beam Counters Info
	//
	_BCData.cer1adc = beamCounters->CK1adc();	
	_BCData.cer2adc = beamCounters->CK2adc();
	_BCData.cer3adc = beamCounters->CK3adc();
	_BCData.s1adc = beamCounters->S1adc();
	_BCData.s2adc = beamCounters->S2adc();
	_BCData.s3adc = beamCounters->S3adc();
	_BCData.s4adc = beamCounters->S4adc();

	s1->Fill(_BCData.s1adc);
	s2->Fill(_BCData.s2adc);
	s3->Fill(_BCData.s3adc);
	s4->Fill(_BCData.s4adc);

	//
	//	Extract Timing Info	
	//
	_timing.s1Count = timing->S1Count();
	_timing.s2Count = timing->S2Count();
	_timing.s3Count = timing->S3Count();
	_timing.s4Count = timing->S4Count();
	_timing.triggerTime = timing->triggerTime();
	_timing.ttcL1Atime = timing->ttcL1Atime();
	
	int numChs = 0;
	if (_verbosity>0)
	{
		cout << "### Before Loop: " << endl
			<< "### HF Digis=" << hfDigiCollection->size() << endl
			<< "### HBHE Digis=" << hbheDigiCollection->size() << endl
			<< "### HO Digis=" << hoDigiCollection->size() << endl;

		cout << "### Triggers: " << endl
			<< "### PED Trigger: " << _triggers.ped << endl
			<< "### Led Trigger: " << _triggers.led << endl
			<< "### Laser Trigger: " << _triggers.laser << endl
			<< "### Beam Trigger: " << _triggers.beam << endl
			<< "### In Spill Trigger: " << _triggers.inSpillTrg << endl;

		cout << "### Wire Chamber A: NHits: " << wcX[0].size() 
			<< "  " << wcY[0].size() << endl
			<< "### Wire Chamber B: NHits=" << wcX[1].size() 
			<< "  " << wcY[1].size() << endl
			<< "### Wire Chamber C: NHits=" << wcX[2].size() 
			<< "  " << wcY[2].size() << endl
			<< "### Wire Chamber D: NHits=" << wcX[3].size() 
			<< "  " << wcY[3].size() << endl
			<< "### Wire Chamber E: NHits=" << wcX[4].size() 
			<< "  " << wcY[4].size() << endl;
//		cout << "HF Table Pos(X, Y, V) = " << tableX << "  " << tableY << "  " <<
//			tableV << endl;
			

		cout << "### Beam Counters:" << endl
			<< "### Cerenkov1: " << _BCData.cer1adc << endl
			<< "### Cerenkov2: " << _BCData.cer2adc << endl
			<< "### Cerenkov3: " << _BCData.cer3adc << endl
			<< "### s1adc=" << _BCData.s1adc << endl
			<< "### s2adc=" << _BCData.s2adc << endl
			<< "### s3adc=" << _BCData.s3adc << endl
			<< "### s4adc=" << _BCData.s4adc << endl;

		cout << "### Beam Timing: " << endl
			<< "### S1Count: " << _timing.s1Count << endl
			<< "### S2Count: " << _timing.s2Count << endl
			<< "### S3Count: " << _timing.s3Count << endl
			<< "### S4Count: " << _timing.s4Count << endl;
	}

	
	for (HFDigiCollection::const_iterator digi=hfDigiCollection->begin();
			digi!=hfDigiCollection->end(); ++digi)
	{
//		cout << "Processing Digi" << endl;

		int iphi = digi->id().iphi();
		int ieta = digi->id().ieta();
		int depth = digi->id().depth();
		int nTS = digi->size();

		int fiberChanId = digi->elecId().fiberChanId();
		int fiberIndex = digi->elecId().fiberIndex();
		int slbChannelIndex = digi->elecId().slbChannelIndex();
		int slbSiteNumber = digi->elecId().slbSiteNumber();
		string slbChannelCode = digi->elecId().slbChannelCode();
		int htrChanId = digi->elecId().htrChanId();
		int spigot = digi->elecId().spigot();
		int dccid = digi->elecId().dccid();
		int htrSlot = digi->elecId().htrSlot();
		int htrTopBottom = digi->elecId().htrTopBottom();
		int readoutVMECrateId = digi->elecId().readoutVMECrateId();
		int linearIndex = digi->elecId().linearIndex();

		if (_verbosity>1)
		{
			cout << "### Digi->elecId:" << endl;
			cout << fiberChanId << "  " << fiberIndex << "  "
				<< slbChannelIndex << "  " << slbSiteNumber << "  "
				<< slbChannelCode << "  " << htrChanId << "  "
				<< spigot << "  " << dccid << "  " << htrSlot << "  "
				<< htrTopBottom << "  " << readoutVMECrateId << "  "
				<< linearIndex << endl;
			cout << "### Digi->detId:" << endl;
			cout << iphi << "  " << ieta << "  " << depth << endl;
		}

		//
		//	Set the Branched arrays
		//
		_hfInfo.iphi[numChs] = iphi;
		_hfInfo.ieta[numChs] = ieta;
		_hfInfo.depth[numChs] = depth;
		_hfInfo.numTS = nTS;
	
		for (int iTS=0; iTS<nTS; iTS++)
			_hfInfo.pulse[numChs][iTS] = adc2fC[digi->sample(iTS).adc()&0xff];

		if (_verbosity>1)
		{
			cout << "### Digi->Data:" << endl;
			for (int iTS=0; iTS<nTS; iTS++)
				cout << _hfInfo.pulse[numChs][iTS] << "  ";
			cout << endl;
		}

		numChs++;
	}
	_hfInfo.numChs = numChs;
	_treeHF->Fill();
	_treeTriggers->Fill();
	_treeWC->Fill();
	_treeBC->Fill();
	_treeTiming->Fill();

	return;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
H2TestBeamAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	getData(iEvent, iSetup);

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
H2TestBeamAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
H2TestBeamAnalyzer::endJob() 
{
//	_file->Write();
//	_file->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
H2TestBeamAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
H2TestBeamAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
//	_file->Write();
//	_file->Close();
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
H2TestBeamAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
H2TestBeamAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
H2TestBeamAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(H2TestBeamAnalyzer);
