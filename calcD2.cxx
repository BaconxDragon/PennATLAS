//std
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <boost/range/combine.hpp>
#include <math.h>

//ROOT
#include "TTree.h"
#include "TRandom.h"
#include "TFile.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Boost.h"
#include "Math/VectorUtil.h"
#include <TH1.h>
#include <TF1.h>
#include "TCanvas.h"

//ATLAS
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include <xAODJet/Jet.h>
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticle.h"



/***
Program takes in a file specified in first command line argument
Program takes in an output ID in second command line argument
Program calculate the energy correlation functions for different jets in the file:
e2 = (1/M^2) * sum{ (Ei * Ej * angle ij) }
e3 = (1/M^3) * sum{ (Ei * Ej * Ek * angle ij * angle ik * angle jk) }
Program calculates D2 discriminant for each jet D2 = e3 / (e2^3) in the rest frame of the jet, and in the lab frame
Program stores the jet 4-vector, D2's in different frames, and decay angle to a tree
Program write the tree to a file [output ID]output.root

Version 1.3
Trying to better match the D2 Paper
-Pt cut at > 400 Gev
-Mass cut changed to < 100 Gev
Shand Seiffert 6/9/2023
***/



/***
class from riley, is for truth matching jets
holds the TruthParticle for a boson, and a TruthParticle vector of its products
holds member functions for indentifying its type of decay
holds function decayAngle which gets the angle between the boson and its nth decay product, in the boson frame
***/
struct BosonWithDecay
{
    const xAOD::TruthParticle* v = nullptr;
    std::vector<const xAOD::TruthParticle*> decay;

    bool isLep() const
    {
        for (auto c : decay)
            if (c->isLepton()) return true;
        return false;
    }

    bool isHad() const
    {
        for (auto c : decay)
            if (c->isQuark()) return true;
        return false;
    }

    bool isOk() const
    {
        if (!v) return false;
        bool hasLep = false;
        bool hasHad = false;
        for (auto c : decay)
        {
            if (c->isLepton()) hasLep = true;
            else if (c->isQuark()) hasHad = true;
        }
        return hasLep ^ hasHad;
    }

    const xAOD::TruthParticle * lepChild() const
    {
        for (auto c : decay)
            if (c->isChLepton()) return c;
        return nullptr;
    }
    const xAOD::TruthParticle * nuChild() const
    {
        for (auto c : decay)
            if (c->isNeutrino()) return c;
        return nullptr;
    }

	double decayAngle(int decayNum) 
	{
		ROOT::Math::PtEtaPhiMVector bosVec(v->pt(), v->eta(), v->phi(), v->m());
		ROOT::Math::PtEtaPhiMVector decayVec(decay[decayNum]->pt(), decay[decayNum]->eta(), decay[decayNum]->phi(), decay[decayNum]->m());

		ROOT::Math::Boost bosFrame(bosVec.BoostToCM());
    	ROOT::Math::PtEtaPhiMVector decayVecBoosted(bosFrame(decayVec));

		return ROOT::Math::VectorUtil::Angle(decayVecBoosted, bosVec);
	}
};



double getD2rest(xAOD::JetConstituentVector, ROOT::Math::PtEtaPhiMVector, long long);
//calculates d2 in the jet rest frame.
//Input: the array of particles given by JetContainer::GetConstituents, 4-vector of the jet, integer with number of particles in the jet
//Output: a double which is the D2 in the rest frame of the jet

double getD2lab(xAOD::JetConstituentVector, ROOT::Math::PtEtaPhiMVector, long long);
//calculates d2 in the lab frame.
//Input: the array of particles given by JetContainer::GetConstituents, 4-vector of the jet, integer with number of particles in the jet
//Output: a double which is the D2 in the lab frame

std::vector<BosonWithDecay> findBosons(const xAOD::TruthParticleContainer*);
//Matches up bosons with their decay products
//Input: truth particle container which holds the truth particles including bosons and their decay products
//Output: vector of the bosons and their corresponding decay products

const xAOD::Jet* findMinRJet(const xAOD::JetContainer*, ROOT::Math::PtEtaPhiMVector);
//Find the jet with the minimum Radius to the vector of interest for truth matching
//Input: the container with the jets from an event, the vector to compare them to (boson vector)
//Output: the jet with the minimum radius to the vectir



//power dependence for angles
const double betaVar = 1.0;
//radius cutoff for truth matching jets
const double deltaRcut = 0.4;



int main(int argc, char* argv[]) {

    //command line args
    if (argc < 2) return 1;
    std::string inputFilePath(argv[1]);

    //input file
    std::unique_ptr<TFile> iFile(TFile::Open(inputFilePath.c_str(), "READ"));
    if (!iFile) {
        std::cout << "Couldn't open " << inputFilePath << std::endl;
        return 1;
    }
	
    //output file
	std::string jobID(argv[2]);
	std::string outPath = jobID + "output.root";
    std::unique_ptr<TFile> oFile(TFile::Open(outPath.c_str(), "RECREATE"));
    if (!oFile) {
        std::cout << "Couldn't open " << "ExampleOutput.root" << std::endl;
        return 1;
    }

	//declaring tree and branches
	TTree * t = new TTree("tree", "D2 Test");
	double D2rest, D2lab, ptj, phij, etaj, mj, weightj, d2refj, qAngle, minR, numJets;
	t->Branch("d2rest", &D2rest);
	t->Branch("d2lab", &D2lab);
	t->Branch("pt", &ptj);
	t->Branch("phi", &phij);
	t->Branch("eta", &etaj);
	t->Branch("m", &mj);
	t->Branch("weight", &weightj);
	t->Branch("d2Reference", &d2refj);
	t->Branch("decayAngle", &qAngle);
	t->Branch("minRadius", &minR);
	t->Branch("jetCount", &numJets);

    //vector to be used for jet boost
    ROOT::Math::PtEtaPhiMVector jetFrame;
	
	//preping file for read
	xAOD::Init();
	xAOD::TEvent event;
	event.readFrom(iFile.get());

	//set number of entries to be read
	long long numEntries = event.getEntries();

	//loop over each jet
	for (long long eventNumber = 0; eventNumber < numEntries; ++eventNumber) {
		
		//getting the events and jet info from the file
		event.getEntry(eventNumber);

		const xAOD::EventInfo* ei = nullptr;
		event.retrieve(ei, "EventInfo");
		const std::vector<float> & weights = ei->mcEventWeights();
		
		const xAOD::JetContainer* jets = nullptr;
		event.retrieve(jets, "AntiKt10TruthJets");

		const xAOD::TruthParticleContainer* bosonsWithDecayParticles = nullptr;
		event.retrieve(bosonsWithDecayParticles, "TruthBosonsWithDecayParticles");

		//matches bosons to decay products
		std::vector<BosonWithDecay> eventBosons = findBosons(bosonsWithDecayParticles);
		
		//hadronically decaying boson is always the first one in the out array in this particular truth simulation
		ROOT::Math::PtEtaPhiMVector bos1(eventBosons[0].v->pt(), eventBosons[0].v->eta(), eventBosons[0].v->phi(), eventBosons[0].v->m());
		numJets = jets->size();

		//gets decay angle between the hadronic boson and its decay quark
		//will only have two decay products, two quarks, which go in opposite directions, thus it takes the smaller angle
		qAngle = eventBosons[0].decayAngle(0);
		if (qAngle > M_PI_2) qAngle = M_PI - qAngle;

		//declaring d2 value to get from the one given in file
		//float because data stored as float in xAod
		float d2 = 0.0;
		ROOT::Math::PtEtaPhiMVector jetMinVec;
		minR = deltaRcut + 1.0;

		//finding mininum distance jet between hadronically decaying boson and jet
		const xAOD::Jet* jetMin = findMinRJet(jets, bos1);
		if (jetMin) {

			ROOT::Math::PtEtaPhiMVector jetVec(jetMin->pt(), jetMin->eta(), jetMin->phi(), jetMin->m());
			jetMinVec = jetVec;
			minR = ROOT::Math::VectorUtil::DeltaR(jetMinVec, bos1);

		}

		//mass cut for only jets below 100 Gev, Pt cut above 400 Gev, radius cut to match with the hadronic boson
		if (jetMinVec.M() < 100000 && jetMinVec.pt() > 400000 && minR < deltaRcut) {
			
			//writing to tree variables
			ptj = jetMinVec.pt();
			phij = jetMinVec.phi();
			etaj = jetMinVec.eta();
			mj = jetMinVec.M();
			mj /= 1000;
			weightj = weights[0];

			//getting d2 in both frames, along with the given d2
			D2rest = getD2rest(jetMin->getConstituents(), jetMinVec, jetMin->numConstituents());
			D2lab = getD2lab(jetMin->getConstituents(), jetMinVec, jetMin->numConstituents());

			//gets D2 given in file, (only applies if taking from some branches)
			if (jetMin->getAttribute< float >("D2", d2)) d2refj = d2;
			else d2refj = -100;
			t->Fill();
			
		}

	}

	//writes and closes the file
    t->Write();
	//RPlot->Write();
    oFile->Close();
    
    return 0;

}



//matches bosons with their decay products
std::vector<BosonWithDecay> findBosons(const xAOD::TruthParticleContainer* particles) {

	std::vector<BosonWithDecay> out;
	BosonWithDecay current;

	for (auto p : *particles)
	{

		//if particle is a boson, store it
		if (p->isW() || p->isZ() || p->isHiggs())
		{

			if (!current.v) {
				current.v = p;
				
			}
			else if (current.v->pdgId() == p->pdgId() && current.v->barcode() < p->barcode()) current.v = p;
			else { // move on to next boson, wrap up old one
				out.push_back(current);
				current = BosonWithDecay();
				current.v = p;
				
			}

		}
		else current.decay.push_back(p);
		//otherwise it is a decay product
		
	}
	if (current.v) out.push_back(current);

	if (out.size() > 2) std::cout << "WARNING! bosonsWithDecay() found more than 2 bosons" << std::endl;

	return out;
}



//Find the jet with the minimum Radius to the vector of interest for truth matching
const xAOD::Jet* findMinRJet(const xAOD::JetContainer* jets, ROOT::Math::PtEtaPhiMVector givenVec) {

	const xAOD::Jet* jetMin = nullptr;
	ROOT::Math::PtEtaPhiMVector jetMinVec(0, 0, 0, 0);

	//looping over jets to find the one with the minimum delta R to the hadronically decaying boson system
	for (const xAOD::Jet* jet : *jets) {
		
		ROOT::Math::PtEtaPhiMVector jetVec(jet->pt(), jet->eta(), jet->phi(), jet->m());

		//if vector is empty
		if (jetMinVec.pt() == 0) {
			
			jetMinVec = jetVec;
			jetMin = jet;

		}

		//finds minimum radius
		if (ROOT::Math::VectorUtil::DeltaR(jetVec, givenVec) < ROOT::Math::VectorUtil::DeltaR(jetMinVec, givenVec)) {

			jetMinVec = jetVec;
			jetMin = jet;

		}

	}

	return jetMin;

}



//calculates d2 in the lab frame
double getD2lab(xAOD::JetConstituentVector particles, ROOT::Math::PtEtaPhiMVector jetFrame, long long p) {
	
	double e2lab = 0;
	double e3lab = 0;
	
	//declaring variables in the rest frame
	double a1lab, a2lab, a3lab;
    ROOT::Math::PtEtaPhiMVector part1, part2, part3;
	
	//loop over particle combinations in jet to sum for energy correlation functions:

	//loop over each individual particle in jet
	for (int i = 0; i < p; i++) {
		
        //get first particle vector
		ROOT::Math::PtEtaPhiMVector part1(particles[i].pt(), particles[i].eta(), particles[i].phi(), particles[i].m());
		
		//loop over every second particle which combine with the first
		for (int j = 0; j < i; j++) { 
			
            //get second particle vector
			ROOT::Math::PtEtaPhiMVector part2(particles[j].pt(), particles[j].eta(), particles[j].phi(), particles[j].m());
			
			//gets angles
			a1lab = ROOT::Math::VectorUtil::DeltaR(part1, part2);
			
			//adds to e2 energy correlation
			e2lab += part1.pt() * part2.pt() * pow(a1lab, betaVar);
			

			//loop over every third particle to combine with first two
			for (int k = 0; k < j; k++) {
				
                //get third particle vector
				ROOT::Math::PtEtaPhiMVector part3(particles[k].pt(), particles[k].eta(), particles[k].phi(), particles[k].m());
				
				//gets angles
			    a2lab = ROOT::Math::VectorUtil::DeltaR(part2, part3);
				a3lab = ROOT::Math::VectorUtil::DeltaR(part1, part3);
				
				//adds to e3 energy correlation
				e3lab += part1.pt() * part2.pt() * part3.pt() * pow(a1lab * a2lab * a3lab, betaVar); //check about e function again
				//change between pt, m, e
				
			}

		}

	}
	
	//normalize e2 and e3
	e2lab /= pow(jetFrame.pt(), 2);
	e3lab /= pow(jetFrame.pt(), 3);
	//std::cout << e2lab << "  " << e3lab << " \n";
	
	//calculate d2
	return (e3lab / pow(e2lab, 3)); 
    //d2 = e3 / e2^3
}



//calculates d2 in the jet rest frame
double getD2rest(xAOD::JetConstituentVector particles, ROOT::Math::PtEtaPhiMVector jetFrame, long long p) {
	
	double e2rest = 0;
	double e3rest = 0;
	
	//declaring variables in the rest frame
	double a1rest, a2rest, a3rest;
    ROOT::Math::PtEtaPhiMVector part1, part2, part3, part1b, part2b, part3b;
	
	ROOT::Math::Boost jetBoost(jetFrame.BoostToCM());
    ROOT::Math::PtEtaPhiMVector jetFrameb(jetBoost(jetFrame));
	
	//loop over particle combinations in jet to sum for energy correlation functions:

	//loop over each individual particle in jet
	for (int i = 0; i < p; i++) {
		
        //get first particle vector
		ROOT::Math::PtEtaPhiMVector part1(particles[i].pt(), particles[i].eta(), particles[i].phi(), particles[i].m());
		
		//changes to jet frame
		ROOT::Math::PtEtaPhiMVector part1b(jetBoost(part1));
		
		//loop over every second particle which combine with the first
		for (int j = 0; j < i; j++) { 
			
            //get second particle vector
			ROOT::Math::PtEtaPhiMVector part2(particles[j].pt(), particles[j].eta(), particles[j].phi(), particles[j].m());
			
			//changes to jet frame
			ROOT::Math::PtEtaPhiMVector part2b(jetBoost(part2));
			
			//gets angles
			a1rest = ROOT::Math::VectorUtil::Angle(part1b, part2b);
			
			//adds to e2 energy correlation
			e2rest += part1b.e() * part2b.e() * pow(a1rest, betaVar);
			
			
			//loop over every third particle to combine with first two
			for (int k = 0; k < j; k++) {
				
                //get third particle vector
				ROOT::Math::PtEtaPhiMVector part3(particles[k].pt(), particles[k].eta(), particles[k].phi(), particles[k].m());
				
				//changes to jet frame
				ROOT::Math::PtEtaPhiMVector part3b(jetBoost(part3));
				
				//gets angles
				a2rest = ROOT::Math::VectorUtil::Angle(part2b, part3b);
				a3rest = ROOT::Math::VectorUtil::Angle(part1b, part3b);
				
				//adds to e3 energy correlation
				e3rest += part1b.e() * part2b.e() * part3b.e() * pow(a1rest * a2rest * a3rest, betaVar); //check about e function again
				
			}
		}
	}
	
	//normalize e2 and e3
	e2rest /= pow(jetFrameb.e(), 2);
	e3rest /= pow(jetFrameb.e(), 3);
	
	//calculate d2
	return (e3rest / pow(e2rest, 3)); 
    //d2 = e3 / e2^3
}