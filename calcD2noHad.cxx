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
#include "TGraph.h"
#include <TGraph2D.h>

//ATLAS
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include <xAODJet/Jet.h>
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticle.h"
#include <xAODJet/JetConstituentVector.h>

//fastjet
#include "fastjet/ClusterSequence.hh"

/***
 * Needs to be updated:
 * 
Program takes in a file specified in first command line argument
Program takes in an output ID in second command line argument
Program calculate the energy correlation functions for different jets in the file:
e2 = (1/M^2) * sum{ (Ei * Ej * angle ij) }
e3 = (1/M^3) * sum{ (Ei * Ej * Ek * angle ij * angle ik * angle jk) }
Program calculates D2 discriminant for each jet D2 = e3 / (e2^3) in the rest frame of the jet, and in the lab frame
Program stores the jet 4-vector, D2's in different frames, decay information, and substructure information to a tree
Program write the tree to a file [output ID]output.root

Version 2.4
-Clusters into subjets with antikt and R=0.4 in jet rest frame
-Calculates the percentage of energy contained in the leading two subjets

Shand Seiffert 6/23/2023
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



const double betaVar = 1.0;
//power dependence for angles in e2, e3 formulas

const double betaVarRest = 2.0;
//power dependence for angles in e2, e3 rest formulas

double getD2rest(xAOD::JetConstituentVector, ROOT::Math::PtEtaPhiMVector, long long, double);
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

const double deltaRcut = 0.4;
//radius cutoff for truth matching jets

const xAOD::Jet* findMinRJet(const xAOD::JetContainer*, ROOT::Math::PtEtaPhiMVector);
//Find the jet with the minimum Radius to the vector of interest for truth matching
//Input: the container with the jets from an event, the vector to compare them to (boson vector)
//Output: the jet with the minimum radius to the vector

const double subjetRadius = 0.6;
//radius used for reconstructing subjets

bool subjetOutput = false;
//boolean for whether the clustering function prints subjet info to the terminal as it processes

const std::vector<fastjet::PseudoJet> restClustering(const xAOD::Jet*, double, ROOT::Math::PtEtaPhiMVector);
//clusters them using antikt to get subjet info in the rest frame
//Input: The jet, R for clustering, the vector of the frame in which to cluster
//Output: a vector of subjets

const double extERadiusCutoff = 0.4;
//the radius cutoff for the external energy calculation below

double externalEnergy(xAOD::JetConstituentVector, std::vector<fastjet::PseudoJet>, ROOT::Math::PtEtaPhiMVector, double);
//boosts subjet axes into jet frame, and finds how much energy of the jet falls outside a certain radius around the subjet axes
//input: the jet constituents, the vector of 2 subjets, the number of constituents.
//output: a double representing the energy percentage.

double subjetAngles(std::vector<fastjet::PseudoJet>, ROOT::Math::PtEtaPhiMVector);
//takes a vector of two subjets, and boosts into the boson frame and check their angle 
//(if these are indeed the two subjets desired, will be roughly opposite)
//Input: Vector with two subjets, vector of the boson frame
//Output: The angle between these two subjet axes




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
	double leadSubjetPt, numSubjets, subjetAngle, extEnergy1, extEnergy2, extEnergy3;
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
	t->Branch("leadSubjetPt", &leadSubjetPt);
	t->Branch("numSubjets", &numSubjets);
	t->Branch("subjetAngle", &subjetAngle);
	t->Branch("d2rest1_3", &extEnergy1);
	t->Branch("d2rest1_7", &extEnergy2);
	t->Branch("d2rest2_0", &extEnergy3);

	int eventCounter = 0;
	
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

		float d2 = 0.0;

		for (const xAOD::Jet* jet : *jets) {
			
			ROOT::Math::PtEtaPhiMVector jetMinVec(jet->pt(), jet->eta(), jet->phi(), jet->m());

			if (jetMinVec.M() < 100000 && jetMinVec.pt() > 400000) {
				
				//writing to tree variables
				ptj = jetMinVec.pt();
				phij = jetMinVec.phi();
				etaj = jetMinVec.eta();
				mj = jetMinVec.M();
				mj /= 1000;
				weightj = weights[0];
				
				//getting d2 in both frames, along with the given d2
				D2rest = getD2rest(jet->getConstituents(), jetMinVec, jet->numConstituents(), 1.0);
				D2lab = getD2lab(jet->getConstituents(), jetMinVec, jet->numConstituents());

				//Perform clustering and store subjets
				std::vector<fastjet::PseudoJet> subjets = restClustering(jet, subjetRadius, jetMinVec);

				//storing variables for tree analysis
				numSubjets = subjets.size();
				leadSubjetPt = subjets[0].perp();
				subjetAngle = subjetAngles(subjets, jetMinVec);
				extEnergy1 = getD2rest(jet->getConstituents(), jetMinVec, jet->numConstituents(), 1.3);
				extEnergy2 = getD2rest(jet->getConstituents(), jetMinVec, jet->numConstituents(), 1.7);
				extEnergy3 = getD2rest(jet->getConstituents(), jetMinVec, jet->numConstituents(), 2.0);

				//std::cout << extEnergy1 << std::endl;

				/******/ /*
				if (eventCounter == 11) {

					xAOD::JetConstituentVector particles = jet->getConstituents();
					int numParts = jet->numConstituents();
					
					ROOT::Math::Boost jetBoost(jetMinVec.BoostToCM());

					double eta[numParts], phi[numParts];
					for (int i = 0; i < numParts; i++) {

						ROOT::Math::PtEtaPhiMVector particle(particles[i].pt(), particles[i].eta(), particles[i].phi(), particles[i].m());
						ROOT::Math::PtEtaPhiMVector particleb(jetBoost(particle));

						eta[i] = particleb.eta();
						phi[i] = particleb.phi();

					}

					TGraph * eventPlot = new TGraph(numParts, eta, phi);
					eventPlot->Write();

				}

				eventCounter++;
				*/ /******/
				
				//gets D2 given in file, (only applies if taking from some branches)
				if (jet->getAttribute< float >("D2", d2)) d2refj = d2;
				else d2refj = -100;
				t->Fill();
				
			}

		}

	}

	//writes and closes the file
    t->Write();
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
double getD2rest(xAOD::JetConstituentVector particles, ROOT::Math::PtEtaPhiMVector jetFrame, long long p, double betaGiven) {
	
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
			a1rest = (1 / ROOT::Math::VectorUtil::Angle(part1b, part2b));
			
			//adds to e2 energy correlation
			e2rest += part1b.e() * part2b.e() * pow(a1rest, betaGiven);
			
			
			//loop over every third particle to combine with first two
			for (int k = 0; k < j; k++) {
				
                //get third particle vector
				ROOT::Math::PtEtaPhiMVector part3(particles[k].pt(), particles[k].eta(), particles[k].phi(), particles[k].m());
				
				//changes to jet frame
				ROOT::Math::PtEtaPhiMVector part3b(jetBoost(part3));
				
				//gets angles
				a2rest = (1 / ROOT::Math::VectorUtil::Angle(part2b, part3b));
				a3rest = (1 / ROOT::Math::VectorUtil::Angle(part1b, part3b));
				
				//adds to e3 energy correlation
				e3rest += part1b.e() * part2b.e() * part3b.e() * pow(a1rest * a2rest * a3rest, betaGiven); //check about e function again
				
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



//clusters to find subjets in rest frame of jet
const std::vector<fastjet::PseudoJet> restClustering(const xAOD::Jet* jet, double subjetAngle, ROOT::Math::PtEtaPhiMVector restVec) {

	//declare constituent vector to loop, and boost
	std::vector<fastjet::PseudoJet> constituents;
	ROOT::Math::Boost restBoost(restVec.BoostToCM());
	
	//looping over constituents, boosting to jet frame and adding them for clustering
	for (const auto* constituent : jet->getConstituents()) {
		
		ROOT::Math::PtEtaPhiMVector particle(constituent->pt(), constituent->eta(), constituent->phi(), constituent->m());
		ROOT::Math::PtEtaPhiMVector particleBoosted(restBoost(particle));
		constituents.push_back(fastjet::PseudoJet( particleBoosted.px(), particleBoosted.py(), particleBoosted.pz(), particleBoosted.E()));
		
	}

	//Perform clustering and store subjets
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, subjetAngle);
	fastjet::ClusterSequence cs(constituents, jet_def);
	std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs.inclusive_jets());

	//if subjet output is toggled on, it will read out subjet information
	if (subjetOutput) {

		std::cout << "     pt y phi" << std::endl;
		for (unsigned i = 0; i < subjets.size(); i++) {
			
			std::cout << "jet " << i << ": "<< subjets[i].perp() << " "
				<< subjets[i].rap() << " " << subjets[i].phi() << std::endl;
			
			std::vector<fastjet::PseudoJet> constituents = subjets[i].constituents();
			
			for (unsigned j = 0; j < constituents.size(); j++) {
				
				std::cout << " constituent " << j << "’s pt: "<< constituents[j].perp() << std::endl;

			}

		}

	}
	
	return subjets;

}



//takes a vector of two subjets, and boosts into the boson frame and check their angle
double subjetAngles(std::vector<fastjet::PseudoJet> subjets, ROOT::Math::PtEtaPhiMVector bosonFrame) {

	//for boosting into boson frame
	ROOT::Math::Boost bosonBoost(bosonFrame.BoostToCM());

	if (subjets.size() < 2) return -100;

	//boosting subjets
	ROOT::Math::PtEtaPhiMVector subjet1(subjets[0].perp(), subjets[0].eta(), subjets[0].phi(), subjets[0].m());
	ROOT::Math::PtEtaPhiMVector subjet1b(bosonBoost(subjet1));

	ROOT::Math::PtEtaPhiMVector subjet2(subjets[1].perp(), subjets[1].eta(), subjets[1].phi(), subjets[1].m());
	ROOT::Math::PtEtaPhiMVector subjet2b(bosonBoost(subjet2));

	//find and return angle
	return ROOT::Math::VectorUtil::Angle(subjet1, subjet2);

}



//boosts subjet axes into jet frame, and finds how much energy of the jet falls outside a certain radius around the subjet axes
double externalEnergy(xAOD::JetConstituentVector constituents, std::vector<fastjet::PseudoJet> subjets, ROOT::Math::PtEtaPhiMVector jetFrame, double angleCuttoff) {

	long long numConstituents = constituents.size();
	ROOT::Math::Boost jetBoost(jetFrame.BoostToCM());

	//energy totals
	double energyInt = 0;
	double energyIntsub1 = 0;
	double energyIntsub2 = 0;
	double energyTotal = 0;

	if (subjets.size() < 2) return -100;

	//subjet vectors
	ROOT::Math::PtEtaPhiMVector subjet1(subjets[0].perp(), subjets[0].eta(), subjets[0].phi(), subjets[0].m());
	ROOT::Math::PtEtaPhiMVector subjet2(subjets[1].perp(), subjets[1].eta(), subjets[1].phi(), subjets[1].m());

	for (int i = 0; i < numConstituents; i++) {

		ROOT::Math::PtEtaPhiMVector constituentVec(constituents[i].pt(), constituents[i].eta(), constituents[i].phi(), constituents[i].m());
		ROOT::Math::PtEtaPhiMVector constituentVecb(jetBoost(constituentVec));
		energyTotal += constituentVecb.E();
		double angle1 = ROOT::Math::VectorUtil::Angle(subjet1, constituentVecb);
		double angle2 = ROOT::Math::VectorUtil::Angle(subjet2, constituentVecb);

		if (angle1 < angleCuttoff) energyIntsub1 += constituentVecb.E();
		if (angle2 < angleCuttoff) energyIntsub2 += constituentVecb.E();

		if (angle1 < angleCuttoff || angle2 < angleCuttoff) energyInt += constituentVecb.E();

	}

	return 1 - (energyIntsub1 / energyTotal) - (energyIntsub2 / energyTotal);
	

}
