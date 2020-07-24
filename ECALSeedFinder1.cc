// -*- C++ -*-
//
// Package:    ECALSeedFinder1
// Class:      ECALSeedFinder1
// 
/**\class ECALSeedFinder1 ECALSeedFinder1.cc RecoElectron/ECALSeedFinder1/src/ECALSeedFinder1.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  James Pivarski
//         Created:  Thu Mar 23 17:28:30 CST 2006
// $Id$
//
//


// system include files
#include <memory>
#include <cassert>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EGammaReco/interface/BasicCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DMatchedLocalPosCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DLocalPosCollection.h"
#include "RecoTracker/RoadMapRecord/interface/Roads.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

//
// class decleration
//

class ECALSeedFinder1 : public edm::EDProducer {
   public:
      explicit ECALSeedFinder1(const edm::ParameterSet&);
      ~ECALSeedFinder1();


      virtual void produce(edm::Event&, const edm::EventSetup&);
   private:
      void findHits(double hypothesisPhi, const Ring ring, std::vector<double>& hits);
      void findSeeds(edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::BasicCluster& cluster, double hypothesisCharge);
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
ECALSeedFinder1::ECALSeedFinder1(const edm::ParameterSet& iConfig)
{
   //register your products
  produces<double>("biggestClusterEnergy");

   //now do what ever other initialization is needed

}


ECALSeedFinder1::~ECALSeedFinder1()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
ECALSeedFinder1::findHits(double hypothesisPhi, const Ring ring, std::vector<double>& hits) {
   using namespace edm;
   using namespace std;

   const double phiRange = 2.*M_PI/20.;  // consider detector components within 1/20th of phi from our hypothesis track

   // make sure that our phi hypothesis is in the right range (0..2pi)
   while (hypothesisPhi < 0.) hypothesisPhi += 2.*M_PI;
   while (hypothesisPhi >= 2.*M_PI) hypothesisPhi -= 2.*M_PI;

   Ring::const_iterator firstAcceptable = ring.lower_bound(hypothesisPhi - phiRange);
   Ring::const_iterator lastAcceptable =  ring.upper_bound(hypothesisPhi + phiRange);
   ++lastAcceptable; // we want to include the last acceptable ring in the upcoming loop

   for (Ring::const_iterator acceptable = firstAcceptable;  acceptable != lastAcceptable;  ++acceptable) {
     // cout << "look in detId " << acceptable->second.rawId() << " (which is at phi " << acceptable->first << ")" << endl;
     hits.push_back(acceptable->first);
   }

   if (hypothesisPhi - phiRange < 0.) {  // go get the ones below zero (loop around)
     firstAcceptable = ring.lower_bound(hypothesisPhi - phiRange + 2.*M_PI);
     lastAcceptable = ring.upper_bound(2.*M_PI);
     ++lastAcceptable; // for inclusion (see above)
     for (Ring::const_iterator acceptable = firstAcceptable;  acceptable != lastAcceptable;  ++acceptable) {
       // cout << "look in detId " << acceptable->second.rawId() << " (which is at phi " << acceptable->first << ")" << endl;
       hits.push_back(acceptable->first);
     }
   }

   if (hypothesisPhi + phiRange > 2.*M_PI) {  // go get the ones above 2pi (loop around)
     firstAcceptable = ring.lower_bound(0.);
     lastAcceptable = ring.upper_bound(hypothesisPhi + phiRange - 2.*M_PI);
     ++lastAcceptable; // for inclusion (see above)
     for (Ring::const_iterator acceptable = firstAcceptable;  acceptable != lastAcceptable;  ++acceptable) {
       // cout << "look in detId " << acceptable->second.rawId() << " (which is at phi " << acceptable->first << ")" << endl;
       hits.push_back(acceptable->first);
     }
   }

   // modified "hits" is the output of this function
}

void
ECALSeedFinder1::findSeeds(edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::BasicCluster& cluster, double hypothesisCharge) {
   using namespace edm;
   using namespace std;

   ESHandle<Roads> roads;
   iSetup.get<TrackerDigiGeometryRecord>().get(roads);

   // correction for B given in T, distances in cm, momenta in GeV/c
   const double speedOfLight = 2.99792458e8;
   const double unitCorrection = speedOfLight * 1e-2 * 1e-9;
   // B in T, right now hardcoded, has to come from magnetic field service
   double B = 4.0;
   
   assert(cluster.position().z() != 0.);      // set of measure zero?
   assert(cluster.position().perp2() != 0.);  // this would be strange

   // calculate track parameters from the biggest-energy cluster
   const double tanDip = cluster.position().z() / sqrt(cluster.position().perp2());
   const double innerZposRange = 20.; // centimeters (how far we allow hypothesis track to miss the detector element)
   const double outerZposRange = 20.; // centimeters
   const double innerRadiusRange = 20.; // centimeters
   const double outerRadiusRange = 20.; // centimeters
   const double trackRadius = cluster.energy() / unitCorrection / B / sqrt(1. + tanDip*tanDip);
   const double phiAtEcal = atan2(cluster.position().y(), cluster.position().x());

   // loop over seed Ring pairs
   for (Roads::const_iterator road = roads->begin();  road != roads->end();  ++road) {
     const Roads::RoadSeed seed = road->first;
     const Ring innerRing = seed.first;
     const Ring outerRing = seed.second;

     vector<double> innerHits;  // "double" is a place-holder for hits
     vector<double> outerHits;

     // collect plausible hits from the inner ring
     if (innerRing.getType() == Ring::TIBRing  ||  innerRing.getType() == Ring::TOBRing) { // detector element has a definite radius, extended in z position
       const double ringRadius = (innerRing.getrmin() + innerRing.getrmax())/2.;

       const double hypothesisZpos = tanDip * ringRadius;
       if ((innerRing.getzmin() < hypothesisZpos  &&  hypothesisZpos < innerRing.getzmax())  ||
	   innerRing.getzmin() - innerZposRange < hypothesisZpos                             ||
	   innerRing.getzmax() + innerZposRange > hypothesisZpos                                ) {

	 const double phiCorrection = atan2(ringRadius/2., trackRadius) * hypothesisCharge;  // correction to phi position due to B field
	 findHits(phiAtEcal + phiCorrection, innerRing, innerHits);
       }
     }
     else if (innerRing.getType() == Ring::TIDRing  ||  innerRing.getType() == Ring::TECRing) { // detector element has a definite z position, extended in radius
       const double ringZpos = (innerRing.getzmin() + innerRing.getzmax())/2.;

       const double hypothesisRadius = fabs(ringZpos) / tanDip;
       if ((innerRing.getrmin() < hypothesisRadius  &&  hypothesisRadius < innerRing.getrmax())  ||
	   innerRing.getrmin() - innerRadiusRange < hypothesisRadius                             ||
	   innerRing.getrmax() + innerRadiusRange > hypothesisRadius                                ) {
	 
	 const double phiCorrection = atan2(hypothesisRadius, trackRadius) * hypothesisCharge;  // correction to phi position due to B field
	 findHits(phiAtEcal + phiCorrection, innerRing, innerHits);
       }
     }

     // collect plausible hits from the outer ring
     if (outerRing.getType() == Ring::TOBRing  ||  outerRing.getType() == Ring::TIBRing) { // detector element has a definite radius, extended in z position
       const double ringRadius = (outerRing.getrmin() + outerRing.getrmax())/2.;

       const double hypothesisZpos = tanDip * ringRadius;
       if ((outerRing.getzmin() < hypothesisZpos  &&  hypothesisZpos < outerRing.getzmax())  ||
	   outerRing.getzmin() - outerZposRange < hypothesisZpos                             ||
	   outerRing.getzmax() + outerZposRange > hypothesisZpos                                ) {

	 const double phiCorrection = atan2(ringRadius/2., trackRadius) * hypothesisCharge;  // correction to phi position due to B field
	 findHits(phiAtEcal + phiCorrection, outerRing, outerHits);
       }
     }
     else if (outerRing.getType() == Ring::TECRing  ||  outerRing.getType() == Ring::TIDRing) { // detector element has a definite z position, extended in radius
       const double ringZpos = (outerRing.getzmin() + outerRing.getzmax())/2.;

       const double hypothesisRadius = fabs(ringZpos) / tanDip;
       if ((outerRing.getrmin() < hypothesisRadius  &&  hypothesisRadius < outerRing.getrmax())  ||
	   outerRing.getrmin() - outerRadiusRange < hypothesisRadius                             ||
	   outerRing.getrmax() + outerRadiusRange > hypothesisRadius                                ) {
	 
	 const double phiCorrection = atan2(hypothesisRadius, trackRadius) * hypothesisCharge;  // correction to phi position due to B field
	 findHits(phiAtEcal + phiCorrection, outerRing, outerHits);
       }
     }

     // loop over all inner/outer hit combinations
     for (vector<double>::const_iterator innerHit = innerHits.begin();  innerHit != innerHits.end();  ++innerHit) {
       for (vector<double>::const_iterator outerHit = outerHits.begin();  outerHit != outerHits.end();  ++outerHit) {

	 // calculate the projection of origin + innerHit + outerHit: how does it compare with the shower position/energy?

       }
     }

   } // end loop over roads
}

// ------------ method called to produce the data  ------------
void
ECALSeedFinder1::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::BasicClusterCollection> clusters;
   iEvent.getByType(clusters);

   // for now, we'll just apply this algorithm to the biggest-energy cluster
   double biggestEnergy = 0.;
   for (reco::BasicClusterCollection::const_iterator cluster = clusters->begin(); cluster != clusters->end(); ++cluster) {
     cout << "cluster energy = " << cluster->energy() << endl;
     if (cluster->energy() > biggestEnergy) {
       biggestEnergy = cluster->energy();
     }

     // this does all of the work
     findSeeds(iEvent, iSetup, *cluster,  1.);
     findSeeds(iEvent, iSetup, *cluster, -1.);
   }

   auto_ptr<double> toPut(new double(biggestEnergy));
   iEvent.put(toPut, "biggestClusterEnergy");
}

//define this as a plug-in
DEFINE_FWK_MODULE(ECALSeedFinder1)
