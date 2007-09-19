/****************************************************************************
 * CrTrappedParticle.cxx:
 ****************************************************************************/

#define MODEL_SAA_FLUX
#undef ALLOW_SAA_SERVER


#include <cmath>
#include <string>
#include <sstream>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrTrappedParticle.hh"

#ifdef ALLOW_SAA_SERVER
#include <sys/types.h>  /* basic system data types */
#include <sys/socket.h> /* basic socket definitions */
#include <sys/time.h>   /* timeval{} for select() */
#include <time.h>       /* timespec{} for pselect() */
#include <netinet/in.h> /* sockaddr_in{} and other Internet defns */
#include <netdb.h>      /* needed by gethostbyname */
#include <arpa/inet.h>  /* needed by inet_ntoa */
#include <unistd.h>
#endif

#ifdef MODEL_SAA_FLUX
#include "astro/IGRField.h"
#include "facilities/Util.h"
#include "psb97/PSB97_model.h"
#endif
// private function definitions.

//
//
//

#define M_LAT_TOLERANCE 0.1
#define M_LON_TOLERANCE 0.1
#define M_ALT_TOLERANCE 1.0

//#######################################################################################

CrTrappedParticle::CrTrappedParticle(const std::string& paramstring="8,psb97,.",const std::string& ptype="p"): CrSpectrum()
{
   std::vector< std::string > tokens;  
   facilities::Util::stringTokenize(paramstring,",",tokens);

   m_spectrumLatitude=-1;
   m_spectrumLongitude=-1;
   m_spectrumAltitude=-1;
   m_modelMinEnergy=0.1;
   m_modelMaxEnergy=1000.;

   m_particleType=invalid;
   m_thresholdEnergy=0.;
  
   
   if(tokens.size()!=3 or facilities::Util::stringToInt(tokens[0])!=8){
        std::cerr<<"Found illegal parameter string in CrTrappedParticle: "<<paramstring<<std::endl;
        std::cerr<<"param for CrTrappedParticle must be of the form \"8,<model>,<dir/serveraddress>\"."<<std::endl;
        std::cerr<<"NO TRAPPED PARTICLE FLUX IS GENERATED."<<std::endl;
	return;
   };
  
   m_model=tokens[1];

   if(tokens[2].find(":") != std::string::npos) {
// we use an external server to retreive SAA particle fluxes.....   
      m_serverAddress=tokens[2];
      m_xmlDirectory="";
   } else {
// we use internal psb97 model to generate fluxes  
      m_serverAddress="";
      m_xmlDirectory=tokens[2];
      if(ptype.find("p")==std::string::npos || tokens[1]!="psb97"){
         std::cerr<<" Only psb97 proton flux implemented in the code. Use trapped particle data server for all other models."<<std::endl; 
         std::cerr<<"Trying to find psb97 tables and to use them."<<std::endl; 
         m_particleType = proton;
         m_model="psb97";
         m_thresholdEnergy=8.;
         m_eStep=2.;
         m_eMax=1000.;
	 return;
      };	 
   };

   if (ptype.find("e")!=std::string::npos) {
     m_particleType = electron;
     m_thresholdEnergy=0.8;
     m_eStep=0.1;
     m_eMax=20.;
   };   
   
   if (ptype.find("p")!=std::string::npos) {
     m_particleType = proton; 
     m_thresholdEnergy=8.;
     m_eStep=2.;
     m_eMax=1000.;
   };
   
   if(!checkModelCompatibility(m_model,ptype)){
        std::cerr<<m_model<<" is not a model for "<<ptype<<". Check configuration. Exit."<<std::endl;
	exit(1);
   };
         
}

//#######################################################################################


CrTrappedParticle::~CrTrappedParticle()
{}

//#######################################################################################


// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrTrappedParticle::dir(double energy, 
					      CLHEP::HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming from north
  // and phi=pi/2 for that comming from east 
{
  double costheta = 2.0*engine->flat()-1.0;
  double phi = engine->flat() * 2 * M_PI;
//  std::cout<<"Particle direction called. Returning: "<<costheta<<", "<<phi<<std::endl;
  return std::pair<double,double>(costheta, phi);

}

//#######################################################################################


// Gives back particle energy (obviously in GeV)
double CrTrappedParticle::energySrc(CLHEP::HepRandomEngine* engine) const						
{
  
//do we need a new spectrum ? if yes, request it from the server....
//...or get it from the PSB97 tables
 
  if(m_serverAddress!="") requestNewSpectrum(m_thresholdEnergy,m_eMax,m_eStep);  
  else  psb97UpdateSpectrum(m_thresholdEnergy,m_eMax,m_eStep);
  
  if(m_maxNonzeroFluxEnergy==0) return 0;

  G4double random = engine->flat();

  std::map<G4double,G4double>::const_iterator spec_it = m_intSpectrum.lower_bound(random);
  if(spec_it == m_intSpectrum.end()) return 0;

  G4double r2=(*spec_it).first;
  G4double e2=(*spec_it).second;
  if(--spec_it == m_intSpectrum.end()) return 0;
  G4double r1=(*spec_it).first;
  G4double e1=(*spec_it).second;
  
  G4double energy=e1+(random-r1)*(e2-e1)/(r2-r1);
  
//  std::cout<<"TrappedParticle "<<particleName()<<" "
//      <<"rand="<<random<<" r1="<<r1<<" r2="<<r2<<" e1="<<e1<<" e2="<<e2<<" energy="<<energy<<std::endl;
  
  return 1e-3*energy;  // in GeV
}

//#######################################################################################


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".

double CrTrappedParticle::flux() const
{
  if(m_serverAddress!="") requestNewSpectrum(m_thresholdEnergy,m_eMax,m_eStep);  
  else  psb97UpdateSpectrum(m_thresholdEnergy,m_eMax,m_eStep);
  

  G4double flux=m_integralFlux*10000./(4.*M_PI);  // [c/s/m^2/sr]
  if(m_particleType==electron) flux*=0.5; // flux from tables is electrons+positrons;
  return flux; 
}


//#######################################################################################

// Gives back solid angle from which particle comes
double CrTrappedParticle::solidAngle() const
{
   // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return  4 * M_PI ;
}

//#######################################################################################


// Gives back particle name
const char* CrTrappedParticle::particleName() const
{
  return "particle";
}

//#######################################################################################


// Gives back the name of the component
std::string CrTrappedParticle::title() const
{
  return  "CrTrappedParticle";
}


//#######################################################################################


bool CrTrappedParticle::coordinatesChanged() const
{
    bool changed=false;
    if(fabs(m_latitude-m_spectrumLatitude)>M_LAT_TOLERANCE) changed=true;
    if(fabs(m_longitude-m_spectrumLongitude)>M_LON_TOLERANCE) changed=true;
    if(fabs(m_altitude-m_spectrumAltitude)>M_ALT_TOLERANCE) changed=true;
    return changed;
}





//#######################################################################################

bool CrTrappedParticle::psb97UpdateSpectrum(G4double minE,G4double maxE,const G4double stepE) {
#ifdef MODEL_SAA_FLUX
    static TrappedParticleModels::PSB97Model psb97(m_xmlDirectory);

// values computed by askGPS implemented in CRSpectrum. We just get them now from the IGRField
    double ll = astro::IGRField::Model().L();
    double bb = astro::IGRField::Model().B();
    
    std::cout<<"update spectrum: "<<m_latitude<<","<<m_longitude;
    m_integralFlux=psb97(ll,bb,minE);
    m_intSpectrum[0.]=minE;
    std::cout<<"integral flux: "<<m_integralFlux<<std::endl;

    if(m_integralFlux>0){
       for(G4double e=minE+stepE;e<=maxE;e+=stepE) {
	  float fluxval=1. - psb97(ll,bb,e)/m_integralFlux;
	  if (fluxval<0.999999){ 
	     m_maxNonzeroFluxEnergy=e;
             m_intSpectrum[fluxval]=e;
	  };
       };        
       m_intSpectrum[1.]=m_maxNonzeroFluxEnergy+stepE;
    } else {
       m_maxNonzeroFluxEnergy=0;   
       m_integralFlux=0;
    };
#else
    m_maxNonzeroFluxEnergy=0;   
    m_integralFlux=0;
#endif    
   return true;
};

//#######################################################################################


bool CrTrappedParticle::requestNewSpectrum(G4double minE,G4double maxE,const G4double stepE) 
{
// communication with the flux server. all the annoying xml parsing is done by hand to 
// avoid requiring an extra library 
// smart way would be to use XML-RPC libraries instead at the cost of an extra dependence....  

#ifdef ALLOW_SAA_SERVER
    std::stringstream xmlRequestStr;
    G4double year= m_time/86400./365. + 2000.;
    char buffer[65536];
    std::string answer;
    std::string spectrum_string;

// do this only after a significant coordinate change. Since it takes a long time to 
// query the server
    if(!coordinatesChanged()) return true;

//  connect to server
    connectToServer();

// set the model we use for the saa flux
    std::string ptype;
    if (m_particleType==electron) ptype="e";
    if (m_particleType==proton) ptype="p";
        
    xmlRequestStr<<"<model model=\""<<m_model<<"\" year=\""<<year<<"\"/>";//<<endl;   

//    std::cout<<"Requesting model: "<<xmlRequestStr.str()<<std::endl;   

    write(m_socketHandle,xmlRequestStr.str().c_str(),xmlRequestStr.str().length());
    int bytes_read = read(m_socketHandle,buffer,65535); buffer[bytes_read]=0;
    answer=buffer;
    if(answer.find("<OK")==std::string::npos) {
       std::cout<<"CrTrappedParticle WARNING: model not set. Will produce bogus results. Answer from server was "<<answer<<std::endl;
       exit(1);
    };   

    std::string search_for="emin=\"";
    unsigned int eminpos = answer.find(search_for) + std::string(search_for).length();
    search_for="\" emax=\"";
    unsigned int emaxpos = answer.find(search_for);
    m_modelMinEnergy=atof(answer.substr(eminpos,emaxpos-eminpos).c_str());
    emaxpos += std::string(search_for).length();
    unsigned int endpos = answer.find("\"/>",emaxpos);
    m_modelMaxEnergy=atof(answer.substr(emaxpos,endpos-emaxpos).c_str());
        
    if(m_modelMinEnergy==0 || m_modelMaxEnergy==0){
       std::cout<<"CrTrappedParticle WARNING: model energy range invalid. Will produce bogus results. Answer from server was "<<answer<<std::endl;
       exit(1);
    };    	
    
    if(minE<m_modelMinEnergy) minE=m_modelMinEnergy;
    if(maxE>m_modelMaxEnergy) maxE=m_modelMaxEnergy;
            
// if model set, request spectrum for the coordinates we have
    
    xmlRequestStr.str("");   
    xmlRequestStr<<"<spectrum lat=\""<<m_latitude<<"\" lon=\""<<m_longitude
                 <<"\" alt=\""<<m_altitude
		 <<"\" energies=\""<<minE;
   
    for(G4double e=minE+stepE;e<=maxE;e+=stepE) {
       xmlRequestStr<<","<<e;
    };
    xmlRequestStr<<"\"/>\n";
    
//    std::cout<<"Requesting new spectrum: "<<xmlRequestStr.str()<<std::endl;   

    write(m_socketHandle,xmlRequestStr.str().c_str(),xmlRequestStr.str().length());
    bytes_read = read(m_socketHandle,buffer,65535); buffer[bytes_read]=0;
    answer=buffer;

//parse answer from server
    search_for="flux=\"";
    unsigned int spectrum_start_pos=answer.find(search_for)+search_for.length();
    unsigned int spectrum_end_pos=answer.find("\"/>");
    if(spectrum_start_pos==std::string::npos || spectrum_end_pos==std::string::npos) {
       std::cout<<"CrTrappedParticle WARNING: spectrum not set. Will produce bogus results. "<<std::endl;
       exit(1);
    };
    
    
   spectrum_string=answer.substr(spectrum_start_pos,spectrum_end_pos-spectrum_start_pos);
    
//    std::cout<<"Got response: "<<answer;
//    std::cout<<"Extracted spectrum: "<<spectrum_string<<std::endl;

// got the spectrum. now we read the values
// first value corresponds to total flux.
// all fluxes are stored normalized to the total flux, in an inverted map 
// optimal for energySrc to sample the spectrum from flat random numbers.... 
    m_intSpectrum.clear();
    
    std::string flux=spectrum_string.substr(0,spectrum_string.find(","));    
    m_integralFlux=atof(flux.c_str());
    m_intSpectrum[0.]=minE;
    
    if(m_integralFlux>0){
       spectrum_string=spectrum_string.substr(spectrum_string.find(",")+1,spectrum_string.length());
       for(G4double e=minE+stepE;e<=maxE;e+=stepE) {
	  flux=spectrum_string.substr(0,spectrum_string.find(","));    
	  float fluxval=1. - atof(flux.c_str())/m_integralFlux;
	  if (fluxval<1.){ 
	     m_maxNonzeroFluxEnergy=e;
             m_intSpectrum[fluxval]=e;
	  };
	  spectrum_string=spectrum_string.substr(spectrum_string.find(",")+1,spectrum_string.length());    
       };        
       m_intSpectrum[1.]=m_maxNonzeroFluxEnergy+stepE;
    } else {
       m_maxNonzeroFluxEnergy=0;   
       m_integralFlux=0;
    };
    
    // disconnect from server
    
    disconnectFromServer();
       
// set new spectrum coordinates      
    m_spectrumLatitude=m_latitude;
    m_spectrumLongitude=m_longitude;
    m_spectrumAltitude=m_altitude;
    
    return true;
#else
    return false;
#endif    
     
}

bool CrTrappedParticle::checkModelCompatibility(const std::string& model,const std::string& particle){
     if(particle.find("e")!=std::string::npos){
       if(model.find("ap8") != std::string::npos) return false;
       if(model.find("crrespro") != std::string::npos) return false;
       if(model.find("psb97") != std::string::npos) return false;
     };
     if(particle.find("p")!= std::string::npos){
       if(model.find("ae8") != std::string::npos) return false;
       if(model.find("crresle") != std::string::npos) return false;
     };
     return true;
};




//#######################################################################################

// open a socket to the server telling us the fluxes within the saa
void CrTrappedParticle::connectToServer(){
#ifdef ALLOW_SAA_SERVER
  unsigned int colon_pos=m_serverAddress.find(":",0);
  std::string server_port_str="";
  std::string server_name="";
  unsigned short server_port;
  sockaddr_in server_c_address;
  hostent *      hp;
  
  if(colon_pos==std::string::npos){
    server_name=m_serverAddress;
    server_port=3555;
  } else {;
    server_name=m_serverAddress.substr(0,colon_pos);
    server_port_str=m_serverAddress.substr(colon_pos+1,m_serverAddress.size());
    server_port=atoi(server_port_str.c_str());
  };

//  std::cout<<paramstring<<" "<<m_serverAddress<<" "<<server_name<<" "<<server_port_str<<std::endl;

  hp = gethostbyname(server_name.c_str());
  if ( !hp) {
    std::cerr<<"FATAL error encountered in resolving hostname: "<<server_name<<". Exit."<<std::endl; exit(1);
  }

  m_socketHandle = socket(AF_INET, SOCK_STREAM, 0);  
  if (m_socketHandle<0) {
    std::cerr<<"FATAL error encountered in creating socket. Exit."<<std::endl;  exit(1);
  }

  memset(&server_c_address,0,sizeof(sockaddr_in));
  server_c_address.sin_family = AF_INET;
  server_c_address.sin_port   = htons(server_port);
  memcpy(&server_c_address.sin_addr, hp->h_addr, hp->h_length);
//  std::cout<<server_c_address.sin_port <<" "<<server_c_address.sin_addr.s_addr<<endl
//      <<" "<<(int)hp->h_addr[0]<<"."<<(int)hp->h_addr[1]<<"."<<(int)hp->h_addr[2]<<"."<<(int)hp->h_addr[3]<<" "<<hp->h_length<<std::endl;
  
  if (connect(m_socketHandle, (const sockaddr*)(&server_c_address), sizeof(server_c_address)) < 0) {
     std::cout<<"FATAL error encountered in connecting to socket. Exit."<<std::endl;  exit(1);
  }
#endif
};

//#######################################################################################


void CrTrappedParticle::disconnectFromServer(){
#ifdef ALLOW_SAA_SERVER
   close(m_socketHandle) ;
#endif
};

