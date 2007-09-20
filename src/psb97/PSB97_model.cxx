#include "PSB97_model.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <xmlBase/XmlParser.h>
#include <xmlBase/Dom.h>

#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMDocument.hpp>

#include <facilities/Util.h>

#ifdef WIN32
#define isnan(x) _isnan(x)
#endif



namespace TrappedParticleModels {


    using namespace std;


    PSB97Model::PSB97Model(const string &xmldir) {
       m_Flux.push_back(new FluxTable(xmldir+"/psb97_20MeV_lbmap.xml"));
       m_Flux.push_back(new FluxTable(xmldir+"/psb97_100MeV_lbmap.xml"));
       m_Flux.push_back(new FluxTable(xmldir+"/psb97_200MeV_lbmap.xml"));
       m_Flux.push_back(new FluxTable(xmldir+"/psb97_400MeV_lbmap.xml"));           
    };
    
    
    PSB97Model::FluxTable::FluxTable(const string &xml_file){
       using namespace xmlBase;
       using namespace facilities;
       
       vector<string> tokens;        
       vector<string> data_tokens;        

       m_Header.resize(8);
       
       string xmlfile(xml_file);
       Util::expandEnvVar(&xmlfile,"$(",")");
       
       try {
	    XmlParser *parser = new XmlParser;

        std::cout<<"Parsing XML file "<<xmlfile<<std::endl; 
  
	    DOMDocument * document = parser->parse(xmlfile.c_str());    
        DOMElement* xmlRootElement = document->getDocumentElement();
	    
        if(string(XMLString::transcode(xmlRootElement->getTagName()))!="PSB97ModelData"){
 	       std::cerr << "Could not find root element PSB97ModelData " <<XMLString::transcode(xmlRootElement->getTagName())<<std::endl;
	       throw ;
	    };
	    
//	    std::cout<<"Trying to find tags in "<<xmlfile<<std::endl;
	    DOMNodeList * headerList = xmlRootElement->getElementsByTagName(XMLString::transcode("header"));
	    DOMNodeList * indexList  = xmlRootElement->getElementsByTagName(XMLString::transcode("index"));
	    DOMNodeList * dataList   = xmlRootElement->getElementsByTagName(XMLString::transcode("data"));
	     
//	    std::cout<<"Checking number of tags. "<<headerList->getLength()<<" "<<indexList->getLength()<<" "<<dataList->getLength()<<std::endl;
	    if(headerList->getLength()!=1 || indexList->getLength()!=1 || dataList->getLength()!=1){
	       std::cerr << "XML file must contain exactly one header, index and data element." << std::endl;
	       throw ;	    
	    } 
       
//	    std::cout<<"Extracting elements. "<<std::endl;
        DOMElement *header=0; 
		DOMElement *index=0;
		DOMElement *data=0;
		DOMNode* headerNode = headerList->item(0);
        DOMNode* indexNode = indexList->item(0);
        DOMNode* dataNode = dataList->item(0);
		
		if(headerNode->getNodeType()==DOMNode::ELEMENT_NODE) header = static_cast<DOMElement*>(headerNode);
		if(indexNode->getNodeType()==DOMNode::ELEMENT_NODE) index = static_cast<DOMElement*>(indexNode);
		if(dataNode->getNodeType()==DOMNode::ELEMENT_NODE) data = static_cast<DOMElement*>(dataNode);
		if(!header || !index || !data) {
	       std::cerr << "Error parsing header. Node type wrong." << std::endl;
		   throw;
		};
		unsigned int header_len = Dom::getIntAttribute(header,"len");
	    unsigned int index_len  = Dom::getIntAttribute(index,"len");
	    unsigned int data_len   = Dom::getIntAttribute(data,"len");
	    
//	    std::cout<<"Extracting element content. "<<std::endl;
	    string header_string= Dom::getTextContent(header);
	    string index_string = Dom::getTextContent(index);
	    string data_string  = Dom::getTextContent(data);
	    
        stringstream header_stream(header_string);

//	    std::cout<<"Filling header. "<<std::endl;

	    vector<float>::iterator h_it=m_Header.begin();
	    m_Header.resize(header_len);
	    for (unsigned int i=0;i<header_len;i++,h_it++) {
	      if (header_stream.eof()) throw;
	      header_stream>> *h_it;
//	      std::cout<< *h_it <<endl ;	      
	    };
	    
        stringstream index_stream(index_string);
	    stringstream data_stream(data_string);	        
	    
//	    std::cout<<"Filling data. "<<std::endl;

	    m_Data.resize(index_len);
	    vector< vector<float> >::iterator  d_it=m_Data.begin();

		unsigned int data_size=0;
        unsigned int index1,index2;

        index_stream>>index1;
	    for (unsigned int i=0;i<index_len; i++,d_it++)  {
	       vector<float>& data_vector = *d_it;
	       if (index_stream.eof()) throw;
	       if (data_stream.eof()) throw;
	       
	       index_stream>>index2;
	       unsigned int length = index2 - index1;
	       
	       data_vector.resize(length);
	       vector<float>::iterator dv_it = data_vector.begin();

//	       std::cout<<"  --> adding vector with length l="<<length<<" "<<index1<<" "<<index2<<std::endl; 

	       for (unsigned int j=0; j<length; j++,dv_it++) {
              data_stream >> *dv_it ; 	
		      data_size++;
	        };
            index1=index2; 	   
	    };    
	    
        if(data_size!=data_len) {
          std::cerr<<"Data length wrong in XML file. "<<data_size<<" <--> "<<data_len<<std::endl;
          throw;
        };
	    
	    tokens.clear(); tokens.resize(1);
	    data_tokens.clear(); data_tokens.resize(1);

        std::cout<<"Done parsing "<<xmlfile<<std::endl; 
	    
        delete parser;
       
      } catch(...) {     
        std::cerr << "An error occurred during parsing psb97 model data xml file "<<xmlfile<<std::endl;
        exit(1);
      };
 
   };   
 
 

    PSB97Model::FluxTable::~FluxTable(){};
    
    
    float PSB97Model::operator()(float ll,float bb,float ee) const {

// using namespace std;

       vector<FluxTable*>::const_iterator ft_it = m_Flux.begin();
       FluxTable* fluxTabLow;
       FluxTable* fluxTabHigh;	
	
        while (ft_it!=m_Flux.end()) {
	  if (ee<=(*ft_it)->Energy()) break;
          ft_it++;
        };
        if(ft_it==m_Flux.begin()) ++ft_it;
        if(ft_it==m_Flux.end()) --ft_it;
        
        fluxTabHigh= *ft_it;
        fluxTabLow=  *(--ft_it);

//       cout<<"ee="<<ee<<", fluxtab_elow="<<fluxTabLow->Energy()<<", fluxtab_ehigh="<<fluxTabHigh->Energy()<<endl;
         	
//   	cout<<"ll,bb="<<ll<<", "<<bb<<endl;
	
	if( (!fluxTabLow->isInLRange(ll)) || (!fluxTabLow->isInBRange(bb)) ) return 0.;

	unsigned int lindex= fluxTabLow->Lindex(ll);
	unsigned int bindex= fluxTabLow->Bindex(bb);
        float dll= fluxTabLow->DeltaL(ll);
	float dbb= fluxTabLow->DeltaB(bb);

//	cout<<"lindex,bindex="<<lindex<<", "<<bindex
//	    <<" dll,dbb="<<dll<<", "<<dbb<<endl;
//	cout<<"lsizes="<<lsizeLow<<", "<<lsizeLow1<<", "<<lsizeHigh<<", "<<lsizeHigh1<<endl;
	
// find flux values and interpolate linear between the lls and bbs	
	float flux11=  fluxTabLow->Flux(lindex,bindex);
	float flux21=  fluxTabLow->Flux(lindex+1,bindex);
	float flux12=  fluxTabLow->Flux(lindex,bindex+1);
	float flux22=  fluxTabLow->Flux(lindex+1,bindex+1);

    if((flux11+flux21+flux12+flux22)<=1.e-5) return 0.;
	
    float log_flux_low = log(linear_interpolation(dll,dbb,flux11,flux21,flux12,flux22));
	
	
//        cout << "fluxes_low="<<flux11<<", "<<flux21<<", "<<flux12<<", "<<flux22
//	     << " log_flux_low="<<log_flux_low<<endl;

	flux11=  fluxTabHigh->Flux(lindex,bindex);
	flux21=  fluxTabHigh->Flux(lindex+1,bindex);
	flux12=  fluxTabHigh->Flux(lindex,bindex+1);
	flux22=  fluxTabHigh->Flux(lindex+1,bindex+1);

    if((flux11+flux21+flux12+flux22)<=1.e-5) return 0.;
        
	float log_flux_high = log(linear_interpolation(dll,dbb,flux11,flux21,flux12,flux22));
 
    if(log_flux_high>=log_flux_low) return 0.;
 
//         cout << "fluxes_high="<<flux11<<", "<<flux21<<", "<<flux12<<", "<<flux22
//	     << " log_flux_high="<<log_flux_high<<endl;
        
    float flux = exp( log_flux_low+ (ee-fluxTabLow->Energy())/(fluxTabHigh->Energy()-fluxTabLow->Energy())*(log_flux_high-log_flux_low) );
	if (flux<0.01 || isnan(flux)) return 0.;

//	cout<<"flux="<<flux<<endl;

	return flux;
	
    };


    float PSB97Model::linear_interpolation(float dx,float dy,float v11, float v21, float v12, float v22) const {

      float value = (1-dx)*(1-dy)*v11
        	  + (1-dx)*dy*v12
        	  + dx*(1-dy)*v21
        	  + dx*dy*v22;
//      cout<<"interpolation value="<<value<<endl;
      return value;

    };
};


