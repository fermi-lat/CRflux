#ifndef PSB97_MODEL_H
#define PSB97_MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

namespace TrappedParticleModels {


    class PSB97Model {
      
      struct FluxTable{
         FluxTable(const std::string &xmlfile);
	 ~FluxTable();
	 
	 bool isInLRange(float ll){ return (ll>=m_Header[1] && ll<=m_Header[2] );}; 
	 bool isInBRange(float bb){ return (bb>=m_Header[4] && bb<=m_Header[5] );}; 
	 unsigned int Lindex(float ll){ return (unsigned int) ((ll-m_Header[1])/m_Header[3]);};
	 unsigned int Bindex(float bb){ return (unsigned int) ((bb-m_Header[4])/m_Header[6]);};
	 float DeltaL(float ll){ return (ll-(m_Header[1]+Lindex(ll)*m_Header[3]))/m_Header[3];};
	 float DeltaB(float bb){ return (bb-(m_Header[4]+Bindex(bb)*m_Header[6]))/m_Header[6];};

	 bool indexExists(unsigned int idxL,unsigned int idxB){ return idxB < m_Data[idxL].size() ;}; 
	 float Flux(unsigned int idxL,unsigned int idxB){ return indexExists(idxL,idxB) ? (m_Data[idxL])[idxB] : 0. ;}; 
         float Energy(){ return m_Header[0];};
	 
         std::vector<float> m_Header;
         std::vector< std::vector<float> > m_Data;

      };	 
	
      std::vector<FluxTable*> m_Flux;	
	 
      public:
         PSB97Model(const std::string& xmldir=".");
	 ~PSB97Model(){};
         	  
         float operator()(float ll,float bb, float ee) const;
      
      private:	 
         float linear_interpolation (float dx,float dy,float v11, float v21, float v12, float v22) const;
    };



};

#endif

