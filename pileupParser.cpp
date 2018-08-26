#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>

using namespace std;

class site
{
public:
	int covA=0;
	int covT=0;
	int covG=0;
	int covC=0;

	int polymorphicbool=0;
	int indelbool=0;
//	int covdel=0;
};

int main(int argc, char *argv[])
{
	string argv1=argv[1];
	//string argv2=argv[2];

	//----------the major construct------------
	map<string, map<int,site>> sites; //here, the string includes both contig and position
	//-----------------------------------------

	ifstream pileup(argv1.c_str());

	int sample_cnt=0;

		string line;

		
		for(;getline(pileup,line);) //each site
		{
			istringstream linestream(line);
			
			string chr;
			string pos;
			char ref;
			int cov;
			string bases;
			
			linestream >> chr >> pos >> ref >> cov >> bases;

			//=======screen for coverage?=======
			//if(cov<2) continue;
			//==================================

			if(sites.find(chr)==sites.end()) 
			{
				map<int,site> tmp;
				sites[chr]=tmp;
			}
			
			//compute the abundance of bases at the site for the ith sample
			//if(sites[chr_pos].is_indel==1) continue; //the site includes indels at some samples, remove

			int count[6]={0,0,0,0,0,0}; //0-ref, 1-A, 2-T, 3-G, 4-C, 5-"del"
			int indel_bool=0;
			for(int si=0;si<bases.size();++si)
			{
				if(bases[si]=='A'||bases[si]=='a') count[1]++;
				else if(bases[si]=='T'||bases[si]=='t') count[2]++;
				else if(bases[si]=='G'||bases[si]=='g') count[3]++;
				else if(bases[si]=='C'||bases[si]=='c') count[4]++;
				else if(bases[si]=='.'||bases[si]==',') count[0]++;
				else if(bases[si]=='^') {++si;continue;}
				else if(bases[si]=='$') continue;
				else if(bases[si]=='+' || bases[si]=='-')
				{
					indel_bool=1;
					++si;
					int skipcnt=bases[si]-'0';
					for(int cnt=0; cnt<skipcnt;cnt++)
						++si;
				}
			}
			
			
				if(ref=='A' || ref=='a') count[1]+=count[0];
				else if(ref=='T' || ref=='t') count[2]+=count[0];
				else if(ref=='G' || ref=='g') count[3]+=count[0];
				else if(ref=='C' || ref=='c') count[4]+=count[0];

				int nonzero=0;

				if(count[1]!=0) nonzero++;
				if(count[2]!=0) nonzero++;
				if(count[3]!=0) nonzero++;
				if(count[4]!=0) nonzero++;
				
				sites[chr][stoi(pos)].covA=/*(int)((double)*/count[1]/**100/(double)cov)*/;
				sites[chr][stoi(pos)].covT=/*(int)((double)*/count[2]/**100/(double)cov)*/;
				sites[chr][stoi(pos)].covG=/*(int)((double)*/count[3]/**100/(double)cov)*/;
				sites[chr][stoi(pos)].covC=/*(int)((double)*/count[4]/**100/(double)cov)*/;

				sites[chr][stoi(pos)].indelbool=indel_bool;
				if(nonzero>=2) sites[chr][stoi(pos)].polymorphicbool=1;
				//cout<<nonzero<<"\t"<<sites[chr][stoi(pos)].covA<<"\t"<<sites[chr][stoi(pos)].covT<<"\t"<<sites[chr][stoi(pos)].covG<<"\t"<<sites[chr][stoi(pos)].covC<<"\t"<<sites[chr][stoi(pos)].indelbool<<sites[chr][stoi(pos)].polymorphicbool<<endl;
				
				
		}


	//polymorphic?biallelic?
	ofstream output("polymorphic.site.coverage");
	output<<"contig\tposition\tAcoverage\tTcoverage\tGcoverage\tCcoverage"<<endl;
	for(auto si : sites)
	{
		string ctg=si.first;
		for(auto sii : si.second)
		{
			//if(sii.second.indelbool==1 && argv2=="1") continue;
			//if(sii.second.polymorphicbool!=1) continue;
			output<<si.first<<"\t"<<sii.first<<"\t"<<sii.second.covA<<"\t"<<sii.second.covT<<"\t"<<sii.second.covG<<"\t"<<sii.second.covC<<"\t"<<endl;
		}
	}	
		
}

