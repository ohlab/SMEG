#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <thread>

using namespace::std;
int sitecnt=0;

class SNP //the SNP at a position
{
public:
	map<char, double> nucleotidecount; // 'A' - counts in the corresponding cluster
	map<char, int> nucleotidecluster; // 'A' - -1 if not observed, -2 if present in more than one clusters, or the cluster number if observed only in that cluster
	int isunique=0;

	SNP()
	{
		nucleotidecount['A']=0;
		nucleotidecount['T']=0;
		nucleotidecount['G']=0;
		nucleotidecount['C']=0;

		nucleotidecluster['A']=-1;
		nucleotidecluster['T']=-1;
		nucleotidecluster['G']=-1;
		nucleotidecluster['C']=-1;
	}
};		

void checkSite(int seqi, vector<SNP> &snps, map<string,string> &sequences, map<string, int> &straincluster)
{
	char ref;
	int isFirst=1;
	int poly=0;
	int gap=0;

	for(auto sequence : sequences)
	{
		if(sequence.second[seqi]=='-') //gap is a no-go anytime
		{
			gap=1;
			break;
		}

		if(isFirst)
		{
			ref=sequence.second[seqi];
			isFirst=0;
			snps[seqi].nucleotidecluster[ref]=straincluster[sequence.first];
			snps[seqi].nucleotidecount[ref]++;
			snps[seqi].isunique++;
		}
		else 
		{
			if(sequence.second[seqi]!=ref) //site being polymorphic
				poly=1;

			char nt=sequence.second[seqi]; //whether or not polymorphic, check if unique

			if(snps[seqi].nucleotidecluster[nt]!=-2) //still ok
			{
				int sc=straincluster[sequence.first];
				if(snps[seqi].nucleotidecluster[nt]==-1) //first occurence of the base
				{
					snps[seqi].nucleotidecluster[nt]=sc;
					snps[seqi].nucleotidecount[nt]++;
					snps[seqi].isunique++;
				}
				else if(snps[seqi].nucleotidecluster[nt]==sc) //matches previous occurences
				{
					snps[seqi].nucleotidecount[nt]++;
				}
				else if(snps[seqi].nucleotidecluster[nt]!=sc)
				{
					snps[seqi].nucleotidecluster[nt]=-2;
					snps[seqi].isunique--; //this won't decrease unique count multiple times since it checks for -2 first
					continue;
				}
			}
		}

	}

	if(gap==1 && poly==0)
	{
		snps[seqi].isunique=0;
	}
	if(sitecnt%5000==0) 
	{
		cout<<sitecnt<<" sites screened."<<endl;
		cout.flush();
	}
	sitecnt++;
}

void checkSites(int starti, int endi, vector<SNP> &snps, map<string,string> &sequences,map<string, int> &straincluster)
{
	for(int seqi=starti; seqi<=endi; seqi++)
	{
		checkSite(seqi, snps, sequences, straincluster);
	}
}

int main(int argc, char *argv[]) //input aln, input clustermap, threads, threshold, output
{
	string argv1=argv[1];
	string argv2=argv[2];
	string argv3=argv[3];
	string argv4=argv[4];
	string argv5=argv[5];

	ifstream aln(argv1.c_str());
	ifstream cluster(argv2.c_str());
	double thres=stod(argv4);
	ofstream output(argv5.c_str());

	string line;


	cout<<"Reading alignment... "<<endl;

	string presentgenome="";
	map<string, string> sequences;

	for(;getline(aln,line);)
	{
		if(line=="") continue;

		if(line[0]=='>')
		{
			istringstream linestream(line);
			string f1,f2;
			getline(linestream,f1,'>');
			getline(linestream,f2,'>');

			presentgenome=f2;
		}
		else
		{
			for(int i=0;i<line.size();i++)
  			{
	    			if(line[i]>='a' && line[i]<='z')
	    			{
					line[i]=line[i]-'a'+'A';
	    			}
  			}
			sequences[presentgenome].append(line);
		}
	}

	cout<<"Reading strain clusters... "<<endl;

	map<string, int> straincluster;
	map<int, double> clustersize;

	for(;getline(cluster,line);)
	{
		istringstream linestream(line);
		string strain,cluster;
		getline(linestream,strain,'\t');
		getline(linestream,cluster,'\t');
		straincluster[strain]=stoi(cluster);
		clustersize[stoi(cluster)]++;
	}


	int numsite=((sequences.begin())->second).size();

	cout<<"Total sites to screen: "<<numsite<<"..."<<endl;
	cout<<"Generating SNP objects..."<<endl;

	vector<SNP> snps;
	snps.resize(numsite);


	int p=stoi(argv[3]);
	
	int numthread;
	if(numsite < p)
		numthread=numsite;
	else
		numthread=p;

	int len=ceil((double)numsite/(double)numthread);
	thread t[numthread];

	vector<int> posvar;
	vector<int> posbeast;
	posvar.resize(numsite);
	posbeast.resize(numsite);

	for(int threadi=0; threadi<numthread; threadi++)
	{
		int starti=threadi*len;
		int endi=(threadi+1)*len-1;
		if(endi>=numsite)
			endi=numsite-1;
		t[threadi] = thread(checkSites, starti, endi, std::ref(snps), std::ref(sequences), std::ref(straincluster));
		//linei++;
	}

	for(int threadi=0; threadi<numthread; threadi++)
	{
		t[threadi].join();
	}
		
	cout<<"outputting..."<<endl;

	output<<"1_base_position\tnucleotide\tcluster_#\tpresent_in_#_strains"<<endl;
	for(int snpi=0; snpi<snps.size(); snpi++)
	{
		if(snps[snpi].isunique>0)
		{
			SNP snp=snps[snpi];
			for(auto ncluster : snp.nucleotidecluster)
			{
				if(ncluster.second>=0)
				{
					if(snp.nucleotidecount[ncluster.first]>=(double)thres*clustersize[ncluster.second])
						output<<snpi+1<<"\t"<<ncluster.first<<"\t"<<ncluster.second<<"\t"<<snp.nucleotidecount[ncluster.first]<<endl;
				}
			}
		}
	}
	
}
